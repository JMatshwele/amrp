#!/usr/bin/env python3
"""Patch 11 — v3.3.1 specs 1 and 3, plus split metal-physics reporting.

  1. Figures are rendered on save for every hit where pose_is_docked is
     true, and land in the session folder, so the ZIP picks them up without
     any change to the archive code.

  2. The metal-coordination term is now reported THREE ways, per request:
       delta_e_metal                   the correction on its own
       vina_affinity_kcal_per_mol      the untouched Vina score
       corrected_score_kcal_per_mol    the two summed
     A reader can take whichever the argument needs and can always recover
     the other two. This is the arrangement that survives review: the term
     is auditable because it is printed separately, and the sum is available
     without anyone having to reconstruct it by hand. What it must not be
     called is a binding free energy — the term is uncalibrated, so the sum
     inherits that and says so in its own units field.

  3. POST /api/vvuq/finalize-audit?code=bamstu-valid-jimmah walks the ledger
     and signs it, or refuses and names the sample that blocked it.
"""
from pathlib import Path

MAIN = Path(__file__).resolve().parent.parent / "backend" / "main.py"
src = MAIN.read_text(encoding="utf-8")
n = 0


def sub(old, new, label):
    global src, n
    if old not in src:
        raise SystemExit(f"ANCHOR MISS [{label}]")
    if src.count(old) != 1:
        raise SystemExit(f"ANCHOR AMBIGUOUS [{label}]: {src.count(old)}")
    src = src.replace(old, new)
    n += 1
    print(f"  ok  {label}")


# ── 1. import the renderer ──────────────────────────────────────────────
sub("import bamstu_vvuq as vvuq\n",
    "import bamstu_vvuq as vvuq\n"
    "import bamstu_render as figures\n",
    "import renderer")

# ── 2. ΔE_metal reported separately and summed ──────────────────────────
sub('''    return {
        "bamstu_composite_score": round(score, 2),
        "units": "arbitrary (screening units) — NOT kcal/mol, NOT a free energy",
        "vina_affinity_kcal_per_mol": round(float(vina_affinity), 2),''',
    '''    # ── Three-way reporting of the metal physics term ───────────────
    #
    # delta_e_metal is the coordination correction ALONE. It is kept
    # separate so a reader can see exactly what was added and can strip it
    # out, and so the term can be argued about on its own merits without
    # having to unpick it from a total.
    #
    # corrected_score_kcal_per_mol is the same term summed onto the Vina
    # score, which is what the E_total = E_vina + dE_metal formulation asks
    # for. It excludes the torsional term, which is a ligand property rather
    # than a metal-interaction term and would muddy the comparison.
    #
    # Neither is a binding free energy. dE_metal is uncalibrated -- its
    # coefficients have never been fitted against measured affinities -- so
    # the sum inherits that status. Run /api/vvuq/calibrate before attaching
    # thermodynamic language to either.
    delta_e_metal = round(-coord_total, 3)
    corrected = round(float(vina_affinity) + delta_e_metal, 2)

    return {
        "bamstu_composite_score": round(score, 2),
        "units": "arbitrary (screening units) — NOT kcal/mol, NOT a free energy",
        "vina_affinity_kcal_per_mol": round(float(vina_affinity), 2),
        "delta_e_metal": delta_e_metal,
        "delta_e_metal_units": (
            "kcal/mol dimensionally, but UNCALIBRATED — no fit against any "
            "measured affinity set"
        ),
        "corrected_score_kcal_per_mol": corrected,
        "corrected_score_definition": "E_total = E_vina + dE_metal",
        "corrected_score_caveat": (
            "Reported because the sum is what the physics-correction "
            "formulation asks for, and because a reader who wants it should "
            "not have to reconstruct it. It is not a measured binding energy "
            "and must not be labelled dG. The Vina score and the correction "
            "are both printed above so either can be recovered."
        ),''',
    "delta_e_metal separate + summed")

sub('''        "vina_affinity_kcal_per_mol": active["delta_g"],
        "bamstu_composite_score":     composite["bamstu_composite_score"],''',
    '''        "vina_affinity_kcal_per_mol": active["delta_g"],
        "delta_e_metal":              composite["delta_e_metal"],
        "corrected_score_kcal_per_mol": composite["corrected_score_kcal_per_mol"],
        "bamstu_composite_score":     composite["bamstu_composite_score"],''',
    "hit carries all three numbers")

# ── 3. render figures when a sample is saved ────────────────────────────
sub('''    manifest_path = sdir / "manifest.json"''',
    '''    # ── Figures ─────────────────────────────────────────────────────
    #
    # Rendered here rather than at ZIP time, so they exist on disk the moment
    # the sample finishes. An overnight batch that is interrupted still has
    # figures for everything that completed, and the archive builder needs no
    # special case because it simply zips the folder.
    figures_dir = sample_dir / "figures"
    fig_report = []
    for h in hits:
        try:
            rec_path = RECEPTOR_STORE / f"{h.get('target_pdb')}.pdbqt"
            if not rec_path.exists():
                rec_path = RECEPTOR_STORE / f"{h.get('target_pdb')}.pdb"
            r = figures.render_hit_figures(sample, h, rec_path, figures_dir)
            fig_report.append(r)
        except Exception as exc:
            fig_report.append({"sample": sample,
                               "target_pdb": h.get("target_pdb"),
                               "errors": [f"{type(exc).__name__}: {exc}"]})
    if fig_report:
        (sample_dir / "figures_manifest.json").write_text(
            json.dumps(fig_report, indent=2), encoding="utf-8")

    manifest_path = sdir / "manifest.json"''',
    "render figures on save")

sub('''        "n_poses_written": n_poses,''',
    '''        "n_poses_written": n_poses,
        "figures": {
            "rendered": sum(1 for f in fig_report if f.get("png") or f.get("svg")),
            "skipped_not_docked": sum(1 for f in fig_report if f.get("skipped")),
            "errors": [e for f in fig_report for e in (f.get("errors") or [])],
        },''',
    "figure summary in the save response")

# ── 4. serve individual figures to the UI ───────────────────────────────
sub('@app.get("/api/vvuq/status")',
    '''@app.get("/api/session/{session_id}/figure")
async def session_figure(session_id: str, sample: str, target_pdb: str,
                         kind: str = "png"):
    """Serve one rendered figure so the visualiser can show it inline."""
    sdir = SESSION_STORE / _safe_id(session_id) / _safe_id(sample) / "figures"
    tag = f"{_safe_id(sample)}_{_safe_id(target_pdb)}"
    name = f"{tag}_3D_Pose.png" if kind == "png" else f"{tag}_LigPlot_2D.svg"
    f = sdir / name
    if not f.is_file():
        raise HTTPException(404, f"No {kind} figure for {sample}/{target_pdb}. "
                                 "Figures are only rendered when "
                                 "pose_is_docked is true.")
    return FileResponse(str(f),
                        media_type="image/png" if kind == "png" else "image/svg+xml")


@app.get("/api/session/{session_id}/pose")
async def session_pose(session_id: str, sample: str, target_pdb: str):
    """Serve one docked pose as PDB text for the 3D canvas."""
    pdir = SESSION_STORE / _safe_id(session_id) / _safe_id(sample) / "poses"
    if not pdir.is_dir():
        raise HTTPException(404, "No poses stored for that sample.")
    pref = f"{_safe_id(sample)}_{_safe_id(target_pdb)}_"
    for f in sorted(pdir.iterdir()):
        if f.name.startswith(pref):
            return PlainTextResponse(f.read_text(encoding="utf-8",
                                                 errors="ignore"))
    raise HTTPException(404, f"No pose for {sample}/{target_pdb}.")


# ════════════════════════════════════════════════════════════════════
#  AUDIT SIGN-OFF  —  /bamstu-valid-jimmah
# ════════════════════════════════════════════════════════════════════

_AUDIT_CODE = "bamstu-valid-jimmah"


def _sha256_file(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _run_audit(session_id: str) -> dict:
    """Walk one session ledger and check every sample.

    WHAT THIS AUDIT DOES AND DOES NOT CERTIFY
    -----------------------------------------
    Everything checked here is RUN INTEGRITY: did each sample dock, is the
    pose a real one, do the figures exist, are the semantic tables kept
    apart, does any narrative claim metal coordination the contact table
    does not show.

    None of it is SCIENTIFIC VALIDATION. A run can pass every check below
    and still be built on a docking setup that cannot reproduce a
    crystallographic pose. Those are separate questions and the ledger keeps
    them separate, because a document headed PASSED that quietly conflates
    the two is precisely how an internal completeness check ends up cited in
    a thesis as evidence the predictions are sound.
    """
    from datetime import datetime as _dt

    sdir = SESSION_STORE / _safe_id(session_id)
    mf = sdir / "manifest.json"
    if not mf.is_file():
        raise HTTPException(404, f"No session '{session_id}' on disk.")
    manifest = json.loads(mf.read_text(encoding="utf-8"))

    checks, failures, samples = [], [], []
    seeds, exhaustivenesses = set(), set()

    for entry in manifest.get("samples", []):
        name = entry["name"]
        rp = sdir / name / f"{name}_result.json"
        if not rp.is_file():
            failures.append({"sample": name, "check": "result_present",
                             "detail": "result JSON missing from disk"})
            continue
        try:
            data = json.loads(rp.read_text(encoding="utf-8"))
        except Exception as exc:
            failures.append({"sample": name, "check": "result_readable",
                             "detail": str(exc)})
            continue

        seeds.add(data.get("seed"))
        exhaustivenesses.add(data.get("exhaustiveness"))

        inhib = data.get("inhibition_targets") or []
        perm = data.get("permeability_targets") or []
        other = data.get("other_semantics") or []
        hits = inhib + perm + other
        srec = {"sample": name, "n_targets": len(hits), "issues": []}

        if not hits:
            failures.append({"sample": name, "check": "has_results",
                             "detail": "no docked targets in the record"})
            srec["issues"].append("no targets")

        for h in hits:
            pdb = h.get("target_pdb", "?")

            # 1. Vina score present and a real pose.
            if h.get("vina_affinity_kcal_per_mol") is None:
                failures.append({"sample": name, "check": "vina_score",
                                 "target": pdb, "detail": "no Vina score"})
                srec["issues"].append(f"{pdb}: no score")
            if not h.get("pose_is_docked"):
                failures.append({
                    "sample": name, "check": "pose_is_docked", "target": pdb,
                    "detail": ("stored geometry is a viewer placement, not a "
                               "Vina pose")})
                srec["issues"].append(f"{pdb}: not a docked pose")
                continue   # figures are intentionally absent for these

            # 2. Both figures exist and are non-trivial.
            fdir = sdir / name / "figures"
            png = fdir / f"{name}_{pdb}_3D_Pose.png"
            svg = fdir / f"{name}_{pdb}_LigPlot_2D.svg"
            for f, kind, floor in ((png, "3D pose PNG", 2048),
                                   (svg, "2D LigPlot SVG", 512)):
                if not f.is_file():
                    failures.append({"sample": name, "check": "figure_present",
                                     "target": pdb,
                                     "detail": f"{kind} missing: {f.name}"})
                    srec["issues"].append(f"{pdb}: {kind} missing")
                elif f.stat().st_size < floor:
                    failures.append({
                        "sample": name, "check": "figure_non_empty",
                        "target": pdb,
                        "detail": (f"{kind} is {f.stat().st_size} B, below the "
                                   f"{floor} B floor — a file-exists check "
                                   "passes on a blank image")})
                    srec["issues"].append(f"{pdb}: {kind} blank")

            # 3. No metal assertion where none was detected.
            inter = h.get("interactions") or {}
            n_metal = inter.get("n_metal_coordination", 0)
            if n_metal == 0:
                dem = h.get("delta_e_metal")
                if dem not in (None, 0, 0.0):
                    failures.append({
                        "sample": name, "check": "metal_term_gated",
                        "target": pdb,
                        "detail": (f"delta_e_metal = {dem} with zero detected "
                                   "coordination contacts")})
                    srec["issues"].append(f"{pdb}: ungated metal term")
                for c in (inter.get("metal_coordination") or []):
                    failures.append({
                        "sample": name, "check": "metal_table_consistent",
                        "target": pdb,
                        "detail": f"contact listed while count is 0: {c}"})

        # 3b. Narrative must not assert coordination the tables deny.
        moa = data.get("MoA") or {}
        text = moa.get("inference") or ""
        top_metal = max((((h.get("interactions") or {})
                          .get("n_metal_coordination", 0)) for h in hits),
                        default=0)
        if text and top_metal == 0:
            chk = main_claim_filter(text, 0)
            if not chk["clean"]:
                failures.append({
                    "sample": name, "check": "moa_claims",
                    "detail": (f"{chk['n_violations']} unsupportable claim(s) "
                               "in released MoA text")})
                srec["issues"].append("MoA claims")

        # 4. Semantic tables must not be pooled.
        wrong = [h.get("target_pdb") for h in inhib
                 if h.get("docking_semantics") != "enzyme_inhibition"]
        if wrong:
            failures.append({
                "sample": name, "check": "semantic_separation",
                "detail": (f"non-inhibition target(s) {wrong} inside the "
                           "inhibition table — a porin score is transit, not "
                           "inhibition, and pooling them compares opposite "
                           "mechanisms")})
            srec["issues"].append("semantics pooled")

        samples.append(srec)

    checks = [
        {"id": 1, "check": "Vina score present and pose_is_docked",
         "failures": sum(1 for f in failures
                         if f["check"] in ("vina_score", "pose_is_docked"))},
        {"id": 2, "check": "3D PNG and 2D SVG exist and exceed size floor",
         "failures": sum(1 for f in failures
                         if f["check"].startswith("figure_"))},
        {"id": 3, "check": "no metal assertion where contacts are zero",
         "failures": sum(1 for f in failures
                         if f["check"] in ("metal_term_gated",
                                           "metal_table_consistent",
                                           "moa_claims"))},
        {"id": 4, "check": "semantic tables not pooled",
         "failures": sum(1 for f in failures
                         if f["check"] == "semantic_separation")},
    ]

    passed = not failures and bool(samples)

    # Scientific validation is a separate axis and is reported as such.
    validation = {
        "redock_rmsd": "NOT RECORDED in this session",
        "decoy_enrichment": "NOT RECORDED in this session",
        "affinity_calibration": "NOT RECORDED in this session",
        "meaning": (
            "This audit certifies RUN INTEGRITY only: every sample docked, "
            "every pose is a real Vina pose, every figure exists, no "
            "narrative claims coordination the tables deny, and the semantic "
            "tables are kept apart. It does NOT certify that the predictions "
            "are correct. That requires redocking RMSD under 2.0 A, decoy "
            "enrichment, and a cross-validated calibration, none of which is "
            "checked here."
        ),
    }

    ledger = {
        "audit_status": ("PASSED_BAMSTU_VALID_JIMMAH" if passed
                         else "FAILED_BAMSTU_VALID_JIMMAH"),
        "audit_scope": "RUN_INTEGRITY_ONLY_NOT_SCIENTIFIC_VALIDATION",
        "session_id": session_id,
        "timestamp": _dt.now().isoformat(timespec="seconds"),
        "platform_version": "3.3.1",
        "total_samples_verified": len(samples),
        "total_targets_verified": sum(s["n_targets"] for s in samples),
        "rng_seeds": sorted(s for s in seeds if s is not None),
        "exhaustiveness": sorted(e for e in exhaustivenesses if e is not None),
        "checks": checks,
        "samples": samples,
        "failures": failures,
        "blocking_sample_ids": sorted({f["sample"] for f in failures
                                       if "sample" in f}),
        "scientific_validation": validation,
        "signature_note": (
            "The digest below is a SHA-256 content hash over this ledger and "
            "every figure it certifies. It detects tampering and file "
            "substitution. It is not a cryptographic signature and proves "
            "nothing about who produced the data — calling a bare hash a "
            "signature would overstate it."
        ),
    }

    digests = {}
    for entry in manifest.get("samples", []):
        fdir = sdir / entry["name"] / "figures"
        if fdir.is_dir():
            for f in sorted(fdir.iterdir()):
                if f.is_file():
                    digests[f"{entry['name']}/figures/{f.name}"] = _sha256_file(f)
    ledger["file_digests_sha256"] = digests

    import hashlib
    payload = json.dumps(ledger, sort_keys=True, separators=(",", ":"))
    ledger["ledger_digest_sha256"] = hashlib.sha256(
        payload.encode("utf-8")).hexdigest()

    out = sdir / "BAMSTU_Audit_Ledger_Final.json"
    out.write_text(json.dumps(ledger, indent=2), encoding="utf-8")
    ledger["ledger_path"] = str(out)
    return ledger


@app.post("/api/vvuq/finalize-audit")
async def finalize_audit(session_id: str, code: str = ""):
    """Sign off a session, or refuse and name what blocked it."""
    if code != _AUDIT_CODE:
        raise HTTPException(
            403, f"Audit requires code={_AUDIT_CODE}.")
    ledger = _run_audit(session_id)
    return JSONResponse(ledger, status_code=200 if
                        ledger["audit_status"].startswith("PASSED") else 409)


@app.get("/bamstu-valid-jimmah")
async def bamstu_valid_jimmah(session_id: str):
    """Convenience alias so the audit can be triggered from a browser or curl."""
    ledger = _run_audit(session_id)
    return JSONResponse(ledger, status_code=200 if
                        ledger["audit_status"].startswith("PASSED") else 409)


@app.get("/api/vvuq/status")''',
    "figure/pose serving + audit endpoints")

MAIN.write_text(src, encoding="utf-8")
import ast
ast.parse(src)
print(f"patch 11 applied: {n} edits")
