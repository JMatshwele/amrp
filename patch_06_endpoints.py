#!/usr/bin/env python3
"""Patch 6 — endpoints.

Rewires both docking endpoints onto the v3.2 signatures and adds the
verification surface the platform never had: target checking against RCSB,
seed-replicate ensembles, an exhaustiveness ladder, redocking RMSD, and a
receptor provenance dump.

The pre-screening fallback is also hardened. In v3.1, when Vina was missing,
base_dg silently became a "calibrated literature midpoint" and the only
thing distinguishing an estimate from a docking run was a Data_Quality
column that disappears the moment anyone copies a table into a thesis. Now
the numeric field is null in that mode and the estimate travels under a
different key entirely, so a transcription cannot turn it into a score.
"""
from pathlib import Path

MAIN = Path(__file__).parent / "backend" / "main.py"
src = MAIN.read_text(encoding="utf-8")

# ── /health ─────────────────────────────────────────────────────────────
old_health_start = '@app.get("/health")\ndef health():'
old_health_end = '@app.get("/api/download_receptors")'
i = src.index(old_health_start)
j = src.index(old_health_end, i)

health = '''@app.get("/health")
def health():
    return {
        "status":    "BAMSTU Digital Pathogen Pro v3.2 running",
        "version":   "3.2", "year": 2026,
        "developer": "Dr. James P.T. Matshwele, BAMSTU Labs, University of Botswana",
        "ai":        "Groq connected" if groq_client else "AI disabled",
        "vina":      get_vina_status(),
        "chemistry": chem.status(),
        "engines":   ["AutoDock Vina 1.2", "BAMSTU composite screening layer v3.2"],
        "targets":   {p: len(t) for p, t in TARGETS.items()},
        "v32_changes": [
            "Descriptors from an RDKit molecular graph, or withheld. The v3.1 "
            "element-count HBA/HBD is gone.",
            "cLogP is real Wildman-Crippen on the organic scaffold, or null. "
            "The one-constant-per-element table that carried the "
            "Wildman-Crippen name is gone.",
            "Ki reporting removed. dg_to_ki() returns None.",
            "Cosmetic per-target offset and the +0.01 display nudge removed. "
            "Ties are reported as ties.",
            "Size bonuses per polar and heavy atom removed.",
            "Permeability decoupled from the score.",
            "--seed passed to Vina; runs are reproducible.",
            "Grid box built from active-site residues rather than a fixed 30 A cube.",
            "Contact analysis runs against the docking receptor, so no water "
            "contact can be reported for a pose docked in a dry pocket.",
            "Metal coordination detected with element-specific windows.",
            "LLM mechanistic claims filtered post-generation.",
            "MRSA panel rebuilt: PBP2a site corrected from the C-terminal tail "
            "to the transpeptidase site, DHPS organism corrected, ClfA renamed "
            "from fibronectin to fibrinogen binding.",
        ],
        "validation_state": (
            "No redocking RMSD, decoy enrichment or affinity correlation has "
            "been run on this installation yet. Until /api/vvuq/* returns "
            "results, treat every score as an untested prediction."
        ),
    }


'''
src = src[:i] + health + src[j:]

# ── docking endpoints ───────────────────────────────────────────────────
i = src.index('@app.post("/api/dock_vina")')
j = src.index('@app.post("/api/chat")')

NEW = '''def _target_record(pathogen: str, target_pdb: str) -> dict:
    for t in TARGETS.get(pathogen, []):
        if t["target_pdb"] == target_pdb:
            return t
    return {}


def _dock_one(atoms, descriptors, coord_fp, filename, pathogen, t,
              exhaustiveness, seed, work_dir, run_blind=True):
    """Dock one ligand against one target and assemble the full record."""
    receptor = prepare_receptor(t["target_pdb"], work_dir)
    if not receptor:
        return None

    box = compute_active_site_box(t["target_pdb"], t.get("active_residues", []))
    if box:
        cx, cy, cz, n_found = box[0], box[1], box[2], box[6]
    else:
        cx = cy = cz = None
        n_found = 0

    ligand = prepare_ligand(atoms, filename, work_dir, coord_fp=coord_fp)
    if ligand is None:
        return None

    active = run_vina(ligand, receptor, t["target_pdb"], cx=cx, cy=cy, cz=cz,
                      exhaustiveness=exhaustiveness, blind=False,
                      seed=seed, box=box)
    if not active.get("success"):
        return {"target_pdb": t["target_pdb"], "target_name": t["target_name"],
                "dock_status": "failed", "error": active.get("error")}

    blind = None
    if run_blind:
        blind = run_vina(ligand, receptor, t["target_pdb"],
                         cx=None, cy=None, cz=None,
                         exhaustiveness=max(8, exhaustiveness // 2),
                         blind=True, seed=seed)

    pose = generate_pose(atoms, t["target_pdb"], t.get("active_residues", []),
                         coord_fp=coord_fp)
    inter = contacts.analyse_pose(pose, receptor)
    n_metal = inter.get("n_metal_coordination", 0) if inter.get("ok") else 0

    composite = compute_composite_score(
        atoms, active["delta_g"], pathogen,
        coord_metal_contacts=n_metal, target_pdb=t["target_pdb"],
        descriptors=descriptors)

    selectivity = None
    if blind and blind.get("success"):
        selectivity = compute_selectivity(active["delta_g"], blind["delta_g"])

    return {
        "target_name":   t["target_name"],
        "target_pdb":    t["target_pdb"],
        "compartment":   t["compartment"],
        "organism":      t.get("organism", "not declared"),
        "provenance":    t.get("provenance", "not declared"),
        "docking_semantics": t.get("docking_semantics", "enzyme_inhibition"),
        "site_note":     t.get("site_note", ""),
        "target_verified": t.get("verified", True),
        "dock_status":   "ok",
        "vina_affinity_kcal_per_mol": active["delta_g"],
        "bamstu_composite_score":     composite["bamstu_composite_score"],
        "composite_units":            composite["units"],
        "composite_terms":            composite["terms"],
        "calibration_status":         composite["calibration_status"],
        "heavy_atoms":   descriptors.get("heavy_atoms"),
        "seed":          seed,
        "box":           active.get("box"),
        "metal_proxy":   active.get("metal_proxy"),
        "charge_model":  active.get("charge_model"),
        "active_site_residues_found": n_found,
        "interactions":  inter,
        "blind_affinity_kcal_per_mol": (blind.get("delta_g")
                                        if blind and blind.get("success") else None),
        "selectivity_analysis": selectivity,
        "docked_pose":   pose,
        "pose_format":   "pdb",
        "data_quality":  "REAL_DOCKING",
        "ki":            None,
        "ki_note":       ("Removed in v3.2. An inhibition constant is measured "
                          "by enzyme kinetics, not derived from a docking score."),
    }


@app.post("/api/dock_vina")
async def dock_vina_single(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    target_pdb: str = Form(...),
    target_name: str = Form("Unknown Target"),
    exhaustiveness: int = Form(16),
    seed: int = Form(42),
    n_replicates: int = Form(1),
):
    """Single-target docking with seed control and optional replicate UQ."""
    if not VINA_AVAILABLE:
        return JSONResponse({
            "vina_available": False,
            "message": "AutoDock Vina is required. Place vina.exe in the backend folder.",
            "download_url": "https://github.com/ccsb-scripps/AutoDock-Vina/releases",
        }, status_code=503)

    try:
        content = (await file.read()).decode("utf-8")
    except Exception as e:
        raise HTTPException(400, f"Cannot read .xyz file: {e}")

    try:
        atoms, descriptors, coord_fp = parse_xyz(content)
    except chem.ChemistryError as e:
        raise HTTPException(422, str(e))
    if not atoms:
        raise HTTPException(422, "No valid atomic coordinates in .xyz file.")

    t = _target_record(pathogen, target_pdb) or {
        "target_name": target_name, "target_pdb": target_pdb,
        "compartment": "Unknown", "active_residues": [], "base_dg": -6.0}

    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    stem = file.filename.split(".")[0]

    hit = _dock_one(atoms, descriptors, coord_fp, stem, pathogen, t,
                    exhaustiveness, seed, work_dir)
    if hit is None:
        return JSONResponse({"success": False,
                             "error": f"Could not prepare receptor {target_pdb}."})
    if hit.get("dock_status") != "ok":
        return JSONResponse({"success": False, "error": hit.get("error")})

    # Replicate ensemble: the band, not a point estimate.
    uq = None
    if n_replicates > 1:
        receptor = prepare_receptor(target_pdb, work_dir)
        ligand = prepare_ligand(atoms, stem, work_dir, coord_fp=coord_fp)
        box = compute_active_site_box(target_pdb, t.get("active_residues", []))

        def _run(seed=42, exhaustiveness=exhaustiveness):
            r = run_vina(ligand, receptor, target_pdb,
                         cx=box[0] if box else None,
                         cy=box[1] if box else None,
                         cz=box[2] if box else None,
                         exhaustiveness=exhaustiveness, seed=seed, box=box)
            if not r.get("success"):
                raise RuntimeError(r.get("error", "vina failed"))
            return r["delta_g"]

        uq = vvuq.ensemble_uq(_run, n_replicates=n_replicates,
                              exhaustiveness=exhaustiveness)

    adme = compute_adme(descriptors, pathogen)
    return JSONResponse({
        "vina_available": True, "success": True,
        "pathogen": pathogen, "exhaustiveness": exhaustiveness, "seed": seed,
        "hit": hit,
        "uncertainty": uq,
        "physicochemistry": adme,
        "engine": "AutoDock Vina 1.2 + BAMSTU composite screening layer v3.2",
        "reporting_rule": (
            "Report vina_affinity_kcal_per_mol as an AutoDock Vina 1.2 score, "
            "with the Zn metal proxy stated in Methods. The composite is a "
            "triage number in arbitrary units and is not a binding energy."
        ),
    })


@app.post("/api/dock_whole_cell")
async def dock_whole_cell(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    exhaustiveness: int = Form(16),
    seed: int = Form(42),
    include_unverified: bool = Form(False),
):
    """Whole-panel sweep, with non-inhibition targets held in separate tables."""
    try:
        content = (await file.read()).decode("utf-8")
    except Exception as e:
        raise HTTPException(400, f"Cannot read .xyz file: {e}")

    try:
        atoms, descriptors, coord_fp = parse_xyz(content)
    except chem.ChemistryError as e:
        raise HTTPException(422, str(e))
    if not atoms:
        raise HTTPException(422, "No valid atomic coordinates found.")

    adme = compute_adme(descriptors, pathogen)
    targets = TARGETS.get(pathogen, TARGETS["Staphylococcus aureus"])
    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    stem = file.filename.split(".")[0]

    skipped = [t["target_name"] for t in targets
               if t.get("verified") is False and not include_unverified]
    active_targets = [t for t in targets
                      if t.get("verified") is not False or include_unverified]

    if not VINA_AVAILABLE:
        # No numeric score is emitted in this mode. v3.1 filled base_dg with a
        # literature midpoint and labelled the row; a labelled number still
        # gets transcribed into a table as a number.
        return JSONResponse({
            "mode": "NO_DOCKING_ENGINE",
            "hits": [],
            "physicochemistry": adme,
            "message": (
                "AutoDock Vina is not installed, so no docking was performed "
                "and no score is reported. v3.1 substituted a literature "
                "midpoint here and returned it in the delta_g field with a "
                "quality label attached. That label does not survive being "
                "copied into a thesis table, so the field is now empty "
                "instead. Install Vina from "
                "https://github.com/ccsb-scripps/AutoDock-Vina/releases"
            ),
            "descriptors_still_valid": True,
        }, status_code=503)

    hits, failures = [], []
    for t in active_targets:
        rec = _dock_one(atoms, descriptors, coord_fp, stem, pathogen, t,
                        exhaustiveness, seed, work_dir)
        if rec is None:
            failures.append({"target": t["target_name"],
                             "reason": "receptor preparation failed"})
        elif rec.get("dock_status") != "ok":
            failures.append({"target": t["target_name"],
                             "reason": rec.get("error", "docking failed")})
        else:
            hits.append(rec)

    if not hits:
        raise HTTPException(500, "Every target failed. " + str(failures[:3]))

    # Separate the tables by what a negative score actually means. Ranking a
    # porin against a carbapenemase compares transit against inhibition.
    inhibition = [h for h in hits
                  if h["docking_semantics"] == "enzyme_inhibition"]
    permeability = [h for h in hits if h["docking_semantics"] == "permeability"]
    other = [h for h in hits
             if h["docking_semantics"] not in ("enzyme_inhibition", "permeability")]

    inhibition = rank_hits(inhibition)
    permeability = rank_hits(permeability)
    other = rank_hits(other)

    top = inhibition[0] if inhibition else hits[0]
    moa = get_ai_moa(
        {**adme, "mw": descriptors.get("mw"), "clogp": descriptors.get("clogp")},
        top, pathogen, atoms=atoms,
        interaction_summary=top.get("interactions", {}))

    return JSONResponse({
        "inhibition_targets":   inhibition,
        "permeability_targets": permeability,
        "other_semantics":      other,
        "Hits":                 inhibition,   # legacy key for the old UI
        "MoA":                  moa,
        "physicochemistry":     adme,
        "ADME":                 adme,          # legacy key
        "failures":             failures,
        "skipped_unverified":   skipped,
        "exhaustiveness":       exhaustiveness,
        "seed":                 seed,
        "n_atoms":              len(atoms),
        "n_metals":             descriptors.get("n_metals"),
        "engine": "AutoDock Vina 1.2 + BAMSTU composite screening layer v3.2",
        "data_quality": "REAL_DOCKING",
        "table_separation_note": (
            "Three tables, not one ranked list. A negative score against a "
            "porin means the compound fits the channel, which is transit, not "
            "inhibition. Against an adhesin it means interface occlusion, "
            "which predicts no MIC shift. Pooling them into one affinity "
            "ranking compares opposite mechanisms."
        ),
        "publication_note": (
            "AutoDock Vina 1.2 scores, rigid ligand and rigid receptor, all "
            "transition metals mapped to the Zn atom type, seed "
            f"{seed}, exhaustiveness {exhaustiveness}. No redocking RMSD or "
            "decoy enrichment has been run on this installation; run "
            "/api/vvuq/redock before citing these as validated predictions."
        ),
    })


# ════════════════════════════════════════════════════════════════════
#  VERIFICATION / VALIDATION / UQ ENDPOINTS
# ════════════════════════════════════════════════════════════════════

@app.get("/api/verify_targets")
async def verify_targets(pathogen: str = "all"):
    """Check every panel entry against the live RCSB record.

    Checks the deposited organism against the pathogen the panel claims, and
    confirms every configured active-site residue exists in the chain. A
    residue list that points at nothing makes the centroid fall back to the
    whole-protein centre, and the docking then quietly targets the middle of
    the structure while the output still names the active site.
    """
    panel = TARGETS if pathogen == "all" else {pathogen: TARGETS.get(pathogen, [])}
    return JSONResponse(vvuq.verify_panel(panel))


@app.get("/api/receptor_provenance")
async def receptor_provenance():
    """What was substituted or dropped during receptor preparation."""
    return JSONResponse({
        "receptors": _RECEPTOR_PREP_LOG,
        "note": ("Empty until receptors are prepared. Metals listed under "
                 "metals_mapped_to_zn were given zinc parameters because Vina "
                 "has none for them."),
    })


@app.post("/api/vvuq/ensemble")
async def vvuq_ensemble(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    target_pdb: str = Form(...),
    n_replicates: int = Form(5),
    exhaustiveness: int = Form(16),
):
    """Replicate the same docking across independent seeds; return the band."""
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required.")
    content = (await file.read()).decode("utf-8")
    atoms, descriptors, coord_fp = parse_xyz(content)
    t = _target_record(pathogen, target_pdb)
    work_dir = Path(__file__).parent / "data"
    receptor = prepare_receptor(target_pdb, work_dir)
    ligand = prepare_ligand(atoms, file.filename.split(".")[0], work_dir,
                            coord_fp=coord_fp)
    box = compute_active_site_box(target_pdb, t.get("active_residues", []))

    def _run(seed=42, exhaustiveness=exhaustiveness):
        r = run_vina(ligand, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exhaustiveness, seed=seed, box=box)
        if not r.get("success"):
            raise RuntimeError(r.get("error", "vina failed"))
        return r["delta_g"]

    return JSONResponse({
        "target": target_pdb,
        "uq": vvuq.ensemble_uq(_run, n_replicates=n_replicates,
                               exhaustiveness=exhaustiveness),
        "determinism": vvuq.check_determinism(_run, n_repeats=2, seed=42),
    })


@app.post("/api/vvuq/convergence")
async def vvuq_convergence(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    target_pdb: str = Form(...),
):
    """Run an exhaustiveness ladder and report whether the score has settled."""
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required.")
    content = (await file.read()).decode("utf-8")
    atoms, descriptors, coord_fp = parse_xyz(content)
    t = _target_record(pathogen, target_pdb)
    work_dir = Path(__file__).parent / "data"
    receptor = prepare_receptor(target_pdb, work_dir)
    ligand = prepare_ligand(atoms, file.filename.split(".")[0], work_dir,
                            coord_fp=coord_fp)
    box = compute_active_site_box(target_pdb, t.get("active_residues", []))

    def _dock(exhaustiveness=16, seed=42):
        r = run_vina(ligand, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exhaustiveness, seed=seed, box=box)
        if not r.get("success"):
            raise RuntimeError(r.get("error", "vina failed"))
        return r["delta_g"]

    return JSONResponse({"target": target_pdb,
                         "convergence": vvuq.convergence_ladder(_dock)})


@app.post("/api/vvuq/redock")
async def vvuq_redock(target_pdb: str = Form(...),
                      ligand_resname: str = Form(...),
                      exhaustiveness: int = Form(32),
                      seed: int = Form(42)):
    """Extract the native co-crystallised ligand, redock it, report RMSD.

    The single most informative validation for a docking setup, and the one
    v3.1 never ran. It tests the box, the receptor preparation and the
    scoring function together. If the software cannot reproduce a pose that
    is already known, its predictions for unknown ligands carry no weight.
    """
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required.")
    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        crystal = _fetch_pdb(target_pdb)
    except Exception as e:
        raise HTTPException(502, f"Could not fetch {target_pdb}: {e}")

    native_lines = [ln for ln in crystal.read_text(errors="ignore").splitlines()
                    if ln.startswith("HETATM")
                    and ln[17:20].strip().upper() == ligand_resname.upper()]
    if not native_lines:
        raise HTTPException(
            404,
            f"Ligand {ligand_resname} not present in {target_pdb}. List the "
            "HETATM residue names in the structure and pick the one that "
            "occupies the pocket you are docking into.")

    native_path = work_dir / f"{target_pdb}_{ligand_resname}_native.pdb"
    native_path.write_text("\\n".join(native_lines) + "\\nEND\\n")

    atoms = []
    for ln in native_lines:
        elem = ln[76:78].strip() or ln[12:16].strip()[:1]
        elem = elem[:1].upper() + elem[1:].lower() if len(elem) > 1 else elem.upper()
        if elem == "H":
            continue
        atoms.append({"elem": elem, "x": float(ln[30:38]),
                      "y": float(ln[38:46]), "z": float(ln[46:54])})

    receptor = prepare_receptor(target_pdb, work_dir)
    ligand = prepare_ligand(atoms, f"{target_pdb}_{ligand_resname}_redock",
                            work_dir, coord_fp=chem.coordinate_fingerprint(atoms))
    t = next((x for p in TARGETS.values() for x in p
              if x["target_pdb"] == target_pdb), {})
    box = compute_active_site_box(target_pdb, t.get("active_residues", []))

    r = run_vina(ligand, receptor, target_pdb,
                 cx=box[0] if box else None, cy=box[1] if box else None,
                 cz=box[2] if box else None,
                 exhaustiveness=exhaustiveness, seed=seed, box=box)
    if not r.get("success"):
        raise HTTPException(500, f"Redock failed: {r.get('error')}")

    out_pose = Path(r.get("out_pdbqt", "")) if r.get("out_pdbqt") else None
    if not out_pose or not out_pose.exists():
        return JSONResponse({
            "status": "NOT_RUN",
            "reason": "Vina output pose file not located; cannot compute RMSD.",
            "vina_score": r.get("delta_g"),
        })

    result = vvuq.redock_validation(target_pdb, ligand_resname,
                                    native_path, out_pose)
    result["vina_score_of_native"] = r.get("delta_g")
    result["box"] = r.get("box")
    return JSONResponse(result)


@app.get("/api/vvuq/status")
async def vvuq_status():
    """What has and has not been validated on this installation."""
    return JSONResponse({
        "pass_bars": vvuq.PASS_BARS,
        "verification": {
            "seed_control": "IMPLEMENTED — --seed passed to every Vina call",
            "determinism_test": "available at /api/vvuq/ensemble",
            "exhaustiveness_convergence": "available at /api/vvuq/convergence",
        },
        "validation": {
            "redock_rmsd": "available at /api/vvuq/redock — NOT YET RUN",
            "decoy_enrichment": "NOT IMPLEMENTED — needs a DUD-E style decoy set",
            "affinity_correlation": (
                "NOT IMPLEMENTED — needs 8-10 congeneric complexes with "
                "measured affinities. This is the test that would let the "
                "composite score be calibrated instead of asserted."
            ),
            "target_panel": "available at /api/verify_targets",
        },
        "uq": {
            "search_uncertainty": "available at /api/vvuq/ensemble",
            "scoring_function_uncertainty": (
                "NOT PROPAGATED — approximately 1.5-2.0 kcal/mol for empirical "
                "functions, epistemic, and larger than the search spread. Any "
                "reported band that omits it understates the real uncertainty."
            ),
        },
        "honest_summary": (
            "The engine is now reproducible and its provenance is visible. It "
            "is not yet validated. Reproducibility and correctness are "
            "different claims and this build has only earned the first."
        ),
    })


'''

src = src[:i] + NEW + src[j:]
MAIN.write_text(src, encoding="utf-8")
import ast
ast.parse(src)
print("patch 06 applied: endpoints rewired, VVUQ API added")
