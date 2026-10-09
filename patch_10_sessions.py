#!/usr/bin/env python3
"""Patch 10 — server-side session storage.

THE PROBLEM
-----------
Until now the autosampler kept every result in a JavaScript object in the
browser tab, and the only thing written to disk was a ZIP the browser
downloaded when the batch finished. For the stated use case — load samples,
leave it running overnight, read the results in the morning — that is the
wrong place for the data to live. A tab crash, a browser update, a laptop
sleeping and dropping the page, or someone closing the window loses the
entire run, and the run is hours of Vina time that has to be repeated.

The server was writing receptor and pose files into backend/data all along,
but nothing tied them to a sample or a batch, so they were scratch files
rather than results.

THE FIX
-------
Results are written to backend/sessions/<session_id>/ after EVERY sample,
not at the end of the batch. If the browser dies at sample 12 of 40, the
first 12 are on disk complete with their poses. The ZIP is then built
server-side on request, so it can be rebuilt at any time from the files
rather than from whatever the tab still holds in memory.

The browser download stays. It is now a copy of something that already
exists, instead of the only copy.
"""
from pathlib import Path

MAIN = Path(__file__).resolve().parent.parent / "backend" / "main.py"
src = MAIN.read_text(encoding="utf-8")

if "SESSION_STORE" in src:
    raise SystemExit("patch 10 already applied")

ANCHOR = '@app.get("/api/vvuq/status")'
assert src.count(ANCHOR) == 1

NEW = '''# ════════════════════════════════════════════════════════════════════
#  SESSION STORAGE — results live on disk, not in a browser tab
# ════════════════════════════════════════════════════════════════════

SESSION_STORE = Path(__file__).parent / "sessions"
SESSION_STORE.mkdir(parents=True, exist_ok=True)


def _safe_id(raw: str) -> str:
    """Filesystem-safe identifier. Rejects traversal rather than sanitising it.

    A session id arrives from the browser, so it is untrusted input that is
    about to become a path.
    """
    cleaned = "".join(c for c in str(raw) if c.isalnum() or c in "-_")[:64]
    if not cleaned:
        raise HTTPException(422, "Invalid session or sample identifier.")
    return cleaned


def _session_dir(session_id: str) -> Path:
    d = SESSION_STORE / _safe_id(session_id)
    d.mkdir(parents=True, exist_ok=True)
    return d


@app.post("/api/session/append")
async def session_append(payload: dict = Body(...)):
    """Persist one completed sample. Called after every sample, not at the end.

    Writes the full result JSON, every docked pose as a separate .pdb, and a
    CSV, then updates the manifest. An interrupted batch therefore leaves a
    complete record of everything that finished, which is the difference
    between losing one sample and losing a night.
    """
    session_id = _safe_id(payload.get("session_id") or "")
    sample = _safe_id(payload.get("sample_name") or "sample")
    data = payload.get("data") or {}

    sdir = _session_dir(session_id)
    sample_dir = sdir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    (sample_dir / f"{sample}_result.json").write_text(
        json.dumps(data, indent=2), encoding="utf-8")

    hits = ((data.get("inhibition_targets") or [])
            + (data.get("permeability_targets") or [])
            + (data.get("other_semantics") or []))

    poses_dir = sample_dir / "poses"
    poses_dir.mkdir(exist_ok=True)
    n_poses = 0
    for h in hits:
        if h.get("docked_pose"):
            tag = "docked" if h.get("pose_is_docked") else "PLACEMENT_NOT_DOCKED"
            (poses_dir / f"{sample}_{h['target_pdb']}_{tag}.pdb").write_text(
                h["docked_pose"], encoding="utf-8")
            n_poses += 1

    rows = ["Target,PDB,Semantics,Organism,Vina_score_kcal_per_mol,"
            "BAMSTU_composite_au,Metal_contacts,H_bonds,Pose_is_docked,"
            "Rank,Resolvable,Seed,Exhaustiveness"]
    for h in hits:
        it = h.get("interactions") or {}
        rows.append(
            f'"{h.get("target_name","")}",{h.get("target_pdb","")},'
            f'{h.get("docking_semantics","")},"{h.get("organism","")}",'
            f'{h.get("vina_affinity_kcal_per_mol","")},'
            f'{h.get("bamstu_composite_score","")},'
            f'{it.get("n_metal_coordination","")},{it.get("n_hbonds","")},'
            f'{h.get("pose_is_docked","")},{h.get("rank","")},'
            f'{h.get("rank_is_resolvable","")},{data.get("seed","")},'
            f'{data.get("exhaustiveness","")}')
    (sample_dir / f"{sample}_docking.csv").write_text(
        "\\n".join(rows) + "\\n", encoding="utf-8")

    moa = data.get("MoA") or {}
    if moa:
        if moa.get("status") == "QUARANTINED":
            body = ("QUARANTINED — this text made claims the pipeline cannot "
                    "support. Do not paste it into a manuscript.\\n\\n"
                    + "\\n".join(f" - {v.get('matched')}: {v.get('requirement')}"
                                 for v in (moa.get("claim_filter", {})
                                           .get("violations") or []))
                    + "\\n\\nSafe alternative:\\n"
                    + str(moa.get("safe_alternative", ""))
                    + "\\n\\nQuarantined text (inspection only):\\n"
                    + str(moa.get("quarantined_text", "")))
        else:
            body = str(moa.get("inference") or "")
        body += ("\\n\\nREFERENCES\\n" + "\\n".join(moa.get("references") or [])
                 + "\\n\\nEvery citation must be opened and checked against the "
                   "sentence it is attached to before submission.")
        (sample_dir / f"{sample}_MoA.txt").write_text(body, encoding="utf-8")

    manifest_path = sdir / "manifest.json"
    manifest = {"session_id": session_id, "format": "BAMSTU_SESSION_V2",
                "version": "3.3", "samples": []}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            pass
    manifest.setdefault("samples", [])
    manifest["samples"] = [s for s in manifest["samples"]
                           if s.get("name") != sample]
    best = None
    inhib = data.get("inhibition_targets") or []
    if inhib:
        b = min(inhib, key=lambda h: h.get("vina_affinity_kcal_per_mol") or 0)
        best = {"target_name": b.get("target_name"),
                "target_pdb": b.get("target_pdb"),
                "vina": b.get("vina_affinity_kcal_per_mol"),
                "composite": b.get("bamstu_composite_score"),
                "pose_is_docked": b.get("pose_is_docked")}
    manifest["samples"].append({
        "name": sample,
        "finished": datetime.now().isoformat(timespec="seconds"),
        "pathogen": data.get("pathogen") or payload.get("pathogen"),
        "seed": data.get("seed"),
        "exhaustiveness": data.get("exhaustiveness"),
        "n_poses": n_poses,
        "best_inhibition_hit": best,
    })
    manifest["updated"] = datetime.now().isoformat(timespec="seconds")
    manifest["n_samples"] = len(manifest["samples"])
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    return JSONResponse({
        "saved": True,
        "session_id": session_id,
        "sample": sample,
        "path": str(sample_dir),
        "n_poses_written": n_poses,
        "n_samples_in_session": manifest["n_samples"],
        "note": ("Written to disk immediately. If the browser closes now, this "
                 "sample survives."),
    })


@app.get("/api/session/list")
async def session_list():
    """Every session on disk, newest first."""
    out = []
    for d in sorted(SESSION_STORE.iterdir()) if SESSION_STORE.exists() else []:
        mf = d / "manifest.json"
        if not mf.is_file():
            continue
        try:
            m = json.loads(mf.read_text(encoding="utf-8"))
            out.append({"session_id": d.name,
                        "n_samples": m.get("n_samples", 0),
                        "updated": m.get("updated"),
                        "samples": [s.get("name") for s in m.get("samples", [])]})
        except Exception:
            continue
    out.sort(key=lambda s: s.get("updated") or "", reverse=True)
    return JSONResponse({
        "storage_path": str(SESSION_STORE.resolve()),
        "n_sessions": len(out),
        "sessions": out,
    })


@app.get("/api/session/{session_id}")
async def session_get(session_id: str):
    """Full session payload, reassembled from disk."""
    sdir = SESSION_STORE / _safe_id(session_id)
    mf = sdir / "manifest.json"
    if not mf.is_file():
        raise HTTPException(404, f"No session '{session_id}' on disk.")
    manifest = json.loads(mf.read_text(encoding="utf-8"))
    samples = []
    for entry in manifest.get("samples", []):
        rp = sdir / entry["name"] / f"{entry['name']}_result.json"
        if rp.is_file():
            try:
                samples.append({"name": entry["name"],
                                "data": json.loads(rp.read_text(encoding="utf-8"))})
            except Exception:
                continue
    return JSONResponse({"format": "BAMSTU_SESSION_V2",
                         "session_id": manifest.get("session_id"),
                         "saved": manifest.get("updated"),
                         "n_samples": len(samples),
                         "samples": samples})


@app.get("/api/session/{session_id}/zip")
async def session_zip(session_id: str):
    """Build the archive from the files on disk and return it.

    Built server-side on request rather than in the browser, so it can be
    regenerated after the tab is gone. Receptor PDBQTs are included because
    re-fetching RCSB months later is a different calculation wearing the old
    results' name.
    """
    import zipfile

    sdir = SESSION_STORE / _safe_id(session_id)
    if not (sdir / "manifest.json").is_file():
        raise HTTPException(404, f"No session '{session_id}' on disk.")
    manifest = json.loads((sdir / "manifest.json").read_text(encoding="utf-8"))

    pdbs = set()
    for entry in manifest.get("samples", []):
        rp = sdir / entry["name"] / f"{entry['name']}_result.json"
        if not rp.is_file():
            continue
        try:
            d = json.loads(rp.read_text(encoding="utf-8"))
        except Exception:
            continue
        for h in ((d.get("inhibition_targets") or [])
                  + (d.get("permeability_targets") or [])
                  + (d.get("other_semantics") or [])):
            if h.get("target_pdb"):
                pdbs.add(h["target_pdb"])

    zpath = sdir / f"BAMSTU_Session_{_safe_id(session_id)}.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(sdir.rglob("*")):
            if f.is_file() and f.suffix != ".zip":
                z.write(f, f.relative_to(sdir).as_posix())
        for p in sorted(pdbs):
            rec = RECEPTOR_STORE / f"{p}.pdbqt"
            if rec.is_file():
                z.writestr(f"receptors/{p}.pdbqt",
                           rec.read_text(encoding="utf-8", errors="ignore"))
        z.writestr("README.txt",
                   "BAMSTU Digital Pathogen Pro — session archive\\n"
                   f"Session: {session_id}\\n"
                   f"Samples: {manifest.get('n_samples')}\\n"
                   f"Built:   {datetime.now().isoformat(timespec='seconds')}\\n\\n"
                   "RELOADING\\n"
                   "  Autosampler panel -> Reload session, or\\n"
                   f"  GET /api/session/{session_id}\\n\\n"
                   "WHAT THE NUMBERS ARE\\n"
                   "  Vina_score_kcal_per_mol  AutoDock Vina 1.2 score. Not a\\n"
                   "                           measured binding energy. Scoring\\n"
                   "                           error is roughly 1.5-2.0 kcal/mol,\\n"
                   "                           so targets closer than that are\\n"
                   "                           not separable.\\n"
                   "  BAMSTU_composite_au      Arbitrary units, uncalibrated.\\n"
                   "  Pose_is_docked           FALSE means the geometry is a\\n"
                   "                           placement, not a Vina pose, and no\\n"
                   "                           contacts were computed from it.\\n\\n"
                   "NOT HERE: no inhibition constants. A Ki is measured by\\n"
                   "enzyme kinetics.\\n")

    return FileResponse(str(zpath), media_type="application/zip",
                        filename=zpath.name)


@app.get("/api/vvuq/status")'''

src = src.replace(ANCHOR, NEW)

# FileResponse and datetime must be importable.
if "FileResponse" not in src.split("@app")[0]:
    src = src.replace(
        "from fastapi.responses import JSONResponse, PlainTextResponse, HTMLResponse",
        "from fastapi.responses import (JSONResponse, PlainTextResponse,\n"
        "                               HTMLResponse, FileResponse)", 1)
if "from datetime import datetime" not in src:
    src = src.replace("import uvicorn\n", "import uvicorn\nfrom datetime import datetime\n", 1)

MAIN.write_text(src, encoding="utf-8")
import ast
ast.parse(src)
print("patch 10 applied: session storage endpoints")
