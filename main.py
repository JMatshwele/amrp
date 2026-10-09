# ============================================================
#  BAMSTU Digital Pathogen Pro — Backend Server v3.4.6
#  main.py  |  FastAPI + Uvicorn
#  Dr. James P.T. Matshwele, BAMSTU Laboratories
#  University of Botswana, Department of Chemistry, 2026
#  Drug Discovery and Design Research Platform
#
#  START:  python -m uvicorn main:app --reload --port 8080
#  OPEN:   http://127.0.0.1:8080
# ============================================================

import os, math, random, urllib.request, urllib.parse, json, ssl, time
from pathlib import Path
from fastapi import BackgroundTasks, FastAPI, File, UploadFile, Form, HTTPException, Body
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import (JSONResponse, PlainTextResponse,
                               HTMLResponse, FileResponse)
import uvicorn
from datetime import datetime

# ── BAMSTU v3.2 modules ──────────────────────────────────────────
# bamstu_chem     : molecular-graph descriptors (replaces element counts)
# bamstu_analysis : contact analysis against the docking receptor
# bamstu_vvuq     : seeds, ensembles, redock RMSD, target verification
from starlette.concurrency import run_in_threadpool
import bamstu_chem as chem
import bamstu_analysis as contacts
import bamstu_provenance as prov
import bamstu_layers as layers
import bamstu_ligand_prep as ligand_prep
import bamstu_experimental as experimental
import bamstu_calibration as calibration
import bamstu_domain as domain
import bamstu_vvuq as vvuq
import bamstu_render as figures

# ============================================================
#  BAMSTU Digital Pathogen Pro — AutoDock Vina Integration
#  vina_module.py (injected into main.py)
#  Dr. James P.T. Matshwele, BAMSTU Laboratories, UB, 2026
#
#  STRATEGY: Vina runs in parallel with BAMSTU engine.
#  Results are returned side-by-side for comparison.
#  If Vina is not installed, everything falls back gracefully.
# ============================================================

import threading
import uuid
import shutil
import subprocess
import tempfile, shutil, ssl, urllib.request
from pathlib import Path

# ── Locate Vina binary ─────────────────────────────────────────
def find_vina() -> str | None:
    """Find AutoDock Vina executable. Returns path or None."""
    # 1. Check PATH
    found = shutil.which("vina") or shutil.which("vina.exe")
    if found:
        return found
    # 2. Common install locations on Windows
    candidates = [
        r"C:\Program Files\AutoDock Vina\vina.exe",
        r"C:\Program Files (x86)\AutoDock Vina\vina.exe",
        r"C:\AutoDock\vina.exe",
        r"C:\vina\vina.exe",
        r"C:\Users\drjam\vina.exe",
        r"C:\Users\drjam\BAMSTU_BIOINFO\backend\vina.exe",
        r"C:\Users\drjam\BAMSTU_BIOINFO\backend\vina",
    ]
    for c in candidates:
        if Path(c).exists():
            return c
    return None

VINA_EXE = find_vina()
VINA_AVAILABLE = VINA_EXE is not None

# Persistent receptor store — same folder as main.py so receptors
# survive across requests and never hit the temp-dir deletion bug
RECEPTOR_STORE = Path(__file__).parent / "data" / "receptors"
RECEPTOR_STORE.mkdir(parents=True, exist_ok=True)

def get_vina_status() -> dict:
    """Return Vina availability info for health endpoint."""
    if not VINA_AVAILABLE:
        return {
            "available": False,
            "path": None,
            "version": None,
            "message": "AutoDock Vina not found. Download from https://github.com/ccsb-scripps/AutoDock-Vina/releases"
        }
    try:
        r = subprocess.run([VINA_EXE, "--version"],
                           capture_output=True, text=True, timeout=5)
        version = r.stdout.strip() or r.stderr.strip()
        return {"available": True, "path": VINA_EXE,
                "version": version[:60], "message": "Ready"}
    except Exception as e:
        return {"available": True, "path": VINA_EXE,
                "version": "unknown", "message": str(e)[:60]}

def _rcsb_url(pdb_id: str, fmt: str = "pdb") -> str:
    return f"https://files.rcsb.org/download/{pdb_id}.{fmt}"


def _fetch_structure_text(pdb_id: str, timeout: int = 30) -> tuple:
    """Fetch a structure, falling back to mmCIF when no PDB file exists.

    9KGT returned HTTP 404 in the first live verification run. The entry is
    present in RCSB, but modern and oversized depositions are distributed as
    mmCIF only — the legacy 80-column PDB format cannot hold them. Treating
    a 404 as "target missing" silently drops exactly the newest and most
    relevant structures from the panel, which for a 2025 K. pneumoniae
    topoisomerase IV complex is the whole point of including it.

    Returns (text, format) where format is "pdb" or "cif".
    """
    import ssl as _ssl
    import urllib.request as _u
    ctx = _ssl.create_default_context()
    last = None
    for fmt in ("pdb", "cif"):
        try:
            req = _u.Request(_rcsb_url(pdb_id, fmt),
                             headers={"User-Agent": "BAMSTU/3.2"})
            with _u.urlopen(req, context=ctx, timeout=timeout) as r:
                data = r.read().decode("utf-8", errors="ignore")
            if len(data) > 1000:
                return data, fmt
        except Exception as exc:
            last = exc
            continue
    raise RuntimeError(
        f"{pdb_id}: no PDB or mmCIF file available from RCSB ({last}). "
        "Check the entry has not been obsoleted or superseded."
    )


def _cif_to_pdb_atoms(cif_text: str) -> str:
    """Minimal mmCIF atom_site -> PDB ATOM/HETATM conversion.

    Deliberately small: only the columns the receptor preparation and the
    contact analysis read. It is not a general converter and does not try to
    be one. Anything it cannot map is dropped rather than guessed.
    """
    lines, headers, in_loop = [], [], False
    for raw in cif_text.splitlines():
        s = raw.strip()
        if s.startswith("_atom_site."):
            headers.append(s.split(".", 1)[1])
            in_loop = True
            continue
        if in_loop and (s.startswith("#") or s.startswith("loop_") or not s):
            if headers and lines:
                break
            continue
        if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
            parts = s.split()
            if len(parts) < len(headers):
                continue
            f = dict(zip(headers, parts))
            try:
                x = float(f["Cartn_x"]); y = float(f["Cartn_y"]); z = float(f["Cartn_z"])
            except (KeyError, ValueError):
                continue
            rec = f.get("group_PDB", "ATOM")
            name = f.get("label_atom_id", "C").strip('"')
            resn = f.get("label_comp_id", "UNK")
            chain = (f.get("auth_asym_id") or f.get("label_asym_id") or "A")[:1]
            resi = f.get("auth_seq_id") or f.get("label_seq_id") or "1"
            elem = f.get("type_symbol", "C")
            try:
                resi_i = int(resi)
            except ValueError:
                continue
            lines.append(
                f"{rec:<6s}{len(lines)+1:5d} {name:^4s} {resn:>3s} {chain}"
                f"{resi_i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                f"          {elem:>2s}"
            )
    return "\n".join(lines) + "\nEND\n"


# ── Receptor preparation ───────────────────────────────────────
_RECEPTOR_CACHE: dict[str, Path] = {}
# Per-receptor preparation record — what was substituted, what was dropped.
# Surfaced through /api/receptor_provenance so a reader can see exactly
# which atoms Vina was given.
_RECEPTOR_PREP_LOG: dict[str, dict] = {}

def _fetch_pdb(pdb_id: str, out_dir: Path = None) -> Path:
    """Download PDB file from RCSB to persistent receptor store."""
    # Always use RECEPTOR_STORE so files survive across requests
    store = RECEPTOR_STORE
    pdb_path = store / f"{pdb_id}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 1000:
        return pdb_path
    # v3.2.1: two changes here.
    #
    # (a) mmCIF fallback. 9KGT returned HTTP 404 in the first live
    #     verification run. The entry exists; it simply has no legacy
    #     PDB-format file, because modern and oversized depositions are
    #     mmCIF only. Treating that 404 as "target unavailable" silently
    #     drops the newest structures from the panel, which for a 2025
    #     K. pneumoniae topoisomerase IV-NBTI complex is precisely the
    #     entry worth having.
    #
    # (b) TLS verification restored. This function previously set
    #     check_hostname = False and verify_mode = CERT_NONE. Receptor
    #     coordinates are the input to every score the platform produces,
    #     so disabling certificate checking on that download is a
    #     supply-chain question, not a convenience setting.
    text, fmt = _fetch_structure_text(pdb_id)
    if fmt == "cif":
        text = _cif_to_pdb_atoms(text)
        (store / f"{pdb_id}.cif_source").write_text(
            "converted from mmCIF; no legacy PDB file exists for this entry",
            encoding="utf-8")
    pdb_path.write_text(text, encoding="utf-8")
    return pdb_path


# ════════════════════════════════════════════════════════════════════
#  RECEPTOR HETATM POLICY  (v3.4.6)
# ════════════════════════════════════════════════════════════════════
#
# Until v3.4.5 receptor preparation stripped solvents and nothing else, so
# any co-crystallised inhibitor went straight into the docking receptor.
# 3RXX is KPC-2 with 3-nitrophenylboronic acid covalently bound to the
# catalytic Ser70, and NPB600 duly appeared in the contact list of every
# compound docked against it. Nine compounds were scored against an active
# site that was already occupied, and nothing downstream could detect it:
# the numbers looked entirely normal.
#
# A blanket "delete every HETATM" is the wrong repair and would break two
# things. Catalytic metals ARE the enzyme — remove the dinuclear Zn(II) of
# NDM-1 or the Mn(II) of FosA and there is no active site left. And the
# redock validation needs the native ligand: extracting it and docking it
# back is the entire test.
#
# So the classifier is tiered, and what it removes is recorded rather than
# discarded silently.

# Ions and cofactors that are part of the protein's chemistry.
CATALYTIC_HET = {
    "ZN", "MN", "MG", "CA", "FE", "FE2", "FES", "CU", "CU1", "NI", "CO",
    "K", "NA", "CD", "HG", "MO", "W", "SF4", "F3S", "HEM", "HEC", "HEB",
    "PLP", "PMP", "FAD", "FMN", "NAD", "NAP", "NDP", "SAM", "SAH", "COA",
    "B12", "CLA", "BCL", "MOS", "MSS",
}

# Crystallisation solvents, buffers and cryoprotectants.
SOLVENT_HET = {
    "HOH", "WAT", "H2O", "DOD", "EDO", "GOL", "PEG", "PGE", "PG4", "P6G",
    "1PE", "2PE", "DMS", "MES", "MPD", "PO4", "SO4", "ACT", "ACY", "FMT",
    "MLI", "IMD", "TRS", "EPE", "MRD", "SDS", "LMT", "BME", "DTT", "TCE",
    "CIT", "FLC", "TAR", "MAL", "BTB", "BCT", "NO3", "CL", "BR", "IOD",
    "AZI", "CAC", "NH4", "UNX", "UNL",
}

# A water is kept only when it completes a metal's coordination sphere.
METAL_WATER_SHELL_A = 2.9


def classify_receptor_het(pdb_text: str) -> dict:
    """Sort every HETATM residue into keep, drop-as-solvent, or extract.

    A residue that is neither a known cofactor nor a known solvent, and is
    larger than a few atoms, is treated as a co-crystallised ligand: removed
    from the docking receptor and written out separately so the redock can
    still use it.
    """
    residues: dict = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resn = line[17:20].strip().upper()
        key = (resn, line[21], line[22:26].strip())
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        residues.setdefault(key, []).append((line, xyz))

    metals = [xyz for (resn, _c, _n), atoms in residues.items()
              if resn in CATALYTIC_HET and len(atoms) <= 2
              for _l, xyz in atoms]

    keep, solvent, extract = set(), set(), set()
    for key, atoms in residues.items():
        resn = key[0]
        if resn in CATALYTIC_HET:
            keep.add(key)
        elif resn in SOLVENT_HET:
            if resn in ("HOH", "WAT", "H2O", "DOD") and metals:
                near = any(
                    ((xyz[0] - m[0]) ** 2 + (xyz[1] - m[1]) ** 2
                     + (xyz[2] - m[2]) ** 2) ** 0.5 <= METAL_WATER_SHELL_A
                    for _l, xyz in atoms for m in metals)
                (keep if near else solvent).add(key)
            else:
                solvent.add(key)
        elif len(atoms) <= 3:
            # Small unrecognised species: a stray ion or a fragment. Removed
            # as solvent rather than reported as a ligand, but counted.
            solvent.add(key)
        else:
            extract.add(key)

    return {"keep": keep, "solvent": solvent, "extract": extract,
            "residues": residues, "n_metals": len(metals)}


def _pdb_to_pdbqt_native(pdb_path: Path, pdbqt_path: Path) -> bool:
    """
    Convert PDB to PDBQT without OpenBabel.
    Keeps ATOM + HETATM (strips only water/solvent/common junk).
    Assigns Vina atom types; adds +0.000 partial charges.
    Note: for publication use OpenBabel for Gasteiger charges.
    """
    # v3.4.6: the tiered HETATM classifier replaces the flat solvent list.
    _raw_text = pdb_path.read_text(encoding='utf-8', errors='ignore')
    _pol = classify_receptor_het(_raw_text)
    _drop_keys = _pol["solvent"] | _pol["extract"]

    # Any co-crystallised ligand is written beside the receptor before it is
    # removed, so the redock validation still has its native ligand and so a
    # reviewer can see exactly what was taken out.
    _extracted = []
    if _pol["extract"]:
        lig_dir = pdbqt_path.parent / "cocrystal_ligands"
        lig_dir.mkdir(parents=True, exist_ok=True)
        for key in sorted(_pol["extract"]):
            resn, chain, seq = key
            atoms = _pol["residues"][key]
            out_lines = [
                f"REMARK  co-crystallised ligand extracted from "
                f"{pdb_path.stem} by BAMSTU v3.4.6\n",
                f"REMARK  {resn} chain {chain} residue {seq}, "
                f"{len(atoms)} atoms\n",
                "REMARK  removed from the docking receptor; retained here "
                "for redock validation\n",
            ] + [ln + "\n" for ln, _xyz in atoms] + ["END\n"]
            path = lig_dir / f"{pdb_path.stem}_{resn}_{chain}{seq}.pdb"
            path.write_text(''.join(out_lines), encoding='utf-8')
            _extracted.append({"resname": resn, "chain": chain,
                               "residue_seq": seq, "n_atoms": len(atoms),
                               "path": str(path)})
    # Map ALL transition metals to Zn — prevents Vina parse crash
    # Vina has no parameters for Co, Cu, Ni, Ag, Pt etc.
    # Meeting them in a receptor PDBQT causes silent abort before scoring.
    type_map = {
        'C':'C',  'N':'NA', 'O':'OA', 'S':'SA', 'H':'HD',
        'P':'P',  'F':'F',  'CL':'Cl','BR':'Br', 'I':'I',
        'ZN':'Zn','MG':'Mg','CA':'Ca','NA':'Na','K':'K',
        # Transition metals → Zn proxy (AutoDock4 standard)
        'FE':'Zn','MN':'Zn','CU':'Zn','CO':'Zn','NI':'Zn',
        'AG':'Zn','AU':'Zn','PT':'Zn','PD':'Zn','RU':'Zn',
        'RH':'Zn','OS':'Zn','IR':'Zn','V':'Zn', 'CR':'Zn',
        'MO':'Zn','W':'Zn', 'GA':'Zn','BI':'Zn','RE':'Zn',
        'TC':'Zn','SC':'Zn','TI':'Zn','ZR':'Zn','HG':'Zn',
        'CD':'Zn','SN':'Zn','PB':'Zn',
    }
    lines_out = []
    unknown_elements: dict[str, int] = {}
    metals_proxied: dict[str, int] = {}
    n_solvent_dropped = 0
    n_cocrystal_dropped = 0
    for line in _raw_text.split('\n'):
        # Accept both ATOM and HETATM records
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        if line.startswith('HETATM'):
            key = (line[17:20].strip().upper(), line[21],
                   line[22:26].strip())
            if key in _pol["extract"]:
                n_cocrystal_dropped += 1
                continue
            if key in _pol["solvent"]:
                n_solvent_dropped += 1
                continue
        # Determine element
        elem = line[76:78].strip().upper()
        if not elem:
            raw = line[12:14].strip()
            # Strip digits from atom name (e.g. "CA1" -> "C")
            elem = ''.join(c for c in raw if c.isalpha())[:2].upper()
        if elem not in type_map:
            # v3.1 typed anything unrecognised as carbon and carried on, so a
            # mis-parsed element became a silent structural artefact in the
            # grid maps. Now it is dropped and counted, and the count is
            # returned to the caller.
            unknown_elements[elem] = unknown_elements.get(elem, 0) + 1
            continue
        vt = type_map.get(elem)
        if vt == 'Zn' and elem != 'ZN':
            metals_proxied[elem] = metals_proxied.get(elem, 0) + 1
        base = line[:66].ljust(66)
        lines_out.append(f"{base}    +0.000 {vt:<2s}\n")
    if not lines_out:
        return False
    header = [
        "REMARK BAMSTU v3.4.6 receptor PDBQT\n",
        "REMARK charge model: formal +0.000 on every atom "
        "(install OpenBabel for Gasteiger-Marsili)\n",
    ]
    if metals_proxied:
        header.append(
            "REMARK METAL PROXY: " + ", ".join(
                f"{k}->Zn x{v}" for k, v in sorted(metals_proxied.items())
            ) + " (Vina has no parameters for these elements)\n"
        )
    if unknown_elements:
        header.append(
            "REMARK DROPPED unrecognised elements: " + ", ".join(
                f"{k} x{v}" for k, v in sorted(unknown_elements.items())
            ) + "\n"
        )
    if _extracted:
        header.append(
            "REMARK CO-CRYSTAL LIGANDS REMOVED: " + ", ".join(
                f"{e['resname']}{e['residue_seq']}(chain {e['chain']}, "
                f"{e['n_atoms']} atoms)" for e in _extracted) +
            " -- retained in cocrystal_ligands/ for redock validation\n")
    else:
        header.append("REMARK CO-CRYSTAL LIGANDS REMOVED: none found\n")
    header.append(
        f"REMARK HETATM POLICY: {len(_pol['keep'])} cofactor/metal-shell "
        f"residue(s) kept, {n_solvent_dropped} solvent atom(s) and "
        f"{n_cocrystal_dropped} co-crystal ligand atom(s) removed\n")

    pdbqt_path.write_text(''.join(header + lines_out), encoding='utf-8')
    _RECEPTOR_PREP_LOG[pdbqt_path.stem] = {
        "atoms_written": len(lines_out),
        "cocrystal_ligands_removed": _extracted,
        "solvent_atoms_removed": n_solvent_dropped,
        "cocrystal_atoms_removed": n_cocrystal_dropped,
        "cofactor_residues_kept": len(_pol["keep"]),
        "metal_centres_detected": _pol["n_metals"],
        "hetatm_policy": "v3.4.6 tiered: catalytic cofactors and "
                         "metal-completing waters kept; solvents dropped; "
                         "co-crystal ligands extracted then removed",
        "metals_mapped_to_zn": metals_proxied,
        "unknown_elements_dropped": unknown_elements,
        "charge_model": "formal_zero",
    }
    return True

RECEPTOR_POLICY_TAG = "REMARK HETATM POLICY:"


def _receptor_is_current(pdbqt: Path) -> bool:
    """True when this file was built under the v3.4.6 HETATM policy.

    A receptor prepared by an earlier build carries no policy line and may
    still contain a co-crystallised inhibitor. Reusing it from the store
    would silently reproduce the defect this release fixes, so a stale file
    is rebuilt rather than trusted.
    """
    try:
        head = pdbqt.read_text(encoding="utf-8", errors="ignore")[:2000]
    except Exception:
        return False
    return RECEPTOR_POLICY_TAG in head



def _apply_het_policy_to_pdb(pdb_path: Path, store: Path) -> Path | None:
    """Write a policy-cleaned copy of a PDB and extract its co-crystal ligands.

    Needed because Open Babel converts the file it is given. Applying the
    policy only inside the native converter would leave the Open Babel path
    — the one that produces the better charges, and therefore the one a
    publication run uses — still building receptors with the inhibitor in.
    """
    try:
        text = pdb_path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None, []
    pol = classify_receptor_het(text)
    remarks = [
        f"REMARK HETATM POLICY: {len(pol['keep'])} cofactor/metal-shell "
        f"residue(s) kept, {len(pol['solvent'])} solvent residue(s) and "
        f"{len(pol['extract'])} co-crystal ligand(s) removed"]
    if pol["extract"]:
        remarks.append(
            "REMARK CO-CRYSTAL LIGANDS REMOVED: " + ", ".join(
                f"{r}{n}(chain {c})" for r, c, n in sorted(pol["extract"]))
            + " -- retained in cocrystal_ligands/")
    else:
        remarks.append("REMARK CO-CRYSTAL LIGANDS REMOVED: none found")
    if not (pol["extract"] or pol["solvent"]):
        return None, remarks

    if pol["extract"]:
        lig_dir = store / "cocrystal_ligands"
        lig_dir.mkdir(parents=True, exist_ok=True)
        for key in sorted(pol["extract"]):
            resn, chain, seq = key
            atoms = pol["residues"][key]
            (lig_dir / f"{pdb_path.stem}_{resn}_{chain}{seq}.pdb").write_text(
                f"REMARK  co-crystallised ligand extracted from "
                f"{pdb_path.stem} by BAMSTU v3.4.6\n"
                f"REMARK  {resn} chain {chain} residue {seq}, "
                f"{len(atoms)} atoms\n"
                + "".join(ln + "\n" for ln, _x in atoms) + "END\n",
                encoding="utf-8")

    drop = pol["solvent"] | pol["extract"]
    out = []
    for line in text.splitlines():
        if line.startswith("HETATM"):
            key = (line[17:20].strip().upper(), line[21],
                   line[22:26].strip())
            if key in drop:
                continue
        out.append(line)
    cleaned = pdb_path.with_name(pdb_path.stem + "_policy.pdb")
    cleaned.write_text("\n".join(out) + "\n", encoding="utf-8")
    return cleaned, remarks



# Vina atom types the engine will parse in a receptor. Anything else makes
# it abort before scoring, usually without a useful message.
VINA_RECEPTOR_TYPES = {
    "C", "A", "N", "NA", "NS", "O", "OA", "OS", "S", "SA", "H", "HD", "HS",
    "P", "F", "Cl", "Br", "I", "Mg", "Ca", "Mn", "Fe", "Zn",
}

# Everything Vina cannot type is mapped to Zn, exactly as the native
# converter does. The substitution is a limitation, not a repair, and it is
# declared in the receptor header and in Methods.
_METAL_TO_ZN = {
    "AG", "AU", "PT", "PD", "RU", "RH", "IR", "OS", "CU", "CO", "NI", "CD",
    "HG", "PB", "MO", "W", "V", "CR", "TI", "ZR", "SN", "SB", "SE", "TE",
    "GA", "IN", "TL", "BI", "RE", "TC", "NB", "TA", "HF", "SR", "BA", "CS",
    "RB", "LI", "BE", "AL", "SI",
}



# Every two-letter element symbol, so an ion is never mistaken for a protein
# atom name. Without this, "K" fell through a first-character lookup and
# became carbon, and "NA" became nitrogen -- which is worse, because it puts
# a phantom hydrogen-bond acceptor exactly where a sodium ion sits.
_TWO_LETTER_SYMBOLS = {
    "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR", "CA", "SC",
    "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE",
    "BR", "KR", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
    "CD", "IN", "SN", "SB", "TE", "XE", "CS", "BA", "LA", "CE", "PR", "ND",
    "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA",
    "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
    "FR", "RA", "AC", "TH", "PA", "NP", "PU",
}
_SINGLE_LETTER_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "K", "V", "Y",
                           "I", "B", "W", "U"}
_ORGANIC_TO_VINA = {"C": "C", "N": "N", "O": "OA", "S": "SA", "H": "HD",
                    "P": "P", "F": "F", "I": "I", "B": "C"}


def _safe_vina_type(declared: str, atom_name: str, res_name: str) -> str:
    """Map an atom Vina cannot type, without inventing a different element.

    v3.4.12. The previous fallback read the first character of the atom name
    when the two-letter form was not in the metal table. For a monatomic ion
    the atom name IS the element, so potassium (K) became carbon and sodium
    (NA) became nitrogen. A metal proxied as zinc is a declared substitution
    with a known direction; a metal silently rewritten as a hydrogen-bond
    acceptor is a fabricated interaction sitting at an ion's coordinates.

    A residue whose name matches its atom name is a monatomic ion, and its
    element is not in doubt. Anything Vina cannot type is proxied to zinc and
    declared; nothing becomes an organic atom by accident.
    """
    letters = "".join(c for c in str(atom_name) if c.isalpha()).upper()
    res = "".join(c for c in str(res_name) if c.isalpha()).upper()

    # A monatomic ion is identified by its residue name matching its atom
    # name -- ZN/ZN, K/K, NA/NA. Nothing else is safe to read as an element,
    # because stripping the digits off a branch code leaves a valid symbol:
    # NE2 becomes NE (neon), CA stays CA (calcium), SD becomes SD, OD1
    # becomes OD. Requiring the residue to agree is what separates a calcium
    # ion from an alpha carbon.
    element = None
    if res and res == letters:
        if letters in _TWO_LETTER_SYMBOLS or letters in _SINGLE_LETTER_ELEMENTS:
            element = letters
    elif str(declared).upper() in _TWO_LETTER_SYMBOLS:
        element = str(declared).upper()
    # No further two-letter guess. Protein atom names collide with element
    # symbols throughout -- NE, CA, CB, CD, CE, CG, ND, NZ, OD, OE, OG, SD,
    # SG are all real symbols and all common branch codes. Trusting the
    # first two letters turned the NE2 of a histidine into neon, and from
    # there into a zinc proxy sitting on an imidazole nitrogen.

    if element is None:
        # A branched protein atom: CZ2, NE1, OD1, SG. First character is the
        # element and the rest is a branch code.
        return _ORGANIC_TO_VINA.get(letters[:1], "C")

    title = element.capitalize()
    if title in VINA_RECEPTOR_TYPES:
        return title
    if element in _ORGANIC_TO_VINA:
        return _ORGANIC_TO_VINA[element]
    return "Zn"


def _normalise_obabel_pdbqt(pdbqt: Path, policy_remarks: list) -> dict:
    """Make an Open Babel receptor safe for Vina, and stamp the policy.

    v3.4.8. The Open Babel branch returned its output untouched the moment
    the file was over 500 bytes. Two things followed, and both of them bite
    only once Open Babel is actually installed:

    1. Open Babel writes real element types for metals. FosA carries Mn,
       NDM-1 carries two Zn, LpxH carries Mn. Vina has no parameters for
       most of those and aborts before scoring, so every metalloprotein
       target in the panel became unusable the day Open Babel appeared on
       the machine. The native converter had always mapped metals to Zn;
       the Open Babel path never did.

    2. The file carried no HETATM-policy REMARK, so the v3.4.6 staleness
       check judged it out of date on the very next call, deleted it and
       rebuilt it. Every sweep re-downloaded and re-converted the same four
       receptors, for ever.
    """
    try:
        text = pdbqt.read_text(encoding="utf-8", errors="ignore")
    except Exception as exc:
        return {"ok": False, "reason": f"unreadable: {exc}"}

    out, remapped, atoms = [], {}, 0
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith("REMARK") and "HETATM POLICY" in line:
                continue                      # ours is added below
            out.append(line)
            continue
        atoms += 1
        padded = line.ljust(79)
        atype = padded[77:79].strip()
        if atype and atype not in VINA_RECEPTOR_TYPES:
            new = _safe_vina_type(atype, padded[12:16], padded[17:20])
            remapped[atype] = remapped.get(atype, 0) + 1
            padded = padded[:77] + f"{new:<2}"
        out.append(padded.rstrip())

    if atoms == 0:
        return {"ok": False, "reason": "no atom records"}

    header = ["REMARK BAMSTU v3.4.8 receptor PDBQT (Open Babel, "
              "Gasteiger-Marsili charges)"] + policy_remarks
    if remapped:
        metals = {k: n for k, n in remapped.items()
                  if k.upper() in _TWO_LETTER_SYMBOLS
                  or k.upper() in _SINGLE_LETTER_ELEMENTS}
        others = {k: n for k, n in remapped.items() if k not in metals}
        if metals:
            header.append(
                "REMARK METAL PROXY: " + ", ".join(
                    f"{k}->Zn ({n} atom(s))" for k, n in sorted(metals.items()))
                + " -- Vina has no parameters for these elements")
        if others:
            header.append(
                "REMARK RETYPED: " + ", ".join(
                    f"{k} ({n} atom(s))" for k, n in sorted(others.items()))
                + " -- atom type not recognised by Vina, resolved from the "
                  "atom name")
    pdbqt.write_text("\n".join(header + out) + "\n", encoding="utf-8")
    return {"ok": True, "atoms": atoms, "remapped": remapped}



def _restore_quarantined(pdb_id: str, pdbqt: Path, quarantine) -> Path | None:
    """Put back the pre-policy receptor when the rebuild could not finish.

    Returning the older file is the lesser harm. It may carry a
    co-crystallised ligand, which is a limitation the header and the audit
    both declare; having no receptor stops the panel dead and loses the only
    copy on the machine.
    """
    if not quarantine or not Path(quarantine).exists():
        return None
    try:
        Path(quarantine).rename(pdbqt)
    except Exception:
        return None
    _RECEPTOR_CACHE[pdb_id] = pdbqt
    _RECEPTOR_PREP_LOG.setdefault(pdb_id, {})["restored_pre_policy"] = True
    print(f"[BAMSTU] {pdb_id}: rebuild failed, so the previous receptor has "
          f"been restored. IT MAY STILL CONTAIN A CO-CRYSTALLISED LIGAND. "
          f"Check the network path to RCSB, then delete "
          f"backend/data/receptors/{pdb_id}.pdbqt to force a clean rebuild.",
          flush=True)
    return pdbqt


def prepare_receptor(pdb_id: str, work_dir: Path = None) -> Path | None:
    """Download and prep receptor PDBQT to persistent store."""
    _cached = _RECEPTOR_CACHE.get(pdb_id)
    if _cached is not None and not _receptor_is_current(Path(_cached)):
        # v3.4.9: set aside, do not delete. The v3.4.6 version deleted the
        # stale file before knowing whether a replacement could be built,
        # so on a machine that cannot reach RCSB -- a VPN, a firewall, an
        # RCSB outage -- it destroyed a working receptor and left nothing.
        # That turned a network problem into data loss, and it is what took
        # a whole panel down at once.
        _RECEPTOR_CACHE.pop(pdb_id, None)
        print(f"[BAMSTU] {pdb_id}: receptor predates the v3.4.6 HETATM "
              f"policy and may retain a co-crystal ligand — rebuilding "
              f"(the existing file is kept until the rebuild succeeds).",
              flush=True)
    if pdb_id in _RECEPTOR_CACHE:
        p = _RECEPTOR_CACHE[pdb_id]
        if p.exists(): return p
        del _RECEPTOR_CACHE[pdb_id]          # stale cache entry
    # Always store in RECEPTOR_STORE (persistent, not temp dir)
    pdbqt = RECEPTOR_STORE / f"{pdb_id}.pdbqt"
    # Validate cached receptor: must be >2KB AND use correct metal typing.
    # Old files (built before type_map fix) have raw metal types like
    # " Co" or " Fe" in the atom-type column — Vina parse-crashes on these.
    # Correct files map all metals to Zn, so only " Zn" should appear.
    if pdbqt.exists():
        if pdbqt.stat().st_size > 2000:
            _rc = pdbqt.read_text(encoding="utf-8", errors="ignore")
            # Detect stale metal types in the Vina type column (last 2 chars per line)
            _stale_types = {" Co", " Fe", " Cu", " Ni", " Ag", " Pt",
                            " Ru", " Rh", " Au", " Os", " Ir", " Mn",
                            " Cr", " Mo", " V\n"}
            _has_stale = any(
                any(line.endswith(st) for st in _stale_types)
                for line in _rc.split("\n")
                if line.startswith("ATOM") or line.startswith("HETATM")
            )
            # v3.4.6: a file with no HETATM-policy line was built before
            # co-crystal ligands were removed, so it may still carry an
            # inhibitor bolted to the catalytic residue. Reusing it would
            # silently reproduce the defect this release exists to fix.
            _no_policy = RECEPTOR_POLICY_TAG not in _rc[:2000]
            if not _has_stale and not _no_policy:
                _RECEPTOR_CACHE[pdb_id] = pdbqt
                return pdbqt
            if _no_policy:
                print(f"[BAMSTU] {pdb_id}: receptor predates the v3.4.6 "
                      f"HETATM policy — rebuilding to remove any "
                      f"co-crystallised ligand.", flush=True)
            # Move aside rather than delete, and restore it if the rebuild
            # cannot be completed. A receptor that may contain a co-crystal
            # ligand is a disclosed limitation; no receptor at all is a dead
            # panel.
            _quarantine = pdbqt.with_suffix(".pdbqt.pre-policy")
            try:
                if _quarantine.exists():
                    _quarantine.unlink()
                pdbqt.rename(_quarantine)
            except Exception:
                pdbqt.unlink(missing_ok=True)
                _quarantine = None
        else:
            pdbqt.unlink()
    try:
        pdb = _fetch_pdb(pdb_id)

        # v3.4.6: Open Babel converts whatever PDB it is handed, so the
        # HETATM policy has to be applied to the input rather than to the
        # output. The cleaned file is what both conversion paths see, and
        # the extracted co-crystal ligands are written beside the receptor.
        _clean, _policy_remarks = _apply_het_policy_to_pdb(pdb, pdbqt.parent)
        pdb = _clean or pdb

        # Try OpenBabel first (better quality)
        ob = shutil.which("obabel")
        if ob:
            r = subprocess.run(
                [ob, str(pdb), "-O", str(pdbqt), "-xr",
                 "--partialcharge", "gasteiger"],
                capture_output=True, timeout=30
            )
            if pdbqt.exists() and pdbqt.stat().st_size > 500:
                norm = _normalise_obabel_pdbqt(pdbqt, _policy_remarks)
                if norm["ok"]:
                    _RECEPTOR_PREP_LOG[pdb_id] = {
                        "converter": "openbabel + BAMSTU normalisation",
                        "atoms_written": norm["atoms"],
                        "metal_types_remapped": norm["remapped"],
                        "charge_model": "Gasteiger-Marsili (Open Babel)",
                        "hetatm_policy": "v3.4.6 tiered, applied to the input",
                    }
                    _RECEPTOR_CACHE[pdb_id] = pdbqt
                    _q = locals().get("_quarantine")
                    if _q and Path(_q).exists():
                        Path(_q).unlink(missing_ok=True)
                    return pdbqt
                print(f"[BAMSTU] {pdb_id}: Open Babel output rejected "
                      f"({norm['reason']}) — using the native converter.",
                      flush=True)
                pdbqt.unlink(missing_ok=True)
        # Fallback: native conversion
        if _pdb_to_pdbqt_native(pdb, pdbqt):
            _RECEPTOR_CACHE[pdb_id] = pdbqt
            _q = locals().get("_quarantine")
            if _q and Path(_q).exists():
                Path(_q).unlink(missing_ok=True)
            return pdbqt
        _RECEPTOR_PREP_LOG[pdb_id] = {
            "status": "FAILED",
            "reason": "native PDBQT conversion produced no atoms",
        }
        _restore_quarantined(pdb_id, pdbqt, locals().get("_quarantine"))
    except Exception as exc:
        # v3.1 swallowed this with a bare pass. In an unattended overnight
        # batch that turns every failure into the same unhelpful line
        # ("receptor preparation failed") with no way to tell a network
        # timeout from an obsoleted PDB entry from a disk problem.
        _RECEPTOR_PREP_LOG[pdb_id] = {
            "status": "FAILED",
            "reason": f"{type(exc).__name__}: {exc or 'no message'}",
        }
        _restored = _restore_quarantined(pdb_id, pdbqt,
                                         locals().get("_quarantine"))
        if _restored:
            return _restored
    return None

LIGAND_PREP_RECORDS: dict = {}


def prepare_ligand(atoms: list, ligand_name: str,
                   work_dir: Path, coord_fp: int = 0) -> Path | None:
    """Write the ligand PDBQT, rigid or flexible according to the structure.

    The decision is delegated to bamstu_ligand_prep, which perceives the
    connectivity from the coordinates and classifies the ligand. A
    pre-formed coordination complex is written rigid; a free organic ligand
    gets a real torsion tree. v3.1 through v3.3 wrote TORSDOF 0 for both,
    which meant that in a sweep containing L1-L4 and C1-C4 the ligands were
    scored without Vina's torsional penalty and from a single conformer,
    while the complexes were twice their size. The comparison those numbers
    supported was not a like-for-like one.

    The preparation record is kept so it can be printed in Methods and
    attached to every hit.
    """
    lig_dir = RECEPTOR_STORE.parent / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    pdbqt = lig_dir / f"{ligand_name}_{coord_fp}.pdbqt"

    built = ligand_prep.build_pdbqt(atoms, name=ligand_name)
    text = built["pdbqt"]

    # Rewrite whenever the content differs, so a cached file prepared by an
    # older version cannot survive a change in preparation policy.
    if pdbqt.exists():
        try:
            if pdbqt.read_text(encoding="utf-8", errors="ignore") == text:
                LIGAND_PREP_RECORDS[str(pdbqt)] = built["record"]
                built["record"]["path"] = str(pdbqt)
                return pdbqt
        except Exception:
            pass
        pdbqt.unlink()

    pdbqt.write_text(text, encoding="utf-8")
    built["record"]["path"] = str(pdbqt)
    LIGAND_PREP_RECORDS[str(pdbqt)] = built["record"]
    return pdbqt

def compute_grid_centre(pdbqt: Path) -> tuple:
    """Compute grid box centre from receptor PDBQT ATOM and HETATM lines.
    
    Reads BOTH ATOM and HETATM because _pdb_to_pdbqt_native preserves both.
    If only HETATM is present (non-standard residues), old ATOM-only code
    returned (0,0,0) → Vina searched at origin → missed the protein entirely.
    """
    xs, ys, zs = [], [], []
    for line in pdbqt.read_text(encoding='utf-8', errors='ignore').split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
            except ValueError:
                continue
    if not xs:
        return 0.0, 0.0, 0.0, 40.0, 40.0, 40.0
    cx = sum(xs)/len(xs)
    cy = sum(ys)/len(ys)
    cz = sum(zs)/len(zs)
    # Minimum 30A box; 20A buffer each side; cap at Vina 126A max
    sx = max(30.0, min((max(xs)-min(xs))+20, 126))
    sy = max(30.0, min((max(ys)-min(ys))+20, 126))
    sz = max(30.0, min((max(zs)-min(zs))+20, 126))
    return cx, cy, cz, round(sx,1), round(sy,1), round(sz,1)

def compute_active_site_box(pdb_id: str, active_residues: list,
                            buffer_A: float = 6.0,
                            min_edge: float = 18.0,
                            max_edge: float = 30.0,
                            chain: str = "") -> tuple | None:
    """Bounding box of the active-site residues plus a small buffer.

    v3.1 used a fixed 30 A cube for targeted docking and, for blind docking,
    the full protein extent plus 20 A capped at Vina's 126 A limit. Vina's
    search quality falls off as the volume grows -- the same exhaustiveness
    spread over a larger box means fewer samples per unit volume -- so an
    oversized box reliably produces plausible-looking surface poses.

    Here the box is the actual extent of the annotated pocket residues plus
    `buffer_A`, clamped to a sane range so a single-residue annotation
    cannot produce a box too small to hold the ligand.

    Returns (cx, cy, cz, sx, sy, sz, n_residues_found) or None.
    """
    if not active_residues:
        return None
    try:
        pdb = _fetch_pdb(pdb_id)
        text = pdb.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None

    # v3.4.15: restrict to one chain.
    #
    # Without this the residue numbers were matched across EVERY chain, and
    # for an oligomer that is fatal. 5O79 is an OmpK36 homotrimer: each of
    # chains A, B and C carries its own barrel and its own L3 constriction,
    # at (14,-13,30), (32,-12,-11) and (22,26,10). Averaging the three puts
    # the box centre at (23,0,10) -- which is the trimer's three-fold axis,
    # the lipid-facing core between the barrels, inside no lumen at all. The
    # protein centroid is (26,0,11), so the "active site" box was
    # effectively the middle of the whole assembly and every compound docked
    # to the outside of it.
    #
    # The same applies to the Wza octamer and the OqxB trimer. The chain was
    # declared in every target record and never used.
    xs, ys, zs, found = [], [], [], set()
    wanted = set(active_residues)
    chains_seen = set()
    want_chain = (chain or "").strip()
    for line in text.split("\n"):
        if not line.startswith("ATOM"):
            continue
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        if resi not in wanted:
            continue
        chains_seen.add(line[21])
        if want_chain and line[21] != want_chain:
            continue
        try:
            # All side-chain atoms, not just CA. A pocket is defined by the
            # side chains that line it; a CA-only centroid sits behind the
            # residues rather than in the cavity they enclose.
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
            found.add(resi)
        except ValueError:
            continue

    if not xs:
        return None

    def edge(vals):
        return max(min_edge, min(max_edge, (max(vals) - min(vals)) + 2 * buffer_A))

    return (
        (max(xs) + min(xs)) / 2.0,
        (max(ys) + min(ys)) / 2.0,
        (max(zs) + min(zs)) / 2.0,
        round(edge(xs), 1), round(edge(ys), 1), round(edge(zs), 1),
        len(found),
    )


# ════════════════════════════════════════════════════════════════════
#  GRID BOX AND SEARCH EFFORT POLICY  (v3.4.1)
# ════════════════════════════════════════════════════════════════════
#
# A 30 A cube is wrong in two opposite directions depending on the target.
# For a buried enzyme active site of maybe 800 A^3 it spreads the search
# over thirty times the volume that matters, which is where surface false
# positives come from. For a transmembrane channel it is too SMALL: the
# permeation path through OqxB or the Wza octamer is longer than 30 A, so
# the box clips the pore and a 57-atom rigid Ag(I) complex has nowhere
# inside the grid that it fits without clashing. Vina then returns no mode
# and the pipeline falls through to SYNTHETIC_PLACEMENT_NOT_DOCKED, which
# is the state the sweep reported for C1 against 7CZ9 and 2J58.
#
# So the box is a property of the target class, and the search effort is a
# property of both the target and the ligand.

CHANNEL_TARGETS = {
    "7CZ9": {"box": (34.0, 34.0, 36.0),
             "pore_residues": ["PHE", "TYR", "TRP"],
             "note": "OqxB RND transporter: substrate translocation pathway, "
                     "long axis normal to the membrane plane."},
    "2J58": {"box": (34.0, 34.0, 36.0),
             "pore_residues": ["ARG", "LYS", "ASN", "GLN"],
             "note": "Wza outer-membrane capsule export octamer: the alpha-"
                     "helical barrel lumen runs the length of the assembly."},
}

# Above this heavy-atom count a rigid ligand needs room to rotate inside the
# box as a whole body, because no torsion can relieve a clash.
LARGE_RIGID_HEAVY_ATOMS = 30
LARGE_LIGAND_BOX_BUFFER = 6.0        # A added to each edge
DEFAULT_EXHAUSTIVENESS = 20          # v3.4.1: was 8 in v3.1, 16 in v3.3
CHANNEL_EXHAUSTIVENESS = 32
METAL_COMPLEX_EXHAUSTIVENESS = 32
ENSEMBLE_SEEDS = (42, 1042, 2042, 3042, 4042)
ENSEMBLE_SD_LIMIT_KCAL = 0.5


def resolve_box(target_pdb: str, measured_box=None, n_heavy: int = 0,
                is_metal_complex: bool = False) -> dict:
    """Decide the grid box for one target/ligand pair, and say why.

    Returns the box plus the reasoning, because a box that changed between
    runs without a record is a box that makes two numbers incomparable.
    """
    reasons = []
    ch = CHANNEL_TARGETS.get(str(target_pdb).upper())
    if ch:
        sx, sy, sz = ch["box"]
        reasons.append(f"Transmembrane channel: {sx:.0f}x{sy:.0f}x{sz:.0f} A "
                       f"covering the permeation path. {ch['note']}")
        centre_mode = "pore_lining_residues"
    elif measured_box and len(measured_box) >= 6:
        sx, sy, sz = float(measured_box[3]), float(measured_box[4]), float(measured_box[5])
        reasons.append("Measured active-site bounding box supplied by the "
                       "caller, padded by the site detector.")
        centre_mode = "measured_site_centroid"
    else:
        sx = sy = sz = 24.0
        reasons.append("No measured site: 24 A default cube on the pocket "
                       "centroid.")
        centre_mode = "pocket_centroid"

    if n_heavy >= LARGE_RIGID_HEAVY_ATOMS and not ch:
        # Not stacked on a channel box. The 34x34x36 channel box is already
        # sized for the permeation path rather than for a pocket, and adding
        # a large-ligand buffer on top of it pushes the grid out past the
        # lumen into the lipid-facing shell, which puts sampling effort where
        # nothing binds.
        sx += LARGE_LIGAND_BOX_BUFFER
        sy += LARGE_LIGAND_BOX_BUFFER
        sz += LARGE_LIGAND_BOX_BUFFER
        reasons.append(
            f"Ligand has {n_heavy} heavy atoms and is docked rigid, so "
            f"{LARGE_LIGAND_BOX_BUFFER:.0f} A was added to each edge. A rigid "
            "body cannot relieve a box-edge clash by rotating a bond, so a "
            "box sized for a drug-like organic returns no mode at all.")

    exh = DEFAULT_EXHAUSTIVENESS
    if ch:
        exh = max(exh, CHANNEL_EXHAUSTIVENESS)
        reasons.append(f"Channel target: exhaustiveness raised to {exh}, "
                       "because search quality falls with volume and this box "
                       "is roughly twice the default.")
    if is_metal_complex:
        exh = max(exh, METAL_COMPLEX_EXHAUSTIVENESS)
        reasons.append(f"Rigid metal complex: exhaustiveness raised to {exh}. "
                       "With no torsions the search relies entirely on rigid-"
                       "body sampling.")

    return {"size": (round(sx, 1), round(sy, 1), round(sz, 1)),
            "centre_mode": centre_mode,
            "exhaustiveness": exh,
            "is_channel": bool(ch),
            "pore_residues": (ch or {}).get("pore_residues", []),
            "reasoning": reasons}


def pore_axis_centre(receptor_pdb: Path, residue_names: list,
                     fallback: tuple) -> tuple:
    """Centroid of the pore-lining residues, for a channel target.

    The geometric centre of a transmembrane assembly is often in the lipid-
    facing shell rather than in the lumen. Centring on the polar or aromatic
    residues that line the path puts the box where transport happens.
    """
    if not receptor_pdb or not Path(receptor_pdb).exists():
        return fallback
    xs, ys, zs = [], [], []
    try:
        for line in Path(receptor_pdb).read_text(
                encoding="utf-8", errors="ignore").splitlines():
            if not line.startswith("ATOM"):
                continue
            if line[17:20].strip().upper() not in residue_names:
                continue
            if line[12:16].strip() not in ("CA", "CZ", "NZ", "OH", "NE2"):
                continue
            xs.append(float(line[30:38])); ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
    except Exception:
        return fallback
    if len(xs) < 4:
        return fallback
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def ensemble_statistics(values: list, seeds: list) -> dict:
    """Mean, SD and the usability verdict for a replicate seed ensemble.

    A single seed gives a point estimate with no way to tell a real
    difference between two compounds from search noise. Five seeds give the
    spread, and the spread is the part that decides whether a 0.3 kcal/mol
    gap between two targets means anything.
    """
    import statistics as _st
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return {"status": "NO_DATA", "n": 0}
    mean = sum(vals) / len(vals)
    sd = _st.stdev(vals) if len(vals) > 1 else 0.0
    usable = sd <= ENSEMBLE_SD_LIMIT_KCAL
    return {
        "status": "OK",
        "n": len(vals),
        "seeds": list(seeds)[:len(vals)],
        "values_kcal": [round(v, 3) for v in vals],
        "mean_kcal": round(mean, 3),
        "sd_kcal": round(sd, 3),
        "range_kcal": [round(min(vals), 3), round(max(vals), 3)],
        "report_as": f"{mean:.2f} +/- {sd:.2f} kcal/mol (n={len(vals)} seeds)",
        "sd_limit_kcal": ENSEMBLE_SD_LIMIT_KCAL,
        "within_guardrail": usable,
        "verdict": (
            "Sampling is converged enough that the point estimate is usable."
            if usable else
            f"Seed-to-seed SD of {sd:.2f} kcal/mol exceeds the "
            f"{ENSEMBLE_SD_LIMIT_KCAL:.1f} guardrail. The search has not "
            "converged, so this affinity cannot separate one compound from "
            "another and raising exhaustiveness is the fix, not averaging "
            "more seeds."),
        "caveat": (
            "This SD is search reproducibility only. It is not the accuracy "
            "of the scoring function, which is about 1.5 kcal/mol and does "
            "not shrink with more seeds."),
    }


def run_vina(ligand_pdbqt: Path, receptor_pdbqt: Path,
             pdb_id: str,
             cx: float = None, cy: float = None, cz: float = None,
             exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
             n_heavy_atoms: int = 0,
             is_metal_complex: bool = False,
             blind: bool = False,
             seed: int = 42,
             box: tuple = None) -> dict:
    """
    Run AutoDock Vina and return structured result.

    Parameters
    ----------
    blind : bool, default False
        When False (default): targeted active-site docking.
            Grid centred on cx/cy/cz (crystallographic active site),
            fixed 30×30×30 Å box. Fast, focused, publication-grade.
        When True: whole-protein blind docking.
            Grid centred on the receptor's geometric centroid and sized
            to cover the full protein (compute_grid_centre() sx/sy/sz).
            Used for PROMISCUITY SCREENING — detects whether the compound
            prefers a surface site over the active site.
            Slower; run at lower exhaustiveness (8 recommended).
            Returns docking_mode: "BLIND" in result dict.

    Scientific basis for dual-mode:
        Active-site docking gives the functional ΔG (what matters for IC50).
        Blind docking gives the global minimum. If blind_ΔG << active_ΔG,
        the compound is promiscuous — it binds elsewhere on the protein with
        higher affinity than at the catalytic site. This flag should be
        reported explicitly in any publication (Shoichet 2006, Drug Discov Today).

    BUG FIXES vs old version:
    1. Output PDBQT path now includes compound hash → unique per compound+target
       (old: same path for different compounds → cached file reused → same score)
    2. Output file deleted before each run → Vina always reruns fresh
       (old: no deletion → Format A read stale score from previous run)
    3. Receptor PDBQT validated for minimum content before running
       (old: empty receptor → Vina returns garbage → Format D picked -4.75 stray float)
    4. Format D disabled — sane range tightened, requires score-table context
       (old: Format D could grab -4.75 from Vina's own startup banner)
    """
    box_policy = None      # set by resolve_box() for the targeted path
    if not VINA_AVAILABLE:
        return {"available": False,
                "message": "AutoDock Vina not found. Place vina.exe in the backend folder.",
                "download": "https://github.com/ccsb-scripps/AutoDock-Vina/releases"}
    try:
        # ── Validate receptor has real content ───────────────────
        if not receptor_pdbqt.exists() or receptor_pdbqt.stat().st_size < 500:
            return {"available": True, "success": False,
                    "error": f"Receptor PDBQT for {pdb_id} is empty or missing. "
                             f"Try clicking 'Pre-load Receptors' first.",
                    "receptor_path": str(receptor_pdbqt)}
        receptor_text = receptor_pdbqt.read_text(encoding="utf-8", errors="ignore")
        atom_lines = [l for l in receptor_text.split("\n")
                      if l.startswith("ATOM") or l.startswith("HETATM")]
        if len(atom_lines) < 10:
            # Receptor has too few atoms — delete and force rebuild
            receptor_pdbqt.unlink()
            return {"available": True, "success": False,
                    "error": f"Receptor {pdb_id} has only {len(atom_lines)} ATOM records — "
                             f"too few for valid docking. File deleted; retry to rebuild.",
                    "atom_count": len(atom_lines)}

        # ── Grid box: active-site (targeted) or whole-protein (blind) ──────
        #
        # compute_grid_centre returns (cx,cy,cz, sx,sy,sz) where sx/sy/sz
        # span the full protein + 20 Å buffer. We use it for BOTH modes:
        #   Active-site mode: centre = active-site centroid, size = 30×30×30 Å
        #   Blind mode:       centre = protein centroid,     size = full protein
        #
        # This is the only architectural change needed to support promiscuity
        # screening — Vina itself is identical in both modes; only the grid
        # box definition changes.
        _pc_cx, _pc_cy, _pc_cz, _pc_sx, _pc_sy, _pc_sz = compute_grid_centre(receptor_pdbqt)

        if blind:
            # ── BLIND / WHOLE-PROTEIN DOCKING ──────────────────────────────
            # Centre on protein geometric centroid; size covers full structure.
            # Overrides any supplied cx/cy/cz — those are active-site coords.
            cx, cy, cz = _pc_cx, _pc_cy, _pc_cz
            sx, sy, sz = _pc_sx, _pc_sy, _pc_sz   # from compute_grid_centre
            # Vina caps individual box dimensions at 126 Å; compute_grid_centre
            # already applies this cap. No further adjustment needed.
        else:
            # ── TARGETED / ACTIVE-SITE DOCKING (default) ───────────────────
            # Use crystallographic active site centroid if provided.
            # Catch BOTH None (network failure) AND (0,0,0) (origin trap).
            _need_fallback = (
                cx is None or cy is None or cz is None or
                (cx == 0.0 and cy == 0.0 and cz == 0.0)
            )
            if _need_fallback:
                cx = _pc_cx if (cx is None or cx == 0.0) else cx
                cy = _pc_cy if (cy is None or cy == 0.0) else cy
                cz = _pc_cz if (cz is None or cz == 0.0) else cz
            # v3.2: use the measured active-site bounding box when the
            # caller supplied one. The old fixed 30 A cube was roughly 27000
            # A^3 for a pocket that is typically under 1000 A^3, which spreads
            # the same exhaustiveness over 27x the volume and is where the
            # surface false positives came from.
            _pol = resolve_box(pdb_id, measured_box=box,
                               n_heavy=int(n_heavy_atoms or 0),
                               is_metal_complex=bool(is_metal_complex))
            sx, sy, sz = _pol["size"]
            box_policy = _pol
            if _pol["is_channel"] and _pol["pore_residues"]:
                cx, cy, cz = pore_axis_centre(
                    receptor_pdbqt.with_suffix(".pdb")
                    if receptor_pdbqt else None,
                    _pol["pore_residues"], (cx, cy, cz))
            if _pol["exhaustiveness"] > exhaustiveness:
                exhaustiveness = _pol["exhaustiveness"]

        # ── FIX 1: Unique output path per compound+target ────────
        # Hash of ligand file content → different compounds → different output files
        # Prevents stale score from previous compound being read as current result
        import hashlib as _hl
        lig_hash  = int(_hl.md5(ligand_pdbqt.read_text(encoding="utf-8", errors="ignore")[:500].encode()).hexdigest(), 16) % 9999999
        out_pdbqt = ligand_pdbqt.parent / f"{pdb_id}_{lig_hash}_out.pdbqt"
        log_file  = out_pdbqt.with_suffix(".log")

        # Delete stale output files before each run
        for _stale in [out_pdbqt, log_file]:
            if _stale.exists():
                try: _stale.unlink()
                except Exception: pass
        cmd = [
            VINA_EXE,
            "--receptor",       str(receptor_pdbqt),
            "--ligand",         str(ligand_pdbqt),
            "--center_x",       f"{cx:.3f}",
            "--center_y",       f"{cy:.3f}",
            "--center_z",       f"{cz:.3f}",
            "--size_x",         f"{sx:.1f}",
            "--size_y",         f"{sy:.1f}",
            "--size_z",         f"{sz:.1f}",
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes",      "5",
            # v3.2: without --seed no run was bit-reproducible, so no result
            # in v3.1 could be regenerated and no regression could be
            # distinguished from search noise.
            "--seed",           str(int(seed)),
            "--out",            str(out_pdbqt),
            "--cpu",            "1",
            # --log removed: not supported by Vina 1.2.5 on this system
        ]

        proc   = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        stdout = (proc.stdout or "")
        stderr = (proc.stderr or "")
        # Read log file (most reliable on Windows — bypasses stdout/stderr ambiguity)
        log_content = ""
        if log_file.exists():
            try: log_content = log_file.read_text(encoding="utf-8", errors="ignore")
            except Exception: pass
        all_output = stdout + stderr + log_content

        # ── Parse score: 4 independent strategies ───────────────
        best_dg = None

        # Strategy A: PDBQT REMARK line (most reliable — Vina always writes this)
        # Format: "REMARK VINA RESULT:    -8.147      0.000      0.000"
        if out_pdbqt.exists() and out_pdbqt.stat().st_size > 10:
            for line in out_pdbqt.read_text(encoding="utf-8", errors="ignore").split("\n"):
                if "VINA RESULT" in line.upper():
                    for tok in line.split():
                        try:
                            v = float(tok)
                            if -30.0 < v < -0.5:
                                best_dg = round(v, 2)
                                break
                        except ValueError:
                            continue
                    if best_dg is not None:
                        break

        # Strategy B: score table — Vina 1.1.x, 1.2.x, GPU-Vina
        # Scans all_output (stdout+stderr) because Vina 1.2.x on Windows
        # writes its docking table to STDERR, not stdout.
        # Header:    "mode |   affinity | dist from best mode"  (1.2.x)
        #            "mode |   affinity | rmsd l.b.| rmsd u.b." (1.1.x)
        # Separator: "-----+------------+----------+----------"
        # Score row: "   1       -8.147       0.000       0.000"
        if best_dg is None:
            in_table = False
            for line in all_output.split("\n"):
                low = line.lower()
                if "affinity" in low and "|" in line:
                    in_table = True
                    continue
                if in_table and ("|" in line or "---" in line):
                    continue
                if in_table:
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[0].isdigit():
                        try:
                            v = float(parts[1])
                            if -30.0 < v < -0.5:
                                best_dg = round(v, 2)
                                break
                        except ValueError:
                            pass
                    elif parts:
                        in_table = False

        # Strategy C: affinity/kcal keyword lines in all_output
        if best_dg is None:
            for line in all_output.split("\n"):
                low = line.lower()
                if "affinity" in low or "kcal/mol" in low:
                    for tok in line.split():
                        try:
                            v = float(tok)
                            if -30.0 < v < -0.5:
                                best_dg = round(v, 2)
                                break
                        except ValueError:
                            continue
                    if best_dg is not None:
                        break

        # Strategy D: all tokens in all_output — last resort
        # Gated: only if Vina printed an affinity table at all
        if best_dg is None:
            if "affinity" in all_output.lower():
                for tok in all_output.split():
                    try:
                        v = float(tok)
                        if -20.0 < v < -2.0:
                            best_dg = round(v, 2)
                            break
                    except ValueError:
                        continue

        # Strategy E: scan ANY REMARK line in output PDBQT for negative float
        # Catches Vina builds that omit "VINA RESULT" text from REMARK
        if best_dg is None and out_pdbqt.exists() and out_pdbqt.stat().st_size > 10:
            for line in out_pdbqt.read_text(encoding="utf-8", errors="ignore").split("\n"):
                if line.startswith("REMARK"):
                    for tok in line.split():
                        try:
                            v = float(tok)
                            if -30.0 < v < -0.5:
                                best_dg = round(v, 2)
                                break
                        except ValueError:
                            continue
                    if best_dg is not None:
                        break

        if best_dg is None:
            # Return detailed debug info to help diagnose
            return {
                "available": True,
                "success":   False,
                "error": (
                    f"Vina ran but no score found. "
                    f"out_pdbqt: {out_pdbqt.exists()} "
                    f"({out_pdbqt.stat().st_size if out_pdbqt.exists() else 0}B). "
                    f"Vina output ({len(all_output)} chars): "
                    + (all_output[:300] if all_output else "[empty — check Vina installation]")
                ),
                "debug": {
                    "receptor_atoms":   len(atom_lines),
                    "receptor_size_kb": round(receptor_pdbqt.stat().st_size/1024, 1),
                    "stdout_lines":     stdout.count("\n"),
                    "stderr_lines":     stderr.count("\n"),
                    "all_output_lines": all_output.count("\n"),
                    "out_pdbqt_exists": out_pdbqt.exists(),
                    "out_pdbqt_size":   out_pdbqt.stat().st_size if out_pdbqt.exists() else 0,
                    "stdout_preview":   stdout[:300],
                    "stderr_preview":   stderr[:300],
                    "all_output_preview": all_output[:500],
                    "grid":             [round(cx,2),round(cy,2),round(cz,2),round(sx,1),round(sy,1),round(sz,1)],
                    "cmd":              " ".join(cmd),
                },
                "hint": (
                    "Most common causes: (1) receptor PDBQT has no usable ATOM records "
                    f"— this receptor has {len(atom_lines)} ATOM lines; "
                    "(2) ligand PDBQT is malformed; (3) grid box is wrong. "
                    "Check the debug.stdout_preview field for Vina's actual output. "
                    "Install OpenBabel for more robust file preparation."
                )
            }

        # ── Estimated docking Ki ──────────────────────────────────
        # One implementation, in bamstu_layers, so the constant cannot drift
        # between call sites. This copy used RT = 0.5922 while the Layer 1
        # function used 0.59248, which is small but means two places in one
        # report could disagree in the third digit for the same score.
        _kiq = layers.estimated_docking_ki(best_dg)
        ki_str = _kiq["display"] if _kiq.get("available") else "not available"

        best_pose = ""
        if out_pdbqt.exists():
            try:
                best_pose = out_pdbqt.read_text(encoding="utf-8", errors="ignore")[:3000]
            except Exception:
                pass

        return {
            "available":       True,
            "success":         True,
            "delta_g":         best_dg,
            "ki":              ki_str,
            "pose_pdbqt":      best_pose,
            "grid_centre":     [round(cx,2), round(cy,2), round(cz,2)],
            "grid_size":       [round(sx,1), round(sy,1), round(sz,1)],
            "exhaustiveness":  exhaustiveness,
            "box_policy":      box_policy,
            "estimated_docking_ki": _kiq,
            "receptor_atoms":  len(atom_lines),
            "engine":          "AutoDock Vina",
            "docking_mode":    "BLIND_WHOLE_PROTEIN" if blind else "TARGETED_ACTIVE_SITE",
            "seed":            int(seed),
            "reproducible":    True,
            # Path to the docked pose. Needed by /api/vvuq/redock to compute
            # RMSD against the crystallographic pose; v3.1 never exposed it,
            # which is part of why redocking validation was never run.
            "out_pdbqt":       str(out_pdbqt),
            "box": {
                "center": [round(cx, 3), round(cy, 3), round(cz, 3)],
                "size":   [sx, sy, sz],
                "volume_A3": round(sx * sy * sz, 1),
                "source": ("whole protein (blind)" if blind else
                           ("active-site residue bounding box + 6 A buffer"
                            if box else "default 24 A cube (no residue "
                                        "annotation available)")),
            },
            "metal_proxy": {
                "applied": True,
                "mapping": "all transition metals -> Zn atom type",
                "reason": ("AutoDock Vina 1.2 has no force-field parameters for "
                           "Ag, Au, Pt, Ru, Co, Cu and the rest, and aborts on "
                           "them during parsing."),
                "consequence": (
                    "Every metal contribution in this score is a zinc "
                    "contribution. Zn(II) has an ionic radius of 0.74 A against "
                    "1.15 A for Ag(I), prefers tetrahedral geometry where Ag(I) "
                    "is linear or trigonal, and is a borderline acid where Ag(I) "
                    "is soft and thiophilic. For a thiophilic centre the "
                    "substitution is systematically wrong in a known direction."
                ),
                "disclosure_required": (
                    "State this substitution in the Methods section. Disclosed, "
                    "it is a recognised limitation; undisclosed, it is not."
                ),
            },
            "charge_model": ("Gasteiger-Marsili (OpenBabel)"
                             if shutil.which("obabel") else
                             "formal +0.000 on all atoms — install OpenBabel"),
            "scoring_function_error_kcal": 1.5,
            "score_label": "AutoDock Vina 1.2 score (NOT a measured binding energy)",
        }

    except subprocess.TimeoutExpired:
        return {"available": True, "success": False,
                "error": "Vina timed out (>5 min). Try exhaustiveness=8."}
    except FileNotFoundError:
        return {"available": False,
                "error": f"Vina executable not found at: {VINA_EXE}"}
    except Exception as e:
        return {"available": True, "success": False, "error": str(e)[:300]}




def compute_selectivity(active_dg: float, blind_dg: float) -> dict:
    """
    Compare active-site vs blind docking ΔG to classify binding selectivity.

    Parameters
    ----------
    active_dg : float
        ΔG from targeted active-site docking (kcal/mol, negative = binding).
    blind_dg : float
        ΔG from blind whole-protein docking — the global energy minimum
        anywhere on the receptor surface (kcal/mol).

    Returns
    -------
    dict with:
        selectivity_ratio   : active_dg / blind_dg
                              Values near 1.0 → compound prefers active site.
                              Values > 1.0 → active site is WORSE than surface.
        binding_mode        : classification string
        selectivity_score   : 0–100 percentage scale for reporting
        promiscuous         : bool — True if compound binds surface >> active site
        interpretation      : human-readable sentence for thesis Discussion
        publication_note    : what to state in Methods

    Scientific basis
    ----------------
    If blind_ΔG << active_ΔG (i.e. the compound binds elsewhere more strongly),
    the compound is promiscuous — it associates with non-functional surface sites
    preferentially. This does not necessarily mean it has no activity, but:
    (a) the active-site ΔG over-estimates functional affinity
    (b) the compound would need to outcompete its own surface binding to reach Kcat
    (c) MD simulations and experimental IC50 may diverge from docking prediction

    Thresholds (consensus from Shoichet 2006 Drug Discov Today;
    Mysinger & Shoichet 2010; Baber et al. 2006):
        ratio ≥ 0.90 : ACTIVE_SITE_SELECTIVE   — within 10% of blind min
        ratio ≥ 0.75 : MILDLY_PROMISCUOUS       — binds ~1–2 kcal better elsewhere
        ratio < 0.75 : HIGHLY_PROMISCUOUS        — strong off-target preference
    """
    if blind_dg == 0 or blind_dg is None:
        return {"error": "blind_dg is zero or missing — cannot compute selectivity"}

    ratio = round(active_dg / blind_dg, 4)
    score = round(min(100.0, max(0.0, ratio * 100)), 1)

    if ratio >= 0.90:
        mode   = "ACTIVE_SITE_SELECTIVE"
        flag   = False
        interp = (
            f"The compound docks at the active site (ΔG = {active_dg:.2f} kcal/mol) "
            f"within 10% of its global minimum energy on this receptor "
            f"(blind ΔG = {blind_dg:.2f} kcal/mol, selectivity ratio = {ratio:.2f}). "
            "This indicates genuine active-site preference and supports the "
            "mechanistic interpretation reported in the docking section."
        )
        pub = (
            "Blind whole-protein docking confirmed active-site preference "
            f"(selectivity ratio = {ratio:.2f}; active-site ΔG = {active_dg:.2f}, "
            f"global minimum ΔG = {blind_dg:.2f} kcal/mol)."
        )
    elif ratio >= 0.75:
        mode   = "MILDLY_PROMISCUOUS"
        flag   = True
        interp = (
            f"The compound shows mild surface promiscuity: blind docking "
            f"(ΔG = {blind_dg:.2f} kcal/mol) is moderately better than active-site "
            f"docking (ΔG = {active_dg:.2f} kcal/mol), selectivity ratio = {ratio:.2f}. "
            "The active-site affinity is still plausible, but this should be "
            "acknowledged in the Discussion and supported by experimental IC50 data."
        )
        pub = (
            "Mild surface promiscuity observed (selectivity ratio = "
            f"{ratio:.2f}; active-site ΔG = {active_dg:.2f}, global min "
            f"ΔG = {blind_dg:.2f} kcal/mol). Active-site binding is consistent "
            "but experimental validation is required."
        )
    else:
        mode   = "HIGHLY_PROMISCUOUS"
        flag   = True
        interp = (
            f"The compound is highly promiscuous: its best binding site on this "
            f"receptor is NOT the active site. Blind docking ΔG = {blind_dg:.2f} "
            f"kcal/mol vs active-site ΔG = {active_dg:.2f} kcal/mol "
            f"(selectivity ratio = {ratio:.2f}). The active-site ΔG likely "
            "overestimates functional affinity. MD simulation or experimental "
            "binding assays are required before citing this docking score."
        )
        pub = (
            "HIGH PROMISCUITY FLAG: compound preferentially binds outside the "
            f"active site (selectivity ratio = {ratio:.2f}; active-site "
            f"ΔG = {active_dg:.2f} kcal/mol, global minimum ΔG = "
            f"{blind_dg:.2f} kcal/mol). This result should not be cited as "
            "evidence of active-site inhibition without further experimental support."
        )

    return {
        "active_dg":         active_dg,
        "blind_dg":          blind_dg,
        "selectivity_ratio": ratio,
        "selectivity_score": score,
        "binding_mode":      mode,
        "promiscuous":       flag,
        "interpretation":    interp,
        "publication_note":  pub,
    }


# ── FastAPI app + Groq client ────────────────────────────────
app = FastAPI(title="BAMSTU Digital Pathogen Pro", version="3.1")
app.add_middleware(
    CORSMiddleware,
    # v3.4.5: a wildcard origin on a server that reads and deletes local
    # session files means any page the user visits can drive this API on
    # 127.0.0.1. The frontend is a local file, so a null/file origin plus
    # localhost is all that is needed.
    allow_origins=["null", "http://localhost:8080", "http://127.0.0.1:8080",
                   "http://localhost:5500", "http://127.0.0.1:5500"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"], allow_credentials=True
)

# ── API key resolution ───────────────────────────────────────────
# The argument to os.environ.get is the NAME of a variable, never the secret.
# A key pasted in as the name does two things at once: it commits the secret
# to a file that goes into thesis appendices, and it looks up a variable that
# will never exist, so the client silently stays None with no error to
# explain the silence. Both have happened in this codebase, which is why the
# tripwire below exists rather than a comment asking people not to do it.

GROQ_KEY_NAMES = ("BAMSTU_GROQ_API_KEY", "GROQ_API_KEY")
GROQ_KEY_SOURCE = "not resolved"
SECRET_PREFIXES = ("gsk_", "sk-", "AKIA", "ghp_", "xoxb-")


def _resolve_groq_key():
    """Read the key from the environment, and refuse a pasted literal."""
    global GROQ_KEY_SOURCE
    bad = [n for n in GROQ_KEY_NAMES if n.startswith(SECRET_PREFIXES)]
    if bad:
        GROQ_KEY_SOURCE = (
            "REFUSED: a literal API key has been pasted into GROQ_KEY_NAMES "
            "in main.py. That is a variable NAME, so the lookup can never "
            "succeed, and the key is now in the source file. Revoke that key "
            "at the provider, restore the name, and set the value in the "
            "environment instead.")
        print("[BAMSTU] " + GROQ_KEY_SOURCE, flush=True)
        return None
    for name in GROQ_KEY_NAMES:
        val = (os.environ.get(name) or "").strip()
        if val:
            GROQ_KEY_SOURCE = f"environment variable {name}"
            return val
    GROQ_KEY_SOURCE = (
        "absent — set BAMSTU_GROQ_API_KEY in the environment before starting "
        "the server. The narrative layer falls back to a template that makes "
        "no mechanistic claims; docking is unaffected.")
    return None


try:
    from groq import Groq
    # Read from the environment. A key committed to a source file is a key
    # published to everyone the file is ever sent to, and this file is bound
    # into thesis appendices. Absent, the narrative layer degrades to a
    # template that makes no mechanistic claims; docking is unaffected.
    # The argument to os.environ.get is the NAME of a variable, never the
    # secret itself. A key pasted in here does two things at once: it is
    # committed to source, and it is looked up as a variable name that will
    # never exist, so the client silently stays None and the narrative layer
    # goes quiet with no error to explain why. Both have happened.
    GROQ_API_KEY = _resolve_groq_key()
    groq_client = Groq(api_key=GROQ_API_KEY) if GROQ_API_KEY else None
except Exception as _groq_exc:
    groq_client = None
    GROQ_KEY_SOURCE = f"groq package unavailable ({_groq_exc})"

# ── BAMSTU Literature Module ─────────────────────────────────
# Queries PubMed, Semantic Scholar, Europe PMC, OpenAlex
# Returns real abstracts + DOIs for Groq MoA grounding
# Dr. James P.T. Matshwele, BAMSTU Labs, UB, 2026

import urllib.request, urllib.parse, json, ssl, time

_SSL = ssl.create_default_context()

# See bamstu_vvuq: a VPN sets HTTP_PROXY machine-wide and urllib obeys it,
# so RCSB becomes unreachable from this process while the browser is fine.
_DIRECT_OPENER = urllib.request.build_opener(
    urllib.request.ProxyHandler({}),
    urllib.request.HTTPSHandler(context=_SSL))
_SSL.check_hostname = False
_SSL.verify_mode = ssl.CERT_NONE

def _get(url: str, timeout=8) -> dict | None:
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "BAMSTU-DigitalPathogen/3.0 (bamstu@ub.bw)"
        })
        with urllib.request.urlopen(req, context=_SSL, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8", errors="ignore"))
    except Exception:
        return None

def _build_query(target_name: str, pathogen: str, compound_type: str) -> str:
    """Build a targeted PubMed/scholar search query."""
    # Strip long parenthetical from target name
    short = target_name.split("(")[0].strip()
    # Pathogen short name
    bug = "MRSA" if "aureus" in pathogen else "Klebsiella"
    # Compound
    ctype = "metal complex" if compound_type == "metal" else "antibacterial"
    return f"{short} {bug} {ctype} inhibition"

def search_pubmed(query: str, max_results=4) -> list[dict]:
    """Search PubMed and fetch abstracts."""
    results = []
    q = urllib.parse.quote(query)
    # Step 1: search for PMIDs
    search_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&term={q}&retmax={max_results}&retmode=json&sort=relevance"
    )
    data = _get(search_url)
    if not data:
        return results
    pmids = data.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return results

    # Step 2: fetch summaries (title + abstract)
    ids = ",".join(pmids)
    fetch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={ids}&rettype=abstract&retmode=json"
    )
    # Use esummary for structured JSON
    summary_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=pubmed&id={ids}&retmode=json"
    )
    summ = _get(summary_url)
    if summ and "result" in summ:
        for pmid in pmids:
            art = summ["result"].get(pmid, {})
            title = art.get("title", "")
            authors = art.get("authors", [])
            author_str = authors[0].get("name","") if authors else "Unknown"
            year = art.get("pubdate","")[:4]
            doi = ""
            for artid in art.get("articleids",[]):
                if artid.get("idtype") == "doi":
                    doi = artid.get("value","")
                    break
            if title:
                results.append({
                    "source": "PubMed",
                    "title": title,
                    "authors": author_str,
                    "year": year,
                    "doi": doi,
                    "pmid": pmid,
                    "abstract": "",  # filled below
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
    return results

def search_semantic_scholar(query: str, max_results=4) -> list[dict]:
    """Search Semantic Scholar for open-access papers with abstracts."""
    q = urllib.parse.quote(query)
    url = (
        f"https://api.semanticscholar.org/graph/v1/paper/search"
        f"?query={q}&limit={max_results}"
        f"&fields=title,authors,year,abstract,openAccessPdf,externalIds,publicationVenue"
    )
    data = _get(url)
    if not data:
        return []
    results = []
    for p in data.get("data", []):
        doi = p.get("externalIds", {}).get("DOI", "")
        authors = p.get("authors", [])
        author_str = authors[0].get("name","") + " et al." if len(authors) > 1 else (authors[0].get("name","") if authors else "Unknown")
        pdf_url = (p.get("openAccessPdf") or {}).get("url","")
        abstract = (p.get("abstract") or "")[:600]
        if p.get("title") and abstract:  # only include if we have an abstract
            results.append({
                "source": "Semantic Scholar",
                "title": p.get("title",""),
                "authors": author_str,
                "year": str(p.get("year","")),
                "doi": doi,
                "abstract": abstract,
                "pdf_url": pdf_url,
                "url": f"https://doi.org/{doi}" if doi else ""
            })
    return results

def search_europe_pmc(query: str, max_results=3) -> list[dict]:
    """Search Europe PMC — strong for open-access biomedical literature."""
    q = urllib.parse.quote(query)
    url = (
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        f"?query={q}&format=json&pageSize={max_results}"
        f"&resultType=core&sort=CITED"
    )
    data = _get(url)
    if not data:
        return []
    results = []
    for p in data.get("resultList", {}).get("result", []):
        abstract = (p.get("abstractText") or "")[:600]
        doi = p.get("doi","")
        if p.get("title") and abstract:
            results.append({
                "source": "Europe PMC",
                "title": p.get("title",""),
                "authors": p.get("authorString",""),
                "year": str(p.get("pubYear","")),
                "doi": doi,
                "abstract": abstract,
                "url": f"https://doi.org/{doi}" if doi else f"https://europepmc.org/article/MED/{p.get('pmid','')}"
            })
    return results

def search_openalex(query: str, max_results=3) -> list[dict]:
    """Search OpenAlex — 250M+ scholarly works."""
    q = urllib.parse.quote(query)
    url = (
        f"https://api.openalex.org/works"
        f"?search={q}&per-page={max_results}"
        f"&select=title,authorships,publication_year,abstract_inverted_index,doi,open_access"
        f"&mailto=bamstu@ub.bw"
    )
    data = _get(url)
    if not data:
        return []
    results = []
    for p in data.get("results", []):
        doi = (p.get("doi") or "").replace("https://doi.org/","")
        authors = p.get("authorships",[])
        author_str = ""
        if authors:
            a0 = authors[0].get("author",{}).get("display_name","")
            author_str = a0 + (" et al." if len(authors)>1 else "")
        # Reconstruct abstract from inverted index
        inv = p.get("abstract_inverted_index") or {}
        if inv:
            words = [""] * (max(max(pos) for pos in inv.values()) + 1)
            for word, positions in inv.items():
                for pos in positions:
                    if pos < len(words):
                        words[pos] = word
            abstract = " ".join(words)[:600]
        else:
            abstract = ""
        oa_url = (p.get("open_access") or {}).get("oa_url","")
        if p.get("title") and abstract:
            results.append({
                "source": "OpenAlex",
                "title": p.get("title",""),
                "authors": author_str,
                "year": str(p.get("publication_year","")),
                "doi": doi,
                "abstract": abstract,
                "pdf_url": oa_url,
                "url": f"https://doi.org/{doi}" if doi else oa_url
            })
    return results

def fetch_literature(target_name: str, pathogen: str,
                     has_metals: bool, max_papers=6) -> list[dict]:
    """
    Master function: queries all sources, deduplicates by title,
    returns up to max_papers results with abstracts.
    """
    compound_type = "metal" if has_metals else "organic"
    query = _build_query(target_name, pathogen, compound_type)

    # Query all sources in parallel-ish (sequential, fast enough)
    all_papers = []
    all_papers.extend(search_semantic_scholar(query, 4))
    time.sleep(0.15)  # polite delay
    all_papers.extend(search_europe_pmc(query, 3))
    time.sleep(0.15)
    all_papers.extend(search_pubmed(query, 3))
    time.sleep(0.15)
    all_papers.extend(search_openalex(query, 3))

    # Deduplicate by title similarity (first 40 chars)
    seen = set()
    unique = []
    for p in all_papers:
        key = p["title"][:40].lower().strip()
        if key not in seen and p.get("abstract"):
            seen.add(key)
            unique.append(p)

    # Return top results with abstracts, sorted: those with DOI first
    unique.sort(key=lambda x: (0 if x.get("doi") else 1))
    return unique[:max_papers]

def format_citations(papers: list[dict]) -> list[str]:
    """Format papers as citation strings for MoA references."""
    refs = []
    for p in papers:
        author = p.get("authors","Unknown")
        year   = p.get("year","")
        title  = p.get("title","")[:80]
        doi    = p.get("doi","")
        source = p.get("source","")
        url    = p.get("url","") or (f"https://doi.org/{doi}" if doi else "")
        if doi:
            refs.append(f"{author} ({year}). {title}. DOI:{doi} [{source}]")
        else:
            refs.append(f"{author} ({year}). {title}. {url} [{source}]")
    return refs

def build_literature_context(papers: list[dict]) -> str:
    """Build a compact context string for Groq from paper abstracts."""
    if not papers:
        return ""
    ctx_parts = []
    for i, p in enumerate(papers[:5], 1):
        ctx_parts.append(
            f"[Paper {i}] {p['authors']} ({p['year']}). {p['title']}.\n"
            f"Abstract: {p['abstract'][:350]}"
        )
    return "\n\n".join(ctx_parts)


# ── Biological target database (active site residue lists) ───
# active_residues: known catalytic/binding residues for centroid calculation
TARGETS = {
    "Staphylococcus aureus": [
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # BAMSTU Digital Pathogen Pro v3.2 — MRSA / MDR S. aureus Panel
        # Rebuilt August 2026. Every entry now carries the same four fields
        # the K. pneumoniae panel already had, because the Gram-positive
        # panel was the weaker half and was being read as though it were not:
        #   organism           what the PDB actually contains
        #   provenance         the deposition citation
        #   docking_semantics  what a negative score here means
        #   site_note          what the residue list points at
        #
        # Verify at runtime with GET /api/verify_targets, which checks the
        # organism string and confirms every configured residue exists in
        # the chain. Residue numbering changes between depositions; a list
        # that was right two years ago can point at nothing today, and the
        # centroid then falls back silently to the whole-protein centre.
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

        # ── METHICILLIN RESISTANCE DETERMINANT: PBP2a (mecA) ─────────────
        # The protein that defines MRSA. Acquired via SCCmec, low affinity
        # for beta-lactams, completes cross-linking when PBP1-4 are acylated.
        #
        # CORRECTION: v3.1 listed active_residues [645,646,647,648,649,650].
        # PBP2a is 668 residues; 645-650 is C-terminal tail, roughly 40 A
        # from the catalytic serine. Every PBP2a docking run in v3.1
        # therefore targeted a surface loop, not the transpeptidase site.
        # Catalytic machinery: Ser403 (nucleophile), Lys406, Ser462-Asn464
        # (SXN motif), Lys597-Thr598-Gly599 (KTG motif), Thr600.
        # Bouley et al. 2015, J Am Chem Soc 137:1738 (4CJN, quinazolinone).
        {"target_name": "PBP2a — transpeptidase active site (mecA, MRSA determinant)",
         "target_pdb": "4CJN", "compartment": "Cell Wall", "base_dg": -6.75,
         "active_residues": [403, 406, 462, 464, 597, 598, 599, 600],
         "chain": "A",
         "organism": "Staphylococcus aureus subsp. aureus Mu50 (native MRSA)",
         "uniprot": "Q53707 (mecA / PBP2a)",
         "provenance": "Bouley et al. 2015, J Am Chem Soc 137:1738; PDB 4CJN",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "Transpeptidase site around catalytic Ser403. PBP2a also has "
                      "an allosteric site ~60 A away (Otero 2013 PNAS, PDB 3ZG0); "
                      "a compound docked here is NOT an allosteric hit and the two "
                      "must not be conflated in the Discussion."},

        # ── BETA-LACTAM RESISTANCE ARM: PBP2a allosteric domain ──────────
        # Same structure, different site. Occupancy of the allosteric domain
        # opens the active site; this is the mechanism ceftaroline exploits.
        # Separated deliberately so a hit at one site is never reported as a
        # hit at the other. Allosteric ligand contacts: Lys273, Tyr297,
        # Arg241, Glu239, His293 (muramic acid / ceftaroline site).
        {"target_name": "PBP2a — allosteric domain (ceftaroline-type site)",
         "target_pdb": "4CJN", "compartment": "Cell Wall", "base_dg": -6.0,
         "active_residues": [239, 241, 273, 293, 297],
         "chain": "A",
         "organism": "Staphylococcus aureus subsp. aureus Mu50 (native MRSA)",
         "uniprot": "Q53707 (mecA / PBP2a)",
         "provenance": "Otero et al. 2013, PNAS 110:16808; PDB 3ZG0/4CJN",
         "docking_semantics": "allosteric_modulation",
         "site_note": "Binding here does NOT inhibit. It opens the active site. "
                      "A favourable score is a modulator hypothesis, not a kill "
                      "mechanism, and requires the active-site entry to be scored "
                      "as well before any claim is made."},

        # ── DNA REPLICATION: GyrB ATPase domain ──────────────────────────
        # 4URO is the S. aureus GyrB 24 kDa N-terminal ATPase domain with
        # novobiocin bound (Lu et al. 2014, ACS Chem Biol 9:2023).
        # CORRECTION: v3.1 residues [42,43,44,85,86,136,165] were a mix of
        # numbering conventions. The novobiocin/ATP pocket in this construct
        # is Asn54, Glu58, Asp81, Arg84, Gly85, Ile86, Arg144.
        {"target_name": "DNA Gyrase B — ATPase pocket (novobiocin site)",
         "target_pdb": "4URO", "compartment": "Cytoplasm", "base_dg": -6.25,
         "active_residues": [54, 58, 81, 84, 85, 86, 144],
         "chain": "A",
         "organism": "Staphylococcus aureus (native)",
         "uniprot": "P0A0K8 (gyrB)",
         "provenance": "Lu et al. 2014, ACS Chem Biol 9:2023; PDB 4URO",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "ATP-competitive pocket of the 24 kDa N-terminal construct. "
                      "This is NOT the fluoroquinolone site — that sits on GyrA at "
                      "the DNA gate and is a separate target entry."},

        # ── FOLATE PATHWAY: DHPS ─────────────────────────────────────────
        # CORRECTION: v3.1 used 3TYE. That entry belongs to the St Jude DHPS
        # catalytic-mechanism series and is not an S. aureus deposition, so
        # the panel was docking against the wrong organism while labelling
        # the result S. aureus. 1AD4 is the S. aureus DHPS structure
        # (Hampele et al. 1997, J Mol Biol 268:21).
        # Pterin pocket: Asn11, Asp84, Asp167, Lys203, Arg239.
        # Sulfa-resistance positions F17, S18, T51 included so a resistant
        # background can be modelled explicitly.
        {"target_name": "Dihydropteroate Synthase (DHPS) — pterin/pABA pocket",
         "target_pdb": "1AD4", "compartment": "Cytoplasm", "base_dg": -5.75,
         "active_residues": [11, 51, 84, 167, 203, 239],
         "chain": "A",
         "organism": "Staphylococcus aureus (native)",
         "uniprot": "P0A017 (folP / sulA)",
         "provenance": "Hampele et al. 1997, J Mol Biol 268:21; PDB 1AD4",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "Pterin/pABA pocket. F17 and S18 are recurrent sulfonamide-resistance "
                      "positions (Griffith et al. 2018, Front Microbiol 9:1369) but are NOT "
                      "resolved in chain A of this deposition, so they were removed from the "
                      "centroid in v3.2.1 after target verification flagged them. T51 is "
                      "resolved and is retained. To model an F17L/S18L resistant background, "
                      "a structure that resolves the N-terminal region is needed."},

        # ── FATTY ACID SYNTHESIS: FabI ───────────────────────────────────
        # 4ALI: S. aureus FabI with NADP and triclosan (Schiebel et al. 2012,
        # Structure 20:802). Catalytic YX6K motif: Tyr157, Lys164.
        {"target_name": "Enoyl-ACP Reductase (FabI) — NADPH/triclosan site",
         "target_pdb": "4ALI", "compartment": "Cytoplasm", "base_dg": -5.5,
         "active_residues": [93, 95, 147, 157, 164, 204],
         "chain": "A",
         "organism": "Staphylococcus aureus (native)",
         "uniprot": "Q9FDG0 (fabI)",
         "provenance": "Schiebel et al. 2012, Structure 20:802; PDB 4ALI",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "Ternary site: inhibitor stacks on the nicotinamide ring of "
                      "bound NADP+. The cofactor is retained during receptor prep "
                      "for this target; docking into the apo pocket gives a "
                      "different and much weaker answer."},

        # ── VIRULENCE / SURFACE ANCHORING: Sortase A ─────────────────────
        # 1T2W: S. aureus SrtA solution structure. Catalytic triad
        # His120, Cys184, Arg197 — the Cys is a soft thiolate and is the
        # most credible Ag(I) target in the whole Gram-positive panel.
        # CORRECTION: v3.1 residues [119..122, 166..168] missed Cys184 and
        # Arg197 entirely, i.e. missed two thirds of the catalytic triad.
        {"target_name": "Sortase A (SrtA) — Cys184 catalytic triad",
         "target_pdb": "1T2W", "compartment": "Cell Wall", "base_dg": -5.0,
         "active_residues": [120, 180, 184, 197],
         "chain": "A",
         "organism": "Staphylococcus aureus (native)",
         "uniprot": "Q9RI31 (srtA)",
         "provenance": "Zong et al. 2004, J Biol Chem 279:31383; PDB 1T2W",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "Cys184 thiolate is a soft HSAB donor. For Ag(I) and Au(I) "
                      "this is the mechanistically strongest hypothesis in the "
                      "panel. Note the caveat that Vina samples with the metal "
                      "typed as Zn, which is a borderline acid and systematically "
                      "under-weights the Ag-S interaction it is standing in for."},

        # ── ADHESION: Clumping Factor A ──────────────────────────────────
        # CORRECTION: v3.1 called 2WTP the "ClfA Fibronectin-Binding
        # Protein". ClfA binds FIBRINOGEN. Fibronectin binding is done by
        # FnBPA/FnBPB, different genes and different proteins. Replaced with
        # 1N67 (Deivanayagam et al. 2002, EMBO J 21:6660) and renamed.
        # Latching trench: Tyr256, Pro336, Tyr338, Lys389.
        {"target_name": "Clumping Factor A (ClfA) — fibrinogen latching trench",
         "target_pdb": "1N67", "compartment": "Cell Surface", "base_dg": -4.5,
         "active_residues": [256, 310, 336, 338, 389],
         "chain": "A",
         "organism": "Staphylococcus aureus (native)",
         "uniprot": "Q53653 (clfA)",
         "provenance": "Deivanayagam et al. 2002, EMBO J 21:6660; PDB 1N67",
         "docking_semantics": "protein_protein_interface_blockade",
         "site_note": "This is an adhesin, not an enzyme. A favourable score means "
                      "the compound may occlude fibrinogen gamma-chain docking. It "
                      "is an anti-virulence hypothesis and predicts no growth "
                      "inhibition, so it must never be pooled with MIC-linked "
                      "targets in a single ranked table. P336 and Y338 are the "
                      "residues whose mutation abolishes binding in vivo "
                      "(Josefsson et al. 2008, PLoS One 3:e2206)."},

        # ── CENTRAL METABOLISM: Pyruvate Kinase ──────────────────────────
        # Retained from v3.1 but flagged: the organism of 3T05 has not been
        # confirmed against RCSB in this revision. /api/verify_targets will
        # settle it. Until it returns PASS this entry is excluded from
        # ranked output by default (see include_unverified).
        {"target_name": "Pyruvate Kinase (PykA) — effector/active site",
         "target_pdb": "3T05", "compartment": "Cytoplasm", "base_dg": -5.0,
         "active_residues": [72, 73, 74, 220, 221, 222],
         "chain": "A",
         "organism": "UNVERIFIED — confirm via /api/verify_targets before use",
         "uniprot": "unconfirmed",
         "provenance": "Zoraghi et al. 2011, J Biol Chem 286:44716 (association "
                       "not independently confirmed in this revision)",
         "docking_semantics": "enzyme_inhibition",
         "verified": False,
         "site_note": "Held back from ranked output until the RCSB organism check "
                      "passes. An unverified entry that silently ranks is exactly "
                      "how the DHPS organism error survived three versions."},

        # ── MEMBRANE / CELL ENVELOPE: LtaS ───────────────────────────────
        # CORRECTION: v3.1 used 3L7H, which could not be confirmed as the
        # S. aureus LtaS extracellular domain. Replaced with 2W5T, the
        # published eLtaS structure (Lu et al. 2009, PNAS 106:1584).
        # Mn(II) site: Glu255, Thr300 (catalytic), Asp475, His476.
        {"target_name": "Lipoteichoic Acid Synthase (eLtaS) — Mn(II) catalytic site",
         "target_pdb": "2W5T", "compartment": "Inner Membrane", "base_dg": -4.5,
         "active_residues": [255, 300, 475, 476],
         "chain": "A",
         "organism": "Staphylococcus aureus (native, extracellular domain)",
         "uniprot": "Q2G1S6 (ltaS)",
         "provenance": "Lu et al. 2009, PNAS 106:1584; PDB 2W5T",
         "docking_semantics": "enzyme_inhibition",
         "site_note": "Mn(II)-dependent phosphoglycerol transferase. A divalent "
                      "metal sits in this pocket natively, so a metal complex may "
                      "score well by displacing it — that is a metal-substitution "
                      "hypothesis, distinct from occupying the substrate channel, "
                      "and the two cannot be told apart from the score alone."},

        # ── PORE-FORMING TOXIN: alpha-Hemolysin ──────────────────────────
        # Anti-virulence target, present in build_receptors.py but missing
        # from the main panel, so the two files disagreed. Now reconciled.
        {"target_name": "alpha-Hemolysin (Hla) — cap/rim domain",
         "target_pdb": "7AHL", "compartment": "Secreted", "base_dg": -4.75,
         "active_residues": [103, 109, 130, 132, 213],
         "chain": "A",
         "organism": "Staphylococcus aureus (native, heptameric pore)",
         "uniprot": "P09616 (hla)",
         "provenance": "Song et al. 1996, Science 274:1859; PDB 7AHL",
         "docking_semantics": "assembly_blockade",
         "site_note": "7AHL is the assembled heptamer. Docking here models "
                      "interference with pore assembly or lumen occlusion, not "
                      "enzyme inhibition. Anti-virulence readout: expect no MIC "
                      "shift, so validate by haemolysis assay rather than by "
                      "broth microdilution."},
    ],
        "Klebsiella pneumoniae": [
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # BAMSTU Digital Pathogen Pro v3.3 — MDR K. pneumoniae Panel
        # Panel reviewed June 2026 by multi-professor examination panel.
        # Covers: ALL THREE carbapenemase classes (A/B/D), dominant efflux
        # pump, essential biosynthesis, intrinsic resistance, and capsule
        # export. Together these represent the clinically complete MDR-Kp
        # resistome as defined by WHO/ESKAPE 2024 priority pathogen criteria.
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

        # ── CLASS A CARBAPENEMASE: KPC-2 ─────────────────────────────────
        # Plasmid-borne (pKpQIL/IncFIIK), ST258/ST512/ST11 backgrounds.
        # Ser70 catalytic nucleophile; Lys73 general base; oxyanion hole
        # (Ser70-NH, Thr237-OH). KPC-2 is the dominant US/Asian CRKP enzyme.
        # Ke et al. 2012, Antimicrob Agents Chemother 56:2295. [K. pneumoniae]
        # NOTE: Interpret hits as relevant only for KPC+ genotyped strains.
        {"target_name": "KPC-2 Carbapenemase — Class A (KPC+ strains)",
         "target_pdb": "3RXX", "compartment": "Periplasm", "base_dg": -6.25,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native, plasmid-borne KPC-2)",
         "uniprot": "Q9F663 (blaKPC-2)",
         "provenance": "Ke et al. 2012, Antimicrob Agents Chemother 56:2295; PDB 3RXX",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [69, 70, 71, 130, 166, 234, 235]},

        # ── CLASS B CARBAPENEMASE: NDM-1 ─────────────────────────────────
        # Plasmid-borne (IncA/C, IncF), ST14/ST147/ST231/ST383 backgrounds.
        # Di-Zn active site: Zn1(His120/His122/His189), Zn2(Asp124/Cys208/His250).
        # Cys208 is a soft donor — Ag(I) and Cu(II) complexes displace Zn(II)
        # at Cys via HSAB competition. Mechanistic basis for Ag(I) work.
        # Zhang et al. 2012, J Med Chem; PDB 4EXS, L-captopril inhibitor complex.
        # Covers KPC-negative NDM+ MDR CRKP lineages. [K. pneumoniae origin]
        {"target_name": "NDM-1 Metallo-beta-Lactamase — Class B (NDM+ strains)",
         "target_pdb": "4EXS", "compartment": "Periplasm", "base_dg": -6.75,
         "chain": "A",
         "organism": "Klebsiella pneumoniae origin, NDM-1 (native enzyme)",
         "uniprot": "C7C4X0 (blaNDM-1)",
         "provenance": "Zhang & Hao 2011, FASEB J 25:2574; L-captopril complex; PDB 4EXS",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [120, 122, 124, 189, 208, 250]},
        # ── OXA-48 (3HBR) REMOVED in v3.2.1 ──────────────────────────
        # ATCC BAA-1705 is the KPC-POSITIVE Klebsiella pneumoniae
        # reference strain; BAA-1706 is the KPC-negative counterpart.
        # Neither carries blaOXA-48. Ranking a compound against a Class D
        # carbapenemase the validation organisms do not express produced
        # a top hit against an enzyme that is not present, and in the
        # first live run that hit ranked first at -9.2 kcal/mol.
        #
        # The entry also failed target verification independently:
        # configured residues 73 and 150 are absent from chain A of 3HBR.
        #
        # Restore ONLY if the panel is repointed at an OXA-48-producing
        # isolate, and then re-derive the residue list against the
        # deposition rather than reusing the numbers removed here.

        # ── LIPID A BIOSYNTHESIS: LpxH ───────────────────────────────────
        # Universal in Enterobacteriaceae (kpn:KPN_00117, conserved).
        # Mn(II)-dependent pyrophosphohydrolase. No mammalian homolog.
        # Not plasmid-borne, not subject to resistance by porin loss.
        # AZ1 sulfonyl piperazine active in CRKP in vitro (Cho et al. 2020).
        # Cho et al. 2020, PNAS 117:4109. [K. pneumoniae, Mn(II) confirmed]
        {"target_name": "LpxH — Lipid A Biosynthesis (Mn-dependent, universal)",
         "target_pdb": "6PJ3", "compartment": "Inner Membrane", "base_dg": -6.5,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native, Mn(II)-dependent LpxH)",
         "uniprot": "A6T4Z5 (lpxH)",
         "provenance": "Cho et al. 2020, PNAS 117:4109; PDB 6PJ3",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [9, 10, 33, 73, 117, 159]},

        # ── DNA REPLICATION: Topoisomerase IV ParC (NBTI site) ───────────
        # ParC is GyrA-homologous subunit of Topo IV (NOT GyrA itself).
        # QRDR mutations (ParC S80I, GyrA S83L) confer fluoroquinolone
        # resistance but do NOT affect the NBTI compound binding site at
        # the DNA-gate interface (Bax et al. 2010, Nat Chem Biol 6:716).
        # 9KGT: K. pneumoniae Topo IV + DNA + BWC0977 NBTI inhibitor, 2025.
        # Defensible for NBTI-class compounds; NOT for quinolone comparison.
        {"target_name": "Topoisomerase IV ParC — NBTI binding site (not QRDR)",
         "target_pdb": "9KGT", "compartment": "Cytoplasm", "base_dg": -6.0,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native, Topo IV + DNA + BWC0977)",
         "uniprot": "ParC subunit",
         "provenance": "Bax et al. 2010, Nat Chem Biol 6:716 for the NBTI site; PDB 9KGT",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [83, 84, 87, 426, 447, 448]},

        # ── EFFLUX PUMP: OqxB (RND family) ──────────────────────────────
        # OqxAB is chromosomally encoded in K. pneumoniae and detected in
        # >90% of clinical isolates (Bialek-Davenet et al. 2015 JAC;
        # clinical survey data: >95% oqxA/oqxB positive in CRKP).
        # Confers resistance to quinolones, nitrofurantoin, tigecycline,
        # chloramphenicol. Contributes to pan-MDR phenotype independently
        # of carbapenemase carriage — present in KPC, NDM AND OXA isolates.
        # PDB 7CZ9: OqxB + DDM detergent at 1.85 Å. K. pneumoniae native.
        # Bharatham et al. 2021, Nat Commun 12:5484.
        # Ag(I) mechanism at OqxB: substrate-channel blockade at hydrophobic
        # binding pocket; distinct from HSAB metalloenzyme mechanism at NDM-1.
        # This must be stated explicitly in Discussion.
        {"target_name": "OqxB Efflux Pump — RND family (chromosomal, >90% Kp)",
         "target_pdb": "7CZ9", "compartment": "Inner Membrane", "base_dg": -5.5,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native OqxB, 1.85 A)",
         "uniprot": "A0A0E1RZP6 (oqxB)",
         "provenance": "Bharatham et al. 2021, Nat Commun 12:5484; PDB 7CZ9",
         "docking_semantics": "efflux_substrate_channel",
         "active_residues": [134, 176, 178, 179, 616, 618, 621]},

        # ── INTRINSIC RESISTANCE: FosA ───────────────────────────────────
        # Chromosomal fosA5/fosA6 alleles intrinsic to all K. pneumoniae
        # (Ito et al. 2017). Mn(II)-dependent glutathione transferase.
        # Inhibiting FosA restores fosfomycin susceptibility — clinically
        # relevant for last-resort fosfomycin + meropenem combinations.
        # Klontz et al. 2017, Antimicrob Agents Chemother 61:e01572-17.
        {"target_name": "FosA — Fosfomycin Resistance Enzyme (intrinsic, Mn-dep.)",
         "target_pdb": "5V3D", "compartment": "Cytoplasm", "base_dg": -5.25,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native FosA, Mn(II)-dependent)",
         "uniprot": "fosA5/fosA6",
         "provenance": "Klontz et al. 2017, Antimicrob Agents Chemother 61:e01572-17; PDB 5V3D",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [7, 8, 9, 58, 65, 110]},

        # ── CAPSULE EXPORT: Wza (E. coli homolog, >70% Kp identity) ─────
        # Wza forms the outer membrane trans-envelope channel for capsular
        # polysaccharide export. Present in all encapsulated K. pneumoniae.
        # No native K. pneumoniae Wza crystal structure exists in PDB (June
        # 2026 exhaustive search). PDB 2J58 (E. coli K-30, Dong et al. 2006,
        # PNAS 103:2066) is the only reference structure; >70% sequence
        # identity to Kp Wza. Explicitly labeled as E. coli homolog.
        # PREVIOUS ENTRY (7W1D, K1 phage lyase): removed — wrong biology,
        # wrong strain background. Wza is the correct capsule system target.
        {"target_name": "Wza Capsule Export Channel (E. coli homolog, >70% Kp id.)",
         "target_pdb": "2J58", "compartment": "Outer Membrane", "base_dg": -4.25,
         "chain": "A",
         "organism": "Escherichia coli K-30 - HOMOLOG, >70% identity to Kp Wza",
         "uniprot": "wza",
         "provenance": "Dong et al. 2006, PNAS 103:2066; PDB 2J58. No native Kp Wza structure exists.",
         "active_residues": [169, 173, 175, 176, 211, 227, 240],
         "docking_semantics": "channel_blockade",
         "site_note": "Active_residues are amphipathic helix residues lining the central "
                      "Wza pore lumen (Dong 2006, PNAS). Docking here models compounds "
                      "that physically occlude capsule export — a channel-blockade mechanism, "
                      "not classical enzyme active-site binding."},

        # ── OUTER MEMBRANE PERMEABILITY: OmpK36 porin ────────────────────
        # RESTORATION NOTE: OmpK36 was incorrectly removed in v3.2 under the
        # critique that "ompK36 is inactivated in CRKP". That critique was
        # correct about using OmpK36 as an antibacterial KILL target. But the
        # original objection was about the active-site residues, not the
        # protein itself. OmpK36 is the correct entry here — NOT as a target
        # to inhibit — but as the structural channel whose constriction zone
        # residues define which compounds can physically enter K. pneumoniae.
        #
        # Docking to the L3 constriction zone (residues ~111-117 in OmpK36)
        # models permeability, not inhibition. Compounds docking in the
        # channel passage with low ΔG are MORE permeable (they fit the channel
        # geometry and are NOT sterically excluded), which is opposite to an
        # enzyme inhibitor interpretation. This must be stated explicitly in
        # the Discussion and in any publication.
        #
        # 5O79: Acosta-Gutiérrez et al. 2018, J Mol Biol. K. pneumoniae native.
        # L3 constriction loop residues: Tyr106, Arg132, Phe133, Asn136, Tyr226.
        # The GD insertion (Gly115-Asp116) in CRKP ST258 narrows this zone
        # (Begg et al. 2019, Nat Commun). Wild-type 5O79 models the open channel.
        {"target_name": "OmpK36 Porin — L3 constriction zone (permeability, not inhibition)",
         "target_pdb": "5O79", "compartment": "Outer Membrane", "base_dg": -4.0,
         "chain": "A",
         "organism": "Klebsiella pneumoniae (native OmpK36)",
         "uniprot": "A6T5P9 (ompK36)",
         "provenance": "Acosta-Gutierrez et al. 2018, ACS Infect Dis; Begg et al. 2019, Nat Commun; PDB 5O79",
         "active_residues": [106, 111, 113, 116, 132, 226],
         "docking_semantics": "permeability",
         "site_note": "ΔG here models compound transit through the L3 constriction zone, "
                      "NOT inhibitor binding. More negative = compound fits the pore = MORE permeable. "
                      "Active_residues are the L3 loop lining residues (Begg 2019, Nat Commun)."},

        # ── DNA GYRASE: GyrA subunit (breakage-reunion domain, gate site) ─
        # RESTORATION NOTE: GyrA was removed because 1AB4 is E. coli, not
        # K. pneumoniae native. That organism error was real and 1AB4 should
        # not have been labeled as "K. pneumoniae GyrA". However, removing the
        # protein entirely was wrong — GyrA is a validated, conserved target in
        # ALL K. pneumoniae. No native Kp GyrA-DNA complex exists in PDB
        # (exhaustive search, June 2026). 1AB4 (E. coli, 59-kDa GyrA fragment,
        # >85% aa identity to Kp GyrA) is the legitimate reference used in all
        # Enterobacteriaceae GyrA drug discovery papers. It is explicitly labeled
        # as E. coli homolog here — honest, defensible, used consistently across
        # the literature.
        #
        # CRITICAL SITE NOTE: active_residues MUST point to the GyrA-DNA
        # breakage-reunion GATE interface (Tyr122, Arg121, Lys374, Lys432),
        # NOT the QRDR (residues 67-106). QRDR is mutated in MDR strains;
        # the catalytic gate is conserved regardless of quinolone resistance.
        # This means any compound docked at these residues has a binding
        # argument that does NOT depend on whether the strain is FQ-resistant.
        #
        # 1AB4: E. coli 59-kDa GyrA fragment, PDB 1AB4 (Morais Cabral et al.
        # 1997). The only available high-resolution GyrA breakage-reunion
        # domain structure in Enterobacteriaceae. Stated as E. coli homolog.
        {"target_name": "GyrA — DNA breakage-reunion gate (E. coli homolog, >85% Kp id.)",
         "target_pdb": "1AB4", "compartment": "Cytoplasm", "base_dg": -5.75,
         "chain": "A",
         "organism": "Escherichia coli - HOMOLOG, >85% identity to Kp GyrA",
         "uniprot": "gyrA",
         "provenance": "Morais Cabral et al. 1997, Nature 388:903; PDB 1AB4. No native Kp GyrA-DNA complex exists.",
         "docking_semantics": "enzyme_inhibition",
         "active_residues": [118, 119, 120, 121, 122, 374]},
    ],
}

# ── PDB active-site centroid fetcher ─────────────────────────
_pdb_centroid_cache = {}

def get_active_site_centroid(pdb_id: str, active_residues: list) -> tuple:
    """
    Fetch PDB, extract Cα atoms of active site residues,
    return (cx, cy, cz) centroid in crystallographic frame.
    Falls back to protein geometric centre if residues not found.
    """
    if pdb_id in _pdb_centroid_cache:
        return _pdb_centroid_cache[pdb_id]

    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        req = urllib.request.Request(url, headers={"User-Agent": "BAMSTU/3.0"})
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(req, context=ctx, timeout=12) as r:
            pdb_text = r.read().decode("utf-8", errors="ignore")

        xs, ys, zs = [], [], []
        fallback_xs, fallback_ys, fallback_zs = [], [], []

        for line in pdb_text.split("\n"):
            if not line.startswith("ATOM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resi = int(line[22:26].strip())
                aname = line[12:16].strip()
                fallback_xs.append(x)
                fallback_ys.append(y)
                fallback_zs.append(z)
                if resi in active_residues and aname == "CA":
                    xs.append(x); ys.append(y); zs.append(z)
            except (ValueError, IndexError):
                continue

        if xs:
            cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
        elif fallback_xs:
            # Geometric centre of whole protein
            cx = sum(fallback_xs)/len(fallback_xs)
            cy = sum(fallback_ys)/len(fallback_ys)
            cz = sum(fallback_zs)/len(fallback_zs)
        else:
            # No atoms found at all — signal failure with None
            return None, None, None

        _pdb_centroid_cache[pdb_id] = (cx, cy, cz)
        return cx, cy, cz

    except Exception:
        return None, None, None

# ── XYZ parser (v3.2 — coordinates only; descriptors come from the graph) ──
def parse_xyz(content: str):
    """Read coordinates. Descriptors are NOT computed here any more.

    v3.1 computed MW, HBA and HBD inside this function by counting element
    symbols. An .xyz file has no bonds, so those counts could not be right:
    HBD incremented for every O and N whether or not it carried a hydrogen,
    which is why every complex came back with HBD equal to HBA. Descriptors
    now come from bamstu_chem.compute_descriptors(), which rebuilds the
    molecular graph first.

    Returns (atoms, descriptors, coord_fp). `descriptors` may contain nulls;
    callers must check descriptors["provenance"] before using a value.
    """
    atoms = []
    lines = content.strip().split("\n")
    # Standard XYZ has an atom count on line 0 and a comment on line 1.
    # Bare XYZ starts straight into "elem x y z".
    try:
        int(lines[0].split()[0])
        _start = 2
    except (ValueError, IndexError):
        _start = 0

    bad_lines = []
    for ln_no, line in enumerate(lines[_start:], start=_start + 1):
        parts = line.split()
        if not parts:
            continue
        if len(parts) < 4:
            bad_lines.append(ln_no)
            continue
        try:
            elem = parts[0].strip()
            # Normalise 'AG' / 'ag' -> 'Ag' so the mass table and the metal
            # sets match regardless of how the QM package wrote the file.
            elem = elem[:1].upper() + elem[1:].lower() if len(elem) > 1 else elem.upper()
            atoms.append({"elem": elem,
                          "x": float(parts[1]),
                          "y": float(parts[2]),
                          "z": float(parts[3])})
        except ValueError:
            bad_lines.append(ln_no)

    if bad_lines:
        raise HTTPException(
            422,
            f"Malformed coordinate line(s) at {bad_lines[:10]}"
            + (" ..." if len(bad_lines) > 10 else "")
            + ". v3.1 skipped these silently, so a truncated file docked as a "
              "fragment and reported a score for a molecule that was never "
              "fully read."
        )

    descriptors = chem.compute_descriptors(atoms)
    coord_fp = chem.coordinate_fingerprint(atoms)
    return atoms, descriptors, coord_fp

# ── cLogP (v3.2 — real Wildman-Crippen, or nothing) ──────────────────
def compute_clogp(atoms: list):
    """Return the Crippen logP of the organic scaffold, or None.

    REMOVED IN v3.2: the element-additive table that lived here. It assigned
    one constant per element -- a carbonyl, ether and hydroxyl oxygen all
    scored -0.7779 -- and the docstring credited it to Wildman & Crippen
    (1999). Wildman-Crippen defines roughly 68 atom types by chemical
    environment; one constant per element is not a simplified version of it,
    it is a different quantity with someone else's name on it. Misattributing
    a named method is the specific thing a referee checks.

    The metal values ("scaled from octanol-water partition of metal
    complexes", no source) are gone with it. No additive scheme covers a
    coordinated metal centre; the honest answer is to measure logD.
    """
    d = chem.compute_descriptors(atoms)
    return d["clogp"]

# ── Upgraded docking score engine ────────────────────────────
# ── Per-metal Vina-equivalent coordination parameters ─────────────────────
#
# Parameters sourced from:
#   AutoDock4 atom type tables (Morris et al. J Comput Chem 2009)
#   CCDC CSD survey of metal-ligand bond lengths (2023)
#   Pearson HSAB theory for donor atom preference
#   Shannon ionic radii (Acta Cryst A32:751, 1976)
#
# Format: (coord_E, vdW_r, geom_factor, preferred_donors)
#   coord_E      : coordination energy contribution (kcal/mol)
#                  derived from AutoDock4 well-depth parameters
#   vdW_r        : van der Waals radius (Angstrom) for steric term
#   geom_factor  : geometry-dependent occupancy multiplier
#                  octahedral=1.00 > square-planar=0.92 >
#                  tetrahedral=0.85 > trigonal=0.80 > linear=0.75
#   preferred_donors : Pearson HSAB donor atom preference list
#                      Hard metal → O > N; Soft metal → S > N > O
#
# Oxidation-state differentiation:
#   Higher oxidation state → stronger Lewis acid → greater coord_E
#   Higher OS → smaller ionic radius → smaller vdW_r
#   Higher OS → harder acid → prefers harder donors (O > N > S)
# ─────────────────────────────────────────────────────────────────────────
_METAL_PARAMS: dict = {

    # ── COBALT ──────────────────────────────────────────────
    # Co(II): d7, octahedral or tetrahedral, borderline hard/soft
    # Co(III): d6, strongly prefers octahedral, hard Lewis acid
    "Co":   (1.05, 1.25, 1.00, ["N","O","S"]),  # Co(II) octahedral
    "CoII": (1.05, 1.25, 1.00, ["N","O","S"]),  # Co(II) explicit
    "CoIII":(1.18, 1.19, 1.00, ["N","O"]),       # Co(III) hard acid, smaller radius

    # ── COPPER ──────────────────────────────────────────────
    # Cu(I): d10, linear or tetrahedral, soft Lewis acid
    # Cu(II): d9, square-planar Jahn-Teller, borderline acid
    "Cu":   (0.98, 1.28, 0.92, ["N","S","O"]),  # Cu(II) square-planar (default)
    "CuI":  (0.88, 1.46, 0.80, ["S","N"]),       # Cu(I) tetrahedral/linear, soft
    "CuII": (0.98, 1.28, 0.92, ["N","S","O"]),  # Cu(II) square-planar

    # ── RUTHENIUM ───────────────────────────────────────────
    # Ru(II): d6, octahedral, 18-electron, extremely stable
    # Ru(III): d5, octahedral, kinetically inert, harder acid
    "Ru":    (1.08, 1.46, 1.00, ["N","O"]),      # Ru(II) octahedral (default)
    "RuII":  (1.08, 1.46, 1.00, ["N","O"]),      # Ru(II) e.g. [Ru(bpy)3]2+
    "RuIII": (1.14, 1.40, 1.00, ["N","O"]),      # Ru(III) harder, slightly smaller

    # ── IRON ────────────────────────────────────────────────
    # Fe(II): d6 LS/HS, octahedral, borderline Lewis acid
    # Fe(III): d5, octahedral, hard Lewis acid, prefers O donors
    "Fe":    (0.95, 1.26, 1.00, ["N","O","S"]),  # Fe(II) default
    "FeII":  (0.95, 1.26, 1.00, ["N","O","S"]),  # Fe(II) LS octahedral
    "FeIII": (1.05, 1.20, 1.00, ["O","N"]),       # Fe(III) hard acid, O preferred

    # ── SILVER ──────────────────────────────────────────────
    # Ag(I): d10, linear or tetrahedral, soft Lewis acid
    # Ag(I) antibacterial mechanism: disrupts cell membrane,
    # coordinates with Cys/His/Met in bacterial proteins
    "Ag":   (0.82, 1.72, 0.75, ["N","S"]),       # Ag(I) linear default
    "AgI":  (0.82, 1.72, 0.75, ["S","N"]),        # Ag(I) soft, prefers S

    # ── GOLD ────────────────────────────────────────────────
    # Au(I): d10, linear, very soft, strong trans-influence
    # Au(III): d8, square-planar, harder, anticancer
    "Au":    (1.12, 1.66, 0.75, ["S","N"]),       # Au(I) linear
    "AuI":   (1.12, 1.66, 0.75, ["S","N"]),       # Au(I) linear, very soft
    "AuIII": (1.25, 1.50, 0.92, ["N","O","S"]),  # Au(III) square-planar

    # ── NICKEL ──────────────────────────────────────────────
    # Ni(II): d8, square-planar preferred (strong field ligands)
    # or octahedral (weak field). Borderline Lewis acid.
    "Ni":   (1.02, 1.24, 0.92, ["N","S"]),       # Ni(II) square-planar
    "NiII": (1.02, 1.24, 0.92, ["N","S"]),        # Ni(II) square-planar

    # ── ZINC ────────────────────────────────────────────────
    # Zn(II): d10, tetrahedral, hard Lewis acid
    # Present in many enzyme active sites (metalloenzymes)
    "Zn":   (0.90, 1.20, 0.85, ["N","S","O"]),  # Zn(II) tetrahedral

    # ── MANGANESE ───────────────────────────────────────────
    # Mn(II): d5, octahedral, hard Lewis acid, O donors preferred
    # Mn(III): d4, octahedral, Jahn-Teller active
    "Mn":   (0.88, 1.61, 1.00, ["O","N"]),       # Mn(II) octahedral
    "MnII": (0.88, 1.61, 1.00, ["O","N"]),        # Mn(II)
    "MnIII":(0.95, 1.45, 1.00, ["O","N"]),        # Mn(III) Jahn-Teller

    # ── PLATINUM ────────────────────────────────────────────
    # Pt(II): d8, square-planar, kinetically inert, anticancer
    # Pt(IV): d6, octahedral, prodrug forms
    "Pt":   (1.15, 1.39, 0.92, ["N","S"]),       # Pt(II) square-planar
    "PtII": (1.15, 1.39, 0.92, ["N","S"]),        # Pt(II) cisplatin-type
    "PtIV": (1.22, 1.30, 1.00, ["N","O","S"]),   # Pt(IV) octahedral prodrug

    # ── PALLADIUM ───────────────────────────────────────────
    # Pd(II): d8, square-planar, softer than Pt, similar chemistry
    "Pd":   (1.10, 1.63, 0.92, ["N","S"]),       # Pd(II) square-planar
    "PdII": (1.10, 1.63, 0.92, ["N","S"]),        # Pd(II)

    # ── RHODIUM ─────────────────────────────────────────────
    # Rh(I): d8, square-planar, soft
    # Rh(III): d6, octahedral, kinetically inert
    "Rh":    (1.06, 1.42, 1.00, ["N","O","S"]),  # Rh(III) octahedral default
    "RhI":   (1.00, 1.50, 0.92, ["N","S","P"]),  # Rh(I) square-planar, soft
    "RhIII": (1.06, 1.42, 1.00, ["N","O","S"]),  # Rh(III) octahedral

    # ── OSMIUM ──────────────────────────────────────────────
    # Os(II): d6, octahedral, very inert (slower ligand exchange than Ru)
    # Os(III): d5, octahedral, harder
    "Os":    (1.06, 1.44, 1.00, ["N","O"]),      # Os(II) octahedral
    "OsII":  (1.06, 1.44, 1.00, ["N","O"]),       # Os(II)
    "OsIII": (1.12, 1.38, 1.00, ["N","O"]),       # Os(III)

    # ── IRIDIUM ─────────────────────────────────────────────
    # Ir(I): d8, square-planar, soft
    # Ir(III): d6, octahedral, most inert of platinum group
    "Ir":    (1.10, 1.41, 1.00, ["N","O","S"]),  # Ir(III) default
    "IrI":   (1.04, 1.50, 0.92, ["N","S","P"]),  # Ir(I) square-planar
    "IrIII": (1.10, 1.41, 1.00, ["N","O","S"]),  # Ir(III) octahedral

    # ── VANADIUM ────────────────────────────────────────────
    # V(IV): d1, VO2+ vanadyl, hard acid, antidiabetic properties
    # V(V): d0, hard acid, phosphatase inhibitor
    "V":    (0.82, 1.58, 1.00, ["O","N"]),       # V(IV) vanadyl default
    "VIV":  (0.82, 1.58, 1.00, ["O","N"]),        # V(IV) VO2+
    "VV":   (0.88, 1.42, 0.85, ["O"]),            # V(V) hard, O only

    # ── CHROMIUM ────────────────────────────────────────────
    # Cr(III): d3, octahedral, kinetically inert, hard Lewis acid
    "Cr":    (0.85, 1.52, 1.00, ["O","N"]),      # Cr(III) octahedral
    "CrIII": (0.85, 1.52, 1.00, ["O","N"]),       # Cr(III)

    # ── MOLYBDENUM / TUNGSTEN ────────────────────────────────
    # Mo(IV/VI): enzyme cofactors, oxo-transfer chemistry
    "Mo":   (0.90, 1.54, 1.00, ["O","S","N"]),   # Mo(IV) default
    "W":    (0.92, 1.54, 1.00, ["O","S","N"]),   # W(IV) similar to Mo

    # ── GALLIUM / BISMUTH (antimicrobial metals) ─────────────
    # Ga(III): d10, octahedral, Fe(III) mimic — Trojan horse strategy
    # Bi(III): 6s2 lone pair, used in H. pylori treatment
    "Ga":   (0.88, 1.22, 1.00, ["O","N"]),       # Ga(III) Fe mimic
    "GaIII":(0.88, 1.22, 1.00, ["O","N"]),        # Ga(III) explicit
    "Bi":   (0.75, 1.54, 0.80, ["O","S"]),        # Bi(III) antimicrobial
}

# Target-specific metal coordination affinity
# Based on known active site residue composition in TARGETS database
_TARGET_METAL_AFFINITY: dict = {
    # PDB : {metal: bonus kcal/mol}  based on His/Cys/Asp density at active site
    "4CJN": {"Ag":0.15,"Cu":0.20,"Ni":0.18,"Co":0.10},  # PBP2a — Ser/Lys/Cys active site
    "4URO": {"Co":0.25,"Ni":0.22,"Zn":0.20},  # GyrB — Asn/Asp chelation
    "3TYE": {"Zn":0.18,"Mg":0.15,"Co":0.12},  # DHPS — Mg-dependent
    "4ALI": {"Co":0.20,"Ni":0.18,"Fe":0.15},  # FabI — NAD+ binding
    "1T2W": {"Cu":0.15,"Zn":0.12,"Ag":0.10},  # SrtA — Cys active site
    "2WTP": {"Zn":0.10,"Co":0.08},             # ClfA — surface protein
    "3T05": {"Mn":0.25,"Mg":0.20,"Co":0.18},  # PykA — Mg/K dependent
    "3L7H": {"Mn":0.20,"Zn":0.15},            # LtaS — Mn-dependent
    "3RXX": {"Zn":0.15,"Co":0.12,"Ni":0.10},  # KPC-2 beta-lactamase
    "4EXS": {"Zn":0.30,"Ag":0.28,"Cu":0.25,"Co":0.20,"Ni":0.18},  # NDM-1 di-Zn (His120/122/189, Asp124, Cys208/His250)
                                                # Ag(I) and Cu(II) coordinate Cys208 and His ligands directly
                                                # — strongest metal affinity in the panel. Zhang et al. 2012.
    "6PJ3": {"Mn":0.28,"Co":0.20,"Zn":0.15},  # LpxH — native Mn(II)-dependent pyrophosphohydrolase
                                                # (PDB 6PJ3 ligand record: MANGANESE (II) ION; Cho et al. 2020 PNAS)
    "9KGT": {"Mg":0.20,"Co":0.15},             # Topoisomerase IV ParC — Mg-dependent (K. pneumoniae)
    "2J58": {"Cu":0.08,"Ag":0.08},              # Wza capsule export channel (E. coli homolog)
    "5V3D": {"Mn":0.25,"Mg":0.15,"Co":0.18},  # FosA — Mn(II)-dependent glutathione transferase
                                                # (K. pneumoniae fosfomycin resistance enzyme)
    "3HBR": {"Zn":0.10,"Co":0.08,"Cu":0.12},  # OXA-48: Ser70/Lys73/Arg250 active site;
                                                # no Cys — metal via Ser-OH/Arg coordination
    "7CZ9": {"Cu":0.10,"Ag":0.08,"Co":0.08},  # OqxB: hydrophobic substrate-binding pocket;
                                                # aromatic stacking rather than HSAB coordination
    "5O79": {"Cu":0.06,"Ag":0.06},             # OmpK36: nonspecific channel — very low affinity;
                                                # docking here models permeability, not binding
    "1AB4": {"Mg":0.20,"Co":0.15,"Mn":0.12},  # GyrA: Mg²⁺ catalytic gate (Tyr122 site)
                                                # metal coordination at phosphodiester cleavage site
}

# ── Deterministic target tie-breaker (NOT noise — see compute_docking_score) ──
# Explicit, permanent, unique micro-offset per target (±0.21 kcal/mol max,
# far below the AD4 metal correction's 0.5-1.5 kcal/mol scale). This alone
# reduces — but does NOT fully eliminate — coincidental score collisions
# in PRE-SCREENING mode (several targets share identical literature-
# calibrated base_dg values: 4URO/3RXX both -6.25, 1T2W/3T05 both -5.0,
# 2WTP/3L7H/7W1D all three -4.5, so even unique per-target offsets can
# combine with a metal-specific bonus to land back on the same rounded
# total for a DIFFERENT pair of targets). The mathematically complete
# fix is disambiguate_hits() below, which operates on the full sweep
# after all targets are scored and guarantees zero collisions (verified
# 0/2000 trials) regardless of how base_dg values happen to overlap.
# This table remains in place because it still does the useful first-pass
# work of separating same-base_dg targets in the common case, reducing
# how often disambiguate_hits() needs to nudge a value at all.
# _TARGET_TIEBREAK was removed in v3.4. It was a per-target cosmetic
# offset of up to 0.33 kcal/mol added so that two targets could not
# print the same rounded score. Because it entered the score before the
# ranking, it reordered targets whose real separation was smaller than
# the offset -- for the free ligands, four to six of ten. A display
# collision is a display problem and is fixed by printing more digits.


def rank_hits(hits: list) -> list:
    """Sort hits and report ties as ties.

    REMOVED IN v3.2: the +0.01 kcal/mol nudge that guaranteed no two hits
    could ever display the same value. The stated justification was that
    identical displayed scores were a rounding artefact. They were not an
    artefact; they were the correct answer. The empirical error of an
    empirical docking scoring function is roughly 1.5-2.0 kcal/mol, so two
    targets within 0.2 kcal/mol are indistinguishable and a table that shows
    them as distinct is claiming a resolution the method does not have.

    Ties are now surfaced. Where a secondary key is genuinely needed for
    ordering, ligand efficiency is used, because it is a defined quantity
    with a physical meaning rather than a constant chosen to break a tie.
    """
    if not hits:
        return hits

    for h in hits:
        raw = h.get("vina_affinity_kcal_per_mol")
        heavy = h.get("heavy_atoms") or 0
        h["ligand_efficiency"] = (round(raw / heavy, 4)
                                  if raw is not None and heavy else None)

    ordered = sorted(
        hits,
        key=lambda h: (
            h.get("bamstu_composite_score")
            if h.get("bamstu_composite_score") is not None else 0.0,
            -(h.get("ligand_efficiency") or 0.0),
        ),
    )

    # Flag statistical ties rather than hiding them.
    TIE_WINDOW = 0.5   # kcal/mol; conservative relative to the ~1.5-2.0
                       # kcal/mol error of empirical scoring functions
    for idx, h in enumerate(ordered):
        h["rank"] = idx + 1
        tied = []
        for other in ordered:
            if other is h:
                continue
            a = h.get("bamstu_composite_score")
            b = other.get("bamstu_composite_score")
            if a is not None and b is not None and abs(a - b) <= TIE_WINDOW:
                tied.append(other["target_pdb"])
        h["statistically_tied_with"] = tied
        h["rank_is_resolvable"] = not tied

    if any(h["statistically_tied_with"] for h in ordered):
        for h in ordered:
            h["rank_note"] = (
                f"Targets within {TIE_WINDOW} kcal/mol of each other are not "
                "separable by this method. Where 'statistically_tied_with' is "
                "non-empty, the rank order is presentational and must not be "
                "reported as a target preference."
            )
    return ordered


# Backwards-compatible alias so older callers do not break. The nudging
# behaviour is gone; only the sort remains.
disambiguate_hits = rank_hits


def compute_composite_score(atoms: list, vina_affinity: float,
                            pathogen: str,
                            coord_metal_contacts: int = 0,
                            target_pdb: str = "",
                            descriptors: dict | None = None) -> dict:
    """BAMSTU composite screening score — arbitrary units, not a free energy.

    WHAT THIS IS
    ------------
    A triage number for ordering hypotheses within one compound series. It
    starts from the AutoDock Vina affinity and applies a metal-coordination
    term that Vina cannot express, because Vina types every transition metal
    as Zn.

    WHAT THIS IS NOT
    ----------------
    A binding free energy. The coefficients below have never been fitted
    against a measured affinity dataset. There is no PDBbind benchmark, no
    calibration set, no reported correlation with experiment. Until there
    is, the units are arbitrary and the value carries no thermodynamic
    meaning. It is returned under a name that says so, and dg_to_ki() will
    refuse to convert it.

    REMOVED FROM v3.1
    -----------------
      * -0.08 kcal/mol per polar atom, -0.020 per heavy atom (capped 0.90).
        These rewarded molecular bulk whether or not the atom contacted the
        receptor, on top of Vina's own documented bias toward large ligands.
        Between AgL4 and AgSD this bonus alone was worth up to 0.83 kcal/mol
        and it is where the ranking inversion came from.
      * cLogP membrane and lipophilicity-deviation terms. Permeability is
        now reported separately and never summed into a score.
      * The per-target cosmetic offset (_TARGET_TIEBREAK).
      * The rotatable-bond proxy (C_count//6 for metal complexes, //4 for
        organics) that fed the torsional entropy term. Rotatable bonds now
        come from the perceived graph, or the term is not applied at all.

    RETAINED, WITH ITS UNCERTAINTY STATED
    -------------------------------------
    The metal coordination term. It is the one correction with a defensible
    rationale, and it is now gated: it applies only in proportion to metal
    contacts actually detected in the pose by bamstu_analysis. A metal that
    coordinates nothing gets no bonus. In v3.1 every metal centre collected
    the full term whether or not it was anywhere near a donor atom.
    """
    descriptors = descriptors or {}
    terms = []
    score = float(vina_affinity)

    _BASE_METALS = {e.rstrip("IVX0123456789") for e in _METAL_PARAMS}
    metal_atoms = [a for a in atoms
                   if a["elem"] in _METAL_PARAMS or a["elem"] in _BASE_METALS]
    donor_elems = {a["elem"] for a in atoms if a["elem"] in ("O", "N", "S", "F")}

    # ── Metal coordination, gated on observed contacts ──────────────────
    coord_total = 0.0
    if metal_atoms:
        if coord_metal_contacts <= 0:
            terms.append({
                "term": "metal_coordination",
                "value": 0.0,
                "basis": (
                    f"{len(metal_atoms)} metal centre(s) present but zero "
                    "metal-donor contacts detected in the pose within "
                    "CSD-derived coordination windows. No coordination bonus "
                    "applied. v3.1 applied the full term regardless, which is "
                    "why outputs claimed metal coordination while the contact "
                    "table listed none."
                ),
            })
        else:
            # Scale by the fraction of metal centres that actually coordinate.
            occupancy = min(1.0, coord_metal_contacts / max(len(metal_atoms), 1))
            for ma in metal_atoms:
                params = _METAL_PARAMS.get(ma["elem"])
                if not params:
                    continue
                coord_E, _vdw_r, geom_f, pref_donors = params
                match = (sum(1 for d in pref_donors if d in donor_elems)
                         / max(len(pref_donors), 1))
                coord_total += coord_E * geom_f * (0.7 + 0.3 * match) * occupancy
            score -= coord_total
            terms.append({
                "term": "metal_coordination",
                "value": round(-coord_total, 3),
                "basis": (
                    f"{coord_metal_contacts} observed metal-donor contact(s), "
                    f"occupancy {occupancy:.2f}. Coordination energies from "
                    "Rulisek & Vondrasek 1998 and Dudev & Lim 2003; HSAB donor "
                    "preference from Pearson 1963. UNCALIBRATED against any "
                    "measured affinity set."
                ),
            })

    # ── Torsional entropy, only with a real rotatable-bond count ────────
    n_rot = descriptors.get("n_rotatable_bonds")
    if n_rot is not None:
        penalty = n_rot * 0.048   # Vina coefficient, Trott & Olson 2010
        score += penalty
        terms.append({
            "term": "torsional_entropy",
            "value": round(penalty, 3),
            "basis": (f"{n_rot} rotatable bonds from the perceived graph "
                      "x 0.048 kcal/mol (Vina coefficient). Note this "
                      "double-counts slightly: Vina already applies its own "
                      "torsional term internally."),
        })
    else:
        terms.append({
            "term": "torsional_entropy",
            "value": 0.0,
            "basis": ("Not applied. Rotatable-bond count unavailable, and the "
                      "v3.1 proxy (carbon count divided by 6) was not a "
                      "torsion count."),
        })

    # ── Three-way reporting of the metal physics term ───────────────
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
        ),
        "terms": terms,
        "total_adjustment": round(score - float(vina_affinity), 3),
        "calibration_status": "UNCALIBRATED — no fit against measured affinities",
        "reporting_rule": (
            "Report vina_affinity_kcal_per_mol as 'AutoDock Vina 1.2 score'. "
            "Report this composite, if at all, as 'BAMSTU composite screening "
            "score, arbitrary units' with the term table above. Do not attach "
            "thermodynamic language to either."
        ),
    }


def compute_docking_score(atoms: list, base_dg: float, pathogen: str,
                          clogp: float = None, exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
                          coord_fp: int = 0, target_pdb: str = "",
                          coord_metal_contacts: int = 0,
                          descriptors: dict | None = None) -> float:
    """Legacy signature, kept so existing callers do not break.

    Returns only the composite number. New code should call
    compute_composite_score() and keep the term breakdown, because the
    breakdown is what makes the number defensible.
    """
    return compute_composite_score(
        atoms, base_dg, pathogen,
        coord_metal_contacts=coord_metal_contacts,
        target_pdb=target_pdb,
        descriptors=descriptors,
    )["bamstu_composite_score"]


def dg_to_ki(dg: float) -> None:
    """Refuses, on purpose, and names the function that does not.

    v3.1 applied Ki = exp(dG/RT) to the CORRECTED score -- a Vina number
    with uncalibrated offsets and a permeability penalty inside it -- and
    printed the result in nanomolar to two decimals. Two of four compounds
    landed on the identical 4.44 nM.

    An estimated Ki from the RAW Vina affinity is a legitimate Layer 1
    quantity and v3.4 restores it, in bamstu_layers.estimated_docking_ki().
    It is a separate function because the input matters more than the
    arithmetic: the same transform is defensible on one number and
    meaningless on another, and keeping them in one function is how the
    meaningless one got printed for a year.
    """
    return None


# ── KEY FIX: Pose generator with active site translation ──────
# AutoDock atom type -> element. Vina writes its own types in the last
# column of a PDBQT, and several of them are not element symbols: A is
# aromatic carbon, NA/OA/SA are hydrogen-bond acceptors, HD is a polar
# hydrogen. Reading those two characters as an element turns aromatic
# carbon into argon and acceptor nitrogen into sodium, which would put
# spurious metals into the contact analysis.
_AD_TYPE_TO_ELEMENT = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O",
    "OS": "O", "O": "O", "SA": "S", "S": "S", "HD": "H", "H": "H",
    "P": "P", "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br",
    "I": "I", "Mg": "Mg", "MG": "Mg", "Ca": "Ca", "CA": "Ca",
    "Mn": "Mn", "MN": "Mn", "Fe": "Fe", "FE": "Fe", "Zn": "Zn",
    "ZN": "Zn", "Si": "Si", "B": "B",
}


def extract_vina_pose(out_pdbqt: Path, model: int = 1,
                      metal_map: dict = None) -> str | None:
    """Return one binding mode from a Vina output PDBQT as PDB text.

    This is the pose Vina actually produced: translated, rotated and, where
    torsions were allowed, re-dialled. It is the only ligand geometry in the
    system with any claim to physical meaning, and it is what the contact
    analysis must read.

    `metal_map` reverses the Zn proxy for display and analysis. Vina was
    given a Zn atom where a silver sits, so the output PDBQT says Zn. The
    element is restored here so the contact analysis applies Ag-S windows
    rather than Zn-S windows. The SAMPLING was still done with zinc
    parameters -- restoring the label does not restore the physics, and the
    metal_proxy block on the result says so.
    """
    try:
        text = out_pdbqt.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None
    if not text.strip():
        return None

    blocks, current, in_model = [], [], False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            in_model, current = True, []
        elif line.startswith("ENDMDL"):
            if current:
                blocks.append(current)
            in_model, current = False, []
        elif in_model:
            current.append(line)
        elif line.startswith(("ATOM", "HETATM")):
            current.append(line)
    if current:
        blocks.append(current)
    if not blocks or model < 1 or model > len(blocks):
        return None

    metal_map = metal_map or {}
    out = [f"REMARK BAMSTU v3.2.2 — AutoDock Vina docked pose, "
           f"MODEL {model} of {len(blocks)}",
           "REMARK This is Vina's output geometry, not a placed copy of the "
           "input."]
    if metal_map:
        out.append("REMARK Metal element labels restored for analysis: "
                   + ", ".join(f"Zn->{v}" for v in sorted(set(metal_map.values())))
                   + " (sampling used Zn parameters regardless)")

    kept = 0
    for line in blocks[model - 1]:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ad_type = line[77:79].strip() if len(line) > 77 else ""
        elem = _AD_TYPE_TO_ELEMENT.get(ad_type, ad_type[:2] or "C")
        if elem == "Zn" and metal_map:
            # Only one metal species per ligand is handled; a mixed-metal
            # complex is left as Zn and flagged, rather than guessed at.
            species = sorted(set(metal_map.values()))
            if len(species) == 1:
                elem = species[0]
        base = line[:66].ljust(66)
        out.append(f"{base}          {elem:>2s}")
        kept += 1

    if kept == 0:
        return None
    out.append("END")
    return "\n".join(out) + "\n"


def generate_pose(atoms: list, pdb_code: str,
                  active_residues: list,
                  coord_fp: int = 0) -> str:
    """
    VIEWER PLACEMENT ONLY — this is NOT a docked pose.

    Translates the input geometry so its centroid sits on the pocket
    centroid, with a small offset. There is no rotation, no torsional
    sampling and no energy evaluation, so the result has no physical
    meaning and must not be used for contact analysis, interaction tables,
    or anything that feeds a score.

    Use extract_vina_pose() for the real pose. This function exists so the
    3D viewer has something to show when Vina has not run, and results
    built on it are labelled SYNTHETIC_PLACEMENT_NOT_DOCKED.
    """
    # 1. Get protein active site centroid (crystallographic frame)
    _centroid = get_active_site_centroid(pdb_code, active_residues)
    cx = _centroid[0] if _centroid[0] is not None else 0.0
    cy = _centroid[1] if _centroid[1] is not None else 0.0
    cz = _centroid[2] if _centroid[2] is not None else 0.0

    # 2. Compute ligand centroid (local XYZ frame)
    if atoms:
        lx = sum(a["x"] for a in atoms) / len(atoms)
        ly = sum(a["y"] for a in atoms) / len(atoms)
        lz = sum(a["z"] for a in atoms) / len(atoms)
    else:
        lx = ly = lz = 0.0

    # 3. Translation vector: move ligand centroid → active site centroid
    dx, dy, dz = cx - lx, cy - ly, cz - lz

    # 4. Viewer offset.
    #
    # v3.1 called this jitter that "simulates pose sampling". It does not
    # simulate sampling -- it is a rigid translation with no rotation and no
    # torsional change, so it cannot represent a docked pose and must never
    # be analysed as one. See extract_vina_pose() for the real pose.
    #
    # The seeding was also not reproducible: hash() of a str is randomised
    # per process unless PYTHONHASHSEED is set, so the same compound and
    # target gave a different offset after every server restart. Replaced
    # with a stable checksum. A viewer offset does not need to be
    # reproducible, but nothing in a scientific pipeline should be silently
    # non-deterministic -- that is how a difference between two runs stops
    # being diagnosable.
    _stable = sum((i + 1) * ord(c) for i, c in enumerate(pdb_code))
    random.seed((coord_fp + _stable + len(atoms) * 31) % 2147483647)
    jx = random.uniform(-2.5, 2.5)
    jy = random.uniform(-2.5, 2.5)
    jz = random.uniform(-2.5, 2.5)

    lines = [f"REMARK BAMSTU Digital Pathogen Pro v3.1 | Docked vs {pdb_code}",
             f"REMARK ACTIVE SITE CENTROID {cx:.3f} {cy:.3f} {cz:.3f}",
             f"REMARK LIGAND TRANSLATION   {dx:.3f} {dy:.3f} {dz:.3f}"]

    for i, a in enumerate(atoms):     # no truncation
        nx = a["x"] + dx + jx
        ny = a["y"] + dy + jy
        nz = a["z"] + dz + jz
        elem = a["elem"]
        # Correct PDB column alignment for 3D viewer and distance calculation
        if len(elem) == 1:
            aname = f" {elem}  "
        else:
            aname = f"{elem:<4s}"
        lines.append(
            f"HETATM{i+1:5d} {aname} LIG A   1    "
            f"{nx:8.3f}{ny:8.3f}{nz:8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    return "\n".join(lines)

# ── Physicochemical / permeability assessment (v3.2) ─────────────────
def compute_adme(descriptors: dict, pathogen: str) -> dict:
    """Envelope-permeability assessment, kept strictly out of the score.

    Two changes from v3.1:

    1. The Lipinski and Veber verdicts are withheld unless the underlying
       descriptors came from a molecular graph. v3.1 emitted PASS/FAIL from a
       wrong HBD count, so the verdict was wrong every time it was printed.

    2. Nothing computed here is added to a binding score. In v3.1 the cLogP
       permeability penalty went straight into the reported delta-G. Binding
       affinity and bilayer transit are different physical quantities;
       summing them destroys the meaning of both. They are now reported side
       by side and never combined.
    """
    verdicts = chem.lipinski_veber(descriptors)
    mw = descriptors.get("mw")
    hba = descriptors.get("hba")
    clogp = descriptors.get("clogp")
    graph_ok = descriptors.get("graph_ok", False)

    reasons = []
    if "Klebsiella" in pathogen or "coli" in pathogen.lower():
        barrier = ("Gram-negative double membrane: LPS outer membrane, "
                   "periplasm and peptidoglycan, inner membrane")
        if mw is not None and mw > 600:
            reasons.append(f"MW {mw} Da is above the general porin exclusion "
                           "limit (~600 Da) for OmpK36/OmpF-type channels")
        if hba is not None and hba > 8:
            reasons.append(f"HBA {hba}: high polar-atom count raises the "
                           "desolvation cost of inner-membrane crossing")
        if clogp is not None and clogp < 0:
            reasons.append(f"cLogP {clogp}: too hydrophilic for passive "
                           "inner-membrane permeation")
        entry_note = (
            "For Gram-negatives the eNTRy rules (Richter et al. 2017, Nature "
            "545:299) predict accumulation better than Lipinski: ionisable "
            "nitrogen, low globularity, at least three rotatable bonds. "
            "Lipinski was derived from an oral human-absorption set and was "
            "never intended to describe bacterial envelope entry."
        )
    else:
        barrier = "Gram-positive thick peptidoglycan wall and lipid bilayer"
        if mw is not None and mw > 900:
            reasons.append(f"MW {mw} Da is large even for a Gram-positive "
                           "envelope")
        if hba is not None and hba > 10:
            reasons.append(f"HBA {hba} may limit passive diffusion")
        entry_note = (
            "Gram-positive entry is far less size-restricted than "
            "Gram-negative; the peptidoglycan mesh admits molecules well "
            "above the porin limit. Efflux (NorA, MepA) is usually the "
            "governing factor, and docking says nothing about it."
        )

    if not graph_ok:
        permeability = "NOT_ASSESSED"
        perm_basis = (
            "Withheld: descriptors unavailable "
            f"({descriptors.get('provenance')}). A permeability call built on "
            "element counts is not a permeability call."
        )
    elif reasons:
        permeability = "LIKELY_LIMITED"
        perm_basis = "; ".join(reasons)
    else:
        permeability = "NO_OBVIOUS_BARRIER_FLAG"
        perm_basis = ("No descriptor exceeds its threshold. This is the "
                      "absence of a red flag, not evidence of uptake.")

    return {
        "mw": mw,
        "hba": hba,
        "hbd": descriptors.get("hbd"),
        "clogp": clogp,
        "tpsa": descriptors.get("tpsa"),
        "n_rotatable_bonds": descriptors.get("n_rotatable_bonds"),
        "clogp_method": descriptors.get("clogp_method"),
        "descriptor_provenance": descriptors.get("provenance"),
        "descriptor_warnings": descriptors.get("warnings", []),
        "barrier": barrier,
        "permeability": permeability,
        "permeability_basis": perm_basis,
        "entry_rules_note": entry_note,
        "ro5": verdicts["ro5"],
        "ro5_violations": verdicts.get("ro5_violations"),
        "veber": verdicts["veber"],
        "verdict_reason": verdicts.get("reason"),
        "scoring_note": (
            "None of these values enters the binding score. v3.1 added a "
            "cLogP permeability penalty directly into the reported delta-G, "
            "which made the published number a composite of affinity and "
            "transit rather than a binding energy."
        ),
    }

# ── Mechanistic claim filter ─────────────────────────────────────────
# Patterns that a docking score cannot support. Each entry is
# (regex, what would be needed to earn the claim).
_FORBIDDEN_CLAIMS = [
    (r"\bcompetitive(ly)?\s+(inhibit|inhibition|inhibitor|binds)",
     "inhibition class requires enzyme kinetics (Km/Vmax at varying substrate)"),
    (r"\bnon-?competitive|\buncompetitive",
     "inhibition class requires enzyme kinetics"),
    (r"\bspares?\s+mammalian|\bmammalian\s+cells?\s+(are\s+)?spared",
     "selectivity over host requires comparative docking against human "
     "orthologues plus a cytotoxicity assay"),
    (r"\bselectiv(e|ity)\s+for\s+bacterial\b",
     "host-vs-pathogen selectivity requires a comparative calculation"),
    (r"\bnon-?toxic\b|\bsafe\s+for\s+human",
     "toxicity requires an assay"),
    (r"\bIC50\b|\bKi\s*=|\bKd\s*=",
     "a potency constant requires a measurement"),
    (r"\bwill\s+(inhibit|kill|eradicate)|\bdemonstrates?\s+that\b",
     "docking generates hypotheses; it does not demonstrate"),
    (r"\bconfirms?\b|\bproves?\b|\bestablishes?\b",
     "a score confirms nothing"),
    (r"\bbactericidal\b|\bbacteriostatic\b",
     "cidal/static distinction requires MBC:MIC and time-kill data"),
]


def _filter_mechanistic_claims(text: str, metal_contacts: int) -> dict:
    """Scan generated prose for claims the pipeline cannot support."""
    import re
    violations = []
    for pattern, requirement in _FORBIDDEN_CLAIMS:
        for m in re.finditer(pattern, text, flags=re.IGNORECASE):
            snippet = text[max(0, m.start() - 60):m.end() + 60].replace("\n", " ")
            violations.append({
                "matched": m.group(0),
                "requirement": requirement,
                "context": "..." + snippet.strip() + "...",
            })

    # The specific contradiction that appeared in every v3.1 output package.
    if metal_contacts == 0 and re.search(
            r"metal\s+coordinat|coordinat\w*\s+(to|with)\s+"
            r"(the\s+)?(active[- ]site|catalytic|Cys|His|Ser|Asp|Glu)",
            text, flags=re.IGNORECASE):
        violations.append({
            "matched": "metal coordination assertion",
            "requirement": (
                "The contact analysis detected zero metal-donor contacts "
                "within CSD coordination windows for this pose. Text claiming "
                "coordination directly contradicts the pipeline's own table."
            ),
            "context": "whole-text check",
        })

    return {
        "clean": not violations,
        "violations": violations,
        "n_violations": len(violations),
    }


# ── AI MoA — hypothesis generation only, claims filtered ─────────────
def get_ai_moa(cinfo: dict, top_hit: dict, pathogen: str,
               atoms: list = None, interaction_summary: dict = None) -> dict:
    """Generate a hypothesis paragraph, then check it against what was measured.

    Output is labelled a hypothesis, not a mechanism. If the filter finds an
    unsupportable claim the generated text is withheld and the violations are
    returned instead, so a quarantined paragraph cannot be pasted into a
    manuscript by accident.
    """
    has_metals = bool(atoms and any(
        a["elem"] in {"Co", "Cu", "Zn", "Ni", "Fe", "Mn", "Pt", "Pd", "Ru",
                      "Rh", "Ag", "Au", "Os", "Ir", "V", "Cr", "Mo", "W",
                      "Ga", "Bi"}
        for a in (atoms or [])
    ))
    interaction_summary = interaction_summary or {}
    metal_contacts = interaction_summary.get("n_metal_coordination", 0)
    n_hbonds = interaction_summary.get("n_hbonds", 0)
    contact_residues = interaction_summary.get("contact_residues", [])

    papers = []
    try:
        papers = fetch_literature(
            target_name=top_hit["target_name"], pathogen=pathogen,
            has_metals=has_metals, max_papers=6)
    except Exception:
        pass

    lit_context = build_literature_context(papers)
    citations = format_citations(papers)
    landmark_refs = [
        f"RCSB PDB: {top_hit['target_pdb']} — "
        f"https://www.rcsb.org/structure/{top_hit['target_pdb']}",
        "Trott O. & Olson A.J. (2010). AutoDock Vina. J Comput Chem 31:455. "
        "DOI:10.1002/jcc.21334",
        "Berman H.M. et al. (2000). The Protein Data Bank. "
        "Nucleic Acids Res 28:235. DOI:10.1093/nar/28.1.235",
    ]
    all_refs = citations + landmark_refs

    # The measured facts. Anything not on this list is not available to the
    # model, and the model is told so.
    measured = (
        f"MEASURED BY THIS PIPELINE (these are the only facts available):\n"
        f"  Target: {top_hit['target_name']} (PDB {top_hit['target_pdb']})\n"
        f"  AutoDock Vina score: "
        f"{top_hit.get('vina_affinity_kcal_per_mol', 'not available')} kcal/mol "
        f"(scoring-function error approx +/-1.5 kcal/mol)\n"
        f"  Metal-donor contacts detected in the pose: {metal_contacts}\n"
        f"  Putative hydrogen bonds (heavy-atom distance only, no angle "
        f"check): {n_hbonds}\n"
        f"  Contact residues: "
        f"{', '.join(contact_residues[:12]) if contact_residues else 'none recorded'}\n"
        f"  Molecular weight: {cinfo.get('mw')} Da\n"
        f"  cLogP: {cinfo.get('clogp')} "
        f"[{cinfo.get('clogp_method', 'method not recorded')}]\n"
        f"  Descriptor provenance: {cinfo.get('descriptor_provenance')}\n"
        f"  Metal typed as Zn during sampling: {'yes' if has_metals else 'n/a'}\n"
    )

    fallback_inference = (
        f"Hypothesis, not a mechanism. The compound (MW {cinfo.get('mw')} Da) "
        f"scored {top_hit.get('vina_affinity_kcal_per_mol')} kcal/mol against "
        f"{top_hit['target_name']}, which ranks it first in this panel. The "
        f"pose analysis records {metal_contacts} metal-donor contact(s) and "
        f"{n_hbonds} putative hydrogen bond(s). Because the scoring function "
        "carries roughly +/-1.5 kcal/mol of error and the metal centre was "
        "sampled with zinc parameters, this result identifies a target worth "
        "testing rather than a target that has been identified. Enzyme "
        "kinetics would be needed to say anything about inhibition mode, and "
        "comparative work against human orthologues would be needed to say "
        "anything about selectivity."
    )

    if not groq_client:
        return {"inference": fallback_inference, "references": all_refs,
                "literature_count": len(papers), "claim_filter": {"clean": True},
                "generated_by": "template (no LLM configured)"}

    lit_block = (f"\n\nRETRIEVED LITERATURE:\n{lit_context}" if lit_context
                 else "\n\n[No literature retrieved. Do not cite anything.]")

    prompt = (
        "You are drafting a cautious in-silico hypothesis paragraph for a "
        "medicinal-chemistry manuscript. A referee will read it.\n\n"
        + measured + lit_block + "\n\n"
        "HARD CONSTRAINTS. Violating any of these makes the paragraph "
        "unusable:\n"
        "1. Never state an inhibition class. You cannot know whether binding "
        "is competitive, non-competitive or allosteric without enzyme "
        "kinetics, and none were run.\n"
        "2. Never claim selectivity over mammalian or host cells. No "
        "comparative calculation against human proteins was performed.\n"
        "3. If the metal-donor contact count above is zero, do not describe "
        "metal coordination to any residue. Say that no coordination was "
        "detected in the pose.\n"
        "4. Use hypothesis language throughout: 'is consistent with', "
        "'suggests', 'would be tested by'. Never 'demonstrates', 'confirms', "
        "'proves', 'establishes'.\n"
        "5. Cite only the numbered papers above, and only for statements "
        "those papers actually support. If unsure, cite nothing. A reference "
        "attached to a claim it does not support is a fabricated citation.\n"
        "6. State the Zn-proxy limitation once if a metal is present.\n"
        "7. Do not invent residue names, distances, geometries or constants. "
        "Use only the numbers given above.\n\n"
        "Write four short paragraphs, under 200 words total, no markdown:\n"
        "(a) What was computed and what its error bar is.\n"
        "(b) What the contacts observed are consistent with.\n"
        "(c) What the metal centre does or does not appear to do, given the "
        "contact count.\n"
        "(d) The single experiment that would test this hypothesis."
    )

    try:
        _models = ["llama-3.3-70b-versatile", "llama-3.1-8b-instant", "gemma2-9b-it"]
        r = None
        for _m in _models:
            try:
                r = groq_client.chat.completions.create(
                    model=_m, messages=[{"role": "user", "content": prompt}],
                    max_tokens=450, temperature=0.1)
                break
            except Exception as _me:
                if any(x in str(_me).lower()
                       for x in ["decommissioned", "not found", "deprecated"]):
                    continue
                raise
        if r is None:
            raise Exception("All Groq models unavailable")

        text = r.choices[0].message.content.strip()
        check = _filter_mechanistic_claims(text, metal_contacts)

        if not check["clean"]:
            return {
                "inference": None,
                "quarantined_text": text,
                "status": "QUARANTINED",
                "claim_filter": check,
                "references": all_refs,
                "literature_count": len(papers),
                "message": (
                    f"The generated paragraph made {check['n_violations']} "
                    "claim(s) this pipeline cannot support and has been "
                    "withheld. The text is kept under 'quarantined_text' for "
                    "inspection only. Do not paste it into a manuscript. Use "
                    "the template hypothesis below instead."
                ),
                "safe_alternative": fallback_inference,
            }

        return {
            "inference": text,
            "status": "PASSED_CLAIM_FILTER",
            "claim_filter": check,
            "references": all_refs,
            "literature_count": len(papers),
            "sources_queried": ["PubMed", "Semantic Scholar", "Europe PMC",
                                "OpenAlex"],
            "reader_warning": (
                "LLM-generated hypothesis text. Every citation must be opened "
                "and checked against the sentence it is attached to before "
                "submission. Automated citation attachment has produced "
                "mismatches in this pipeline before."
            ),
        }
    except Exception as e:
        return {"inference": fallback_inference, "references": all_refs,
                "literature_count": len(papers),
                "claim_filter": {"clean": True},
                "generated_by": f"template (LLM error: {e})"}


# LEGACY SIGNATURE ALIAS — kept for compatibility
def get_ai_moa_legacy(cinfo, top_hit, pathogen):
    """Legacy alias — dead fallback code removed."""
    return get_ai_moa(cinfo, top_hit, pathogen)


# ════════════════════════════════════════════════════════════
#  API ENDPOINTS
# ════════════════════════════════════════════════════════════

@app.get("/", response_class=HTMLResponse)
def serve_ui():
    """Serve frontend index.html."""
    paths = [
        Path(__file__).parent.parent / "frontend" / "index.html",
        Path(__file__).parent / "index.html",
        Path("../frontend/index.html"),
        Path("index.html"),
    ]
    for p in paths:
        if p.exists():
            return HTMLResponse(content=p.read_text(encoding="utf-8"))
    return HTMLResponse(
        content="<h2 style='font-family:sans-serif;padding:40px;color:#0072b8'>"
                "BAMSTU backend running. Place index.html in the frontend folder, "
                "then visit <a href='http://127.0.0.1:8080'>http://127.0.0.1:8080</a></h2>"
    )


@app.get("/health")
def health():
    return {
        "status":    "BAMSTU Digital Pathogen Pro v3.4.15 running",
        "version":   "3.4.15", "year": 2026,
        # Defaults published so the client cannot silently disagree with the
        # server about what a standard run is. Before v3.4.1 the dropdown
        # said 16 while the server used its own value, and the dropdown was
        # not being sent at all.
        "docking": {
            "default_exhaustiveness": DEFAULT_EXHAUSTIVENESS,
            "channel_exhaustiveness": CHANNEL_EXHAUSTIVENESS,
            "metal_complex_exhaustiveness": METAL_COMPLEX_EXHAUSTIVENESS,
            "ensemble_seeds": list(ENSEMBLE_SEEDS),
            "ensemble_sd_limit_kcal": ENSEMBLE_SD_LIMIT_KCAL,
            "scoring_function_error_kcal": layers.VINA_SCORING_ERROR_KCAL,
            "channel_targets": sorted(CHANNEL_TARGETS),
        },
        "api_key_source": GROQ_KEY_SOURCE,
        "developer": "Dr. James P.T. Matshwele, BAMSTU Labs, University of Botswana",
        "ai":        "Groq connected" if groq_client else "AI disabled",
        "vina":      get_vina_status(),
        "chemistry": chem.status(),
        "engines":   ["AutoDock Vina 1.2", "BAMSTU physics-informed coordination layer v3.4"],
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


@app.get("/api/download_receptors")
async def download_receptors(pathogen: str = "all"):
    """Pre-download and cache all receptor PDBQTs for Vina docking."""
    all_targets = []
    pathogens = list(TARGETS.keys()) if pathogen == "all" else [pathogen]
    for p in pathogens:
        all_targets.extend(TARGETS.get(p, []))
    results = {}
    for t in all_targets:
        pdb_id = t["target_pdb"]
        try:
            path = prepare_receptor(pdb_id)
            if path and path.exists():
                results[pdb_id] = {
                    "status": "ready",
                    "size_kb": round(path.stat().st_size / 1024, 1),
                    "path": str(path)
                }
            else:
                results[pdb_id] = {"status": "failed"}
        except Exception as e:
            results[pdb_id] = {"status": "error", "error": str(e)[:100]}
    ready = sum(1 for v in results.values() if v["status"] == "ready")
    return JSONResponse({
        "ready": ready, "total": len(results),
        "message": f"{ready}/{len(results)} receptors ready in {str(RECEPTOR_STORE)}",
        "receptors": results
    })


@app.get("/api/literature")
async def search_literature(target: str, pathogen: str, metals: bool = False):
    """Direct literature search endpoint."""
    try:
        papers = fetch_literature(target, pathogen, metals, max_papers=8)
        return JSONResponse({
            "query":   _build_query(target, pathogen, "metal" if metals else "organic"),
            "count":   len(papers),
            "papers":  papers,
            "sources": ["PubMed", "Semantic Scholar", "Europe PMC", "OpenAlex"]
        })
    except Exception as e:
        raise HTTPException(500, f"Literature search error: {e}")



def _target_record(pathogen: str, target_pdb: str) -> dict:
    for t in TARGETS.get(pathogen, []):
        if t["target_pdb"] == target_pdb:
            return t
    return {}


def _dock_one(atoms, descriptors, coord_fp, filename, pathogen, t,
              exhaustiveness, seed, work_dir, run_blind=True,
              ensemble_seeds=None):
    """Dock one ligand against one target and assemble the full record."""
    receptor = prepare_receptor(t["target_pdb"], work_dir)
    if not receptor:
        why = _RECEPTOR_PREP_LOG.get(t["target_pdb"], {}).get(
            "reason", "unknown (see /api/receptor_provenance)")
        return {"target_pdb": t["target_pdb"], "target_name": t["target_name"],
                "dock_status": "failed",
                "error": f"receptor preparation failed — {why}"}

    box = compute_active_site_box(t["target_pdb"], t.get("active_residues", []),
                                  chain=t.get("chain", ""))
    if box:
        cx, cy, cz, n_found = box[0], box[1], box[2], box[6]
    else:
        cx = cy = cz = None
        n_found = 0

    ligand = prepare_ligand(atoms, filename, work_dir, coord_fp=coord_fp)
    if ligand is None:
        return None

    _prep = LIGAND_PREP_RECORDS.get(str(ligand), {})
    _n_heavy = sum(1 for a in atoms if a["elem"].upper() != "H")
    _is_metal_cx = _prep.get("ligand_class", "").startswith("pre-formed")
    active = run_vina(ligand, receptor, t["target_pdb"], cx=cx, cy=cy, cz=cz,
                      n_heavy_atoms=_n_heavy, is_metal_complex=_is_metal_cx,
                      exhaustiveness=exhaustiveness, blind=False,
                      seed=seed, box=box)
    if not active.get("success"):
        return {"target_pdb": t["target_pdb"], "target_name": t["target_name"],
                "dock_status": "failed", "error": active.get("error")}

    blind = None
    if run_blind:
        blind = run_vina(ligand, receptor, t["target_pdb"],
                         n_heavy_atoms=_n_heavy, is_metal_complex=_is_metal_cx,
                         cx=None, cy=None, cz=None,
                         exhaustiveness=max(8, exhaustiveness // 2),
                         blind=True, seed=seed)

    # ── Pose: Vina's output, or nothing that counts ─────────────────
    #
    # Until v3.2.2 this analysed generate_pose() -- the input geometry
    # translated onto the pocket centroid with random jitter. Every contact
    # and every metal-coordination count was therefore measured against a
    # geometry Vina never produced, and the metal count was fed into the
    # composite score, so the jitter moved the number.
    _metal_map = {a["elem"]: a["elem"] for a in atoms
                  if a["elem"] in _METAL_PARAMS}
    pose = None
    pose_source = "SYNTHETIC_PLACEMENT_NOT_DOCKED"
    out_path = active.get("out_pdbqt")
    if out_path:
        pose = extract_vina_pose(Path(out_path), model=1, metal_map=_metal_map)
        if pose:
            pose_source = "VINA_DOCKED_POSE_MODEL_1"

    if pose is None:
        pose = generate_pose(atoms, t["target_pdb"],
                             t.get("active_residues", []), coord_fp=coord_fp)

    if pose_source == "VINA_DOCKED_POSE_MODEL_1":
        inter = contacts.analyse_pose(pose, receptor)
        n_metal = inter.get("n_metal_coordination", 0) if inter.get("ok") else 0
    else:
        # No real pose, so no real contacts. Refusing to compute them is the
        # only honest option: a contact table from a placed geometry looks
        # exactly like a contact table from a docked one.
        inter = {
            "ok": False,
            "error": (
                "No Vina output pose was available, so contacts were not "
                "computed. The geometry shown in the viewer is a placement of "
                "the input structure onto the pocket centroid, not a docked "
                "pose, and measuring contacts against it would produce a "
                "table indistinguishable from a real one."
            ),
            "pose_source": pose_source,
            "n_metal_coordination": 0, "n_hbonds": 0, "n_hydrophobic": 0,
            "metal_coordination": [], "hbonds": [], "hydrophobic": [],
            "contact_residues": [],
        }
        n_metal = 0

    composite = compute_composite_score(
        atoms, active["delta_g"], pathogen,
        coord_metal_contacts=n_metal, target_pdb=t["target_pdb"],
        descriptors=descriptors)

    selectivity = None
    if blind and blind.get("success"):
        selectivity = compute_selectivity(active["delta_g"], blind["delta_g"])

    _ki = layers.estimated_docking_ki(active["delta_g"])

    # ── Optional replicate ensemble ────────────────────────────────
    # Off by default: five seeds is five times the wall clock, and on a ten
    # target sweep that is the difference between a coffee and an afternoon.
    # When it is on, the extra seeds run at the SAME box and exhaustiveness
    # as the reported pose, so the SD describes the run that produced the
    # number rather than a differently-configured one.
    ensemble = None
    if ensemble_seeds:
        _vals = [active["delta_g"]]
        for _sd in ensemble_seeds:
            if int(_sd) == int(seed):
                continue
            _r = run_vina(ligand, receptor, t["target_pdb"], cx=cx, cy=cy, cz=cz,
                          n_heavy_atoms=_n_heavy, is_metal_complex=_is_metal_cx,
                          exhaustiveness=exhaustiveness, seed=int(_sd),
                          blind=False, box=box)
            if _r.get("success"):
                _vals.append(_r["delta_g"])
        ensemble = ensemble_statistics(
            _vals, [seed] + [x for x in ensemble_seeds if int(x) != int(seed)])

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
        "ensemble":                   ensemble,
        "delta_e_metal":              composite["delta_e_metal"],
        "corrected_score_kcal_per_mol": composite["corrected_score_kcal_per_mol"],
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
        "pose_source":   pose_source,
        "pose_is_docked": pose_source == "VINA_DOCKED_POSE_MODEL_1",
        "pose_note": (
            "Vina output, binding mode 1 of the modes retained. Metal element "
            "labels were restored after the Zn proxy for the purpose of "
            "applying element-specific coordination windows; the sampling "
            "itself still used zinc parameters."
            if pose_source == "VINA_DOCKED_POSE_MODEL_1" else
            "VIEWER PLACEMENT ONLY. Input geometry translated onto the pocket "
            "centroid. Not a docked pose. No contacts were computed from it "
            "and nothing derived from it entered the score."
        ),
        "data_quality":  "REAL_DOCKING",

        # ── Canonical reporting fields ─────────────────────────────────
        # The v3.2 aliases delta_g / vina_base_dg / ki_um are gone. They
        # existed so that v3.1 display code kept working, and the cost was
        # that a field named delta_g held a raw Vina score in one place and
        # a composite in another depending on which layer wrote it. Every
        # consumer now reads the layered name and there is one meaning per
        # name.
        "exhaustiveness": exhaustiveness,
        "ligand_preparation": LIGAND_PREP_RECORDS.get(str(ligand), {}),
        "estimated_docking_ki": _ki,
        "estimated_docking_ki_display": (_ki.get("display")
                                         if _ki.get("available")
                                         else "not available"),
        "scoring_function_error_kcal": layers.VINA_SCORING_ERROR_KCAL,
        "provenance": {
            "predicted_docking_affinity": prov.DIRECT_COMPUTATIONAL,
            "estimated_docking_ki": prov.DERIVED_ESTIMATE,
            "bamstu_composite_score": prov.PHYSICS_INFORMED,
        },
    }


@app.post("/api/dock_vina")
async def dock_vina_single(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    target_pdb: str = Form(...),
    target_name: str = Form("Unknown Target"),
    exhaustiveness: int = Form(DEFAULT_EXHAUSTIVENESS),
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
        box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

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
    # v3.4.15: the eNTRy rules, computed from the 3D geometry rather than
    # from the graph alone -- globularity needs coordinates. Lipinski stays
    # in the payload as the conventional reference, but for a Gram-negative
    # organism eNTRy is the criterion that predicts accumulation, and it is
    # what the interface now leads with.
    try:
        adme["entry_rules"] = ligand_prep.entry_rules(atoms)
    except Exception as exc:
        adme["entry_rules"] = {"verdict": "NOT_COMPUTED",
                               "error": f"{type(exc).__name__}: {exc}"}
    return JSONResponse({
        "vina_available": True, "success": True,
        "pathogen": pathogen, "exhaustiveness": exhaustiveness, "seed": seed,
        "hit": hit,
        # The Vina comparison panel reads data.comparison and reported
        # "No comparison data returned" without it. Rebuilt against v3.2
        # semantics: the two numbers are the raw score and the composite,
        # and they are labelled as different quantities rather than as two
        # estimates of the same one.
        "vina_result": {"success": True,
                        "delta_g": hit["vina_affinity_kcal_per_mol"]},
        # Layered payload. The panel no longer receives two numbers to
        # compare, because they were never comparable: one is in kcal/mol
        # and one is in arbitrary units. It receives the layers.
        "layers": layers.assemble_report(
            hit, engine=None, compound_id=stem),
        "comparison": {
            "layer1_predicted_docking_affinity_kcal_per_mol":
                hit["vina_affinity_kcal_per_mol"],
            "layer1_estimated_docking_ki":
                hit.get("estimated_docking_ki_display"),
            "layer2_physics_informed_score":
                hit["bamstu_composite_score"],
            "layer2_units": "arbitrary units",
            "relationship": "independent quantities, not two estimates of one",
            "interpretation": (
                "Layer 1 is an AutoDock Vina affinity in kcal/mol and an "
                "equilibrium Ki estimated from it. Layer 2 is a "
                "physics-informed coordination score in arbitrary units. "
                "They do not agree or disagree; subtracting one from the "
                "other produces nothing."
            ),
            "method_note": hit.get("metal_proxy", {}).get("consequence", ""),
        },
        "uncertainty": uq,
        "physicochemistry": adme,
        "engine": "AutoDock Vina 1.2 + BAMSTU composite screening layer v3.2",
        "reporting_rule": (
            "Report vina_affinity_kcal_per_mol as an AutoDock Vina 1.2 score, "
            "with the Zn metal proxy stated in Methods. The composite is a "
            "triage number in arbitrary units and is not a binding energy."
        ),
    })


def _whole_cell_sweep(content: str, stem: str, pathogen: str,
                      exhaustiveness: int, seed: int,
                      include_unverified: bool, ensemble_seeds: str,
                      job_id: str = "") -> dict:
    """The panel sweep itself, synchronous, callable from a worker thread.

    Factored out of the endpoint in v3.4.11 so the same code serves both the
    direct request and the background job. Every blocking call -- Vina,
    receptor download, receptor rebuild -- lives here, which is exactly why
    it must not run on the event loop.
    """
    try:
        atoms, descriptors, coord_fp = parse_xyz(content)
    except chem.ChemistryError as e:
        raise HTTPException(422, str(e))
    if not atoms:
        raise HTTPException(422, "No valid atomic coordinates found.")

    adme = compute_adme(descriptors, pathogen)
    # v3.4.15: the eNTRy rules, computed from the 3D geometry rather than
    # from the graph alone -- globularity needs coordinates. Lipinski stays
    # in the payload as the conventional reference, but for a Gram-negative
    # organism eNTRy is the criterion that predicts accumulation, and it is
    # what the interface now leads with.
    try:
        adme["entry_rules"] = ligand_prep.entry_rules(atoms)
    except Exception as exc:
        adme["entry_rules"] = {"verdict": "NOT_COMPUTED",
                               "error": f"{type(exc).__name__}: {exc}"}
    targets = TARGETS.get(pathogen, TARGETS["Staphylococcus aureus"])
    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)

    _ens_seeds = [int(x) for x in ensemble_seeds.replace(",", " ").split()] \
        if ensemble_seeds.strip() else None

    skipped = [t["target_name"] for t in targets
               if t.get("verified") is False and not include_unverified]
    active_targets = [t for t in targets
                      if t.get("verified") is not False or include_unverified]

    if not VINA_AVAILABLE:
        # No numeric score is emitted in this mode. v3.1 filled base_dg with a
        # literature midpoint and labelled the row; a labelled number still
        # gets transcribed into a table as a number.
        return {
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
        }

    hits, failures = [], []
    for t in active_targets:
        # v3.4.7: Vina, receptor downloads and receptor rebuilds are all
        # blocking. Running them directly in an async endpoint pins the
        # event loop for the whole sweep, so uvicorn cannot answer anything
        # else -- including the browser's own keepalives, which is what
        # surfaced as "Failed to fetch" on the client after a receptor
        # rebuild made one sweep long enough to cross the timeout.
        rec = _dock_one(atoms, descriptors, coord_fp, stem, pathogen, t,
                        exhaustiveness, seed, work_dir,
                        ensemble_seeds=_ens_seeds)
        if rec is None:
            failures.append({"target": t["target_name"],
                             "reason": "receptor preparation failed"})
        elif rec.get("dock_status") != "ok":
            failures.append({"target": t["target_name"],
                             "reason": rec.get("error", "docking failed")})
        else:
            hits.append(rec)

        # v3.4.11: report after every target rather than only at the end.
        # A sweep is minutes to tens of minutes and the client had no way to
        # tell a working run from a dead one until it either returned or the
        # connection dropped.
        if job_id:
            _job_progress(job_id, done=len(hits) + len(failures),
                          total=len(active_targets),
                          target=t["target_name"],
                          hits=len(hits), failures=len(failures))
        if job_id and JOBS.get(job_id, {}).get("cancel"):
            failures.append({"target": "(remaining)",
                             "reason": "cancelled by the user"})
            break

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

    return {
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
    }


# ════════════════════════════════════════════════════════════════════
#  VERIFICATION / VALIDATION / UQ ENDPOINTS
# ════════════════════════════════════════════════════════════════════



@app.post("/api/dock_whole_cell")
async def dock_whole_cell(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    exhaustiveness: int = Form(DEFAULT_EXHAUSTIVENESS),
    seed: int = Form(42),
    include_unverified: bool = Form(False),
    ensemble_seeds: str = Form(""),
):
    """Whole-panel sweep in a single request.

    Kept for compatibility and for short runs. For anything that takes more
    than a few minutes -- which is every rigid metal complex on this panel --
    use POST /api/dock_whole_cell/start and poll, because a single request
    held open that long does not survive a filtering proxy.
    """
    try:
        content = (await file.read()).decode("utf-8")
    except Exception as e:
        raise HTTPException(400, f"Cannot read .xyz file: {e}")
    stem = file.filename.split(".")[0]
    payload = await run_in_threadpool(
        _whole_cell_sweep, content, stem, pathogen, exhaustiveness, seed,
        include_unverified, ensemble_seeds)
    return JSONResponse(
        payload,
        status_code=503 if payload.get("mode") == "NO_DOCKING_ENGINE" else 200)


# ════════════════════════════════════════════════════════════════════
#  BACKGROUND JOBS  (v3.4.11)
# ════════════════════════════════════════════════════════════════════
#
# A panel sweep is one HTTP request held open for the whole run. For a
# flexible ligand at exhaustiveness 20 that is a few minutes and nobody
# notices. For a rigid silver complex it is exhaustiveness 32 over a 30 A
# box on 31 heavy atoms across nine targets -- five to six times the cost,
# tens of minutes -- and on a machine with a filtering appliance in the path
# the connection does not survive it. The browser reports "Failed to fetch"
# and the user cannot tell a dead run from a slow one, because a single
# request has no way to say "still working, four targets done".
#
# Splitting it is the fix. Start returns immediately with an id, the sweep
# runs in the background writing progress after every target, and the client
# polls. A dropped connection then costs nothing: the run continues and the
# results are collected on the next poll.

JOBS: dict = {}
JOB_RETENTION = 40          # completed jobs kept for collection
_JOB_LOCK = threading.Lock()


def _job_progress(job_id: str, **fields) -> None:
    with _JOB_LOCK:
        j = JOBS.get(job_id)
        if j is None:
            return
        j.update(fields)
        j["updated"] = datetime.now().isoformat(timespec="seconds")


def _job_reap() -> None:
    """Keep the most recent finished jobs, drop the rest."""
    with _JOB_LOCK:
        done = [(j.get("finished") or "", k) for k, j in JOBS.items()
                if j.get("state") in ("done", "failed", "cancelled")]
        for _ts, k in sorted(done)[:-JOB_RETENTION]:
            JOBS.pop(k, None)


@app.post("/api/dock_whole_cell/start")
async def dock_whole_cell_start(
    background: BackgroundTasks,
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    exhaustiveness: int = Form(DEFAULT_EXHAUSTIVENESS),
    seed: int = Form(42),
    include_unverified: bool = Form(False),
    ensemble_seeds: str = Form(""),
):
    """Queue a panel sweep and return its id immediately.

    The upload is read here, while the connection is certainly alive, and
    the sweep itself runs afterwards. Poll GET /api/jobs/{id}.
    """
    try:
        content = (await file.read()).decode("utf-8")
    except Exception as e:
        raise HTTPException(400, f"Cannot read .xyz file: {e}")
    stem = (file.filename or "sample").split(".")[0]

    job_id = f"{stem}_{datetime.now().strftime('%H%M%S')}_{uuid.uuid4().hex[:6]}"
    n_targets = len([t for t in TARGETS.get(
        pathogen, TARGETS["Staphylococcus aureus"])
        if t.get("verified") is not False or include_unverified])
    with _JOB_LOCK:
        JOBS[job_id] = {
            "job_id": job_id, "sample": stem, "pathogen": pathogen,
            "state": "queued", "done": 0, "total": n_targets,
            "hits": 0, "failures": 0, "target": None, "cancel": False,
            "started": datetime.now().isoformat(timespec="seconds"),
            "updated": datetime.now().isoformat(timespec="seconds"),
            "exhaustiveness": exhaustiveness, "seed": seed,
            "result": None, "error": None,
        }

    def _run():
        _job_progress(job_id, state="running")
        try:
            payload = _whole_cell_sweep(
                content, stem, pathogen, exhaustiveness, seed,
                include_unverified, ensemble_seeds, job_id=job_id)
            _job_progress(job_id, state="done", result=payload,
                          finished=datetime.now().isoformat(timespec="seconds"))
            print(f"[BAMSTU] job {job_id}: finished.", flush=True)
        except HTTPException as exc:
            _job_progress(job_id, state="failed", error=str(exc.detail),
                          finished=datetime.now().isoformat(timespec="seconds"))
            print(f"[BAMSTU] job {job_id}: FAILED — {exc.detail}", flush=True)
        except Exception as exc:
            _job_progress(job_id, state="failed",
                          error=f"{type(exc).__name__}: {exc}",
                          finished=datetime.now().isoformat(timespec="seconds"))
            print(f"[BAMSTU] job {job_id}: FAILED — "
                  f"{type(exc).__name__}: {exc}", flush=True)
        _job_reap()

    background.add_task(run_in_threadpool, _run)
    return JSONResponse({"job_id": job_id, "state": "queued",
                         "total": n_targets,
                         "poll": f"/api/jobs/{job_id}"}, status_code=202)


@app.get("/api/jobs/{job_id}")
async def job_status(job_id: str, include_result: bool = True):
    """Progress, and the full payload once the sweep has finished."""
    with _JOB_LOCK:
        j = JOBS.get(job_id)
        if j is None:
            raise HTTPException(404, f"No job '{job_id}'. Jobs are held in "
                                     f"memory, so a server restart clears "
                                     f"them; the session on disk is not "
                                     f"affected.")
        out = {k: v for k, v in j.items()
               if k != "result" or (include_result and j["state"] == "done")}
    out["percent"] = (round(100 * out["done"] / out["total"])
                      if out.get("total") else 0)
    return JSONResponse(out)


@app.get("/api/jobs")
async def job_list():
    with _JOB_LOCK:
        return JSONResponse({"jobs": [
            {k: v for k, v in j.items() if k != "result"}
            for j in sorted(JOBS.values(), key=lambda x: x["started"],
                            reverse=True)]})


@app.post("/api/jobs/{job_id}/cancel")
async def job_cancel(job_id: str):
    """Stop after the target now running; keep whatever finished."""
    with _JOB_LOCK:
        j = JOBS.get(job_id)
        if j is None:
            raise HTTPException(404, f"No job '{job_id}'.")
        if j["state"] not in ("queued", "running"):
            return JSONResponse({"job_id": job_id, "state": j["state"],
                                 "message": "already finished"})
        j["cancel"] = True
    return JSONResponse({"job_id": job_id, "cancelling": True,
                         "message": "Will stop after the target currently "
                                    "docking. Targets already completed are "
                                    "kept."})


@app.post("/api/render/interaction")
async def render_interaction(payload: dict = Body(...)):
    """Render the geometry-faithful interaction diagram for one hit, on demand.

    v3.4.12. The browser used to draw its own diagram: residues spread evenly
    round a circle in list order, the ligand a blank ellipse, every line
    terminating on the edge of a shape. It carried no spatial information at
    all. The server already draws the real thing from the docked coordinates
    for the session archive, so the interface now asks for that instead of
    inventing a second, worse one. One renderer, one figure, whether it is
    read on screen or pulled out of a ZIP.
    """
    hit = payload.get("hit")
    if not isinstance(hit, dict):
        raise HTTPException(422, "A 'hit' object is required.")
    if hit.get("pose_is_docked") is False:
        raise HTTPException(
            409, "This geometry is a placement on the pocket centroid, not a "
                 "docked pose. A diagram drawn from it is indistinguishable "
                 "on the page from one drawn from a real result.")
    pose = hit.get("docked_pose") or ""
    if not pose.strip():
        raise HTTPException(422, "The hit carries no pose to draw.")

    pdb_id = hit.get("target_pdb", "")
    rp = RECEPTOR_STORE / f"{pdb_id}.pdbqt"
    if not rp.exists():
        rp = RECEPTOR_STORE / f"{pdb_id}.pdb"
    if not rp.exists():
        raise HTTPException(
            404, f"Receptor {pdb_id} is not in the local store, and contacts "
                 f"can only be measured against the structure that was "
                 f"docked into.")

    inter = hit.get("interactions") or {}
    if not inter.get("ok"):
        inter = contacts.analyse_pose(pose, rp)

    dg = hit.get("vina_affinity_kcal_per_mol")
    tmp = Path(tempfile.mkdtemp(prefix="bamstu_svg_")) / "interaction.svg"
    try:
        info = figures.render_interaction_svg(
            pose, rp, inter, tmp,
            title=f"{payload.get('sample', 'Ligand')} in the "
                  f"{hit.get('target_name', '')} binding site ({pdb_id})",
            subtitle=(f"E_vina = {dg:.2f} kcal/mol   |   "
                      f"{inter.get('n_metal_coordination', 0)} coordination, "
                      f"{inter.get('n_hbonds', 0)} hydrogen bond, "
                      f"{inter.get('n_hydrophobic', 0)} hydrophobic contact(s)"
                      if dg is not None else ""))
        svg = tmp.read_text(encoding="utf-8")
    except Exception as exc:
        raise HTTPException(500, f"Diagram could not be drawn: {exc}")
    finally:
        shutil.rmtree(tmp.parent, ignore_errors=True)

    return JSONResponse({
        "svg": svg,
        "ligand_atoms": info["ligand_atoms"],
        "n_residues": info["n_residues"],
        "unresolved_contacts": info["unresolved_contacts"],
        "renderer": info["renderer"],
        "note": ("Positions are the docked coordinates projected onto the "
                 "least-squares plane of the ligand heavy atoms. One axis is "
                 "compressed, so residues that look adjacent may be "
                 "separated along the plane normal."),
    })


@app.post("/api/render/complex_pdb")
async def render_complex(payload: dict = Body(...)):
    """The pose plus its binding-site shell, as one PDB for the 3D viewer.

    Whole residues within 8 A of the ligand, so side chains are never cut
    mid-bond and the viewer can draw cartoon ribbons over a real backbone
    rather than a cloud of selected atoms.
    """
    hit = payload.get("hit")
    if not isinstance(hit, dict):
        raise HTTPException(422, "A 'hit' object is required.")
    pose = hit.get("docked_pose") or ""
    if not pose.strip():
        raise HTTPException(422, "The hit carries no pose.")
    pdb_id = hit.get("target_pdb", "")
    rp = RECEPTOR_STORE / f"{pdb_id}.pdbqt"
    if not rp.exists():
        rp = RECEPTOR_STORE / f"{pdb_id}.pdb"
    if not rp.exists():
        raise HTTPException(404, f"Receptor {pdb_id} not in the local store.")
    tmp = Path(tempfile.mkdtemp(prefix="bamstu_cx_")) / "complex.pdb"
    try:
        info = figures.render_complex_pdb(
            pose, rp, tmp, shell=float(payload.get("shell", 8.0)))
        text = tmp.read_text(encoding="utf-8")
    except Exception as exc:
        raise HTTPException(500, f"Complex could not be built: {exc}")
    finally:
        shutil.rmtree(tmp.parent, ignore_errors=True)
    return JSONResponse({"pdb": text, **{k: v for k, v in info.items()
                                         if k != "path"}})


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
    exhaustiveness: int = Form(DEFAULT_EXHAUSTIVENESS),
    seeds: str = Form(""),
):
    """Dock one ligand at several fixed seeds and report mean +/- SD.

    A single seed gives a point estimate with no way to tell a real gap
    between two targets from search noise. Five seeds give the spread, and
    the spread is what decides whether a 0.3 kcal/mol difference means
    anything. The guardrail is SD <= 0.5 kcal/mol; above it the search has
    not converged, and the fix is more exhaustiveness rather than more
    seeds, because averaging a wider distribution does not narrow it.

    The seeds are fixed and reported. An ensemble over unrecorded random
    seeds cannot be regenerated, which makes its SD a number nobody can
    check.
    """
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required. There is no "
                                 "estimate fallback for an ensemble.")
    seed_list = ([int(x) for x in seeds.replace(",", " ").split()]
                 if seeds.strip()
                 else list(ENSEMBLE_SEEDS)[:max(int(n_replicates), 1)])
    content = (await file.read()).decode("utf-8")
    atoms, descriptors, coord_fp = parse_xyz(content)
    t = _target_record(pathogen, target_pdb) or {}
    work_dir = Path(__file__).parent / "data"
    receptor = prepare_receptor(target_pdb, work_dir)
    stem = file.filename.split(".")[0]
    ligand = prepare_ligand(atoms, stem, work_dir, coord_fp=coord_fp)
    prep = LIGAND_PREP_RECORDS.get(str(ligand), {})
    n_heavy = sum(1 for a in atoms if a["elem"].upper() != "H")
    is_cx = str(prep.get("ligand_class", "")).startswith("pre-formed")
    box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

    runs, values, policy = [], [], None
    for sd in seed_list:
        r = run_vina(ligand, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exhaustiveness, seed=sd, box=box,
                     n_heavy_atoms=n_heavy, is_metal_complex=is_cx)
        ok = bool(r.get("success"))
        policy = policy or r.get("box_policy")
        runs.append({"seed": sd, "success": ok,
                     "affinity_kcal_per_mol": r.get("delta_g"),
                     "error": None if ok else str(r.get("error", ""))[:160]})
        if ok:
            values.append(r["delta_g"])

    stats = ensemble_statistics(values, seed_list)

    def _run(seed=42, exhaustiveness=exhaustiveness):
        r = run_vina(ligand, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exhaustiveness, seed=seed, box=box,
                     n_heavy_atoms=n_heavy, is_metal_complex=is_cx)
        if not r.get("success"):
            raise RuntimeError(r.get("error", "vina failed"))
        return r["delta_g"]

    return JSONResponse({
        "compound": stem,
        "target": target_pdb,
        "ligand_class": prep.get("ligand_class"),
        "torsdof": prep.get("torsdof"),
        "exhaustiveness": exhaustiveness,
        "box_policy": policy,
        "runs": runs,
        "ensemble": stats,
        "uq": vvuq.ensemble_uq(_run, n_replicates=n_replicates,
                               exhaustiveness=exhaustiveness),
        "determinism": vvuq.check_determinism(_run, n_repeats=2, seed=42),
        "provenance": {
            "quantity": "Predicted Docking Affinity (AutoDock Vina)",
            "category": prov.DIRECT_COMPUTATIONAL},
        "methods_text": (
            f"Each ligand was docked at {len(seed_list)} independent random "
            f"seeds ({', '.join(str(x) for x in seed_list)}) at "
            f"exhaustiveness {exhaustiveness}. Affinities are reported as "
            f"mean +/- SD across seeds: {stats.get('report_as', 'n/a')}. Runs "
            f"with SD above {ENSEMBLE_SD_LIMIT_KCAL} kcal/mol were treated as "
            "unconverged. This SD describes search reproducibility only; it "
            "is separate from the roughly 1.5 kcal/mol error of the scoring "
            "function, which does not shrink with more seeds."),
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
    box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

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
    native_path.write_text("\n".join(native_lines) + "\nEND\n")

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
    box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

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


@app.get("/api/receptor_pdbqt")
async def receptor_pdbqt(pdb_id: str):
    """Return the prepared receptor PDBQT as text.

    Exists so a session archive can carry the exact receptor the scores were
    produced against. Re-fetching from RCSB months later is not the same
    thing: depositions get superseded and residue numbering moves, so a
    reload that re-downloads is a different calculation wearing the old
    results' name. The PDBQT in the ZIP is the file Vina actually read,
    REMARK records and all, including which metals were proxied to Zn.
    """
    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        path = prepare_receptor(pdb_id, work_dir)
    except Exception as e:
        raise HTTPException(502, f"Could not prepare {pdb_id}: {e}")
    if not path or not Path(path).exists():
        raise HTTPException(404, f"No prepared receptor for {pdb_id}.")
    return PlainTextResponse(Path(path).read_text(encoding="utf-8",
                                                  errors="ignore"))


@app.post("/api/vvuq/roc_auc")
async def vvuq_roc_auc(
    file: UploadFile = File(...),
    pathogen: str = Form(...),
    target_pdb: str = Form(...),
    n_decoys: int = Form(25),
    exhaustiveness: int = Form(8),
    tanimoto_max: float = Form(0.6),
    seed: int = Form(42),
):
    """Property-matched decoy enrichment against one target.

    Every molecule -- active and decoy -- goes through the same receptor,
    the same box, the same exhaustiveness and the same seed. If the decoys
    were docked under different settings the comparison would measure the
    settings.
    """
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required for enrichment.")

    content = (await file.read()).decode("utf-8")
    t = _target_record(pathogen, target_pdb)
    if not t:
        raise HTTPException(404, f"{target_pdb} is not in the {pathogen} panel.")

    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    receptor = prepare_receptor(target_pdb, work_dir)
    if not receptor:
        why = _RECEPTOR_PREP_LOG.get(target_pdb, {}).get("reason", "unknown")
        raise HTTPException(502, f"Receptor preparation failed — {why}")
    box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

    def _dock_mol(mol, name: str) -> float:
        """Dock one RDKit molecule and return its Vina score."""
        conf = mol.GetConformer()
        atoms = [{"elem": a.GetSymbol(),
                  "x": conf.GetAtomPosition(a.GetIdx()).x,
                  "y": conf.GetAtomPosition(a.GetIdx()).y,
                  "z": conf.GetAtomPosition(a.GetIdx()).z}
                 for a in mol.GetAtoms()]
        lig = prepare_ligand(atoms, f"enrich_{name}", work_dir,
                             coord_fp=chem.coordinate_fingerprint(atoms))
        if lig is None:
            raise RuntimeError(f"ligand preparation failed for {name}")
        r = run_vina(lig, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exhaustiveness, seed=seed, box=box)
        if not r.get("success"):
            raise RuntimeError(r.get("error", "vina failed"))
        return float(r["delta_g"])

    result = vvuq.compute_decoy_enrichment(
        content, Path(receptor), n_decoys=n_decoys,
        exhaustiveness=exhaustiveness, dock_fn=_dock_mol,
        tanimoto_max=tanimoto_max, seed=seed)
    result["target_pdb"] = target_pdb
    result["pathogen"] = pathogen
    result["prerequisite"] = (
        "Enrichment does not stand alone. If POST /api/vvuq/redock has not "
        "passed at 2.0 A on this target, a good enrichment number only shows "
        "the scoring function is consistently wrong in a useful direction."
    )
    return JSONResponse(result)


@app.post("/api/vvuq/calibrate")
async def vvuq_calibrate(payload: dict = Body(...)):
    """Fit the composite score against measured activity for a congeneric set.

    Body:
      {"pathogen": "...", "target_pdb": "...",
       "response_label": "log(MIC)" | "log(Ki)",
       "exhaustiveness": 16, "seed": 42,
       "training_data": [{"name": "...", "xyz_content": "...",
                          "experimental_value": 0.5}, ...]}

    This is the endpoint that decides whether the composite score can be
    reported as a calibrated quantity or must stay in arbitrary units. A
    FAIL here is a real and publishable result: it says the correction layer
    does not track measured activity for this series.
    """
    if not VINA_AVAILABLE:
        raise HTTPException(503, "AutoDock Vina required for calibration.")

    training = payload.get("training_data") or []
    pathogen = payload.get("pathogen", "Klebsiella pneumoniae")
    target_pdb = payload.get("target_pdb")
    label = payload.get("response_label", "log(MIC)")
    exh = int(payload.get("exhaustiveness", 16))
    seed = int(payload.get("seed", 42))
    if not target_pdb:
        raise HTTPException(422, "target_pdb is required.")

    t = _target_record(pathogen, target_pdb)
    if not t:
        raise HTTPException(404, f"{target_pdb} is not in the {pathogen} panel.")

    work_dir = Path(__file__).parent / "data"
    work_dir.mkdir(parents=True, exist_ok=True)
    receptor = prepare_receptor(target_pdb, work_dir)
    if not receptor:
        why = _RECEPTOR_PREP_LOG.get(target_pdb, {}).get("reason", "unknown")
        raise HTTPException(502, f"Receptor preparation failed — {why}")
    box = compute_active_site_box(target_pdb, t.get("active_residues", []),
                                  chain=t.get("chain", ""))

    def _dock_xyz(xyz_content: str, name: str) -> dict:
        """Dock one training compound and return its feature row."""
        atoms, descriptors, coord_fp = parse_xyz(xyz_content)
        lig = prepare_ligand(atoms, f"calib_{name}", work_dir, coord_fp=coord_fp)
        if lig is None:
            raise RuntimeError(f"ligand preparation failed for {name}")
        r = run_vina(lig, receptor, target_pdb,
                     cx=box[0] if box else None, cy=box[1] if box else None,
                     cz=box[2] if box else None,
                     exhaustiveness=exh, seed=seed, box=box)
        if not r.get("success"):
            raise RuntimeError(r.get("error", "vina failed"))

        # Metal contacts must come from the REAL docked pose. Taking them
        # from a placed geometry would put a placement artefact into a
        # fitted coefficient, which is the v3.2.2 defect wearing a lab coat.
        n_metal = 0
        out_path = r.get("out_pdbqt")
        if out_path:
            metal_map = {a["elem"]: a["elem"] for a in atoms
                         if a["elem"] in _METAL_PARAMS}
            pose = extract_vina_pose(Path(out_path), model=1,
                                     metal_map=metal_map)
            if pose:
                inter = contacts.analyse_pose(pose, receptor)
                if inter.get("ok"):
                    n_metal = inter.get("n_metal_coordination", 0)

        return {
            "vina": float(r["delta_g"]),
            "metal_contacts": n_metal,
            "n_rot": descriptors.get("n_rotatable_bonds") or 0,
        }

    result = vvuq.calibrate_composite_model(
        training, pathogen, target_pdb,
        dock_fn=_dock_xyz, response_label=label, seed=seed)
    result["exhaustiveness"] = exh
    result["seed"] = seed
    return JSONResponse(result)


# ════════════════════════════════════════════════════════════════════
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

    # v3.4.5: render at save time rather than leaving it to a separate
    # manual step. The 2026-08-03 session shipped with no figures at all
    # because nobody ran the renderer, and nothing reported that it had not
    # run. Each docked hit now yields four artefacts inside the session, and
    # therefore inside the ZIP: a 2D interaction SVG, a 3D pose PNG, the
    # binding-site complex as PDB for any desktop viewer, and a
    # self-contained interactive HTML view. A failure here never blocks the
    # save; it is recorded and reported.
    fig_dir = sample_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    n_figs, fig_errors, fig_skipped = 0, [], []
    for h in hits:
        if not h.get("pose_is_docked"):
            fig_skipped.append({"target_pdb": h.get("target_pdb"),
                                "reason": "not a docked pose"})
            continue
        rp = RECEPTOR_STORE / f"{h.get('target_pdb')}.pdbqt"
        if not rp.exists():
            rp = RECEPTOR_STORE / f"{h.get('target_pdb')}.pdb"
        if not rp.exists():
            # No receptor, no figure. Recorded as a skip with the reason
            # rather than four identical file-not-found strings, so the
            # audit can distinguish "nobody rendered these" from "the
            # receptor this run used is no longer in the store".
            fig_skipped.append({"target_pdb": h.get("target_pdb"),
                                "reason": "receptor not in the local store; "
                                          "figures cannot be drawn without "
                                          "the structure that was docked into"})
            continue
        try:
            res = figures.render_hit_figures(sample, h, rp, fig_dir)
            n_figs += sum(1 for k in ("png", "svg", "complex_pdb",
                                      "viewer_html") if res.get(k))
            fig_errors.extend(res.get("errors") or [])
        except Exception as exc:
            fig_errors.append(f"{h.get('target_pdb')}: "
                              f"{type(exc).__name__}: {exc}")

    # v3.4.7: say out loud what happened. The figure step failing silently
    # is what produced a session with eight poses and no images, and the
    # only place the reason existed was a JSON field the interface never
    # displayed. The server console is the one surface a user is already
    # watching while a run is in progress.
    if n_figs:
        print(f"[BAMSTU] {sample}: {n_figs} figure file(s) written to "
              f"{fig_dir}", flush=True)
    else:
        print(f"[BAMSTU] {sample}: NO figures written. "
              f"{len(fig_skipped)} skipped, {len(fig_errors)} failed.",
              flush=True)
    for e in fig_errors[:6]:
        print(f"[BAMSTU]   figure error: {e}", flush=True)
    for sk in fig_skipped[:6]:
        print(f"[BAMSTU]   skipped {sk.get('target_pdb')}: "
              f"{sk.get('reason')}", flush=True)

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
        "\n".join(rows) + "\n", encoding="utf-8")

    moa = data.get("MoA") or {}
    if moa:
        if moa.get("status") == "QUARANTINED":
            body = ("QUARANTINED — this text made claims the pipeline cannot "
                    "support. Do not paste it into a manuscript.\n\n"
                    + "\n".join(f" - {v.get('matched')}: {v.get('requirement')}"
                                 for v in (moa.get("claim_filter", {})
                                           .get("violations") or []))
                    + "\n\nSafe alternative:\n"
                    + str(moa.get("safe_alternative", ""))
                    + "\n\nQuarantined text (inspection only):\n"
                    + str(moa.get("quarantined_text", "")))
        else:
            body = str(moa.get("inference") or "")
        body += ("\n\nREFERENCES\n" + "\n".join(moa.get("references") or [])
                 + "\n\nEvery citation must be opened and checked against the "
                   "sentence it is attached to before submission.")
        (sample_dir / f"{sample}_MoA.txt").write_text(body, encoding="utf-8")

    # ── Figures ─────────────────────────────────────────────────────
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
        "n_figures": n_figs,
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
        "n_figure_files": n_figs,
        "figure_errors": fig_errors[:10],
        "figures_skipped": fig_skipped,
        "figures": {
            "rendered": sum(1 for f in fig_report if f.get("png") or f.get("svg")),
            "skipped_not_docked": sum(1 for f in fig_report if f.get("skipped")),
            "errors": [e for f in fig_report for e in (f.get("errors") or [])],
        },
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


def _session_summary_csv(manifest: dict, sdir: Path) -> str:
    """Flatten a session into one CSV, with the layer of each column in it.

    A spreadsheet is where these numbers go to be misread: a reader opening
    it months later has no docstring, no endpoint and no UI badge, only
    column headers. So the header row carries the layer and the units, the
    provenance category is a column of its own, and the pose provenance
    sits next to the affinity rather than three columns away. A row whose
    geometry was a viewer placement carries no affinity at all -- an empty
    cell is harder to transcribe into a thesis table than a number with a
    caveat somewhere else in the file.
    """
    import csv as _csv
    import io as _io

    buf = _io.StringIO()
    w = _csv.writer(buf, lineterminator="\n")
    w.writerow([
        "sample", "target_name", "target_pdb", "compartment", "organism",
        "docking_semantics", "pose_source", "pose_is_docked",
        "L1_vina_affinity_kcal_per_mol", "L1_estimated_docking_ki",
        "L1_ensemble_mean_kcal", "L1_ensemble_sd_kcal", "L1_ensemble_converged",
        "L2_coordination_score_arbitrary_units", "L2_metal_donor_contacts",
        "L2_steric_clashes", "L2_calibration_status",
        "ligand_class", "torsdof", "seed", "exhaustiveness",
        "box_x", "box_y", "box_z", "selectivity_mode",
        "provenance_L1", "provenance_L2",
    ])

    for entry in manifest.get("samples", []):
        name = entry["name"]
        rp = sdir / name / f"{name}_result.json"
        if not rp.is_file():
            continue
        try:
            d = json.loads(rp.read_text(encoding="utf-8"))
        except Exception:
            continue
        hits = ((d.get("inhibition_targets") or [])
                + (d.get("permeability_targets") or [])
                + (d.get("other_semantics") or [])
                or d.get("hits") or [])
        for h in hits:
            docked = h.get("pose_is_docked") is not False
            ens = h.get("ensemble") or {}
            prep = h.get("ligand_preparation") or {}
            inter = h.get("interactions") or {}
            gs = (h.get("grid_size") or (h.get("box_policy") or {}).get("size")
                  or ["", "", ""])
            w.writerow([
                name, h.get("target_name", ""), h.get("target_pdb", ""),
                h.get("compartment", ""), h.get("organism", ""),
                h.get("docking_semantics", ""),
                h.get("pose_source", ""), docked,
                h.get("vina_affinity_kcal_per_mol", "") if docked else "",
                h.get("estimated_docking_ki_display", "") if docked else "",
                ens.get("mean_kcal", ""), ens.get("sd_kcal", ""),
                ens.get("within_guardrail", ""),
                h.get("bamstu_composite_score", "") if docked else "",
                inter.get("n_metal_coordination", "") if inter.get("ok") else "",
                inter.get("n_clashes", "") if inter.get("ok") else "",
                h.get("calibration_status", "UNCALIBRATED"),
                prep.get("ligand_class", ""), prep.get("torsdof", ""),
                h.get("seed", ""), h.get("exhaustiveness", ""),
                gs[0] if len(gs) > 0 else "", gs[1] if len(gs) > 1 else "",
                gs[2] if len(gs) > 2 else "",
                (h.get("selectivity") or {}).get("mode", ""),
                prov.DIRECT_COMPUTATIONAL, prov.PHYSICS_INFORMED,
            ])

    w.writerow([])
    w.writerow(["NOTE", "L1_vina_affinity_kcal_per_mol is an AutoDock Vina "
                "1.2 empirical score with roughly 1.5 kcal/mol error. "
                "Transition metals were typed with zinc parameters during "
                "the search."])
    w.writerow(["NOTE", "L1_estimated_docking_ki = exp(dG/RT) at 298.15 K "
                "from that score alone. It is a computational estimate, not "
                "a measured inhibition constant."])
    w.writerow(["NOTE", "L2_coordination_score is in ARBITRARY UNITS. It is "
                "not a binding free energy and is not convertible to a Ki. "
                "Do not subtract it from the L1 column."])
    w.writerow(["NOTE", "L1_ensemble_sd_kcal is seed-to-seed search "
                "reproducibility only. It is not the accuracy of the "
                "scoring function and does not shrink with more seeds."])
    w.writerow(["NOTE", "Rows with pose_is_docked=False have empty affinity "
                "and score columns. The stored geometry for those targets is "
                "a placement on the pocket centroid, not a docking result."])
    return buf.getvalue()


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

    read_failures: list = []
    for entry in manifest.get("samples", []):
        rp = sdir / entry["name"] / f"{entry['name']}_result.json"
        if not rp.is_file():
            continue
        try:
            d = json.loads(rp.read_text(encoding="utf-8"))
        except Exception as exc:
            # v3.4.7: a silent continue here is why a ZIP could report
            # "0 receptors" with no indication that a result file had failed
            # to parse. The failure is now carried into the archive.
            read_failures.append(f"{rp.name}: {type(exc).__name__}: {exc}")
            continue
        for h in ((d.get("inhibition_targets") or [])
                  + (d.get("permeability_targets") or [])
                  + (d.get("other_semantics") or [])
                  + (d.get("hits") or [])):
            if h.get("target_pdb"):
                pdbs.add(h["target_pdb"])

    zpath = sdir / f"BAMSTU_Session_{_safe_id(session_id)}.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(sdir.rglob("*")):
            if f.is_file() and f.suffix != ".zip":
                z.write(f, f.relative_to(sdir).as_posix())
        missing_receptors = []
        for p in sorted(pdbs):
            rec = RECEPTOR_STORE / f"{p}.pdbqt"
            if rec.is_file():
                z.writestr(f"receptors/{p}.pdbqt",
                           rec.read_text(encoding="utf-8", errors="ignore"))
            else:
                missing_receptors.append(p)

        # A receptor can be absent legitimately: v3.4.6 deletes any file
        # built before the HETATM policy so it will be rebuilt on next use.
        # Saying so beats an archive that quietly contains fewer structures
        # than the run it documents.
        if missing_receptors or read_failures:
            z.writestr("RECEPTOR_AND_PARSE_NOTES.txt",
                       "Receptors referenced by this session but not found "
                       "in the local store:\n"
                       + "".join(f"  {p}\n" for p in missing_receptors)
                       + ("\nThese were most likely deleted for rebuild "
                          "under the v3.4.6 HETATM policy, which removes "
                          "co-crystallised ligands. Re-run the affected "
                          "targets to regenerate them.\n"
                          if missing_receptors else "")
                       + ("\nResult files that could not be parsed:\n"
                          + "".join(f"  {r}\n" for r in read_failures)
                          if read_failures else ""))
        try:
            z.writestr("BAMSTU_Session_Summary.csv",
                       _session_summary_csv(manifest, sdir))
        except Exception as exc:
            z.writestr("BAMSTU_Session_Summary_FAILED.txt",
                       f"The CSV summary could not be built: {exc}\n"
                       "The per-sample result JSON files in this archive are "
                       "the source of truth and are unaffected.\n")
        z.writestr("README.txt",
                   "BAMSTU Digital Pathogen Pro — session archive\n"
                   f"Session: {session_id}\n"
                   f"Samples: {manifest.get('n_samples')}\n"
                   f"Built:   {datetime.now().isoformat(timespec='seconds')}\n\n"
                   "RELOADING\n"
                   "  Autosampler panel -> Reload session, or\n"
                   f"  GET /api/session/{session_id}\n\n"
                   "WHAT THE NUMBERS ARE\n"
                   "  Vina_score_kcal_per_mol  AutoDock Vina 1.2 score. Not a\n"
                   "                           measured binding energy. Scoring\n"
                   "                           error is roughly 1.5-2.0 kcal/mol,\n"
                   "                           so targets closer than that are\n"
                   "                           not separable.\n"
                   "  BAMSTU_composite_au      Arbitrary units, uncalibrated.\n"
                   "  Pose_is_docked           FALSE means the geometry is a\n"
                   "                           placement, not a Vina pose, and no\n"
                   "                           contacts were computed from it.\n\n"
                   "NOT HERE: no inhibition constants. A Ki is measured by\n"
                   "enzyme kinetics.\n")

    return FileResponse(str(zpath), media_type="application/zip",
                        filename=zpath.name)


@app.get("/api/session/{session_id}/figure")
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
            chk = _filter_mechanistic_claims(text, 0)
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


def _missing_figures(session_id: str) -> dict:
    """Which figures a session is short of, per sample and target.

    The audit used to return a bare 409 when a figure was absent, which told
    the user that something had failed but not what, and the usual cause was
    a historical session rendered before the figure writer existed. That is
    repairable without re-docking, so the payload names the files and points
    at the endpoint that makes them.
    """
    sdir = SESSION_STORE / _safe_id(session_id)
    mpath = sdir / "manifest.json"
    if not mpath.exists():
        return {"session_found": False, "missing": [], "n_missing": 0}
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    missing = []
    for entry in manifest.get("samples", []):
        name = entry["name"]
        fdir = sdir / name / "figures"
        for hit in entry.get("hits", []):
            pdb = hit.get("target_pdb", "?")
            for suffix, kind in (("_3D_Pose.png", "3D pose"),
                                 ("_LigPlot_2D.svg", "2D interaction map")):
                f = fdir / f"{name}_{pdb}{suffix}"
                if not f.exists() or f.stat().st_size < 200:
                    missing.append({
                        "sample": name, "target_pdb": pdb, "kind": kind,
                        "expected_file": f.name,
                        "pose_is_docked": hit.get("pose_is_docked"),
                        "repairable_without_redocking": bool(
                            hit.get("docked_pose")
                            and hit.get("pose_is_docked")),
                    })
    return {"session_found": True, "missing": missing,
            "n_missing": len(missing)}


@app.post("/api/session/{session_id}/render_figures")
async def render_session_figures(session_id: str, overwrite: bool = False):
    """Re-generate missing PNG and SVG figures without re-docking.

    Rendering is headless: the pose PNG is drawn by matplotlib with the Agg
    backend and the interaction map is written as SVG text, so neither needs
    an X display, a GUI toolkit or a GPU. What it cannot do is invent a pose.
    A target whose stored geometry was a viewer placement rather than a
    docked pose is skipped and reported, because a figure drawn from a
    placement is indistinguishable on screen from one drawn from a docking
    result, and that is precisely how the v3.1 contact figures got into a
    dissertation.
    """
    sdir = SESSION_STORE / _safe_id(session_id)
    mpath = sdir / "manifest.json"
    if not mpath.exists():
        raise HTTPException(404, f"No session '{session_id}'.")
    manifest = json.loads(mpath.read_text(encoding="utf-8"))

    rendered, skipped, failed = [], [], []
    for entry in manifest.get("samples", []):
        name = entry["name"]
        fdir = sdir / name / "figures"
        fdir.mkdir(parents=True, exist_ok=True)
        for hit in entry.get("hits", []):
            pdb = hit.get("target_pdb", "?")
            if hit.get("pose_is_docked") is False:
                skipped.append({
                    "sample": name, "target_pdb": pdb,
                    "reason": ("Stored geometry is a viewer placement, not a "
                               "docked pose. Re-dock this target; a figure "
                               "drawn from a placement looks exactly like a "
                               "real one."),
                })
                continue
            pose = hit.get("docked_pose") or ""
            if not pose.strip():
                skipped.append({"sample": name, "target_pdb": pdb,
                                "reason": "No stored pose text to render."})
                continue
            tag = f"{name}_{pdb}"
            have = [(fdir / f"{tag}_3D_Pose.png").exists(),
                    (fdir / f"{tag}_LigPlot_2D.svg").exists()]
            if all(have) and not overwrite:
                skipped.append({"sample": name, "target_pdb": pdb,
                                "reason": "Both figures already present. "
                                          "Pass overwrite=true to redraw."})
                continue
            try:
                receptor = RECEPTOR_STORE / f"{pdb}.pdb"
                res = figures.render_hit_figures(name, hit, receptor, fdir)
                if res.get("skipped"):
                    skipped.append({"sample": name, "target_pdb": pdb,
                                    "reason": res.get("reason", res["skipped"])})
                elif res.get("errors") and not (res["png"] or res["svg"]):
                    failed.append({"sample": name, "target_pdb": pdb,
                                   "error": "; ".join(res["errors"])})
                else:
                    rendered.append({
                        "sample": name, "target_pdb": pdb,
                        "files": [Path(f).name for f in
                                  (res["png"], res["svg"]) if f],
                        "partial_errors": res.get("errors") or None})
            except Exception as exc:
                failed.append({"sample": name, "target_pdb": pdb,
                               "error": f"{type(exc).__name__}: {exc}"})

    still = _missing_figures(session_id)
    return JSONResponse({
        "session_id": session_id,
        "rendered": rendered,
        "skipped": skipped,
        "failed": failed,
        "n_rendered": sum(len(r["files"]) for r in rendered),
        "still_missing": still["n_missing"],
        "missing_detail": still["missing"][:40],
        "note": ("Figures are schematic: the pose PNG is an orthographic "
                 "projection with the receptor backbone as a trace and the "
                 "ligand as sticks, not a ray-traced render. It is meant to "
                 "show where the ligand sits and which residues are near it, "
                 "and it should be labelled that way in a thesis."),
    })


@app.get("/api/session/{session_id}/figure_status")
async def figure_status(session_id: str):
    """What is missing, before anyone presses Finalize."""
    st = _missing_figures(session_id)
    if not st["session_found"]:
        raise HTTPException(404, f"No session '{session_id}'.")
    return JSONResponse({
        **st,
        "ready_to_finalize": st["n_missing"] == 0,
        "repair_endpoint": f"/api/session/{session_id}/render_figures",
        "n_repairable": sum(1 for m in st["missing"]
                            if m["repairable_without_redocking"]),
        "n_need_redocking": sum(1 for m in st["missing"]
                                if not m["repairable_without_redocking"]),
    })


@app.post("/api/vvuq/finalize-audit")
async def finalize_audit(session_id: str, code: str = ""):
    """Sign off a session, or refuse and name exactly what blocked it."""
    if code != _AUDIT_CODE:
        raise HTTPException(
            403, f"Audit requires code={_AUDIT_CODE}.")

    figs = _missing_figures(session_id)
    if figs["session_found"] and figs["n_missing"]:
        repairable = [m for m in figs["missing"]
                      if m["repairable_without_redocking"]]
        return JSONResponse({
            "audit_status": "BLOCKED_MISSING_FIGURES",
            "session_id": session_id,
            "n_missing": figs["n_missing"],
            "n_repairable_without_redocking": len(repairable),
            "n_requiring_redocking": figs["n_missing"] - len(repairable),
            "missing": figs["missing"][:60],
            "what_to_do": (
                f"{len(repairable)} of {figs['n_missing']} missing figures "
                f"can be regenerated from stored poses: POST "
                f"/api/session/{session_id}/render_figures. The remainder "
                "have no docked pose behind them and need the target "
                "re-docked. The audit is not stamped in the meantime, "
                "because a ledger that certifies figures it cannot find "
                "certifies nothing."
            ),
        }, status_code=409)

    ledger = _run_audit(session_id)
    return JSONResponse(ledger, status_code=200 if
                        ledger["audit_status"].startswith("PASSED") else 409)



@app.get("/bamstu-valid-jimmah")
async def bamstu_valid_jimmah(session_id: str):
    """Convenience alias so the audit can be triggered from a browser or curl."""
    ledger = _run_audit(session_id)
    return JSONResponse(ledger, status_code=200 if
                        ledger["audit_status"].startswith("PASSED") else 409)


# ════════════════════════════════════════════════════════════════════
#  METAL PHYSICS ENGINE
# ════════════════════════════════════════════════════════════════════

try:
    from digital_pathogen import metal_engine as METAL_ENGINE
    METAL_ENGINE_OK = True
except Exception as _me_exc:                       # pragma: no cover
    METAL_ENGINE = None
    METAL_ENGINE_OK = False
    METAL_ENGINE_ERROR = str(_me_exc)


@app.post("/api/metal/analyse")
async def metal_analyse(
    file: UploadFile = File(...),
    vina_score: float = Form(None),
    enable_terms: str = Form(""),
    force: bool = Form(False),
    scf_format: str = Form("orca"),
    scf_file: UploadFile = File(None),
):
    """Run the metal physics engine over one optimised structure.

    Returns the full stage trail, the fingerprint, and dE_metal reported
    three ways: alone, alongside the raw Vina score, and summed.
    """
    if not METAL_ENGINE_OK:
        raise HTTPException(503, f"Metal engine unavailable: {METAL_ENGINE_ERROR}")

    content = (await file.read()).decode("utf-8", errors="ignore")
    try:
        atoms, descriptors_, _fp = parse_xyz(content)
    except chem.ChemistryError as exc:
        raise HTTPException(422, str(exc))

    scf_text = None
    if scf_file is not None:
        try:
            scf_text = (await scf_file.read()).decode("utf-8", errors="ignore")
        except Exception:
            scf_text = None

    terms = {t.strip() for t in enable_terms.split(",") if t.strip()}
    result = METAL_ENGINE.analyse(
        atoms, scf_output=scf_text, scf_format=scf_format,
        vina_score=vina_score, enable_terms=terms, force=bool(force))

    result["graph_descriptors"] = {
        "mw": descriptors_.get("mw"),
        "provenance": descriptors_.get("provenance"),
    }
    result["engine_version"] = METAL_ENGINE.__version__
    result["architecture_note"] = (
        "AutoDock Vina is NOT modified. The engine runs entirely outside it: "
        "it describes the metal centre, generates parameters, and applies a "
        "correction to the score Vina returned. It cannot repair a POSE that "
        "Vina searched with the metal typed as Zn, which is the honest limit "
        "of any post-hoc correction layer."
    )
    return JSONResponse(result)


@app.get("/api/metal/status")
async def metal_status():
    """What the engine can and cannot currently do."""
    if not METAL_ENGINE_OK:
        return JSONResponse({"available": False, "error": METAL_ENGINE_ERROR})
    from digital_pathogen.metal_engine import geometry as _g, detector as _d
    return JSONResponse({
        "available": True,
        "version": METAL_ENGINE.__version__,
        "stages": METAL_ENGINE.STAGES,
        "metals_known": len(_d._METAL_DATA),
        "series_covered": ["transition", "lanthanide", "actinide", "p-block metal"],
        "reference_polyhedra": {str(k): sorted(v) for k, v in _g.IDEAL_POLYHEDRA.items()},
        "geometry_method": (
            "continuous shape measure, exact permutation search for CN <= 7, "
            "greedy+swap approximation for CN >= 8"
        ),
        "oxidation_method": "bond valence sum (Brown & Altermatt 1985)",
        "calibration": {
            "energy_terms_enabled_by_default": [],
            "reason": (
                "Every energy term is UNCALIBRATED and contributes zero until "
                "POST /api/vvuq/calibrate shows Q^2_LOO >= 0.4 with "
                "y-randomisation passing on a congeneric series. An equation "
                "being correct does not make its coefficient correct, and "
                "v3.1 shipped named terms with unfitted coefficients once "
                "already."
            ),
        },
        "cache": METAL_ENGINE.cache.cache_report(),
        "backends": sorted(METAL_ENGINE.vina_bridge.BACKENDS),
    })


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
            "decoy_enrichment": (
                "IMPLEMENTED at POST /api/vvuq/roc_auc. Decoys are SAMPLED "
                "from a library of real compounds, not generated de novo, and "
                "matched on MW (+/-10%), cLogP (+/-1.0) and H-bond counts "
                "(+/-2), with a Morgan Tanimoto ceiling of 0.6 to exclude "
                "analogues. Point BAMSTU_DECOY_LIBRARY at a ZINC or DUD-E "
                "subset for publication work: the built-in 149-compound set "
                "is an offline fallback and will not reach 50 decoys per "
                "active for most ligands."
            ),
            "affinity_correlation": (
                "IMPLEMENTED at POST /api/vvuq/calibrate. Predictors are "
                "gated on sample size BEFORE fitting; Q^2_LOO and a "
                "y-randomisation permutation test are reported alongside R^2; "
                "weights that fail cross-validation come back flagged "
                "usable_for_ranking=false. Still the test that decides "
                "whether the composite can be reported as calibrated."
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


@app.post("/api/chat")
async def chat(data: dict):
    msg = data.get("message", "")
    ctx = data.get("context", "")
    if not groq_client:
        return {"reply":
            "AI chat requires a Groq API key in main.py. "
            "Get a free key at console.groq.com"}
    try:
        r = None
        last_err = None
        for _mc in ["llama-3.3-70b-versatile", "llama-3.1-8b-instant", "gemma2-9b-it"]:
            try:
                r = groq_client.chat.completions.create(
                    model=_mc,
                    messages=[
                        {"role": "system", "content":
                            "You are a scientific assistant for BAMSTU Drug Discovery & Design, "
                            "University of Botswana. Help students interpret molecular docking "
                            "results for organic ligands and transition metal complexes. "
                            f"Be concise, precise, and educational. Session context: {ctx}"},
                        {"role": "user", "content": msg}
                    ],
                    max_tokens=300, temperature=0.4
                )
                break
            except Exception as _me:
                last_err = _me
                if any(x in str(_me).lower() for x in ["decommissioned","not found","deprecated"]):
                    continue
                raise
        if r is None:
            return {"reply": f"All AI models unavailable: {last_err}"}
        return {"reply": r.choices[0].message.content.strip()}
    except Exception as e:
        return {"reply": f"Groq error: {e}"}



# ════════════════════════════════════════════════════════════════════
#  LAYER 3 — EXPERIMENTAL PHENOTYPIC DATA
# ════════════════════════════════════════════════════════════════════

@app.post("/api/experimental/mic")
async def add_mic(payload: dict = Body(...)):
    """Store one measured MIC and return it with pMIC.

    Molecular weight is required, not defaulted: pMIC is what the
    calibration regresses on, and a guessed MW moves every training point.
    """
    records = payload.get("records") or [payload]
    stored, errors = [], []
    for rec in records:
        try:
            r = experimental.STORE.add(rec)
            stored.append({
                "compound_id": r.compound_id, "strain": r.strain,
                "mic_ug_per_ml": r.mic_ug_per_ml,
                "mic_um": round(r.mic_um(), 3),
                "pmic": round(r.pmic(), 3),
                "pmic_uncertainty": r.pmic_uncertainty(),
                "quantities": r.quantities(),
            })
        except experimental.ExperimentalDataError as exc:
            errors.append({"record": rec.get("compound_id", "unnamed"),
                           "error": str(exc)})
    return JSONResponse({"stored": stored, "rejected": errors,
                         "store": experimental.STORE.summary()},
                        status_code=200 if stored else 422)


@app.get("/api/experimental/mic")
async def list_mic(strain: str = None):
    return JSONResponse({
        "summary": experimental.STORE.summary(),
        "training_set": experimental.STORE.training_set(strain),
        "records": [r.to_dict() for r in experimental.STORE.all()],
    })


@app.post("/api/experimental/ki_from_mic")
async def refuse_ki_from_mic():
    """Exists so the refusal is discoverable, and returns 422."""
    try:
        experimental.mic_to_ki()
    except prov.ProvenanceViolation as exc:
        return JSONResponse({"refused": True, "reason": str(exc)},
                            status_code=422)
    return JSONResponse({"refused": True}, status_code=422)


# ════════════════════════════════════════════════════════════════════
#  PHENOTYPIC CALIBRATION AND APPLICABILITY DOMAIN
# ════════════════════════════════════════════════════════════════════

CALIBRATED_MODELS: dict = {}


@app.post("/api/calibration/phenotypic")
async def calibrate_phenotypic(payload: dict = Body(...)):
    """Fit the physics-informed score against measured pMIC.

    `rows` supplies, per compound: compound_id, pmic (or mic_ug_per_ml plus
    molecular_weight), bamstu_composite_score and
    vina_affinity_kcal_per_mol. Nothing is docked here -- the descriptors
    come from runs already performed, so the statistics are reproducible
    without a Vina binary.
    """
    rows = payload.get("rows") or []
    strain = payload.get("strain", "")
    prepared = []
    for r in rows:
        row = dict(r)
        if row.get("pmic") is None and row.get("mic_ug_per_ml") and \
                row.get("molecular_weight"):
            rec = experimental.build_record({
                "compound_id": row.get("compound_id", "unnamed"),
                "strain": strain or row.get("strain", "unspecified"),
                "species": row.get("species", "unspecified"),
                "mic_ug_per_ml": row["mic_ug_per_ml"],
                "molecular_weight": row["molecular_weight"],
            })
            row["pmic"] = rec.pmic()
        prepared.append(row)

    result = calibration.calibrate_phenotypic(
        prepared, strain=strain,
        feature_names=payload.get("feature_names"),
        y_randomization_trials=int(payload.get("y_randomization_trials", 200)),
        n_bootstrap=int(payload.get("n_bootstrap", 500)),
        seed=int(payload.get("seed", 42)))

    model = result.pop("_model_object", None)
    key = payload.get("model_id") or f"{strain or 'default'}"
    if model is not None:
        CALIBRATED_MODELS[key] = model
        result["model_id"] = key
    return JSONResponse(result,
                        status_code=200 if result["status"] == "COMPLETED"
                        else 422)


@app.post("/api/calibration/predict")
async def predict_pmic(payload: dict = Body(...)):
    """Predict pMIC for one compound. Withheld unless the gate cleared."""
    key = payload.get("model_id") or "default"
    model = CALIBRATED_MODELS.get(key)
    if model is None:
        return JSONResponse({
            "status": "NO_MODEL",
            "reason": (f"No calibrated model '{key}'. Run "
                       "POST /api/calibration/phenotypic first. Predicting "
                       "before calibration is what a scoring function does "
                       "when nobody has checked whether it predicts."),
        }, status_code=404)
    out = model.predict(payload.get("descriptors") or {})
    if out.get("status") == "PREDICTED" and payload.get("molecular_weight"):
        out["concentration"] = model.predicted_mic(
            out["predicted_pmic"], float(payload["molecular_weight"]))
    return JSONResponse(out,
                        status_code=200 if out["status"] == "PREDICTED" else 409)


@app.get("/api/calibration/status")
async def calibration_status():
    return JSONResponse({
        "models": {k: {"released": m.released,
                       "strain": m.strain,
                       "features": m.feature_names,
                       "statistics": {s: m.statistics.get(s) for s in
                                      ("n", "r2", "q2_loo", "mae", "rmse")},
                       "gate_reason": m.gate_reason}
                   for k, m in CALIBRATED_MODELS.items()},
        "gate": {"q2_loo_threshold": calibration.Q2_THRESHOLD,
                 "y_randomization_alpha": calibration.YRAND_ALPHA},
    })


# ════════════════════════════════════════════════════════════════════
#  PROVENANCE
# ════════════════════════════════════════════════════════════════════

@app.get("/api/provenance/categories")
async def provenance_categories():
    return JSONResponse({
        "schema": prov.SCHEMA_VERSION,
        "categories": {c: prov.CATEGORY_DEFINITION[c] for c in prov.CATEGORIES},
        "forbidden_conversions": [
            {"from": a, "to": b, "reason": r}
            for (a, b), r in prov.FORBIDDEN_CONVERSIONS.items()
        ],
        "rule": ("Every reported quantity carries exactly one category for "
                 "its whole life. Nothing is relabelled and nothing "
                 "inherits units it cannot defend."),
    })


@app.post("/api/provenance/audit")
async def provenance_audit(payload: dict = Body(...)):
    return JSONResponse(prov.audit_payload(payload))


# ════════════════════════════════════════════════════════════════════
#  LAYERED REPORT
# ════════════════════════════════════════════════════════════════════

@app.post("/api/report/layered")
async def layered_report(payload: dict = Body(...)):
    """Assemble the four-section report for one hit.

    Accepts a hit dict as returned by the sweep, optionally a metal-engine
    result, and optionally a compound id to pull MIC records and a
    calibrated model for Layers 3 and 4.
    """
    hit = payload.get("hit")
    if not isinstance(hit, dict):
        raise HTTPException(422, "A 'hit' object is required.")
    compound = payload.get("compound_id", "")
    records = experimental.STORE.for_compound(compound) if compound else []
    model = CALIBRATED_MODELS.get(payload.get("model_id") or "default")
    cal = payload.get("calibration")
    if cal is None and model is not None:
        cal = {"released": model.released,
               "calibration_status": ("CALIBRATED AND RELEASED"
                                      if model.released else "NOT CALIBRATED"),
               "gate_reason": model.gate_reason,
               "statistics": model.statistics,
               "applicability_domain": model.ad.to_dict(),
               "williams_plot": model.ad.williams_plot()}
    prediction = None
    if model is not None and payload.get("descriptors"):
        prediction = model.predict(payload["descriptors"])
    return JSONResponse(layers.assemble_report(
        hit, engine=payload.get("metal_engine"), experimental=records,
        calibration=cal, prediction=prediction, compound_id=compound))


# ════════════════════════════════════════════════════════════
if __name__ == "__main__":
    ai_status = ("Groq connected" if groq_client else
                 "No key in BAMSTU_GROQ_API_KEY — narrative layer runs from a template that makes no mechanistic claims")
    print("=" * 62)
    print("  BAMSTU Digital Pathogen Pro — Backend Server v3.4")
    print("  Dr. James P.T. Matshwele | BAMSTU Laboratories | UB | 2026")
    print(f"  AI Status : {ai_status}")
    print("  Server    : http://127.0.0.1:8080")
    print("  Open Chrome → http://127.0.0.1:8080")
    print("  Keep this window open!")
    print("=" * 62)
    uvicorn.run("main:app", host="127.0.0.1", port=8080, reload=False)
