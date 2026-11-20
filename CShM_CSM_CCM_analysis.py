# -*- coding: utf-8 -*-
import sys
import os
import json
import html
import subprocess
import numpy as np

# ✅ SOLUTION : Force l'utilisation du Python de votre conda
PYTHON_EXTERNE = "/Users/yrousselin/CCDC/ccdc-software/csd-python-api/miniconda/bin/python"
MAX_CN = 12

# Recherche des voisins dans le cristal (MOF, polymères, métaux sur position spéciale)
BOND_TOL = 0.45      # tolérance (Å) ajoutée à la somme des rayons covalents
R_SEARCH = 3.3       # rayon max (Å) de recherche autour du métal
MIN_DIST = 0.6       # en dessous : atome superposé (désordre) -> ignoré

# Rayons covalents (Cordero 2008) pour les cas où l'API ne les fournit pas
COV_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20, "Br": 1.20, "Te": 1.38, "I": 1.39,
    "Li": 1.28, "Na": 1.66, "K": 2.03, "Rb": 2.20, "Cs": 2.44,
    "Be": 0.96, "Mg": 1.41, "Ca": 1.76, "Sr": 1.95, "Ba": 2.15, "Al": 1.21, "Ga": 1.22,
    "In": 1.42, "Sn": 1.39, "Pb": 1.46, "Bi": 1.48,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Sm": 1.98,
    "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90,
    "Yb": 1.87, "Lu": 1.87, "Th": 2.06, "U": 1.96,
}

# ------------------------------------------------------------------ COEUR DE CALCUL
# Chaque fonction accepte un filtre `only` (liste de codes) : cela permet de relancer
# les mesures une par une si cosymlib plante (segfault du code Fortran) sur l'une d'elles.
SYMMETRY_LABELS = {
    2: ['C1', 'C2', 'Cs', 'Ci'],
    3: ['C1', 'C2', 'C3', 'Cs', 'Ci'],
    4: ['C1', 'C2', 'C4', 'S4', 'Cs', 'Ci'],
    5: ['C1', 'C5', 'S4', 'Cs', 'Ci'],
    6: ['C1', 'C2', 'C3', 'C4', 'C6', 'S4', 'S6', 'Cs', 'Ci'],
    7: ['C1', 'C7', 'S4', 'Cs', 'Ci'],
    8: ['C1', 'C2', 'C4', 'C8', 'S4', 'S8', 'Cs', 'Ci'],
    9: ['C1', 'C3', 'C9', 'S4', 'S6', 'Cs', 'Ci'],
    10: ['C1', 'C2', 'C5', 'C10', 'S4', 'S10', 'Cs', 'Ci'],
    11: ['C1', 'C11', 'S4', 'Cs', 'Ci'],
    12: ['C1', 'C2', 'C3', 'C4', 'C6', 'C12', 'S4', 'S6', 'S8', 'S12', 'Cs', 'Ci']
}

def _cn(symbols):
    return len(symbols) - 1

def shape_codes(cn):
    from cosymlib.shape.tools import shape_structure_references as refs
    return [r[0] for r in refs.get(f"{cn} Vertices", [])]

def symmetry_codes(cn):
    return SYMMETRY_LABELS.get(cn, ['C1', 'C2', 'Cs', 'Ci'])

def chirality_codes(cn):
    return [g for g in ("Cs", "Ci", "S4", "S8", "S12") if g in ("Cs", "Ci") or int(g[1:]) <= cn]

def shape_measures(symbols, positions, central_atom=1, only=None):
    from cosymlib.molecule.geometry import Geometry
    from cosymlib.shape.tools import shape_structure_references as refs
    key = f"{_cn(symbols)} Vertices"
    if key not in refs:
        return []
    geom = Geometry(symbols=symbols, positions=positions)
    rows = []
    for code, _num, sym, name in refs[key]:
        if only is not None and code not in only:
            continue
        try:
            cshm = geom.get_shape_measure(code, central_atom=central_atom)
            rows.append((code, name, sym, float(cshm)))
        except Exception:
            continue
    return sorted(rows, key=lambda r: r[3])

def symmetry_measures(symbols, positions, central_atom=1, only=None):
    from cosymlib.molecule.geometry import Geometry
    geom = Geometry(symbols=symbols, positions=positions)
    rows = []
    for label in symmetry_codes(_cn(symbols)):
        if only is not None and label not in only:
            continue
        try:
            csm = geom.get_symmetry_measure(label, central_atom=central_atom)
            rows.append((label, label, label, float(csm)))
        except Exception:
            continue
    return sorted(rows, key=lambda r: r[3])

def _ccm_summary(detail):
    out = {"detail": detail}
    for key in ("geom", "typed"):
        vals = [(r[key], r["group"]) for r in detail if r.get(key) is not None]
        out["ccm_" + key] = min(vals) if vals else None   # (valeur, groupe)
    return out

def chirality_measures(symbols, positions, central_atom=1, only=None):
    """Continuous Chirality Measure (Avnir) du site de coordination.

    CCM = min des CSM calculées pour les groupes ponctuels IMPROPRES indépendants :
      Cs (miroir), Ci (inversion), S4, S8, S12 ... (S_4k ; les S_n d'ordre impair ou 4k+2 contiennent
      déjà un miroir ou un centre d'inversion, donc ne peuvent pas donner un minimum plus bas).
    CCM = 0  -> site achiral ;  CCM > 0 -> chiral (d'autant plus que la valeur est grande).
    Deux versions : 'geom' (tous les ligands traités comme identiques : chiralité du polyèdre seul)
    et 'typed' (la nature des atomes donneurs compte : N != O, etc.).
    """
    from cosymlib.molecule.geometry import Geometry
    cn = _cn(symbols)
    typed = Geometry(symbols=symbols, positions=positions)
    geom_only = Geometry(symbols=[symbols[0]] + ["N"] * cn, positions=positions)
    detail = []
    for g in chirality_codes(cn):
        if only is not None and g not in only:
            continue
        row = {"group": g}
        for key, geo in (("geom", geom_only), ("typed", typed)):
            try:
                row[key] = max(0.0, float(geo.get_symmetry_measure(g, central_atom=central_atom)))  # évite -0.000
            except Exception:
                row[key] = None
        detail.append(row)
    return _ccm_summary(detail)

TASKS = {
    "shape": (shape_measures, shape_codes),
    "symmetry": (symmetry_measures, symmetry_codes),
    "chirality": (chirality_measures, chirality_codes),
}

def _worker():
    import traceback
    if len(sys.argv) < 3:
        json.dump({"error": "Arguments manquants"}, sys.stdout)
        sys.exit(1)
    try:
        d = json.loads(sys.argv[2])
        task = d.get("task", "all")
        cn = _cn(d["symbols"])
        if task == "codes":
            result = {k: TASKS[k][1](cn) for k in TASKS}
        elif task in TASKS:
            result = TASKS[task][0](d["symbols"], d["positions"], d.get("central_atom", 1), d.get("only"))
        else:
            result = {k: TASKS[k][0](d["symbols"], d["positions"], d.get("central_atom", 1)) for k in TASKS}
        json.dump({"result": result}, sys.stdout)
    except ImportError:
        info = ("\nPython utilisé : " + sys.executable
                + "\nsys.path :\n  " + "\n  ".join(sys.path)
                + "\n-> installez cosymlib avec : " + sys.executable + " -m pip install cosymlib")
        json.dump({"error": traceback.format_exc() + info}, sys.stdout)
        sys.exit(1)
    except Exception:
        json.dump({"error": traceback.format_exc()}, sys.stdout)
        sys.exit(1)

if len(sys.argv) > 1 and sys.argv[1] == "--cosymlib-worker":
    _worker()
    sys.exit(0)

def _run_external(payload):
    """Lance le worker dans le Python conda. Renvoie (résultat, message_erreur)."""
    SCRIPT_PATH = os.path.abspath(__file__)
    # Environnement propre : Mercury exporte des variables Python (PYTHONHOME, PYTHONPATH...)
    # qui, héritées par le sous-processus, peuvent lui masquer les paquets installés (cosymlib).
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONHOME", "PYTHONPATH", "PYTHONNOUSERSITE", "PYTHONSTARTUP", "PYTHONEXECUTABLE")}
    p = subprocess.Popen(
        [PYTHON_EXTERNE, SCRIPT_PATH, "--cosymlib-worker", json.dumps(payload)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        cwd=os.path.dirname(SCRIPT_PATH)
    )
    stdout, stderr = p.communicate()
    data = None
    try:
        data = json.loads(stdout) if stdout.strip() else None
    except json.JSONDecodeError:
        pass
    if p.returncode == 0 and data and "result" in data:
        return data["result"], None
    msg = f"code retour {p.returncode}"
    if p.returncode < 0:
        msg += f" (processus tué par le signal {-p.returncode} : plantage de cosymlib)"
    if data and data.get("error"):
        msg += "\n" + data["error"][-1500:]
    elif stdout.strip():
        msg += "\nstdout : " + stdout[-800:]
    if stderr.strip():
        msg += "\nstderr : " + stderr[-1500:]
    return None, msg

def compute(symbols, positions):
    try:
        import cosymlib  # noqa: F401  (cosymlib disponible dans le Python de Mercury)
        return {k: TASKS[k][0](symbols, positions) for k in TASKS}
    except ImportError:
        pass

    base = {"symbols": symbols, "positions": positions, "central_atom": 1}
    codes, err = _run_external(dict(base, task="codes"))
    if codes is None:
        raise RuntimeError(f"Échec du Python externe ({PYTHON_EXTERNE}) :\n{err}")

    out, warnings = {}, []
    for task in TASKS:
        res, err = _run_external(dict(base, task=task))
        if res is None:
            # Une mesure fait planter cosymlib : on relance code par code pour isoler le coupable
            warnings.append(f"{task} : échec global, relance mesure par mesure.\n{err}")
            parts = []
            for code in codes[task]:
                r, e = _run_external(dict(base, task=task, only=[code]))
                if r is None:
                    warnings.append(f"{task} / {code} : ignoré ({e.splitlines()[0]})")
                    continue
                parts.append(r)
            if task == "chirality":
                res = _ccm_summary([row for r in parts for row in r["detail"]])
            else:
                res = sorted((row for r in parts for row in r), key=lambda x: x[3])
        out[task] = res
    out["warnings"] = warnings
    return out

# ------------------------------------------------------------------ RECHERCHE DES SITES
def _xyz(atom):
    c = atom.coordinates
    return np.array([float(c.x), float(c.y), float(c.z)])

def _cov_radius(atom):
    r = getattr(atom, "covalent_radius", None)
    try:
        r = float(r)
        if r > 0:
            return r
    except (TypeError, ValueError):
        pass
    return COV_RADII.get(atom.atomic_symbol, 1.5)

def crystal_neighbours(crystal, metal):
    """Sphère de coordination complète du métal, reconstruite par symétrie.

    Indispensable pour les MOF / polymères de coordination et les métaux sur
    position spéciale : l'unité asymétrique ne contient alors qu'une partie des
    ligands (ex. ZIF-8 : Zn1 sur site -4, un seul N1 indépendant -> CN = 1).
    On remplit une petite boîte autour du métal (atomes individuels) et on
    garde les atomes à distance de liaison.
    """
    if metal.coordinates is None:
        return []
    p_m = _xyz(metal)
    a, b, c = (float(x) for x in crystal.cell_lengths)
    f = getattr(metal, "fractional_coordinates", None)
    if f is not None:
        fx, fy, fz = float(f.x), float(f.y), float(f.z)
        # marge un peu large pour les mailles très obliques
        dx, dy, dz = 1.5 * R_SEARCH / a, 1.5 * R_SEARCH / b, 1.5 * R_SEARCH / c
        box = ((fx - dx, fy - dy, fz - dz), (fx + dx, fy + dy, fz + dz))
    else:
        box = ((-1, -1, -1), (2, 2, 2))
    shell = crystal.packing(box, inclusion="OnlyAtomsIncluded")

    r_m = _cov_radius(metal)
    found = []
    for at in shell.atoms:
        if at.coordinates is None or at.atomic_symbol == "H" or at.is_metal:
            continue
        d = float(np.linalg.norm(_xyz(at) - p_m))
        if d < MIN_DIST or d > R_SEARCH:
            continue
        if d <= r_m + _cov_radius(at) + BOND_TOL:
            # élimine les doublons éventuels (même position générée deux fois)
            if all(np.linalg.norm(_xyz(at) - _xyz(o)) > 0.1 for o, _ in found):
                found.append((at, d))
    found.sort(key=lambda t: t[1])
    return [at for at, _ in found]

def build_sites(entry, selected_atoms):
    sel = list(selected_atoms)
    if sel:
        metals = [a for a in sel if a.is_metal]
        if len(metals) != 1:
            raise ValueError(f"Sélectionnez exactement un métal ({len(metals)} trouvé(s)).")
        m = metals[0]
        return [(m.label, m, [a for a in sel if a is not m and a.index != m.index])]

    crystal = None
    try:
        crystal = entry.crystal
    except Exception:
        crystal = None

    sites = []
    if crystal is not None:
        # Une seule entrée par métal indépendant (unité asymétrique)
        try:
            asu = crystal.asymmetric_unit_molecule
        except Exception:
            asu = entry.molecule
        seen = set()
        for m in (a for a in asu.atoms if a.is_metal):
            if m.label in seen:
                continue
            seen.add(m.label)
            nb = crystal_neighbours(crystal, m)
            if nb:
                sites.append((m.label, m, nb))
    if not sites:
        # Repli : liaisons de la molécule affichée (complexes discrets sans maille)
        for m in (a for a in entry.molecule.atoms if a.is_metal):
            nb = list(m.neighbours)
            if nb:
                sites.append((m.label, m, nb))
    if not sites:
        raise ValueError("Aucun métal lié trouvé.")
    return sites

def xyz_lists(central, ligands):
    atoms = [central] + list(ligands)
    symbols = [a.atomic_symbol for a in atoms]
    positions = np.array([_xyz(a) for a in atoms])
    positions -= positions[0]
    return symbols, positions.tolist()

def html_report(identifier, results, errors):
    h = [
        "<html><head><meta charset='utf-8'><style>",
        "body{font-family:sans-serif; margin:20px;}",
        "h2{color:#2c3e50;} h3{color:#3498db; margin-top:1em;}",
        "table{border-collapse:collapse; margin:1em 0; width:100%;}",
        "td,th{border:1px solid #ddd; padding:8px; text-align:left;}",
        "th{background:#f2f2f2;}",
        "td.n{text-align:right; font-family:monospace;}",
        "tr.best{background:#d8f0d0; font-weight:bold;}",
        ".error{color:#e74c3c;}",
        "</style></head><body>",
        f"<h2>Analyse structurale (cosymlib) - {identifier}</h2>",
        "<p><small>CShM = 0 : géométrie idéale | CSM = 0 : symétrie idéale | CCM = 0 : achiral, CCM &gt; 0 : chiral</small></p>"
    ]
    for label, cn, data, env in results:
        h.append(f"<h3>{label} - Nombre de coordination {cn}</h3>")
        if env:
            h.append(f"<p><small>Sphère : {env}</small></p>")
        for w in data.get("warnings", []):
            h.append(f"<pre class='error'>⚠️ {html.escape(w)}</pre>")
        if data["shape"]:
            h.append("<h4>🔺 Continuous Shape Measures (CShM)</h4><table><tr><th>Code</th><th>Polyèdre</th><th>Symétrie</th><th>CShM</th></tr>")
            for i, (code, name, sym, val) in enumerate(data["shape"]):
                h.append(f"<tr{' class=best' if i == 0 else ''}><td>{code}</td><td>{name}</td><td>{sym}</td><td class='n'>{val:.3f}</td></tr>")
            h.append("</table>")
        if data["symmetry"]:
            h.append("<h4>🔷 Continuous Symmetry Measures (CSM)</h4><table><tr><th>Code</th><th>Symétrie</th><th>Groupe</th><th>CSM</th></tr>")
            for i, (code, name, sym, val) in enumerate(data["symmetry"]):
                h.append(f"<tr{' class=best' if i == 0 else ''}><td>{code}</td><td>{name}</td><td>{sym}</td><td class='n'>{val:.3f}</td></tr>")
            h.append("</table>")
        ch = data.get("chirality")
        if ch and ch.get("detail"):
            h.append("<h4>🌀 Continuous Chirality Measure (CCM)</h4>")
            for key, titre in (("ccm_geom", "polyèdre seul (ligands identiques)"), ("ccm_typed", "avec la nature des atomes")):
                best = ch.get(key)
                if best:
                    etat = "achiral" if best[0] < 0.01 else "chiral"
                    h.append(f"<p><b>CCM {titre} = {best[0]:.3f}</b> (groupe le plus proche : {best[1]}) &rarr; {etat}</p>")
            h.append("<table><tr><th>Groupe impropre</th><th>CSM (polyèdre seul)</th><th>CSM (avec atomes)</th></tr>")
            fmt = lambda v: "-" if v is None else f"{v:.3f}"
            for r in ch["detail"]:
                h.append(f"<tr><td>{r['group']}</td><td class='n'>{fmt(r['geom'])}</td><td class='n'>{fmt(r['typed'])}</td></tr>")
            h.append("</table>")
    for e in errors:
        h.append(f"<pre class='error'>⚠️ {html.escape(e)}</pre>")
    h.append("</body></html>")
    return "\n".join(h)

def main():
    from ccdc.utilities import ApplicationInterface
    helper = ApplicationInterface()
    results, errors = [], []
    try:
        entry = helper.current_entry
        if entry is None or not entry.has_3d_structure:
            raise ValueError("Cette structure n'a pas de coordonnées 3D.")
        sites = build_sites(
            entry,
            helper.selected_atoms if helper.options.get("selected_atoms") else []
        )
        for label, central, ligands in sites:
            cn = len(ligands)
            if cn < 2 or cn > MAX_CN:
                errors.append(f"{label} : CN {cn} ignoré (2-{MAX_CN} autorisé).")
                continue
            sym, pos = xyz_lists(central, ligands)
            env = ", ".join(f"{a.atomic_symbol} {np.linalg.norm(p):.3f} Å"
                            for a, p in zip(ligands, np.array(pos)[1:]))
            try:
                analysis = compute(sym, pos)
                results.append((label, cn, analysis, env))
            except Exception as exc:
                errors.append(f"{label} : {str(exc)}")
    except Exception as exc:
        errors.append(str(exc))
    with open(helper.output_html_file, "w", encoding="utf-8") as f:
        f.write(html_report(getattr(helper, "identifier", "") or "", results, errors))

if __name__ == "__main__":
    main()
