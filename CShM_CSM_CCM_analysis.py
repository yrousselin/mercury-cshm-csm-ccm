# -*- coding: utf-8 -*-
import sys
import os
import json
import subprocess
import numpy as np

# ✅ SOLUTION : Force l'utilisation du Python de votre conda
PYTHON_EXTERNE = "/Users/yrousselin/CCDC/ccdc-software/csd-python-api/miniconda/bin/python"
MAX_CN = 12

# ------------------------------------------------------------------ COEUR DE CALCUL
def shape_measures(symbols, positions, central_atom=1):
    from cosymlib.molecule.geometry import Geometry
    from cosymlib.shape.tools import shape_structure_references as refs
    cn = len(symbols) - 1
    key = f"{cn} Vertices"
    if key not in refs:
        return []
    geom = Geometry(symbols=symbols, positions=positions)
    rows = []
    for code, _num, sym, name in refs[key]:
        try:
            cshm = geom.get_shape_measure(code, central_atom=central_atom)
            rows.append((code, name, sym, float(cshm)))
        except Exception:
            continue
    return sorted(rows, key=lambda r: r[3])

def symmetry_measures(symbols, positions, central_atom=1):
    from cosymlib.molecule.geometry import Geometry
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
    cn = len(symbols) - 1
    labels = SYMMETRY_LABELS.get(cn, ['C1', 'C2', 'Cs', 'Ci'])
    geom = Geometry(symbols=symbols, positions=positions)
    rows = []
    for label in labels:
        try:
            csm = geom.get_symmetry_measure(label, central_atom=central_atom)
            rows.append((label, label, label, float(csm)))
        except Exception:
            continue
    return sorted(rows, key=lambda r: r[3])

def chirality_measures(symbols, positions, central_atom=1):
    """Continuous Chirality Measure (Avnir) du site de coordination.

    CCM = min des CSM calculées pour les groupes ponctuels IMPROPRES indépendants :
      Cs (miroir), Ci (inversion), S4, S8, S12 ... (S_4k ; les S_n d'ordre impair ou 4k+2 contiennent
      déjà un miroir ou un centre d'inversion, donc ne peuvent pas donner un minimum plus bas).
    CCM = 0  -> site achiral ;  CCM > 0 -> chiral (d'autant plus que la valeur est grande).
    Deux versions : 'geom' (tous les ligands traités comme identiques : chiralité du polyèdre seul)
    et 'typed' (la nature des atomes donneurs compte : N != O, etc.).
    """
    from cosymlib.molecule.geometry import Geometry
    cn = len(symbols) - 1
    groups = [g for g in ("Cs", "Ci", "S4", "S8", "S12") if g in ("Cs", "Ci") or int(g[1:]) <= cn]
    typed = Geometry(symbols=symbols, positions=positions)
    geom_only = Geometry(symbols=[symbols[0]] + ["N"] * cn, positions=positions)
    detail = []
    for g in groups:
        row = {"group": g}
        for key, geo in (("geom", geom_only), ("typed", typed)):
            try:
                row[key] = max(0.0, float(geo.get_symmetry_measure(g, central_atom=central_atom)))  # évite -0.000
            except Exception:
                row[key] = None
        detail.append(row)
    out = {"detail": detail}
    for key in ("geom", "typed"):
        vals = [(r[key], r["group"]) for r in detail if r[key] is not None]
        out["ccm_" + key] = min(vals) if vals else None   # (valeur, groupe)
    return out

def _worker():
    if len(sys.argv) < 3:
        json.dump({"error": "Arguments manquants"}, sys.stdout)
        sys.exit(1)
    try:
        d = json.loads(sys.argv[2])
        result = {
            "shape": shape_measures(d["symbols"], d["positions"], d.get("central_atom", 1)),
            "symmetry": symmetry_measures(d["symbols"], d["positions"], d.get("central_atom", 1)),
            "chirality": chirality_measures(d["symbols"], d["positions"], d.get("central_atom", 1))
        }
        json.dump(result, sys.stdout)
    except Exception as e:
        json.dump({"error": str(e)}, sys.stdout)
        sys.exit(1)

if len(sys.argv) > 1 and sys.argv[1] == "--cosymlib-worker":
    _worker()
    sys.exit(0)

def compute(symbols, positions):
    try:
        import cosymlib
        return {
            "shape": shape_measures(symbols, positions),
            "symmetry": symmetry_measures(symbols, positions),
            "chirality": chirality_measures(symbols, positions)
        }
    except ImportError:
        SCRIPT_PATH = os.path.abspath(__file__)
        json_data = json.dumps({"symbols": symbols, "positions": positions, "central_atom": 1})
        p = subprocess.Popen(
            [PYTHON_EXTERNE, SCRIPT_PATH, "--cosymlib-worker", json_data],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise RuntimeError(f"Échec du Python externe :\n{stderr[:500]}")
        try:
            return json.loads(stdout)
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Erreur JSON : {e}\nSortie : {stdout[:500]}")

def build_sites(molecule, selected_atoms):
    sel = list(selected_atoms)
    if sel:
        metals = [a for a in sel if a.is_metal]
        if len(metals) != 1:
            raise ValueError(f"Sélectionnez exactement un métal ({len(metals)} trouvé(s)).")
        m = metals[0]
        return [(m.label, m, [a for a in sel if a is not m and a.index != m.index])]
    sites = []
    for m in (a for a in molecule.atoms if a.is_metal):
        nb = list(m.neighbours)
        if nb:
            sites.append((m.label, m, nb))
    if not sites:
        raise ValueError("Aucun métal lié trouvé.")
    return sites

def xyz_lists(central, ligands):
    atoms = [central] + list(ligands)
    symbols = [a.atomic_symbol for a in atoms]
    positions = np.array([[float(a.coordinates.x), float(a.coordinates.y), float(a.coordinates.z)]
                          for a in atoms])
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
    for label, cn, data in results:
        h.append(f"<h3>{label} - Nombre de coordination {cn}</h3>")
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
        h.append(f"<p class='error'>⚠️ {e}</p>")
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
            entry.molecule,
            helper.selected_atoms if helper.options.get("selected_atoms") else []
        )
        for label, central, ligands in sites:
            cn = len(ligands)
            if cn < 2 or cn > MAX_CN:
                errors.append(f"{label} : CN {cn} ignoré (2-{MAX_CN} autorisé).")
                continue
            sym, pos = xyz_lists(central, ligands)
            try:
                analysis = compute(sym, pos)
                results.append((label, cn, analysis))
            except Exception as exc:
                errors.append(f"{label} : {str(exc)}")
    except Exception as exc:
        errors.append(str(exc))
    with open(helper.output_html_file, "w", encoding="utf-8") as f:
        f.write(html_report(getattr(helper, "identifier", "") or "", results, errors))

if __name__ == "__main__":
    main()