# Coordination geometry analysis (CShM, CSM, CCM) for Mercury

*[Version française : LISEZMOI.md](LISEZMOI.md)*

A Python script for **Mercury** (CCDC) that describes the shape, symmetry and chirality of the **coordination sphere of metal atoms** in the structure you have open. The calculations are done by the **cosymlib** library and the results appear as a report in Mercury.

## What it does

For each metal centre with 2 to 12 donor atoms, the report gives:

- **Continuous Shape Measures (CShM):** the distance between the real coordination polyhedron and each ideal reference polyhedron for that coordination number (octahedron, trigonal prism, square antiprism…). The reference with the lowest value, highlighted in green, is the closest geometry. **CShM = 0** means a perfect shape.
- **Continuous Symmetry Measures (CSM):** how far the coordination sphere is from point groups such as C₂, C₃, Cₛ or Cᵢ. **CSM = 0** means perfect symmetry.
- **Continuous Chirality Measure (CCM):** the smallest CSM over the improper point groups (Cₛ, Cᵢ, S₄, S₈, S₁₂). **CCM = 0** means achiral; **CCM > 0** means chiral. It is given twice:
  - *polyhedron only*: all donor atoms treated as identical (chirality of the shape);
  - *with atom types*: N, O, etc. treated as different atoms.

The report is in French.

**Which atoms are used**

- If you **select atoms** in Mercury before running the script (exactly one metal plus its donor atoms), only that site is analysed.
- Otherwise, every symmetry-independent metal is analysed. Its full coordination sphere is **rebuilt from the unit cell and the symmetry operators**, so the script also works for **MOFs, coordination polymers and metals on special positions** (for example Zn in ZIF-8, where the asymmetric unit holds only one of the four N atoms). Donor atoms are the non-hydrogen, non-metal atoms closer than the sum of the covalent radii + 0.45 Å. If the structure has no unit cell, the bonds defined in Mercury are used instead.
- The report lists the donor atoms and their distances for each metal, so you can check the coordination sphere at a glance.

## Requirements

| What | Why | Where to get it |
|---|---|---|
| **Mercury with the CSD Python API** | runs the script | [CCDC](https://www.ccdc.cam.ac.uk/) (requires a CSD licence) |
| **cosymlib** | computes CShM, CSM and CCM | [cosymlib documentation](https://cosymlib.readthedocs.io/en/latest/) |
| numpy | geometry | included in the CSD Python API |

## Installation

### Step 1: install cosymlib

The script needs a Python in which `import cosymlib` works. When cosymlib is not available inside Mercury, the script calls the Python given by `PYTHON_EXTERNE` (step 2). Two options:

#### Option A (simplest): install cosymlib in the Python of the CSD

Call pip **through the full path** of the CSD Python, so that cosymlib goes into that Python whatever conda environment is active in your terminal:

```
/Users/<name>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python -m pip install cosymlib
/Users/<name>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python -c "import cosymlib; print('cosymlib OK')"
```

(Windows: `"C:\Program Files\CCDC\ccdc-software\csd-python-api\miniconda\python.exe" -m pip install cosymlib`, adapting the path to your installation.)

Precompiled versions of cosymlib exist for Python 3.10 to 3.12 on Linux, Windows and Apple Silicon Macs. On an Intel Mac, pip has to compile cosymlib and needs a Fortran compiler.

> ⚠️ **A CSD update reinstalls this Python and removes cosymlib.** After each update, run the `pip install` command again.

#### Option B: install cosymlib in its own environment

As the [cosymlib installation page](https://cosymlib.readthedocs.io/en/latest/install.html) recommends; this environment survives CSD updates.

1. Install **Miniconda** or **Anaconda** if you do not have it yet.
2. Open a terminal:
   - macOS: the **Terminal** application;
   - Windows: **Anaconda Prompt** from the Start menu.
3. Type these commands one by one, pressing Enter after each:

   ```
   conda create -n cosym_env python=3.11 -y
   conda activate cosym_env
   pip install cosymlib
   python -c "import cosymlib; print('cosymlib OK')"
   ```

   The last line should print `cosymlib OK`. If installation fails, see the platform notes on the cosymlib installation page (compilers on Windows, Apple Silicon Macs).
4. Still in the same terminal, note the full path of this Python:
   - macOS: `which python`, which gives something like `/Users/<name>/miniconda3/envs/cosym_env/bin/python`
   - Windows: `where python`, which gives something like `C:\Users\<name>\miniconda3\envs\cosym_env\python.exe`

### Step 2: point the script to that Python

Open `CShM_CSM_CCM_analysis.py` in a text editor (TextEdit in plain-text mode, Notepad, or any code editor). Replace the path of the `PYTHON_EXTERNE = ...` line (**line 10**) with the Python from step 1:

```python
PYTHON_EXTERNE = "/Users/<name>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python"   # option A, macOS
PYTHON_EXTERNE = "/Users/<name>/miniconda3/envs/cosym_env/bin/python"        # option B, macOS
PYTHON_EXTERNE = "C:/Users/<name>/miniconda3/envs/cosym_env/python.exe"      # option B, Windows
```

On Windows, write the path with forward slashes `/`, as above.

### Step 3: add the script to Mercury

1. Put the script in a folder that you keep for your Mercury scripts.
2. In Mercury, open the **CSD Python API** menu → **Options…** → **Add Location**, and choose that folder.
3. The script now appears in the **CSD Python API** menu.

## Usage

1. Open a structure containing at least one metal.
2. Optional: select one metal and its donor atoms.
3. Run **CSD Python API → CShM_CSM_CCM_analysis**.

## Limitations

- Coordination numbers from 2 to 12 only. Other sites are listed as skipped.
- Donor atoms are found with a distance criterion (sum of covalent radii + `BOND_TOL` = 0.45 Å, at most `R_SEARCH` = 3.3 Å). For long or ambiguous metal–ligand contacts (chelating carboxylates, weak axial contacts), adjust these two values at the top of the script, or select the donor atoms yourself. Metal–metal contacts are not counted as donors.
- With disordered structures, check the donor list in the report: atoms of several disorder components may be counted.

## Troubleshooting

The report shows the full error message. The most common cases:

| Message in the report | Cause | Fix |
|---|---|---|
| `Zn1 : CN 1 ignoré` (or a coordination number that is too low) | old version of the script, which only read the asymmetric unit | use the current version |
| `ModuleNotFoundError: No module named 'cosymlib'` | cosymlib is not installed in the Python given by `PYTHON_EXTERNE` (often after a CSD update) | run the `pip install` command shown at the end of the message |
| `processus tué par le signal 11` | cosymlib crashed on one measure | the other measures are still given; the failed one is listed in red |

## Credits and references

- **cosymlib**, developed by the Electronic Structure and Symmetry group (Universitat de Barcelona), does all the shape, symmetry and chirality calculations: <https://cosymlib.readthedocs.io/>. Please cite cosymlib and the original continuous shape and symmetry measures papers when you publish results.
- The script uses the CSD Python API from the Cambridge Crystallographic Data Centre (CCDC).
- Written by Yoann Rousselin with the help of Claude (Anthropic) as a coding assistant. This script is an interface between Mercury and cosymlib; it does not implement the measures itself.

## Licence

MIT. See [LICENSE](LICENSE). cosymlib is distributed under its own licence.

## Author

**Yoann Rousselin**, CNRS research engineer (X-ray crystallography), ICMUB, Dijon (France).
Questions or bug reports: please open an *Issue* on this page.
