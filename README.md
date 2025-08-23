# Temperature-Dependent Pourbaix Diagram Tool

A Python utility for generating **temperature-aware Pourbaix diagrams** with concentration control and signature-aware **ion free-energy updates**. The script extends pymatgen’s standard Pourbaix workflow by adjusting:

- **Nernst slope** `PREFAC(T)`
- **Water chemical potential** `μ_H2O(T)`
- **Solid phase ΔG** shifts via caloric integrals
- **Ion energies** using either explicit ion thermo (if provided) or **Criss–Cobble**-style fallback for oxy(anions)
- Optional **pH boundary extraction** from the rendered plot

The implementation is split across a *main* script (shown below) and a companion module **`thermo_data.py`** that houses thermochemical registries and utilities.

---

## Table of Contents
1. [Key Features](#key-features)
2. [What This Script Changes vs Pymatgen](#what-this-script-changes-vs-pymatgen)
3. [File Structure](#file-structure)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [How It Works](#how-it-works)
7. [Inputs & CLI Usage](#inputs--cli-usage)
8. [Outputs](#outputs)
9. [Thermochemical Data Contracts (`thermo_data.py`)](#thermochemical-data-contracts-thermo_datapy)
10. [Ion-Energy Update Logic](#ion-energy-update-logic)
11. [Extracting Vertical pH Boundaries](#extracting-vertical-ph-boundaries)
12. [Troubleshooting](#troubleshooting)
13. [Reproducibility & Environments](#reproducibility--environments)
14. [Limitations & Notes](#limitations--notes)
15. [License](#license)

---

## Key Features
- **Temperature-aware Nernst slope**: Sets `pourbaix_diagram.PREFAC = 2.303·R·T/F` at runtime (default 0.0591 V at 298.15 K).
- **Water chemical potential μ_H₂O(T)**: Shifts relative to 298 K by integrating Cp and S(T) for water.
- **Solid ΔG(T)**: Computes formation Gibbs free energy shifts for solids using caloric integrals of Cp(T) and S(T).
- **Signature-aware ion updates**: Per-ion updates keyed by a **signature** `(phase, reduced formula, n_pH, n_Phi, n_H2O)` to distinguish charge states and hydration reactions. Falls back to **Criss–Cobble** for oxyanions/acid oxyanions when explicit ion data are not present; **bare metal cations are left unmodified** unless explicit data exist.
- **Concentration control**: Uniform ion concentrations via `conc_dict` (default 1e-6 M for all solution species present in the entries set).
- **Plot post-processing**: Utility to list **vertical pH boundaries** detected in the diagram.

---

## What This Script Changes vs Pymatgen
Pymatgen’s Pourbaix modules are primarily standard-state (298.15 K) and use a constant Nernst prefactor and water chemical potential. This tool modifies at runtime:

1. `PREFAC(T)`: Nernst slope scaled with **T**.
2. `MU_H2O(T)`: Effective water chemical potential shifted from its 298 K value by ΔG(T)−ΔG(298).
3. Solid phase energies: absolute **per-entry** energy adjusted to reflect **ΔG(T)** of the phase.
4. Ion phase energies: **absolute** per-entry energy re-evaluated using either explicit ion thermo (if present) or **Criss–Cobble** fallback for oxy(anions); **bare cations** are intentionally not adjusted unless explicit ion data exist (prevents incorrect application of CC to cations).

The goal is to produce Pourbaix stability fields that better reflect **non-ambient temperatures** while retaining the Materials Project entry structure and Pymatgen plotting ecosystem.

---

## File Structure
```
project/
├── pourbaixplot.py                 # The main script (your provided code)
├── thermo_data.py          # Thermochem registries & helper functions
├── environment.yml         # (optional) conda environment for reproducibility
└── README.md               # This file
```

---

## Installation
Use either **conda** or **venv + pip**. Conda is recommended for scientific stacks.

### Conda (recommended)
```bash
conda env create -f environment.yml   # or create manually
conda activate pourbaix
```
If creating manually:
```bash
conda create -n pourbaix python=3.11 numpy matplotlib
conda activate pourbaix
pip install pymatgen mp-api
```

---

## Configuration
- **Materials Project API key**: The script will look for an environment variable `MP_API_KEY`. If it is not set, the script will **prompt you to enter it interactively** at runtime.

To set the key permanently:
- Linux/macOS:
```bash
export MP_API_KEY="your_key_here"
```
- Windows (PowerShell):
```powershell
setx MP_API_KEY "your_key_here"
```
- **Matplotlib backend**: If running headless, set `MPLBACKEND=Agg` or save figures without showing.
- **Concentration**: default is `1e-6` M for all solution species present in the entries set; overridden via prompt.

---

## How It Works
### 1) Temperature-Dependent Nernst Factor
```python
pourbaix_diagram.PREFAC = 0.0591 if is_ref(T) else (2.303 * R * T / F)
```

### 2) Water Chemical Potential μ_H₂O(T)
- Baseline at 298 K: `MU_H2O_298_eV = -2.4583`
- Shift by ΔG(T)−ΔG(298) computed from caloric integrals using water’s Cp(T) and S(298):
```python
pourbaix_diagram.MU_H2O = MU_H2O_298_eV if is_ref(T) else mu_H2O_T_eV(T)
```

### 3) Solid ΔG(T)
For a reduced formula `rf`, the script computes per-formula free energy shifts using:
- `delta_f_H_T_kJmol(T)` from integrals of Cp(T)
- `S_at_T_kJmolK(T)` for entropic terms
- Elemental reference states from `standard_species_for_element(...)`

### 4) Ion Energy Updates
- **Explicit path**: If `ion_db` contains a matching **signature** or formula key with Cp(T) and S(298), the script computes ΔG(T) for the full reaction encoded by the Pourbaix entry (accounting for `nH2O`, `npH`, `nPhi`).
- **Fallback (Criss–Cobble)**: For oxyanions/acid oxyanions with known charge **Z** and oxygen count, it estimates S(298) and Cp(T) using CC α(T), β(T) correlations (provided in `thermo_data.py`).
- **Bare cations**: No CC fallback is applied—left at original energy unless explicit ion data are supplied.

### 5) Rebuilding Entries
Each adjusted entry’s **`uncorrected_energy`** is set so that the **total energy at pH=0, E=0** equals the temperature-shifted target energy. This keeps Pymatgen’s plotting logic consistent while reflecting temperature adjustments.

---

## Inputs & CLI Usage
The script is interactive:
1. **Elements** (comma-separated): e.g., `Si,O` or `Al,Si,O`.
2. **Temperature (K)**: e.g., `350`.
3. **Ion concentration (mol/L)**: default `1e-6`.
4. **Output filename** for the figure: default `pourbaix_diagram.png`.

### Example Run
```bash
python pourbaixplot.py
Enter elements (comma-separated): Si, O
Enter temperature in K: 350
Enter concentration in mol/L (default = 1e-6):
Enter output file name (default: pourbaix_diagram.png): si_o_350k.png
```

---

## Outputs
- **Figure**: Saved PNG of the Pourbaix diagram (default `pourbaix_diagram.png`).
- **Console log**: Shows `PREFAC` and `μ_H2O` at T, the applied ion concentrations, and the list of **vertical pH boundaries** detected (see below).

---

## Thermochemical Data Contracts (`thermo_data.py`)
`thermo_data.py` provides data and utilities. Expected contents (examples):

- **Constants & helpers**
  - `T_REF = 298.15`
  - `cp_poly(A, B=0, C=0, D=0, E=0)` → returns a function Cp(T) [J/mol/K]
  - `CAL_TO_KJ = 4.184e-3`

- **Element registry**
  - `element_db: { symbol: {"S298": kJ/mol/K, "cp": callable}, ... }`
  - `standard_species_for_element(Element)` → returns `(sp, atoms_per_sp)` where `sp` is the standard reference (e.g., `O2(g)`, `H2(g)`, `Al(s)`).
  - `get_elem_S298_and_cp(sp)` → `(S298_kJmolK, cp_callable)` for the reference species `sp`.

- **Compound registry**
  - `compound_db[rf] = {"dHf298": kJ/mol, "S298": kJ/mol/K, "cp": callable}`
  - `rf` is reduced formula string (e.g., `SiO2`, `H2O`). Used by solids and by μ_H₂O(T).

- **Aqueous species registry**
  - `species_db["H+(aq)"] = {"S298": kJ/mol/K}`
  - `species_db["e-"]      = {"S298": kJ/mol/K}`
  - `species_db["H2O(l)"]  = {"S298": kJ/mol/K, "cp": callable}`

- **Ion registry (optional, explicit thermo)**
  - `ion_db[key] = {"S298": kJ/mol/K, "cp": callable}` where `key` is either the **ion signature** or the reduced formula string. If present, explicit thermo overrides CC fallback.

- **Ion charge map**
  - `ion_charge[rf] = Z` (integer), e.g., `{ "H4SiO4": -0, "H7SiO6": -1, ... }`

- **Criss–Cobble parameters**
  - `_cc_S298_cal(Z, n_O)` → S(298) estimate in **cal/mol/K** for an oxyanion of charge `Z` with `n_O` oxygens.
  - `_cc_alpha_beta(category, T)` → returns `(α_cal, β)` for the chosen category (`"oxyanion"` or `"acid_oxyanion"`) at temperature `T` (°C). See your `_CC_TABLE`.

> **Important**: Units must be consistent. The script expects S in **kJ/mol/K** and Cp integrals in **kJ/mol** unless otherwise indicated.

---

## Ion-Energy Update Logic
Each `PourbaixEntry` encodes a reaction of the form:
```
Ion + nH2O·H2O + npH·H+ + nPhi·e-  ⇌  Metals (as elements) + nH2O'·H2O + npH'·H+ + nPhi'·e-
```
The code computes **ΔG(T) − ΔG(298)** for that reaction by assembling **reactant** and **product** entropies and heat capacities from either:
- explicit `ion_db` (preferred), or
- **Criss–Cobble** estimates for oxyanions / acid oxyanions when explicit data are missing.

**Bare cations are not modified** (unless explicit data exist) to avoid misapplication of anion correlations.

The final absolute energy at `T` becomes:
```
E_T = E_298 + [ΔG(T) − ΔG(298)] / (96.485 kJ/mol per eV)
```
`uncorrected_energy` is then set so that the total energy used by Pymatgen at pH=0, E=0 equals `E_T` after including concentration and water terms.

---

## Extracting Vertical pH Boundaries
The helper `vertical_ph_all(ax, ...)` scans the Matplotlib **Axes** after plotting and:
- Collects near-vertical segments (within ~1 pixel in x),
- Clusters and merges them to avoid duplicates,
- Optionally removes edges and exact neutral pH,
- Returns a **sorted list of unique pH values** (or full spans with `ymin/ymax`).

**Usage** (already included in the script):
```python
ph_vals = vertical_ph_all(ax)
print("pH with vertical boundaries:", ph_vals)
```
Parameters let you control rounding, pixel tolerances, and whether to return spans.

---

## Troubleshooting
**1) `KeyError`/`NoneType` in thermo lookup**
- Ensure `compound_db[rf]` exists for each **solid** reduced formula present (e.g., `SiO2`, `Al2O3`, `H2O`).
- Ensure `ion_charge[rf]` is defined for ions that rely on CC fallback.
- If an ion is a **bare cation**, either supply explicit `ion_db` data or accept that no T-shift will be applied.

**2) `AttributeError: 'Axes' object has no attribute 'show'`**
- Use `plt.show()` instead of `ax.show()`.

**3) Materials Project authentication**
- 401/403 errors indicate missing/invalid `MP_API_KEY`.

**4) Numeric stability**
- Very low or zero temperatures will raise errors (log terms). Keep `T > 0 K`.

**5) Plot not saving**
- Ensure `plt.tight_layout(); plt.savefig(...); plt.close()` order; make sure the output path is writable.

**6) Unexpected vertical lines**
- Adjust `vert_eps_px` and `cluster_eps_px` to tune what counts as “vertical”.

---

## Reproducibility & Environments
- **Export env**: `conda env export --from-history > environment.yml`
- **Create env**: `conda env create -f environment.yml`
- **Update env**: `conda env update -f environment.yml --prune`
- **venv** alternative**:**
  ```bash
  python -m venv .venv
  source .venv/bin/activate
  pip install -r requirements.txt
  ```

Include `environment.yml` (or `requirements.txt`) alongside the code when transferring.

---

## Limitations & Notes
- The **quality of results depends on your `thermo_data.py`**. Provide accurate Cp(T), S(298), and ΔHf°(298) for solids and water; supply ion data where possible.
- **Criss–Cobble** is a heuristic for oxy(anions) which can be extended to cations.
- The script currently applies **uniform ion concentration**; element-specific activities can be extended by editing `conc_dict`.
- Pymatgen’s Pourbaix diagrams assume ideal dilute solutions and standard electrochemical conventions; interpret high-T diagrams accordingly.

---

## License
This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.


