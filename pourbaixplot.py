#!/usr/bin/env python3
#New version of the Pourbaix diagram script with temperature shifts and ion updates 
import copy
import numpy as np
import matplotlib.pyplot as plt
import os
from mp_api.client import MPRester
from pymatgen.core import Composition, Element
from pymatgen.analysis import pourbaix_diagram
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry
from thermo_data import (
    T_REF, cp_poly, standard_species_for_element,
    element_db, compound_db, species_db, ion_db, ion_charge,
    get_elem_S298_and_cp,
    CAL_TO_KJ, _cc_alpha_beta, _cc_S298_cal
)

# ============================  USER INPUT  =======================================================
elements = input("Enter elements (comma-separated): ").split(",")
T = float(input("Enter temperature in K: "))
conc = float(input("Enter concentration in mol/L (default = 1e-6): ") or "1e-6")
# Get API key from environment or prompt user
MP_API_KEY = os.getenv("MP_API_KEY")
if not MP_API_KEY:
    MP_API_KEY = input("Enter your Materials Project API key: ").strip()

# ============================  CONSTANTS  ========================================================
R = 8.314462618  # J/mol/K
F = 96485.33212  # C/mol
KJMOL_PER_EV = 96.485
T_REF = 298.15

def is_ref(TK: float) -> bool:
    return abs(TK - T_REF) < 5e-4

# Nernst slope
pourbaix_diagram.PREFAC = 0.0591 if is_ref(T) else (2.303 * R * T / F)

# ============================  NUMERIC HELPERS  ==================================================
def cp_poly(A,B=0,C=0,D=0,E=0):
    """Cp(T) [J/mol/K] = A + B T + C T^2 + D T^3 + E/T^2"""
    def _cp(T):
        t = max(float(T), 1e-6)
        return A + B*t + C*t*t + D*t**3 + E/(t*t)
    return _cp

def H_increment_kJmol(T_to, cp_func, Tref=T_REF, n=2001):
    """H(T)-H(298) [kJ/mol] from ∫ Cp dT."""
    if abs(T_to - Tref) < 1e-12:
        return 0.0
    a, b = (Tref, T_to) if T_to > Tref else (T_to, Tref)
    grid = np.linspace(a, b, n)
    cp = np.array([cp_func(float(t)) for t in grid])  # J/mol/K
    dH = np.trapz(cp, grid) / 1000.0                  # kJ/mol
    return dH if T_to > Tref else -dH

def S_at_T_kJmolK(T_to, S298_kJmolK, cp_func, Tref=T_REF, n=2001):
    """S(T) = S(298) + ∫ Cp/T dT  [kJ/mol/K]."""
    if abs(T_to - Tref) < 1e-12:
        return S298_kJmolK
    a, b = (Tref, T_to) if T_to > Tref else (T_to, Tref)
    grid = np.linspace(a, b, n)
    cp_over_T = np.array([cp_func(float(t))/float(t) for t in grid])  # J/mol/K^2
    dS = np.trapz(cp_over_T, grid) / 1000.0                           # kJ/mol/K
    return (S298_kJmolK + dS) if T_to > Tref else (S298_kJmolK - dS)

# Optional explicit ion thermo (HKF, etc.). You can key by signature for charge-specific data.
# Here we show *placeholders* for Cu+ and Cu2+ so they won’t fall back to CC.
def ion_signature(entry: PourbaixEntry):
    """Key distinguishing ions with same reduced formula: phase, formula, npH, nPhi, nH2O."""
    return (
        getattr(entry, "phase_type", None),
        entry.composition.reduced_formula.replace(" ", ""),
        float(getattr(entry, "npH", 0.0)),
        float(getattr(entry, "nPhi", 0.0)),
        float(getattr(entry, "nH2O", 0.0)),
    )

# ============================  ΔfH°(T) and ΔfG°(T)  =============================================
def elemental_terms_from_formula(comp: Composition):
    terms_by_species = {}
    for el, amt in comp.get_el_amt_dict().items():
        sp, atoms_per_sp = standard_species_for_element(Element(el))
        nu = amt / atoms_per_sp
        terms_by_species[sp] = terms_by_species.get(sp, 0.0) + nu
    terms = []
    for sp, nu in terms_by_species.items():
        S298, cp = get_elem_S298_and_cp(sp)
        terms.append({"nu": nu, "cp": cp, "S298": S298})
    return terms

def delta_f_H_T_kJmol(TK, dHf298_kJmol, cp_compound, elem_terms):
    dH_comp = H_increment_kJmol(TK, cp_compound)
    dH_elems = sum(e["nu"] * H_increment_kJmol(TK, e["cp"]) for e in elem_terms)
    return dHf298_kJmol + dH_comp - dH_elems

def delta_f_G_T(TK, dHf298_kJmol, cp_compound, S298_comp_kJmolK, elem_terms):
    dHf_T = delta_f_H_T_kJmol(TK, dHf298_kJmol, cp_compound, elem_terms)
    S_comp_T = S_at_T_kJmolK(TK, S298_comp_kJmolK, cp_compound)
    S_elems_T = sum(e["nu"] * S_at_T_kJmolK(TK, e["S298"], e["cp"]) for e in elem_terms)
    return dHf_T - TK * (S_comp_T - S_elems_T)  # kJ/mol

def solid_dG_kJ_per_fu(TK, rf):
    thermo = compound_db.get(rf)
    if thermo is None:
        return 0.0
    elem_terms = elemental_terms_from_formula(Composition(rf))
    G298 = delta_f_G_T(T_REF, thermo["dHf298"], thermo["cp"], thermo["S298"], elem_terms)
    GT   = delta_f_G_T(TK,    thermo["dHf298"], thermo["cp"], thermo["S298"], elem_terms)
    return GT - G298

# ============================  μ_H2O(T)  =========================================================
MU_H2O_298_eV = -2.4583
def mu_H2O_T_eV(TK: float) -> float:
    if is_ref(TK):
        return MU_H2O_298_eV
    elem_terms = elemental_terms_from_formula(Composition("H2O"))
    th = compound_db["H2O"]
    dG_kJ = delta_f_G_T(TK, th["dHf298"], th["cp"], th["S298"], elem_terms) \
          - delta_f_G_T(T_REF, th["dHf298"], th["cp"], th["S298"], elem_terms)
    return MU_H2O_298_eV + dG_kJ / KJMOL_PER_EV

pourbaix_diagram.MU_H2O = MU_H2O_298_eV if is_ref(T) else mu_H2O_T_eV(T)
print(f"\nSet PREFAC = {pourbaix_diagram.PREFAC:.5f} V and MU_H2O = {pourbaix_diagram.MU_H2O:.6f} eV/H2O at T = {T} K")

# ============================  ION UPDATE (signature-aware, cation-safe)  ========================
DEBUG_IONS = False
_debug_seen_ions = set()
def _debug_ion_once(key, msg):
    if not DEBUG_IONS or key in _debug_seen_ions:
        return
    print(msg); _debug_seen_ions.add(key)

def avg_cp_kJmolK(cp_func, T2, T1=T_REF, n=1501):
    if abs(T2 - T1) < 1e-12:
        return 0.0
    return H_increment_kJmol(T2, cp_func, Tref=T1, n=n) / (T2 - T1)

def metals_from_entry(e):
    metals = {}
    for el, amt in e.composition.reduced_composition.get_el_amt_dict().items():
        if el not in ("H", "O") and amt > 0:
            metals[el] = metals.get(el, 0.0) + amt
    if not metals:
        raise ValueError(f"No metal found in ion {e.composition}")
    return metals

def update_ion_energy(entry: PourbaixEntry, T2: float) -> float:
    """Return NEW absolute energy (eV) at T2 for this ion entry.
       Uses ion signature to separate charge states (e.g., Cu+ vs Cu2+).
       Skips Criss–Cobble for bare cations unless explicit ion_db data exists."""
    sig = ion_signature(entry)  # (phase, rf, npH, nPhi, nH2O)
    rf_key = entry.composition.reduced_formula.replace(" ", "")
    ion_thermo = ion_db.get(sig) or ion_db.get(rf_key)  # prefer signature, fallback formula

    # Partition stoichiometric "reaction" sides from entry coefficients:
    nH2O_reac = max(-entry.nH2O, 0.0); nH2O_prod = max(entry.nH2O, 0.0)
    nHp_reac  = max(-entry.npH,  0.0); nHp_prod  = max(entry.npH,   0.0)
    ne_reac   = max(-entry.nPhi, 0.0); ne_prod   = max(entry.nPhi,  0.0)
    metals = metals_from_entry(entry)

    comp = entry.composition.reduced_composition
    has_H = comp.get_el_amt_dict().get("H", 0.0) > 0
    has_O = comp.get_el_amt_dict().get("O", 0.0) > 0

    if ion_thermo is not None and callable(ion_thermo.get("cp", None)):
        # ----- explicit ion thermo path -----
        S_reac = ion_thermo["S298"] \
               + nHp_reac  * species_db["H+(aq)"]["S298"] \
               + ne_reac   * species_db["e-"]["S298"] \
               + nH2O_reac * species_db["H2O(l)"]["S298"]
        S_prod = sum(metals[m]*get_elem_S298_and_cp(m)[0] for m in metals) \
               + nHp_prod  * species_db["H+(aq)"]["S298"] \
               + ne_prod   * species_db["e-"]["S298"] \
               + nH2O_prod * species_db["H2O(l)"]["S298"]
        dS298 = S_prod - S_reac

        Cp_reac = avg_cp_kJmolK(ion_thermo["cp"], T2) \
                + nHp_reac  * 0.0 \
                + ne_reac   * 0.0 \
                + nH2O_reac * avg_cp_kJmolK(species_db["H2O(l)"]["cp"], T2)
        Cp_prod = sum(metals[m]*avg_cp_kJmolK(get_elem_S298_and_cp(m)[1], T2) for m in metals) \
                + nHp_prod  * 0.0 \
                + ne_prod   * 0.0 \
                + nH2O_prod * avg_cp_kJmolK(species_db["H2O(l)"]["cp"], T2)
        dCp_avg = Cp_prod - Cp_reac

    else:
        # ----- fallback path -----
        if (not has_H) and (not has_O):
            # Bare cation (e.g., Cu+, Cu2+, Al3+) WITHOUT explicit data -> skip correction
            _debug_ion_once(sig, "[ION DEBUG] Bare cation without ion_db — skipping CC fallback; returning original energy.")
            return entry.energy

        # Criss–Cobble for oxyanions / acid oxyanions
        n_O = int(round(comp.get_el_amt_dict().get("O", 0.0)))
        category = "acid_oxyanion" if has_H else "oxyanion"
        Z = ion_charge.get(rf_key)
        if Z is None:
            _debug_ion_once(sig, f"[ION DEBUG] No ion_charge for {rf_key}; skipping correction.")
            return entry.energy

        S298_cal = _cc_S298_cal(Z, n_O); S298_kJ = S298_cal * CAL_TO_KJ
        S_reac = S298_kJ \
               + nHp_reac  * species_db["H+(aq)"]["S298"] \
               + ne_reac   * species_db["e-"]["S298"] \
               + nH2O_reac * species_db["H2O(l)"]["S298"]
        S_prod = sum(metals[m]*get_elem_S298_and_cp(m)[0] for m in metals) \
               + nHp_prod  * species_db["H+(aq)"]["S298"] \
               + ne_prod   * species_db["e-"]["S298"] \
               + nH2O_prod * species_db["H2O(l)"]["S298"]
        dS298 = S_prod - S_reac

        a_cal, b = _cc_alpha_beta(category, T2)
        Cpbar_ion_kJ = (a_cal + b * S298_cal) * CAL_TO_KJ  # kJ/mol/K
        Cp_reac = Cpbar_ion_kJ \
                + nHp_reac  * 0.0 \
                + ne_reac   * 0.0 \
                + nH2O_reac * avg_cp_kJmolK(species_db["H2O(l)"]["cp"], T2)
        Cp_prod = sum(metals[m]*avg_cp_kJmolK(get_elem_S298_and_cp(m)[1], T2) for m in metals) \
                + nHp_prod  * 0.0 \
                + ne_prod   * 0.0 \
                + nH2O_prod * avg_cp_kJmolK(species_db["H2O(l)"]["cp"], T2)
        dCp_avg = Cp_prod - Cp_reac

    dT = T2 - T_REF
    if T2 <= 0.0:
        raise ValueError("Temperature must be > 0 K for the ln(T2/T_REF) term.")
    dG_kJmol = dCp_avg*dT - dS298*dT - T2*dCp_avg*np.log(T2/T_REF)
    return entry.energy + dG_kJmol / KJMOL_PER_EV  # absolute energy at T2 (eV)

# ============================  FETCH ENTRIES  ====================================================
with MPRester(MP_API_KEY) as mpr:
    entries = mpr.get_pourbaix_entries(elements)

# ============================  ION CONCENTRATIONS  ===============================================
solution_elements = set()
for entry in entries:
    if isinstance(entry, PourbaixEntry) and entry.phase_type == "Ion":
        solution_elements.update(el.symbol for el in entry.composition.elements)
conc_dict = {el: conc for el in solution_elements}
print(f"Using [C] = {conc:.1e} mol/L for ions: {', '.join(conc_dict) if conc_dict else '(none)'}")

# ============================  APPLY FULL T-SHIFTS (no fixing to 298 K set)  =====================
if is_ref(T):
    entries_all_fixed = [copy.deepcopy(e) for e in entries]
else:
    # Helper: total energy that PourbaixEntry.energy returns at pH=0, V=0
    def total_E_pH0_V0(entry_obj):
        return (entry_obj.uncorrected_energy
                + pourbaix_diagram.PREFAC * np.log10(entry_obj.concentration)
                - pourbaix_diagram.MU_H2O * entry_obj.nH2O)

    # Solids first
    entries_shifted_total = []
    for e in entries:
        if e.phase_type == "Solid":
            rf = e.composition.reduced_formula
            dE_eV = solid_dG_kJ_per_fu(T, rf) / KJMOL_PER_EV
            E_total_target = e.energy + dE_eV
            e2 = copy.deepcopy(e)
            conc_term = pourbaix_diagram.PREFAC * np.log10(e2.concentration)
            e2.uncorrected_energy = E_total_target - conc_term + pourbaix_diagram.MU_H2O * e2.nH2O
            entries_shifted_total.append(e2)
        else:
            entries_shifted_total.append(copy.deepcopy(e))

    # Then ions — signature-aware update
    entries_all_fixed = []
    for e in entries_shifted_total:
        if e.phase_type == "Ion":
            new_abs_E = update_ion_energy(e, T)   # absolute at T
            dE_eV = new_abs_E - e.energy
            E_total_target = e.energy + dE_eV
            e2 = copy.deepcopy(e)
            conc_term = pourbaix_diagram.PREFAC * np.log10(e2.concentration)
            e2.uncorrected_energy = E_total_target - conc_term + pourbaix_diagram.MU_H2O * e2.nH2O
            entries_all_fixed.append(e2)
        else:
            entries_all_fixed.append(e)

# ============================  BUILD & PLOT  =====================================================
pbx = PourbaixDiagram(entries_all_fixed, conc_dict=conc_dict)
plotter = PourbaixPlotter(pbx)
ax = plotter.get_pourbaix_plot()


mode = "(full T: PREFAC(T), μ_H2O(T), solids ΔG(T), ions w/ signatures)"
ax.set_title(f"Pourbaix Diagram at {T:.2f} K, [C] = {conc:g} mol/L {mode}")
plt.tight_layout()
plt.show()

##################### ============================  SAVE DIAGRAM  =====================================================
output_file = input("Enter output file name (default: pourbaix_diagram.png): ") or "pourbaix_diagram.png"
plt.savefig(output_file, dpi=300)
print(f"Pourbaix diagram saved to {output_file}")
plt.close()
print("Done.")
###### ===================Get pH Values for the Pourbaix diagram========================
import numpy as np

def vertical_ph_all(ax,
                    ndp=4,
                    # visual tolerances (in pixels → mapped to data units)
                    vert_eps_px=0.75,      # treat as vertical if Δx < ~1 px
                    cluster_eps_px=0.75,   # merge segments if x within ~1 px
                    # optional filters
                    exclude_edges=True, edge_eps_px=0.75,
                    exclude_neutral=True, neutral_pH=7.0, neutral_eps=1e-3,
                    # output
                    return_spans=False,
                    debug=False):
    """
    Find ALL vertical Pourbaix boundaries regardless of E=0 crossing.

    Returns:
      - if return_spans=False (default): sorted list of unique pH positions (rounded to ndp)
      - if return_spans=True: list of dicts with {'pH','ymin','ymax'} (pH rounded to ndp)
    """
    # map ~1 pixel to data units
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    bb = ax.bbox
    x_per_px = (x1 - x0) / max(bb.width, 1.0)
    y_per_px = (y1 - y0) / max(bb.height, 1.0)
    x_eps = x_per_px * vert_eps_px
    x_cluster = x_per_px * cluster_eps_px
    edge_eps = x_per_px * edge_eps_px

    verts = []  # (x_mid, y_min, y_max)

    def consider_segment(x1s, y1s, x2s, y2s):
        x1s = float(x1s); x2s = float(x2s); y1s = float(y1s); y2s = float(y2s)
        if not (np.isfinite(x1s) and np.isfinite(x2s) and np.isfinite(y1s) and np.isfinite(y2s)):
            return
        if abs(x1s - x2s) <= x_eps:  # nearly vertical (in data units)
            xm = 0.5 * (x1s + x2s)
            ymin, ymax = (y1s, y2s) if y1s <= y2s else (y2s, y1s)
            verts.append((xm, ymin, ymax))

    # 1) Line2D segments
    for ln in ax.lines:
        x = ln.get_xdata(); y = ln.get_ydata()
        if x is None or y is None or len(x) != len(y) or len(x) < 2:
            continue
        for i in range(len(x) - 1):
            consider_segment(x[i], y[i], x[i+1], y[i+1])

    # 2) LineCollections (typical PB boundaries)
    for coll in ax.collections:
        get_segments = getattr(coll, "get_segments", None)
        if get_segments is None:
            continue
        for (sx1, sy1), (sx2, sy2) in get_segments():
            consider_segment(sx1, sy1, sx2, sy2)

    if debug:
        print(f"[DBG] collected {len(verts)} vertical-ish segments (x_eps={x_eps:.2e})")

    # cluster by x to merge stacked/split segments forming one vertical boundary
    verts.sort(key=lambda t: t[0])
    clusters = []
    cur = []
    for v in verts:
        if not cur or abs(v[0] - cur[-1][0]) <= x_cluster:
            cur.append(v)
        else:
            clusters.append(cur); cur = [v]
    if cur:
        clusters.append(cur)

    # summarize clusters
    results_raw = []
    for cl in clusters:
        xs = [v[0] for v in cl]
        ymin = min(v[1] for v in cl)
        ymax = max(v[2] for v in cl)
        x_rep = float(np.median(xs))
        # optional filters
        if exclude_edges and (abs(x_rep - x0) <= edge_eps or abs(x_rep - x1) <= edge_eps):
            if debug: print(f"[DBG] drop x≈{x_rep:.6f}: near x-limits")
            continue
        if exclude_neutral and abs(x_rep - neutral_pH) <= neutral_eps:
            if debug: print(f"[DBG] drop x≈{x_rep:.6f}: near pH={neutral_pH}")
            continue
        results_raw.append((x_rep, ymin, ymax))

    if not return_spans:
        # round, unique, sort by pH
        ph_list = sorted({round(x, ndp) for x, _, _ in results_raw})
        return ph_list
    else:
        # provide spans too
        out = []
        # merge by rounded pH to avoid duplicates from closely spaced clusters
        by_ph = {}
        for x, ymin, ymax in results_raw:
            xr = round(x, ndp)
            if xr not in by_ph:
                by_ph[xr] = [ymin, ymax]
            else:
                by_ph[xr][0] = min(by_ph[xr][0], ymin)
                by_ph[xr][1] = max(by_ph[xr][1], ymax)
        for xr in sorted(by_ph.keys()):
            ymin, ymax = by_ph[xr]
            out.append({'pH': xr, 'ymin': ymin, 'ymax': ymax})
        return out

ph_vals = vertical_ph_all(ax)
print("pH with vertical boundaries:", ph_vals)

# # pH + y-span for each vertical boundary
# ph_spans = vertical_ph_all(ax, return_spans=True)
# for v in ph_spans:
#     print(f"pH={v['pH']}: spans E in [{v['ymin']:.3f}, {v['ymax']:.3f}] V")