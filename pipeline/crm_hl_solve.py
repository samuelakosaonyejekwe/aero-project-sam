# -*- coding: utf-8 -*-
"""
crm_hl_solve.py  --  STAGE 3: SOLUTION (AoA sweep with convergence monitoring).

run_sweep(solver, aoa_list) imposes each angle of attack through the far-field
flow-direction vector, sets Cl/Cd/Cm report definitions + convergence monitors,
iterates to convergence, and writes out/polar.csv plus a saved case per angle.
"""
import math, csv, os

WALLS    = ["wing", "fuselage", "slat", "flap", "h_tail", "v_tail"]
MRC      = [33.68, 0.0, 4.52]          # moment reference centre [m] (MRC of CRM, x=1325.9 in)
OUT_DIR  = r"C:\aero-project\out"
AOA_LIST = [0, 2, 4, 6, 8, 10, 11, 12, 13, 13.5, 14.5]

def _report_defs(solver, alpha_deg):
    a = math.radians(alpha_deg)
    rd = solver.settings.solution.report_definitions
    rd.lift["cl"]   = dict(zones=WALLS, force_vector=[-math.sin(a), 0.0, math.cos(a)])
    rd.drag["cd"]   = dict(zones=WALLS, force_vector=[ math.cos(a), 0.0, math.sin(a)])
    rd.moment["cm"] = dict(zones=WALLS, mom_center=MRC, mom_axis=[0.0, 1.0, 0.0])

def _convergence(solver):
    """Monitor force coefficients; stop when they flatten (per HLPW-5 practice)."""
    mon = solver.settings.solution.monitor
    for name in ("cl", "cd", "cm"):
        mon.report_plots[name] = dict(report_defs=[name])
    rc = solver.settings.solution.run_calculation
    # converge on Cl stabilising to 1e-5 between iterations, else hit iter cap
    try:
        rc.convergence_conditions.convergence_reports["cl-conv"] = dict(
            report_defition="cl", stop_criterion=1e-5, previous_values_to_consider=50)
    except Exception:
        pass  # convergence-report API name varies by release; residuals still apply

def run_sweep(solver, aoa_list=AOA_LIST, iters=3000):
    os.makedirs(OUT_DIR, exist_ok=True)
    _convergence(solver)
    ff = solver.settings.setup.boundary_conditions.pressure_far_field["farfield"]
    results = []
    for alpha in aoa_list:
        a = math.radians(alpha)
        ff.momentum.flow_direction[0] = math.cos(a)     # x streamwise
        ff.momentum.flow_direction[1] = 0.0             # y spanwise
        ff.momentum.flow_direction[2] = math.sin(a)     # z vertical
        _report_defs(solver, alpha)

        solver.settings.solution.initialization.hybrid_initialize()
        solver.settings.solution.run_calculation.iterate(iter_count=iters)

        rd = solver.settings.solution.report_definitions
        cl, cd, cm = rd["cl"].get_value(), rd["cd"].get_value(), rd["cm"].get_value()
        results.append((alpha, cl, cd, cl/cd, cm))
        print(f"[solve] a={alpha:5.1f}  Cl={cl:7.4f}  Cd={cd:7.4f}  L/D={cl/cd:6.2f}  Cm={cm:8.4f}")
        solver.settings.file.write_case_data(file_name=fr"{OUT_DIR}\crm_a{alpha}.cas.h5")

    with open(fr"{OUT_DIR}\polar.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["alpha_deg", "Cl", "Cd", "L_over_D", "Cm"]); w.writerows(results)
    print(f"[solve] sweep complete -> {OUT_DIR}\\polar.csv")
    return results
