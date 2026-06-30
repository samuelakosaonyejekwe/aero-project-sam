# -*- coding: utf-8 -*-
"""
crm_hl_post.py  --  STAGE 4: POST-PROCESSING.

Two parts:
  (1) export_fields(solver, alpha):  in-Fluent contours & pathlines (pressure,
      Cp, velocity, skin friction) written as images -- call inside the sweep.
  (2) build_polar_report():  reads out/polar.csv, derives L/D, (L/D)max,
      C_L,max and stall angle, writes the engineering charts, and prints a
      validation check against the HLPW-5 reference C_L at 11 deg.
Part (2) needs only NumPy/Matplotlib and runs without Fluent.
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = r"C:\aero-project\out"
FIG_DIR = r"C:\aero-project\figures"
HLPW_CL_AT_11 = 1.07054          # reference value for validation

# ---------- (1) in-solver field export --------------------------------------
def export_fields(solver, alpha):
    g = solver.settings.results.graphics
    defs = {
        "pressure": dict(field="pressure", surfaces=["wing", "fuselage", "slat", "flap"]),
        "cp":       dict(field="pressure-coefficient", surfaces=["wing", "slat", "flap"]),
        "velocity": dict(field="velocity-magnitude", surfaces=["symmetry"]),
        "wall-shear": dict(field="wall-shear-stress", surfaces=["wing"]),
    }
    for name, d in defs.items():
        cname = f"c_{name}"
        g.contour[cname] = dict(field=d["field"], surfaces_list=d["surfaces"])
        g.contour[cname].display()
        solver.settings.results.graphics.views.restore_view(view_name="isometric")
        solver.settings.results.graphics.picture.save_picture(
            file_name=fr"{FIG_DIR}\field_{name}_a{alpha}.png")
    # streamlines on the symmetry plane
    g.pathline["streamlines"] = dict(field="velocity-magnitude", surfaces_list=["symmetry"])
    g.pathline["streamlines"].display()
    solver.settings.results.graphics.picture.save_picture(file_name=fr"{FIG_DIR}\field_pathlines_a{alpha}.png")
    print(f"[post] fields exported for alpha={alpha}")

# ---------- (2) polar report (Fluent-free) ----------------------------------
def build_polar_report(csv_path=fr"{OUT_DIR}\polar.csv"):
    if not os.path.exists(csv_path):
        print(f"[post] {csv_path} not found -- run the solve stage first."); return None
    rows = list(csv.DictReader(open(csv_path)))
    a  = np.array([float(r["alpha_deg"])  for r in rows])
    cl = np.array([float(r["Cl"])         for r in rows])
    cd = np.array([float(r["Cd"])         for r in rows])
    cm = np.array([float(r["Cm"])         for r in rows])
    lod = cl/cd

    iclmax = int(np.argmax(cl)); ilodmax = int(np.argmax(lod))
    summary = dict(CL_max=cl[iclmax], stall_alpha=a[iclmax],
                   LD_max=lod[ilodmax], alpha_LDmax=a[ilodmax], CD_min=cd.min())
    # validation
    if 11 in a:
        clv = cl[list(a).index(11)]
        err = 100*abs(clv-HLPW_CL_AT_11)/HLPW_CL_AT_11
        print(f"[post] VALIDATION  Cl(11 deg)={clv:.4f}  vs HLPW-5 {HLPW_CL_AT_11}  -> {err:.2f}% diff")

    def plot(x, y, xl, yl, t, fn, mark=None):
        fig, ax = plt.subplots(figsize=(6.2,4.2))
        ax.plot(x, y, "o-", color="#1f4e79", lw=1.8, ms=5)
        if mark: ax.axvline(mark, color="#b22222", ls=":", lw=1.2)
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(True, ls=":", alpha=0.6)
        ax.set_title(t, fontsize=11, fontweight="bold", color="#1f4e79")
        fig.tight_layout(); fig.savefig(fn, dpi=150); plt.close(fig)

    plot(a, cl, "alpha [deg]", "C_L", "Lift curve", fr"{FIG_DIR}\post_cl.png", summary["stall_alpha"])
    plot(a, cd, "alpha [deg]", "C_D", "Drag curve", fr"{FIG_DIR}\post_cd.png", summary["stall_alpha"])
    plot(cd, cl, "C_D", "C_L", "Drag polar",       fr"{FIG_DIR}\post_polar.png")
    plot(a, lod, "alpha [deg]", "L/D", "Efficiency",fr"{FIG_DIR}\post_lod.png")
    print(f"[post] summary: C_L,max={summary['CL_max']:.3f} @ {summary['stall_alpha']} deg | "
          f"(L/D)max={summary['LD_max']:.1f} @ {summary['alpha_LDmax']} deg | C_D,min={summary['CD_min']:.4f}")
    return summary

if __name__ == "__main__":
    build_polar_report()
