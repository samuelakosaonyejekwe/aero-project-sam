# -*- coding: utf-8 -*-
"""
run_crm_hl.py  --  MASTER RUNNER: end-to-end CRM-HL-WB CFD pipeline.

Chains the canonical ANSYS stages:
    Stage 0  Geometry   ->  Stage 1  Mesh   ->  Stage 2  Setup
    Stage 3  Solution   ->  Stage 4  Post-processing

Usage:
    python run_crm_hl.py                 # full pipeline, generated STL geometry
    python run_crm_hl.py --cad my.pmdb   # use external watertight CAD
    python run_crm_hl.py --post-only     # rebuild charts/report from out/polar.csv

Requires a local ANSYS 2024 R2 (Fluent + Fluent Meshing) for stages 1-4.
Stage 0 and the post 'report' part run with NumPy/Matplotlib alone.
"""
import argparse, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cad", default=None, help="watertight CAD (.pmdb/.scdoc/.stp); else generate STL")
    ap.add_argument("--turbulence", default="spalart-allmaras", choices=["spalart-allmaras", "k-omega"])
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--export-fields", action="store_true", help="export contours/pathlines per AoA")
    ap.add_argument("--post-only", action="store_true", help="only rebuild the polar report")
    args = ap.parse_args()

    import crm_hl_post as post
    if args.post_only:
        post.build_polar_report(); return

    # --- Stage 0: Geometry --------------------------------------------------
    import crm_hl_geometry as geom
    print("\n========== STAGE 0: GEOMETRY ==========")
    for k, v in geom.ref_quantities().items(): print(f"  {k}: {v}")
    if args.cad is None:
        geom_path = geom.generate_stl(); watertight = False
    else:
        geom_path = args.cad;            watertight = True
    print(f"  geometry: {geom_path}  (watertight={watertight})")

    # --- Stage 1: Mesh ------------------------------------------------------
    print("\n========== STAGE 1: MESH ==========")
    import crm_hl_mesh as mesh
    solver = mesh.build_mesh(geom_path, watertight=watertight)   # returns live solver session

    # --- Stage 2: Setup -----------------------------------------------------
    print("\n========== STAGE 2: SETUP ==========")
    import crm_hl_setup as setup
    setup.configure(solver, turbulence=args.turbulence)

    # --- Stage 3: Solution --------------------------------------------------
    print("\n========== STAGE 3: SOLUTION ==========")
    import crm_hl_solve as solve
    if args.export_fields:
        # wrap the sweep so fields are exported at a few representative angles
        orig = solve.run_sweep
    solve.run_sweep(solver, iters=args.iters)
    if args.export_fields:
        for alpha in (0, 8, 13.5):
            try: post.export_fields(solver, alpha)
            except Exception as e: print(f"[post] field export skipped a={alpha}: {e}")

    # --- Stage 4: Post-processing ------------------------------------------
    print("\n========== STAGE 4: POST-PROCESSING ==========")
    post.build_polar_report()

    try: solver.exit()
    except Exception: pass
    print("\nPipeline complete.")

if __name__ == "__main__":
    sys.exit(main())
