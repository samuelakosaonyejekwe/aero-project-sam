# Aero CFD Case Study — CRM-HL Wing-Body (PyAnsys)

A holistic, end-to-end computational-aerodynamics case study built around the
**NASA High-Lift Common Research Model (CRM-HL) Wing-Body** and implemented as a
fully scripted **PyAnsys** pipeline for a local ANSYS Fluent installation.

The study synthesises six source CFD investigations into a single, reproducible
workflow. The unifying problem is the RANS prediction of the aerodynamic
**force system versus angle of attack** (C_L, C_D, C_M, L/D). The pipeline sweeps
α = 0°–14.5° and reproduces a critical (stall) angle of ≈ 13.5° with
C_L,max ≈ 1.193.

---

## Highlights

- **Physics:** subsonic high-lift wing-body, M∞ = 0.20, Re ≈ 5.6 × 10⁶, RANS (Spalart–Allmaras).
- **Result:** lift-curve slope, drag polar, pitching moment and L/D across a full AoA sweep;
  stall captured at α ≈ 13.5°, C_L,max ≈ 1.193.
- **Live validation:** a boundary-layer-resolved 2.5-D section sweep (y⁺ ≈ 1, ~70 k cells)
  was executed locally on ANSYS Student 2025 R2; lift-curve slope C_lα ≈ 0.103/deg matches
  thin-aerofoil theory.
- **Reproducible:** a 5-stage PyAnsys pipeline (geometry → mesh → setup → solve → post)
  chained by a single master runner.

---

## Deliverable

| File | Description |
|------|-------------|
| [`aero_project_report.pdf`](aero_project_report.pdf) | Full project report — literature synthesis, geometry and reference quantities, the force system, holistic input data, the 5-stage PyAnsys workflow, primary results and the complete collation of output data from the source studies. |

---

## Repository structure

```
aero-project-sam/
├─ README.md
├─ aero_project_report.pdf   # full project report (primary deliverable)
├─ pipeline/                 # PyAnsys CFD pipeline (5 ANSYS stages + master runner)
│   ├─ crm_hl_geometry.py    #   Stage 0  Geometry  (reference quantities, y+ sizing, STL/CAD)
│   ├─ crm_hl_mesh.py        #   Stage 1  Mesh      (Fluent Meshing: watertight / fault-tolerant)
│   ├─ crm_hl_setup.py       #   Stage 2  Setup     (models, materials, BCs, reference values)
│   ├─ crm_hl_solve.py       #   Stage 3  Solution  (AoA sweep + convergence monitors)
│   ├─ crm_hl_post.py        #   Stage 4  Post      (contours, drag polar, validation)
│   ├─ run_crm_hl.py         #   master runner — chains stages 0 → 4
│   ├─ gen_clean_geom.py     #   single watertight wing solid for the live demo
│   ├─ local_mesh_2p5d.py    #   gmsh boundary-layer mesh for the 2.5-D section -> CGNS
│   └─ local_solve_2p5d.py   #   local Fluent RANS sweep -> out/polar_live.csv
├─ geometry/                 # geometry inputs and generated STL solids
├─ mesh/                     # generated meshes (CGNS)
├─ out/                      # solver outputs (polar CSV files)
└─ figures/                  # drawings, result charts and field previews
```

---

## Requirements

- **Python 3.13**, with: `ansys-fluent-core`, `numpy`, `matplotlib`

  ```bash
  pip install ansys-fluent-core numpy matplotlib
  ```

- A local **ANSYS Fluent + Fluent Meshing** for pipeline stages 1–4 (verified on 2025 R2).
- Stage 0 (geometry) runs with NumPy alone — no ANSYS required.

---

## How to run

**Full CFD pipeline (requires local ANSYS):**

```bash
cd pipeline
python run_crm_hl.py                  # full pipeline using the generated STL geometry
python run_crm_hl.py --cad crm.pmdb   # use the official watertight HLPW CAD (recommended)
python run_crm_hl.py --post-only      # rebuild charts from out/polar.csv
python run_crm_hl.py --export-fields  # export contours/pathlines at α = 0, 8, 13.5°
```

**Live 2.5-D section sweep (as executed locally on ANSYS Student 2025 R2):**

```bash
cd pipeline
python local_mesh_2p5d.py     # gmsh boundary-layer mesh -> CGNS (no ANSYS needed)
python local_solve_2p5d.py    # Fluent RANS sweep -> out/polar_live.csv
```

---

## Source studies

Directly relevant — inform the case:

1. Amaya (2025) — HLPW-5 Case 1, CRM-HL-WB (Uniandes thesis) — **primary geometry & dataset**
2. Zore et al. (2018) — ANSYS high-lift configurations, HL-CRM & JSM (AIAA 2018-2844)
3. Hiremath & Malipatil (2014) — Aircraft body vs AoA (IJIRSET)
4. Koç et al. (2024) — UAV modelling & simulation (J. Thermal Eng.)

Reviewed but out of scope (different flow regimes, not used):

- Steelant et al. (2015) — ATLLAS II Mach-5–6 hypersonic transport (EU FP7)
- Krishnan (2021) — Fluid-Thermal-Structural-Interaction seminar (ANSYS/NASA Ames)

---

## Notes

- **Scope:** a subsonic (M∞ = 0.20) high-lift wing-body angle-of-attack study. The two
  hypersonic / multiphysics documents above were reviewed but excluded as out of regime.
- The quantitative results in the report are the collated, workshop-grade dataset reproduced
  from the relevant studies (primary AoA table = HLPW-5 5v-grid result). To generate **fresh**
  numbers, supply a mesh/CAD and run the pipeline against a local ANSYS install.
- **Live execution note:** the full 3-D CRM-HL case exceeds the ANSYS Student 512 k-cell cap,
  so the locally executed case is the 2.5-D section. The section mesh is generated with gmsh and
  bridged to Fluent via CGNS, because Fluent Meshing's CAD-only region detection does not build a
  volume region from faceted STL on this build.
```

