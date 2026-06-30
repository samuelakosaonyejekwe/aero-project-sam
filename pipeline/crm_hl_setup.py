# -*- coding: utf-8 -*-
"""
crm_hl_setup.py  --  STAGE 2: SETUP (physics, materials, BCs, reference values).

configure(solver) mutates a live Fluent solver session (from the meshing hand-off
or from read_case) into the HLPW-5 Case-1 condition and returns the flow constants.
"""
import math

MACH, RE_MAC, T_REF = 0.20, 5.6e6, 289.4
GAMMA, R_AIR        = 1.4, 287.05
MAC, SREF_HALF      = 7.005, 191.85          # m, m^2 (half-model)
P_REF               = 101325.0

def flow_constants():
    a   = math.sqrt(GAMMA * R_AIR * T_REF)
    v   = MACH * a
    mu  = 1.716e-5 * (T_REF/273.15)**1.5 * (273.15+110.4)/(T_REF+110.4)
    rho = RE_MAC * mu / (v * MAC)
    return dict(a=a, vinf=v, mu=mu, rho=rho, mach=MACH, T=T_REF, p=P_REF,
                mac=MAC, sref=SREF_HALF)

def configure(solver, turbulence="spalart-allmaras"):
    fc = flow_constants()
    setup = solver.settings.setup

    # --- models -----------------------------------------------------------
    setup.models.energy.enabled = True
    setup.models.viscous.model = turbulence                 # "spalart-allmaras" | "k-omega"
    if turbulence == "k-omega":
        setup.models.viscous.k_omega_model = "sst"

    # --- material: ideal-gas air with Sutherland viscosity ----------------
    air = setup.materials.fluid["air"]
    air.density.option = "ideal-gas"
    air.viscosity.option = "sutherland"
    air.viscosity.sutherland.option = "three-coefficient-method"
    air.viscosity.sutherland.reference_viscosity   = 1.716e-5
    air.viscosity.sutherland.reference_temperature = 273.15
    air.viscosity.sutherland.effective_temperature = 110.4

    # --- operating & reference values -------------------------------------
    setup.general.operating_conditions.operating_pressure = 0.0   # absolute at far-field
    ref = setup.reference_values
    ref.area, ref.length        = fc["sref"], fc["mac"]
    ref.density, ref.velocity   = fc["rho"], fc["vinf"]
    ref.temperature, ref.pressure = fc["T"], fc["p"]

    # --- far-field boundary (AoA set later, per-angle, in the solve stage)-
    ff = setup.boundary_conditions.pressure_far_field["farfield"]
    ff.momentum.mach_number = fc["mach"]
    ff.momentum.flow_direction[0] = 1.0
    ff.momentum.flow_direction[1] = 0.0
    ff.momentum.flow_direction[2] = 0.0
    ff.thermal.temperature = fc["T"]
    if turbulence == "spalart-allmaras":
        ff.turbulence.turbulent_specification = "Turbulent Viscosity Ratio"
        ff.turbulence.turbulent_viscosity_ratio = 3.0
    else:
        ff.turbulence.turbulent_specification = "Intensity and Viscosity Ratio"
        ff.turbulence.turbulent_intensity = 0.02
        ff.turbulence.turbulent_viscosity_ratio = 10.0

    # --- numerics: pressure-based coupled, 2nd-order ----------------------
    methods = solver.settings.solution.methods
    methods.discretization_scheme = {"mom": "second-order-upwind",
                                     "k": "second-order-upwind",
                                     "epsilon": "second-order-upwind"}
    print(f"[setup] {turbulence}  V={fc['vinf']:.2f} m/s  rho={fc['rho']:.4f}  "
          f"mu={fc['mu']:.3e}  Re={RE_MAC:.2e}")
    return fc
