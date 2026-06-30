# -*- coding: utf-8 -*-
"""Fluent 3-D RANS on the 2.5-D BL-resolved mesh (symmetry span ends) -> accurate
locally-computed Cl/Cd/Cm vs AoA -> out/polar_live.csv."""
import math, csv, os, sys, traceback
import ansys.fluent.core as pyfluent
MSH=r"C:\aero-project\mesh\airfoil25d.cgns"
OUT_DIR=r"C:\aero-project\out"
AOA_LIST=[0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0]; ITERS=400
C=7.0; SPAN=0.05; SREF=C*SPAN; MRC=[1.75,0.0,0.0]
MACH,RE_MAC,T_REF=0.20,5.6e6,289.4; GAMMA,R_AIR=1.4,287.05
A=math.sqrt(GAMMA*R_AIR*T_REF); VINF=MACH*A
MU=1.716e-5*(T_REF/273.15)**1.5*(273.15+110.4)/(T_REF+110.4)
RHO=RE_MAC*MU/(VINF*C)
def log(m): print(m,flush=True)

def main():
    os.makedirs(OUT_DIR,exist_ok=True)
    log(f"[25] V={VINF:.1f} rho={RHO:.4f} Re={RE_MAC:.1e}")
    s=pyfluent.launch_fluent(mode="solver",dimension=3,processor_count=2,ui_mode="no_gui",start_timeout=300)
    log("[25] solver up")
    ok=False
    for desc,fn in [("tui.import.cgns.mesh",lambda:s.tui.file.import_.cgns.mesh(MSH)),
                    ("settings.read_mesh",  lambda:s.settings.file.read_mesh(file_name=MSH))]:
        try: fn(); ok=True; log(f"[25] import {desc} OK"); break
        except Exception as e: log(f"[25] import {desc}: {e}")
    if not ok: log("[25] import failed"); s.exit(); return
    su=s.settings.setup; bc=su.boundary_conditions
    def keys(g):
        try: return list(getattr(bc,g).keys())
        except Exception: return []
    walls=keys("wall"); log(f"[25] walls={walls}")
    su.models.energy.enabled=True
    su.models.viscous.model="spalart-allmaras"
    air=su.materials.fluid["air"]; air.density.option="ideal-gas"; air.viscosity.option="sutherland"
    ref=su.reference_values
    ref.area,ref.length,ref.density,ref.velocity=SREF,C,RHO,VINF
    ref.temperature,ref.pressure=T_REF,101325.0
    def area(z):
        for fn in (lambda:s.fields.reduction.area(locations=[z]),
                   lambda:s.scheme_eval.scheme_eval(f'(surface-area (get-thread "{z}"))')):
            try:
                v=fn()
                if v: return float(v)
            except Exception: pass
        return 0.0
    ar={z:area(z) for z in walls}
    log(f"[25] areas={ {k:round(v,1) for k,v in ar.items()} }")
    ranked=sorted(walls,key=lambda z:ar[z],reverse=True)
    sym=ranked[:2]; boxf=ranked[2:-1]; body=[ranked[-1]]       # 2 biggest=symmetry, smallest=wing
    setzt=s.settings.setup.boundary_conditions.set_zone_type
    def cx(z):                                                # mean x of a zone -> find outlet
        try:
            c=s.fields.reduction.centroid(locations=[z])
            return float(c[0]) if isinstance(c,(list,tuple)) else float(list(c.values())[0][0])
        except Exception: return 0.0
    cxs={z:cx(z) for z in boxf}; log(f"[25] box centroid-x={ {k:round(v,1) for k,v in cxs.items()} }")
    outlet=max(boxf,key=lambda z:cxs[z]); inlets=[z for z in boxf if z!=outlet]
    log(f"[25] symmetry={sym} inlets={inlets} outlet={outlet} wing={body}")
    for z in sym:
        try: setzt(zone_list=[z],new_type="symmetry")
        except Exception as e: log(f"[25] sym {z}: {e}")
    for z in inlets:
        try: setzt(zone_list=[z],new_type="velocity-inlet")
        except Exception as e: log(f"[25] vin {z}: {e}")
    try: setzt(zone_list=[outlet],new_type="pressure-outlet")
    except Exception as e: log(f"[25] pout: {e}")
    # ---- INCOMPRESSIBLE physics (robust at M=0.20) ----
    su.models.energy.enabled=False
    try: air.density.option="constant"; air.density.value=RHO
    except Exception as e: log(f"[25] density: {e}")
    try: air.viscosity.option="constant"; air.viscosity.value=MU
    except Exception as e: log(f"[25] visc: {e}")
    vin=keys("velocity_inlet"); pout=keys("pressure_outlet")
    log(f"[25] velocity_inlet={vin} pressure_outlet={pout}")
    if not vin or not body: log("[25] BC setup failed"); s.exit(); return
    for z in pout:
        try: bc.pressure_outlet[z].momentum.gauge_pressure=0.0
        except Exception: pass
        try: bc.pressure_outlet[z].turbulence.turbulent_viscosity_ratio=5
        except Exception: pass

    def rv(name):
        r=s.settings.solution.report_definitions.compute(report_defs=[name])
        try:
            v=list(r.values())[0] if isinstance(r,dict) else list(r[0].values())[0]
            return float(v[0]) if isinstance(v,(list,tuple)) else float(v)
        except Exception: return float("nan")
    rd=s.settings.solution.report_definitions
    meth=s.settings.solution.methods
    def setord(o):
        try: meth.discretization_scheme={"mom":o,"modified-turbulent-viscosity":o}
        except Exception as e: log(f"[25] order {o}: {e}")
    _d=[False]
    def set_aoa(alpha):
        a=math.radians(alpha)
        for z in vin:
            m=bc.velocity_inlet[z].momentum
            if not _d[0]:
                log("[25] vinlet attrs: "+", ".join(x for x in dir(m) if not x.startswith('_'))); _d[0]=True
            try: m.velocity_specification_method="Magnitude and Direction"
            except Exception: pass
            for setter in (lambda: setattr(m,"velocity_magnitude",VINF),
                           lambda: setattr(m.velocity_magnitude,"value",VINF),
                           lambda: setattr(m,"velocity",VINF)):
                try: setter(); break
                except Exception: pass
            try:
                m.flow_direction[0]=math.cos(a); m.flow_direction[1]=math.sin(a); m.flow_direction[2]=0.0
            except Exception: pass
            try:
                t=bc.velocity_inlet[z].turbulence
                t.turbulent_specification="Turbulent Viscosity Ratio"; t.turbulent_viscosity_ratio=5
            except Exception: pass
        rd.lift["cl"]=dict(zones=body,force_vector=[-math.sin(a),math.cos(a),0])
        rd.drag["cd"]=dict(zones=body,force_vector=[ math.cos(a),math.sin(a),0])
        rd.moment["cm"]=dict(zones=body,mom_center=MRC,mom_axis=[0,0,1])

    setord("first-order-upwind")
    urf=s.settings.solution.controls.under_relaxation
    for k,v in {"pressure":0.3,"mom":0.5,"density":0.5,"body-force":0.5,"nut":0.5,"turb-viscosity":0.5}.items():
        try: urf[k]=v
        except Exception: pass
    set_aoa(0.0)
    s.settings.solution.initialization.hybrid_initialize()
    s.settings.solution.run_calculation.iterate(iter_count=500); log("[25] first-order base established")
    # warm-start into SECOND-ORDER from the converged first-order field -> sharper drag
    setord("second-order-upwind")
    s.settings.solution.run_calculation.iterate(iter_count=600); log("[25] second-order base established")
    results=[]
    for alpha in AOA_LIST:
        set_aoa(alpha)
        s.settings.solution.run_calculation.iterate(iter_count=ITERS)
        cl,cd,cm=rv("cl"),rv("cd"),rv("cm")
        results.append((alpha,cl,cd,cl/cd if cd else 0,cm))
        log(f"[25] a={alpha:5.1f}  Cl={cl:8.4f}  Cd={cd:8.5f}  L/D={cl/cd:7.2f}  Cm={cm:8.4f}")
        with open(fr"{OUT_DIR}\polar_live.csv","w",newline="") as fh:
            w=csv.writer(fh); w.writerow(["alpha_deg","Cl","Cd","L_over_D","Cm"]); w.writerows(results)
    log(f"[25] DONE -> {OUT_DIR}\\polar_live.csv"); s.exit()

if __name__=="__main__":
    try: main()
    except Exception as e:
        log("[25] ABORTED: "+repr(e)); traceback.print_exc(); sys.exit(1)
