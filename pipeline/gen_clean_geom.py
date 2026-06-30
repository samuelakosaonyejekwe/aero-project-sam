# -*- coding: utf-8 -*-
"""GUARANTEED-watertight wing (extruded NACA-0012) + far-field box via trimesh.
Extrusion -> clean closed manifold. Two bodies (wing inside box) for region detection."""
import numpy as np, trimesh
from shapely.geometry import Polygon
OUT=r"C:\aero-project\geometry\wing_box_clean.stl"
CHORD=7.0; SPAN=58.0
BX=(-30.0,70.0); BY=(-70.0,70.0); BZ=(-40.0,40.0)

def airfoil_poly(n=120,chord=CHORD,te_gap=0.004):
    half=n//2
    beta=np.linspace(0,np.pi,half+1); x=(1-np.cos(beta))/2
    yt=5*0.12*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)+te_gap*(x)
    xu,yu=x[::-1],yt[::-1]            # TE->LE upper
    xl,yl=x[1:],-yt[1:][::-1]         # ... build clean ring
    xl=x[1:-1]; yl=-yt[1:-1]
    X=np.concatenate([xu,xl])*chord; Y=np.concatenate([yu,yl])*chord
    return Polygon(np.column_stack([X,Y]))

def main():
    poly=airfoil_poly()
    wing=trimesh.creation.extrude_polygon(poly, height=SPAN)   # extrude along +z
    # rotate so the extrusion (span) lies along y, and centre the span
    wing.apply_transform(trimesh.transformations.rotation_matrix(-np.pi/2,[1,0,0]))
    wing.apply_translation([0,-SPAN/2,0])  # centre span at y=0  (y: -SPAN/2 .. +SPAN/2)
    box=trimesh.creation.box(bounds=[[BX[0],BY[0],BZ[0]],[BX[1],BY[1],BZ[1]]])
    print(f"[geom] wing watertight={wing.is_watertight} faces={len(wing.faces)} vol={wing.volume:.2f} "
          f"bounds={np.round(wing.bounds,2).tolist()}")
    print(f"[geom] box  watertight={box.is_watertight}")
    combined=trimesh.util.concatenate([wing,box]); combined.export(OUT)
    comps=combined.split(only_watertight=False)
    print(f"[geom] exported {OUT}: {len(combined.faces)} faces, {len(comps)} bodies, "
          f"watertight={[bool(c.is_watertight) for c in comps]}")

if __name__=="__main__":
    main()
