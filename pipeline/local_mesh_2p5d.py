# -*- coding: utf-8 -*-
"""2.5-D mesh: 2-D BL-resolved NACA-0012 section extruded 1 cell in span (symmetry ends).
Keeps y+~1 accuracy AND uses the proven 3-D CGNS import path. Runs in-process."""
import numpy as np, gmsh, os, math
C=7.0; SPAN_T=0.05; FAR=40*C    # thin span -> reasonable cell aspect ratios
RE_MAC,T_REF,MACH=5.6e6,289.4,0.20; GAMMA,R=1.4,287.05
VINF=MACH*math.sqrt(GAMMA*R*T_REF)
MU=1.716e-5*(T_REF/273.15)**1.5*(273.15+110.4)/(T_REF+110.4)
RHO=RE_MAC*MU/(VINF*C)
OUT=r"C:\aero-project\mesh\airfoil25d.cgns"
def y_first(yp=1.0):
    cf=0.026/RE_MAC**(1/7); tw=0.5*cf*RHO*VINF**2; ut=math.sqrt(tw/RHO); return yp*MU/(RHO*ut)
def naca(n=160):
    half=n//2; b=np.linspace(0,np.pi,half+1); x=(1-np.cos(b))/2
    yt=5*0.12*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
    xu,yu=x[::-1],yt[::-1]; xl,yl=x[1:-1],-yt[1:-1]
    return np.concatenate([xu,xl])*C, np.concatenate([yu,yl])*C

def main(surf=0.05):
    os.makedirs(os.path.dirname(OUT),exist_ok=True)
    y1=y_first(1.0); print(f"[25d] y1(y+=1)={y1:.2e} m  V={VINF:.1f} rho={RHO:.4f}")
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("af"); occ=gmsh.model.occ
    X,Y=naca(); pts=[occ.addPoint(X[i],Y[i],0,surf) for i in range(len(X))]
    spl=occ.addBSpline(pts+[pts[0]]); afloop=occ.addCurveLoop([spl])
    box=occ.addRectangle(-FAR,-FAR,0,2*FAR,2*FAR)
    base,_=occ.cut([(2,box)],[(2,occ.addPlaneSurface([afloop]))]); occ.synchronize()
    base_surf=base[0][1]
    # airfoil curves (near origin) vs box curves
    af_curves=[]
    for (d,t) in gmsh.model.getEntities(1):
        x0,y0,_,x1,y1b,_=gmsh.model.getBoundingBox(d,t)
        if max(abs(x0),abs(x1),abs(y0),abs(y1b))<FAR-1: af_curves.append(t)
    # boundary layer (2-D, on airfoil)
    bl=gmsh.model.mesh.field.add('BoundaryLayer')
    gmsh.model.mesh.field.setNumbers(bl,'CurvesList',af_curves)
    gmsh.model.mesh.field.setNumber(bl,'Size',y1); gmsh.model.mesh.field.setNumber(bl,'Ratio',1.18)
    gmsh.model.mesh.field.setNumber(bl,'Thickness',0.25*C); gmsh.model.mesh.field.setNumber(bl,'Quads',1)
    try: gmsh.model.mesh.field.setNumbers(bl,'FanPointsList',[pts[0]])
    except Exception: pass
    gmsh.model.mesh.field.setAsBoundaryLayer(bl)
    dd=gmsh.model.mesh.field.add('Distance'); gmsh.model.mesh.field.setNumbers(dd,'CurvesList',af_curves)
    gmsh.model.mesh.field.setNumber(dd,'Sampling',200)
    th=gmsh.model.mesh.field.add('Threshold'); gmsh.model.mesh.field.setNumber(th,'InField',dd)
    gmsh.model.mesh.field.setNumber(th,'SizeMin',surf); gmsh.model.mesh.field.setNumber(th,'SizeMax',FAR/7)
    gmsh.model.mesh.field.setNumber(th,'DistMin',0.5*C); gmsh.model.mesh.field.setNumber(th,'DistMax',8*C)
    gmsh.model.mesh.field.setAsBackgroundMesh(th)
    for o in ("Mesh.MeshSizeExtendFromBoundary","Mesh.MeshSizeFromPoints","Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o,0)
    # extrude 1 layer in z -> 2.5D
    ext=occ.extrude([(2,base_surf)],0,0,SPAN_T,numElements=[1],recombine=True); occ.synchronize()
    vol=[e for e in ext if e[0]==3]
    # classify all surfaces of the 3D body
    eps=1e-4; airfoil=[]; far=[]; sym=[]
    for (d,t) in gmsh.model.getEntities(2):
        x0,y0,z0,x1,y1b,z1=gmsh.model.getBoundingBox(d,t)
        if abs(z1-z0)<eps:                         # constant-z plane -> symmetry end
            sym.append(t)
        elif max(abs(x0),abs(x1),abs(y0),abs(y1b))>FAR-1:
            far.append(t)
        else:
            airfoil.append(t)
    gmsh.model.addPhysicalGroup(2,airfoil,name="wing")
    gmsh.model.addPhysicalGroup(2,far,name="farfield")
    gmsh.model.addPhysicalGroup(2,sym,name="symmetry")
    gmsh.model.addPhysicalGroup(3,[t for (d,t) in vol],name="fluid")
    print(f"[25d] wing surf={len(airfoil)} far={len(far)} sym={len(sym)}")
    gmsh.model.mesh.generate(3)
    nh=len(gmsh.model.mesh.getElementsByType(5)[0]) if gmsh.model.mesh.getElementsByType(5)[0].size else 0
    npri=len(gmsh.model.mesh.getElementsByType(6)[0]) if gmsh.model.mesh.getElementsByType(6)[0].size else 0
    ntet=len(gmsh.model.mesh.getElementsByType(4)[0]) if gmsh.model.mesh.getElementsByType(4)[0].size else 0
    tot=nh+npri+ntet; print(f"[25d] hex={nh} prism={npri} tet={ntet} total={tot}")
    gmsh.write(OUT); gmsh.finalize(); print(f"[25d] wrote {OUT}")
    return tot

if __name__=="__main__":
    n=main(); print("OVER CAP" if n>500000 else "OK (<512k)")
