# -*- coding: utf-8 -*-
"""
crm_hl_geometry.py  --  STAGE 0: GEOMETRY.

Two ways to get geometry into the pipeline:
  (A) Production:  use the official watertight CRM-HL-WB CAD from the AIAA
      High-Lift Prediction Workshop, placed at geometry/crm_hl_wb.pmdb
      (.pmdb / .scdoc / .agdb / .stp).  Mesh with the Watertight workflow.
  (B) Self-contained dry-run (this script):  procedurally generate a faceted
      STL of the parametric CRM-HL-WB (each component a closed manifold) so the
      whole pipeline runs with NO external CAD.  Faceted input is meshed with the
      Fault-Tolerant (wrap) workflow -- see crm_hl_mesh.py.

This module also computes and prints the reference quantities used downstream.
Pure NumPy -- runnable immediately:  python crm_hl_geometry.py
"""
import numpy as np, os

OUT_STL = r"C:\aero-project\geometry\crm_hl_wb.stl"

# ---- parametric definition (full-scale CRM-HL-WB) --------------------------
B, SEMI = 58.76, 29.38
CR, CT  = 10.24, 2.82
MAC     = 7.005
SWEEP, DIH = np.radians(35.0), np.radians(6.0)
LF, DF  = 62.80, 5.90
XLE_R   = 24.50
HTSPAN, HTROOT, HTTIP = 22.0, 5.6, 1.7
VTH, VTROOT, VTTIP    = 9.2, 7.4, 2.6

def ref_quantities():
    SREF = 383.7
    AR   = B**2 / SREF
    taper= CT/CR
    return dict(S_ref=SREF, span=B, MAC=MAC, AR=AR, taper=taper,
                root_chord=CR, tip_chord=CT, sweep_deg=35.0, dihedral_deg=6.0,
                fus_len=LF, fus_dia=DF, slat_deg=30.0, flap_deg=37.0)

# ---- NACA 0012 closed section (loop of points, TE->upper->LE->lower->TE) ----
def naca0012_loop(n=40, chord=1.0):
    x = (1 - np.cos(np.linspace(0, np.pi, n)))/2
    yt = 5*0.12*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                 + 0.2843*x**3 - 0.1015*x**4)
    xu, yu = x[::-1], yt[::-1]            # TE -> LE upper
    xl, yl = x[1:],  -yt[1:]             # LE -> TE lower (skip dup LE)
    X = np.concatenate([xu, xl]) * chord
    Y = np.concatenate([yu, yl]) * chord
    return np.column_stack([X, Y])        # (m,2) closed loop

# ---- STL writer (ASCII, multi-solid) ---------------------------------------
def _tri(f, p1, p2, p3):
    n = np.cross(p2-p1, p3-p1); nn = np.linalg.norm(n)
    n = n/nn if nn > 1e-12 else np.array([0,0,1.0])
    f.write(f"  facet normal {n[0]:.5e} {n[1]:.5e} {n[2]:.5e}\n   outer loop\n")
    for p in (p1, p2, p3):
        f.write(f"    vertex {p[0]:.5e} {p[1]:.5e} {p[2]:.5e}\n")
    f.write("   endloop\n  endfacet\n")

def loft_closed(f, sections):
    """sections: list of (N,3) closed loops (same N). Triangulate side walls + end caps -> watertight."""
    S = [np.asarray(s) for s in sections]; N = len(S[0])
    for k in range(len(S)-1):              # side walls
        A, Bk = S[k], S[k+1]
        for i in range(N):
            j = (i+1) % N
            _tri(f, A[i], A[j], Bk[j]); _tri(f, A[i], Bk[j], Bk[i])
    for cap, flip in ((S[0], False), (S[-1], True)):   # end caps (fan)
        c = cap.mean(axis=0)
        for i in range(N):
            j = (i+1) % N
            if flip: _tri(f, c, cap[j], cap[i])
            else:    _tri(f, c, cap[i], cap[j])

def revolve_closed(f, xs, rs, nseg=24):
    """Body of revolution about x-axis: profile r(x). Watertight tube + apex caps."""
    rings = []
    th = np.linspace(0, 2*np.pi, nseg, endpoint=False)
    for x, r in zip(xs, rs):
        rings.append(np.column_stack([np.full(nseg, x), r*np.cos(th), r*np.sin(th)]))
    for k in range(len(rings)-1):
        A, Bk = rings[k], rings[k+1]
        for i in range(nseg):
            j = (i+1) % nseg
            _tri(f, A[i], A[j], Bk[j]); _tri(f, A[i], Bk[j], Bk[i])
    nose = np.array([xs[0], 0, 0]); tail = np.array([xs[-1], 0, 0])
    for i in range(nseg):                  # apex caps
        j = (i+1) % nseg
        _tri(f, nose, rings[0][j], rings[0][i])
        _tri(f, tail, rings[-1][i], rings[-1][j])

def lifting_surface(f, semi, c_root, c_tip, xle_root, sweep, dih, n_span=14, thick=True):
    secs = []
    for t in np.linspace(0, 1, n_span):
        y = t*semi
        c = c_root + (c_tip-c_root)*t
        xle = xle_root + y*np.tan(sweep)
        z = y*np.tan(dih)
        loop = naca0012_loop(chord=c) if thick else naca0012_loop(chord=c)*[1,0.4]
        pts = np.column_stack([xle + loop[:,0], np.full(len(loop), y), z + loop[:,1]])
        secs.append(pts)
    loft_closed(f, secs)                   # starboard
    secs_p = [s*[1,-1,1] for s in secs]
    loft_closed(f, secs_p)                 # port (mirror)

def generate_stl(path=OUT_STL):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        # fuselage (NACA-like area distribution)
        xs = np.linspace(0, LF, 40)
        rs = (DF/2)*np.sin(np.clip(xs/LF*np.pi, 1e-3, np.pi-1e-3))**0.55
        rs = np.maximum(rs, 0.18*DF); rs[0]=rs[-1]=0.02
        f.write("solid fuselage\n");  revolve_closed(f, xs, rs); f.write("endsolid fuselage\n")
        f.write("solid wing\n");      lifting_surface(f, SEMI, CR, CT, XLE_R, SWEEP, DIH); f.write("endsolid wing\n")
        f.write("solid h_tail\n");    lifting_surface(f, HTSPAN/2, HTROOT, HTTIP, LF-HTROOT-1, np.radians(30), 0); f.write("endsolid h_tail\n")
        # vertical tail: single lofted fin in x-z plane
        with_f = f
        secs=[]
        for t in np.linspace(0,1,8):
            c = VTROOT+(VTTIP-VTROOT)*t; z = DF*0.3 + t*VTH
            xle = LF-VTROOT-1 + t*VTH*np.tan(np.radians(40))
            loop = naca0012_loop(chord=c)
            pts = np.column_stack([xle+loop[:,0], 0.0*loop[:,1], z + 0*loop[:,0]])
            # build fin thickness in y
            pts = np.column_stack([xle+loop[:,0], loop[:,1], np.full(len(loop), z)])
            secs.append(pts)
        f.write("solid v_tail\n"); loft_closed(f, secs); f.write("endsolid v_tail\n")
    return path

if __name__ == "__main__":
    rq = ref_quantities()
    print("=== CRM-HL-WB reference quantities ===")
    for k, v in rq.items():
        print(f"  {k:12s} = {v:.4f}" if isinstance(v, float) else f"  {k:12s} = {v}")
    p = generate_stl()
    size = os.path.getsize(p)/1024
    # count facets
    with open(p) as fh: nf = sum(1 for ln in fh if ln.strip().startswith("facet"))
    print(f"\n[geometry] wrote {p}  ({size:.0f} KB, {nf} facets, 4 named solids)")
    print("[geometry] -> feed to crm_hl_mesh.py (Fault-Tolerant/wrap workflow for faceted STL)")
