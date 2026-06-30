# -*- coding: utf-8 -*-
"""
crm_hl_mesh.py  --  STAGE 1: MESH (PyAnsys / Fluent Meshing).

build_mesh(geom, watertight) launches Fluent Meshing and runs either:
  * the WATERTIGHT GEOMETRY workflow  (clean CAD: .pmdb/.scdoc/.agdb/.stp), or
  * the FAULT-TOLERANT (wrap) workflow (faceted STL from Stage 0, possibly
    non-watertight at junctions),
then builds 15 prism boundary-layer cells sized to y+~=1, a poly-hexcore volume
mesh, checks quality, writes mesh/crm_hl_wb.msh.h5 and returns a live solver
session via switch_to_solver().

Run standalone:  python crm_hl_mesh.py            (uses generated STL, wrap path)
"""
import math, os
import ansys.fluent.core as pyfluent

# ---- reference condition -> principled first-cell height -------------------
MACH, RE_MAC, T_REF = 0.20, 5.6e6, 289.4
GAMMA, R_AIR, MAC   = 1.4, 287.05, 7.005
A_SOUND = math.sqrt(GAMMA * R_AIR * T_REF)
VINF    = MACH * A_SOUND
MU      = 1.716e-5 * (T_REF/273.15)**1.5 * (273.15+110.4)/(T_REF+110.4)
RHO     = RE_MAC * MU / (VINF * MAC)
N_PRISM, GROWTH = 15, 1.2
MESH_OUT = r"C:\aero-project\mesh\crm_hl_wb.msh.h5"
DEF_STL  = r"C:\aero-project\geometry\crm_hl_wb.stl"

def first_layer_height(yplus=1.0):
    cf    = 0.026 / RE_MAC**(1/7)                 # Schlichting turbulent flat-plate Cf
    tau_w = 0.5 * cf * RHO * VINF**2
    u_tau = math.sqrt(tau_w / RHO)
    return yplus * MU / (RHO * u_tau)

def _common_bl_and_volume(wf):
    """Prism boundary layers + poly-hexcore volume mesh (shared by both workflows)."""
    y1 = first_layer_height(1.0)
    print(f"[mesh] first-layer height for y+=1 : {y1:.3e} m")
    wf.add_boundary_layers.add_child_to_task()
    wf.add_boundary_layers.insert_compound_child_task()
    wf.add_boundary_layers.bl_control_name.set_state("prism_bl")
    wf.add_boundary_layers.offset_method_type.set_state("uniform")
    wf.add_boundary_layers.number_of_layers.set_state(N_PRISM)
    wf.add_boundary_layers.growth_rate.set_state(GROWTH)
    wf.add_boundary_layers.first_height.set_state(y1)
    wf.add_boundary_layers()
    wf.create_volume_mesh.volume_fill.set_state("poly-hexcore")
    wf.create_volume_mesh.volume_fill_controls.hex_max_cell_length.set_state(2.0)
    wf.create_volume_mesh()

def _watertight(meshing, geom):
    wf = meshing.watertight()
    wf.import_geometry.file_name.set_state(geom)
    wf.import_geometry.length_unit.set_state("m")
    wf.import_geometry()
    wf.add_local_sizing.add_child_to_task()
    wf.add_local_sizing.insert_compound_child_task()
    wf.add_local_sizing.boi_face_label_list.set_state(["slat", "flap"])
    wf.add_local_sizing.boi_control_name.set_state("highlift_edges")
    wf.add_local_sizing.boi_size.set_state(0.004)
    wf.add_local_sizing()
    smc = wf.create_surface_mesh.cfd_surface_mesh_controls
    smc.min_size.set_state(0.003); smc.max_size.set_state(2.0)
    smc.size_functions.set_state("Curvature & Proximity")
    smc.curvature_normal_angle.set_state(12)
    wf.create_surface_mesh()
    meshing.tui.mesh.check_mesh()
    wf.describe_geometry.setup_type.set_state("The geometry consists of only solid regions")
    wf.describe_geometry()
    wf.update_boundaries.selection_type.set_state("zone"); wf.update_boundaries()
    wf.update_regions()
    _common_bl_and_volume(wf)

def _fault_tolerant(meshing, geom):
    """Wrap workflow for faceted/dirty geometry (the Stage-0 STL)."""
    wf = meshing.fault_tolerant()
    pm = wf.import_cad_and_part_management
    pm.context.set_state(0)
    pm.fmd_file_name.set_state(geom)
    pm.length_unit.set_state("m")
    pm.file_loaded.set_state("yes")
    wf.import_cad_and_part_management()
    # leakage-free wrap around all named solids; create the external flow domain
    wf.create_regions();
    wf.update_regions()
    smc = wf.generate_the_surface_mesh.cfd_surface_mesh_controls
    smc.min_size.set_state(0.004); smc.max_size.set_state(2.0)
    wf.generate_the_surface_mesh()
    meshing.tui.mesh.check_mesh()
    _common_bl_and_volume(wf)

def build_mesh(geom=DEF_STL, watertight=False):
    os.makedirs(os.path.dirname(MESH_OUT), exist_ok=True)
    meshing = pyfluent.launch_fluent(mode="meshing", precision="double",
                                     processor_count=8, dimension=3, ui_mode="no_gui")
    (_watertight if watertight else _fault_tolerant)(meshing, geom)
    meshing.tui.mesh.check_mesh()
    meshing.tui.file.write_mesh(MESH_OUT)
    print(f"[mesh] written -> {MESH_OUT}")
    solver = meshing.switch_to_solver()
    print("[mesh] switched to solver session.")
    return solver

if __name__ == "__main__":
    s = build_mesh(DEF_STL, watertight=False)
    s.settings.file.write_case(file_name=r"C:\aero-project\mesh\crm_hl_wb.cas.h5")
    s.exit()
