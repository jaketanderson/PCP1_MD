# annotate_PCP1.py
#
# This file is a PyMOL script meant to show linear correlations and flexibilites of residues.

from pymol import cmd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def color_res_by_RMSF(R_avg_path="R_avg.npy"):
    R_avg = np.load(R_avg_path, allow_pickle=False)
    
    stored.resis = set()
    cmd.iterate("chain A", "stored.resis.add(int(resi))")

    # 0 to 12 hardcoded RMSF range
    norm = mcolors.Normalize(vmin=0, vmax=12)
    for pair in list(zip(sorted(stored.resis), R_avg)):
        resi, R = pair
        r, g, b, _ = cm.turbo(norm(R))
        rgb = [r, g, b]
        cname = f"c_resi_{resi}"
        cmd.set_color(cname, rgb)
        cmd.color(cname, f"chain A and resi {resi}")

def show_correlations(
    mindist=10.0,
    minmag=0.5,
    avg_DCC_path="avg_DCC_matrix.npy",
    avg_dist_path="avg_dist_matrix.npy",
    chain="A",
    keep_existing=False,
):
    if not keep_existing:
        cmd.delete("dist* lab*")
    avg_DCC_matrix = np.load(avg_DCC_path, allow_pickle=False)
    avg_dist_matrix = np.load(avg_dist_path, allow_pickle=False)

    filtered_DCC_matrix = np.where(avg_dist_matrix >= mindist, avg_DCC_matrix, np.nan)
    i_s, j_s = np.where(np.abs(filtered_DCC_matrix) >= minmag)

    norm = mcolors.Normalize(vmin=-1.0, vmax=1.0)

    # unique i<j pairs
    pairs = {(int(min(i, j)), int(max(i, j))) for i, j in zip(i_s, j_s)}

    for i, j in sorted(pairs):
        val = float(filtered_DCC_matrix[i, j])
        if np.isnan(val):
            continue

        sel1 = f"chain {chain} and resi {i+1} and name CA"
        sel2 = f"chain {chain} and resi {j+1} and name CA"

        # Skip if selections don't resolve to exactly one atom each
        if cmd.count_atoms(sel1) != 1 or cmd.count_atoms(sel2) != 1:
            continue

        distname = f"dist_{i+1}_{j+1}"
        cmd.distance(distname, sel1, sel2)

        # --- robust RGB coloring: define a named PyMOL color ---
        r, g, b, _ = cm.bwr(norm(val))
        rgb = [float(r), float(g), float(b)]
        color_name = f"col_{i+1}_{j+1}"
        cmd.set_color(color_name, rgb)
        cmd.set("dash_color", color_name, distname)
        cmd.set("dash_width", 4, distname)

        # --- custom label: hide distance labels, label a midpoint pseudoatom ---
        cmd.hide("labels", distname)

        c1 = cmd.get_atom_coords(sel1)
        c2 = cmd.get_atom_coords(sel2)
        mid = [(c1[k] + c2[k]) / 2.0 for k in range(3)]

        lab_obj = f"lab_{i+1}_{j+1}"
        cmd.pseudoatom(lab_obj, pos=mid)
        cmd.label(lab_obj, f'"{val:+.2f}"')  # fully custom text
        cmd.set("label_color", color_name, lab_obj)
        cmd.set("label_size", 16, lab_obj)

cmd.extend("color_res_by_RMSF", color_res_by_RMSF)
cmd.extend("show_correlations", show_correlations)
    
