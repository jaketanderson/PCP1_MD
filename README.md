# MD Simulations for PCP1


|Container|Contents|Description|
|----|----|----|
||`annotate_PCP1.py`|A script to be run from PyMOL after running `analysis/analysis.ipynb`. To use, open PyMOL from this directory and use `run annotate_PCP1.py` in the PyMOL command window. Then you can use the commands `color_res_by_RMSF(R_avg_path="R_avg.npy")` and `show_correlations(mindist=10.0, minmag=0.5, avg_DCC_path="avg_DCC_matrix.npy", avg_dist_path="avg_dist_matrix.npy", chain="A", keep_existing=False)`
|`analysis/`|`analysis.ipynb`|Notebook for performing analysis. Should be opened from with the `analysis/` directory as the current working path, where plots will also be deposited
|`input.tar.gz`  |`5u3h.cif`|PDB structure of apo-PCP1
||`5u3h_trimmed.pdb`|apo-PCP1 with the linker lengths matching Nathan Lee's construct
||`5u3h_prepped.pdb`|`tu3h_trimmed` but with proper terminal caps, ready for simulation
||`openmm_env.yml`|The conda environment file to be used as a starting point for installing the needed packages
||`protein.ff19SB.xml`|The AMBER ff19SB force field
||`tip3pfb.xml`|The tip3p-fb force field
||`run.sh`|Script which runs all simulations semi-automatically
||`simulate.py`|Python file for running an individual simulation (this is called by `run.sh`)
|`results_<date&time>.tar.gz`|`<i>/minimized.pdb`|Solvated, energy-minimized structure
||`<i>/minimized_state.xml`|OpenMM state post-minimization
||`<i>/system.xml`|OpenMM system post-solvation
||`<i>/err.txt`|stderr from worker
||`<i>/out.txt`|stdout from worker
||`<i>/production.log`|StateDataReporter output
||`<i>/production.dcd`|DCDReporter output
||`<i>/checkpoint.chk`|OpenMM checkpoint
||`bookkeeping/`|Contains many command outputs and copies of files obtained right before simulation

Note: to get `results_<date&time>.tar.gz`, check Jake's Lab_Reports folder on OneDrive.

