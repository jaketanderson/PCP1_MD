#!/bin/bash
set -euo pipefail

suffix="$(date +'%Y-%m-%d_%H-%M-%S')"

mkdir /workspace/bookkeeping

# Make copies of important files/data for bookkeeping
cp simulate.py "/workspace/bookkeeping/simulate.py_$suffix"
cp run.sh "/workspace/bookkeeping/run.sh_$suffix"
cp protein.ff19SB.xml "/workspace/bookkeeping/protein.ff19SB.xml_$suffix"
cp tip3pfb.xml "/workspace/bookkeeping/tip3pfb.xml_$suffix"
nvidia-smi > "/workspace/bookkeeping/nvidia-smi.txt_$suffix"
conda list > "/workspace/bookkeeping/conda_list.txt_$suffix"
python -m openmm.testInstallation > "/workspace/bookkeeping/testInstallation.txt_$suffix"

run_one () {
  local gpu="$1"
  echo "Launching GPU $gpu..."
  CUDA_VISIBLE_DEVICES="$gpu" python simulate.py "$gpu"
}

pids=()
for gpu in 0 1; do
  mkdir "/workspace/$gpu"
  run_one "$gpu" >"/workspace/${gpu}/out.txt" 2>"/workspace/${gpu}/err.txt" &
  pids+=("$!")
done

# Wait for all; fail fast if any fails
rc=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    rc=1
  fi
done

if (( rc != 0 )); then
  echo "At least one GPU job failed; aborting rest of pipeline."
  exit 1
fi

# Tar up all the results, including intermediate files
tar -I 'pigz -p 6' -cvf results.tar.gz /workspace/0/ /workspace/1/ /workspace/2/ /workspace/bookkeeping/ 

echo "Done!"
echo "End datetime: $(date +'%Y-%m-%d_%H-%M-%S')"
