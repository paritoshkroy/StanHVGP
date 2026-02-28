#!/bin/bash
#SBATCH --array=1-100  # 100 nodes runs this model independently
#SBATCH --account=def-aschmidt  # replace this with your own account
#SBATCH --ntasks=16              # number of processes
#SBATCH --mem-per-cpu=16000M      # memory; default unit is megabytes
#SBATCH --time=12:00:00         # time (HH:MM:SS)
#SBATCH --output=/home/pkroy/StanHVGP/SimulationStudy/%x-%j.out

# Modules
module load StdEnv/2023 r/4.4.0

# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f HVGPSimulationStudy10NeighPriorSet1.R --args $SLURM_ARRAY_TASK_ID
