#!/bin/bash -e
#SBATCH --job-name=XBeach_Nukulaelae # job name (shows up in the queue)
#SBATCH --time=48:00:00 # Walltime (HH:MM:SS)
#SBATCH --mem-per-cpu=2048 # number of tasks (e.g. MPI)
#SBATCH --array=20000-31900%800

BIN=/nesi/project/spc03223/Killo/XBeach/trunk/executables_serial/bin
DATA=/nesi/nobackup/spc03223/TCAP/Forecast/Nukulaelae_forecast

shopt -s dotglob
shopt -s nullglob
subdir_paths=($DATA/t*/Run_*)
cd "${subdir_paths[${SLURM_ARRAY_TASK_ID}]}" && $BIN/xbeach
