#!/bin/bash
#SBATCH -p ampere  #-himem  #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-gpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --gres=gpu:1
#SBATCH --time 06:00:00
#SBATCH --job-name multiome
#SBATCH --output logs/jupyter_log-%j.txt


# Run before installing everything and everytime when using the environment!
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
# module load cuda/11.2
# module load cudnn/8.1_cuda-11.2

### get/set tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)

### print tunneling instructions to jupyter-log-{jobid}.txt
echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport $USER@$SLURM_SUBMIT_HOST.hpc.cam.ac.uk
    -----------------------------------------------------------------

    Then open a browser on your local machine to the following address
    ------------------------------------------------------------------
    localhost:$ipnport
    ------------------------------------------------------------------
    "

### start an ipcluster instance and launch jupyter server from my container
jupyter-lab --no-browser --port=$ipnport --ip=$ipnip

