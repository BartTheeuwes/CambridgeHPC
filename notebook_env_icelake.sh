#!/bin/bash
#SBATCH -p icelake  #-himem  #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 42
#SBATCH --time 10:00:00
#SBATCH --job-name multiome
#SBATCH --output logs/jupyter_log-%j.txt


# Leave this file here UN-MODDIFIED, if you want to make 
# use of it, create a copy for yourself in your own area
# so that you can edit it as required.
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

