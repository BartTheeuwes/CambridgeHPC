#!/bin/bash
#SBATCH -p icelake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --time 4:00:00
#SBATCH --job-name jupyter_job
#SBATCH --output jupyter-log-%J.txt


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
    localhost:$ipnport  (prefix w/ https:// if using password)
    ------------------------------------------------------------------
    "

### start an ipcluster instance and launch jupyter server from my container

singularity exec /rds/project/rds-SDzz0CATGms/containers/rpy_v4_p3_fix1.sif  /usr/local/bin/jupyter lab --no-browser --port=$ipnport --ip=$ipnip

