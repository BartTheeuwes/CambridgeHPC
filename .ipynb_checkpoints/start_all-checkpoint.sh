#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate env_all
#sbatch notebook_env.sh
sbatch notebook_env_cclake.sh
sbatch notebook_env_cclake_himem.sh
sbatch notebook_env_icelake.sh
sbatch notebook_env_icelake_himem.sh
sbatch notebook_env_sapphire.sh

job_id=$(squeue -u bt392 | awk '{print $1}' | sed -n 2p)

while  [ "$(squeue -u bt392 | awk '{print $5}' | sed -n 2p)" != "R" ]
do
  sleep 60
done

sleep 4
cat /rds/project/rds-SDzz0CATGms/users/bt392/logs/jupyter_log.txt | head -n 10

