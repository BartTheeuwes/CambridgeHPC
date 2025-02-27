#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate env_all

# Check if the number of cores was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <num_cores>"
  exit 1
fi

# Capture the number of cores from the first argument
NUM_CORES=$1

#sbatch notebook_env.sh
job1_id=$(sbatch notebook_env_cclake.sh --cpus-per-task=$NUM_CORES | awk '{print $4}')
job2_id=$(sbatch notebook_env_cclake_himem.sh --cpus-per-task=$NUM_CORES | awk '{print $4}')
job3_id=$(sbatch notebook_env_icelake.sh --cpus-per-task=$NUM_CORES | awk '{print $4}')
job4_id=$(sbatch notebook_env_icelake_himem.sh --cpus-per-task=$NUM_CORES | awk '{print $4}')
job5_id=$(sbatch notebook_env_sapphire.sh --cpus-per-task=$NUM_CORES | awk '{print $4}')

# Display job queue
squeue --me
sleep 10

# Sleep as long as there are no jobs running
while [ "$(squeue -u $USER -j $job1_id,$job2_id,$job3_id,$job4_id,$job5_id -t RUNNING --format='%.18i' | grep -v 'JOBID' | wc -l)" -lt 1 ]
do
    echo 'No jobs running yet'
    sleep 60
done

# Determine which jobs are and aren't running
running_jobs=$(squeue -u $USER -j $job1_id,$job2_id,$job3_id,$job4_id,$job5_id -t RUNNING --format='%.18i' | grep -v 'JOBID')
pd_jobs=$(squeue -u $USER -j $job1_id,$job2_id,$job3_id,$job4_id,$job5_id -t PD --format='%.18i' | grep -v 'JOBID')

# Output the log file for each running job
for job_id in $running_jobs; do
  echo $job_id
  cat /rds/project/rds-SDzz0CATGms/users/bt392/logs/jupyter_log-$job_id.txt | head -n 10
done

# scancel jobs that aren't running
for job_id in $pd_jobs; do
  scancel $job_id
done

# Display job queue
squeue --me