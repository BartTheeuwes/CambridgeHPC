#!/bin/bash
#SBATCH -p skylake-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
### Modify this according to your Ray workload.
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
# SBATCH --ntasks=1
#SBATCH --time 22:00:00
### Modify this according to your Ray workload.
#SBATCH --cpus-per-task=30
#SBATCH --mem 370GB
# SBATCH --mem-per-cpu=5GB
#SBATCH --gpus-per-task=0
### outputs
#SBATCH --output logs/02-log-%J.txt
#SBATCH --error logs/02-err-%J.txt

# Getting the node names
nodes=$(scontrol show hostnames "$SLURM_JOB_NODELIST")
nodes_array=($nodes)

head_node=${nodes_array[0]}
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" hostname --ip-address)


port=6379
ip_head=$head_node_ip:$port
export ip_head
echo "IP Head: $ip_head"
echo -e "
        ssh -N -L $port:$head_node_ip:$port $USER@$SLURM_SUBMIT_HOST.hpc.cam.ac.uk
        "

echo "Starting HEAD at $head_node"
srun --nodes=1 --ntasks=1 -w "$head_node" --partition=skylake-himem   --account=gottgens-sl2-cpu --mem 370G \
    ray start --head --node-ip-address="$head_node_ip" --port=$port \
    --num-cpus "${SLURM_CPUS_PER_TASK}" --num-gpus "${SLURM_GPUS_PER_TASK}" --block &
    
# optional, though may be useful in certain versions of Ray < 1.0.
sleep 10

# number of nodes other than the head node
worker_num=$((SLURM_JOB_NUM_NODES - 1))

for ((i = 1; i <= worker_num; i++)); do
    node_i=${nodes_array[$i]}
    echo "Starting WORKER $i at $node_i"
    srun --nodes=1 --ntasks=1 -w "$node_i" --partition=skylake-himem  --account=gottgens-sl2-cpu  --mem 370G \
        ray start --address "$ip_head" \
        --num-cpus "${SLURM_CPUS_PER_TASK}" --num-gpus "${SLURM_GPUS_PER_TASK}" --block &
    sleep 5
done

python 02_TFs_to_genes.py