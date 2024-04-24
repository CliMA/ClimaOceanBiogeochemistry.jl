#!/bin/bash
#SBATCH -N 1
#SBATCH -J CLIMA_TEST
#SBATCH -o CLIMA_TEST_%j.out
#SBATCH -e CLIMA_TEST_%j.err
#SBATCH --mail-user=hengdi95@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00
#SBATCH --mem=1TB

# Upload modules
module purge all
module add spack
module add cuda/11.4-5ncz6p4
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

# MPI specific exports
export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

# Julia specific enviromental variables
export COMMON="/nobackup/users/hengdi95/"
export JULIA_DEPOT_PATH="${COMMON}/depot"
export JULIA_CUDA_MEMORY_POOL=none
export JULIA="${COMMON}/julia/julia"

# Profile specific variable
export JULIA_NVTX_CALLBACKS=gc

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

export RX=1
export RY=4

srun --mpi=pmi2 ./launch.sh $JULIA --check-bounds=no --project hello_world.jl
