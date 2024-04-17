#!/bin/bash

#!/bin/bash
#SBATCH -N 1
#SBATCH -J CLIMA_TEST_1
#SBATCH -o CLIMA_TEST_1_%j.out
#SBATCH -e CLIMA_TEST_1_%j.err
#SBATCH --mail-user=hengdi95@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --mem=1TB
#SBATCH --time=12:00:00
#SBATCH --exclusive

cd /nobackup/users/hengdi95/ClimaOceanBiogeochemistry/examples/
module add spack
module add spack-flat
module load openmpi/3.1.4-gcc-8.3.0-cuda-pmi-ucx
date
/home/software/julia/1.4.1/bin/julia --project ../CLIMA LES_NPZBD.jl
date