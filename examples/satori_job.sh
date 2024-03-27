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
#SBATCH --time=4:00:00
#SBATCH --exclusive

#cd /nobackup/users/jars/projects/runscripts/
module add spack
module add spack-flat

julia --project=. LES_n1p2z2b2d5.jl 
