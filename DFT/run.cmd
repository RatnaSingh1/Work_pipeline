#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 36

# set a partition
#SBATCH --partition normal

# set max wallclock time
#SBATCH --time=24:00:00

# set name of job
#SBATCH --job-name=test123

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output.dat

# send mail to this address
#SBATCH --mail-user=

# In the u0dawin queue, you will need the following line
source /etc/profile.d/modules.sh; source /etc/profile.d/modules_local.sh
module add palma/2019a
module add GCC/8.2.0-2.31.1 
module add OpenMPI/3.1.3
module add ORCA/4.1.2


/Applic.HPC/Easybuild/skylake/2019a/software/ORCA/4.1.2-gompi-2019a/orca penta.inp> penta.out

