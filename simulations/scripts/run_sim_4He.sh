#!/bin/bash
#SBATCH --nodes=1				# 1 node
#SBATCH --mem=8G              			# 8G of Ram for Job
#SBATCH --ntasks-per-node=32     		# 32 tasks per node
#SBATCH --time=24:00:00       			# time limits: 24 hours
#SBATCH --error=logs/myJob_4He.err     		# standard error file
#SBATCH --output=logs/myJob_4He.out	     	# standard output file
#SBATCH --partition=lprod       		# partition name
#SBATCH --mail-type=END       			# type of event notification
#SBATCH --mail-user=luciana.dourado@gssi.it  	# mail address

# Activate Miniconda environment
source ~/miniconda3/bin/activate crpropa3

# Run the Python script
python sim.py 4 2
