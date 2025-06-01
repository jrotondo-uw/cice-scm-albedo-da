#!/bin/bash -l
#PBS -N post_restart
#PBS -A UOKL0045
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=1:00:00
#PBS -l job_priority=economy
#PBS -q main
#PBS -j oe
#PBS -m ae
#PBS -M jrotondo@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

# Load necessary modules and activate your conda environment
# Uncomment this if you're using a module-based Python environment
# module load ncarenv python

conda activate cice-scm-da

# Define variables
CASE=free_SibChuk_final_cat_test
USER=jrotondo

# Run post-processing script with stage 'restart'
python /glade/work/jrotondo/Projects/cice-scm-da/scripts/python/00_post_process.py $CASE $USER restart