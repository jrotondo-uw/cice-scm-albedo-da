#!/bin/bash -l
### Job Name
#PBS -N run_assimilation
### Charging account
#PBS -A UOKL0045
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=4GB
### Allow job to run up to 30 minutes
#PBS -l walltime=11:00:00 
#PBS -l job_priority=economy
### Route the job to the casper queue
#PBS -q main
### Join output and error streams into single file
#PBS -j oe
#PBS -m ae
#PBS -M jrotondo@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
conda activate cice-scm-da

### Run analysis script
python ../python/04a_setup_da_case.py SIAL_CentralArctic_bounded_med_err_2011_em3 jrotondo spinup_CentralArctic_final uniform_ice atm
python ../python/04b_cycle.py SIAL_CentralArctic_bounded_med_err_2011_em3 jrotondo free_CentralArctic_final 3 bounded 2011 6 15 2011 10 1 null SAT_SEAICE_ALBEDODIRVIZ SAT_SEAICE_ALBEDODIRNIR SAT_SEAICE_ALBEDOINDVIZ SAT_SEAICE_ALBEDOINDNIR # SAT_SEAICE_AGREG_THICKNESS
python ../python/05_postprocess.py SIAL_CentralArctic_bounded_med_err_2011_em3 jrotondo 3 all