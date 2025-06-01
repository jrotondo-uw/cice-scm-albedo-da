###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import xarray as xr
import numpy as np
import sys
import shutil

from datetime import datetime, timedelta
from time import sleep, perf_counter
import pandas as pd
###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the assimilation branch case name
case = sys.argv[1] # SIT_f1_NORM_test
user = sys.argv[2]

# set the source free case name
source_case = sys.argv[3] # 'free_test'
truth_member = sys.argv[4]

# set the relevant directories 
# where is this project repo located?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# wherer is the DART installation on your system?
dart_dir = '/glade/work/'+user+'/DART/'
# where are the observations stored on your system?
obs_dir = project_dir + '/data/processed/synthetic_obs_med_err/'+source_case+'/'+truth_member+'/'
# where is the scratch directory for this user?
scratch_dir = '/glade/derecho/scratch/'+user+'/'

# set the assimilation parameters
# filter_kind = int(sys.argv[4]) #1
# regression_kind = sys.argv[5] #'NORM'
qceff_type = sys.argv[5]
postprocessing = 'cice'

# what observations are being assimilated?
if len(sys.argv) > 13:
    obs_type = [sys.argv[i] for i in range(13, len(sys.argv))]
else:
    print('No observations called! Please choose an observation type to assimilate!')
    sys.exit()

if sys.argv[12] == 'cat':
    cat_assim = True
else:
    cat_assim = False
 
# set the assimilation dates
first_assim_date = datetime(int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]))
end_assim_date = datetime(int(sys.argv[9]),int(sys.argv[10]),int(sys.argv[11]) )

###############################################################
## DEFINE HELPER FUNCTIONS                                   ##
###############################################################
def check_completion(file_string): 
    txt = open(file_string).readlines()
    if 'icepack' in file_string:
        check = 'ICEPACK COMPLETED SUCCESSFULLY\n'
    else:
        check = ' Finished ... at YYYY MM DD HH MM SS = \n'

    if  check not in txt:
      print('Process did not finish correctly')
      status = -1
    else:
    #   os.system('rm '+file_string)
      status = 0

    return status

def unpack_categories(file):

    data_vars = ['alvdrn', 'alidrn', 'alvdfn', 'alidfn']
    ncat = 5
    encd = {'alvdr0001': {'_FillValue': None},
            'alvdr0002': {'_FillValue': None},
            'alvdr0003': {'_FillValue': None},
            'alvdr0004': {'_FillValue': None},
            'alvdr0005': {'_FillValue': None},
            'alidr0001': {'_FillValue': None},
            'alidr0002': {'_FillValue': None},
            'alidr0003': {'_FillValue': None},
            'alidr0004': {'_FillValue': None},
            'alidr0005': {'_FillValue': None},
            'alvdf0001': {'_FillValue': None},
            'alvdf0002': {'_FillValue': None},
            'alvdf0003': {'_FillValue': None},
            'alvdf0004': {'_FillValue': None},
            'alvdf0005': {'_FillValue': None},
            'alidf0001': {'_FillValue': None},
            'alidf0002': {'_FillValue': None},
            'alidf0003': {'_FillValue': None},
            'alidf0004': {'_FillValue': None},
            'alidf0005': {'_FillValue': None}}

    ds = xr.open_dataset(file)

    os.rename(file, file+'~')
    #ds = ds[data_vars]

    for n in range(1,ncat+1):
        ds['alvdr000'+'{:01}'.format(n)]=ds.alvdrn.isel(ncat=n-1)
        ds['alidr000'+'{:01}'.format(n)]=ds.alidrn.isel(ncat=n-1)
        ds['alvdf000'+'{:01}'.format(n)]=ds.alvdfn.isel(ncat=n-1)
        ds['alidf000'+'{:01}'.format(n)]=ds.alidfn.isel(ncat=n-1)
        
    #ds=ds.drop(data_vars)

    ds.to_netcdf(file, encoding=encd)  
    
    return len(ds.variables)

def repack_categories(file):

    ncat = 5
    expanded_vars=[]
    for n in range(1,ncat+1):
        expanded_vars.append('alvdr000'+'{:01}'.format(n))
        expanded_vars.append('alidr000'+'{:01}'.format(n))
        expanded_vars.append('alvdf000'+'{:01}'.format(n))
        expanded_vars.append('alidf000'+'{:01}'.format(n))

    encd = {'alvdrn': {'_FillValue': None},
            'alidrn': {'_FillValue': None},
            'alvdfn': {'_FillValue': None},
            'alidfn': {'_FillValue': None}}
            
    ds = xr.open_dataset(file)
    os.rename(file, file+'~')
    alvdrn=np.zeros((ncat,4))
    alidrn=np.zeros((ncat,4))
    alvdfn=np.zeros((ncat,4))
    alidfn=np.zeros((ncat,4))
    
    for n in range(1,ncat+1):
        alvdrn[n-1,:]=ds['alvdr000'+'{:01}'.format(n)].values
        alidrn[n-1,:]=ds['alidr000'+'{:01}'.format(n)].values
        alvdfn[n-1,:]=ds['alvdf000'+'{:01}'.format(n)].values
        alidfn[n-1,:]=ds['alidf000'+'{:01}'.format(n)].values
    
    ds['alvdrn']=xr.DataArray(data=alvdrn, dims=['ncat','ni'])
    ds['alidrn']=xr.DataArray(data=alidrn, dims=['ncat','ni'])
    ds['alvdfn']=xr.DataArray(data=alvdfn, dims=['ncat','ni'])
    ds['alidfn']=xr.DataArray(data=alidfn, dims=['ncat','ni'])

    ds=ds.drop(expanded_vars)
    ds.to_netcdf(file, mode= 'w', encoding=encd)

    os.remove(file+'~')

    return len(ds.variables)

def cycle_da(case, source_case, dart_dir, qceff_type, assim_date, truth_member, obs_type, cat_obs = False):

    case_dir  = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case+'/'
    assim_dir = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/assim_dir/'+case+'/'

    if cat_obs is False:
        obs_subdir = obs_dir + '/aggregate/'
        state_variables = ['alvdr', 'QTY_SEAICE_AGREG_ALBEDODIRVIZ'   , 'UPDATE',
                           'alidr', 'QTY_SEAICE_AGREG_ALBEDODIRNIR'   , 'UPDATE',
                           'alvdf', 'QTY_SEAICE_AGREG_ALBEDOINDVIZ'   , 'UPDATE',
                           'alidf', 'QTY_SEAICE_AGREG_ALBEDOINDNIR'   , 'UPDATE',
                           'aicen', 'QTY_SEAICE_CONCENTR'       , 'NO_COPY_BACK',
                           'vicen', 'QTY_SEAICE_VOLUME'         , 'UPDATE',
                           'vsnon', 'QTY_SEAICE_SNOWVOLUME'     , 'NO_COPY_BACK'] 
    else:
        obs_subdir = obs_dir + '/itd/'
        state_variables = ['aicen', 'QTY_SEAICE_CONCENTR'       , 'UPDATE',
                           'vicen', 'QTY_SEAICE_VOLUME'         , 'UPDATE',
                           'vsnon', 'QTY_SEAICE_SNOWVOLUME'     , 'UPDATE',
                           'alvdr0001', 'QTY_SEAICE_ALVDR0001', 'UPDATE',
                           'alvdr0002', 'QTY_SEAICE_ALVDR0002', 'UPDATE',
                           'alvdr0003', 'QTY_SEAICE_ALVDR0003', 'UPDATE',
                           'alvdr0004', 'QTY_SEAICE_ALVDR0004', 'UPDATE',
                           'alvdr0005', 'QTY_SEAICE_ALVDR0005', 'UPDATE',
                           'alidr0001', 'QTY_SEAICE_ALIDR0001', 'UPDATE',
                           'alidr0002', 'QTY_SEAICE_ALIDR0002', 'UPDATE',
                           'alidr0003', 'QTY_SEAICE_ALIDR0003', 'UPDATE',
                           'alidr0004', 'QTY_SEAICE_ALIDR0004', 'UPDATE',
                           'alidr0005', 'QTY_SEAICE_ALIDR0005', 'UPDATE',
                           'alvdf0001', 'QTY_SEAICE_ALVDF0001', 'UPDATE',
                           'alvdf0002', 'QTY_SEAICE_ALVDF0002', 'UPDATE',
                           'alvdf0003', 'QTY_SEAICE_ALVDF0003', 'UPDATE',
                           'alvdf0004', 'QTY_SEAICE_ALVDF0004', 'UPDATE',
                           'alvdf0005', 'QTY_SEAICE_ALVDF0005', 'UPDATE',
                           'alidf0001', 'QTY_SEAICE_ALIDF0001', 'UPDATE',
                           'alidf0002', 'QTY_SEAICE_ALIDF0002', 'UPDATE',
                           'alidf0003', 'QTY_SEAICE_ALIDF0003', 'UPDATE',
                           'alidf0004', 'QTY_SEAICE_ALIDF0004', 'UPDATE',
                           'alidf0005', 'QTY_SEAICE_ALIDF0005', 'UPDATE']

    # Get date information
    year = assim_date.year
    mon = assim_date.month
    day = assim_date.day
    date_str = '{0}'.format('%04d'%year) + '-{0}'.format('%02d'%mon) + '-{0}'.format('%02d'%day)
    next_day = assim_date + timedelta(days = 1)
    next_day_str = '{0}'.format('%04d'%next_day.year) + '-{0}'.format('%02d'%next_day.month) + '-{0}'.format('%02d'%next_day.day) 

    # Move to assimilation directory
    if os.path.exists(assim_dir) is False:
        os.makedirs(assim_dir)
    else: 
        print('Assimilation directory already exists!')      
    os.chdir(assim_dir)
    shutil.copy(dart_dir + '/models/cice-scm/work/filter', '.')
    shutil.copy(dart_dir+'/models/cice-scm/scripts/templates/DART_input.nml.template', 'input.nml')
    shutil.copy(dart_dir+'/models/cice-scm/scripts/templates/cross_correlations.txt', 'cross_correlations.txt')
    # shutil.copy(dart_dir + '/models/cice-scm2/control_impact_runtime.txt', 'control_impact_runtime.txt')  
    if os.path.exists(dart_dir + '/models/cice-scm/work/cice.r.nc') is False:   
        os.symlink(scratch_dir + '/ICEPACK_RUNS/'+source_case+'/mem0001/restart/iced.2011-01-05-00000.nc', dart_dir + '/models/cice-scm/work/cice.r.nc')
    else:
        print('dummy cice restart file already exists! Not linking another!')
    shutil.copy(dart_dir + '/models/cice-scm/work/cice.r.nc', '.')
    shutil.copy(dart_dir + '/models/cice-scm/work/dart_to_cice', '.')
    shutil.copy(dart_dir + '/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc', '.')
    restarts_in_list = open(assim_dir + 'cice_restarts_in.txt', 'w')
    restarts_out_list = open(assim_dir + 'cice_restarts_out.txt', 'w')

    # if necessary, unpack cice.r.nc
    if cat_obs is True:
        null = unpack_categories('cice.r.nc')

    # Set up directory
    ensemble_size = len(glob.glob('/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+source_case+'/mem*'))
    mem = 1
    while mem <= ensemble_size:
    
        # make a subdirectory for each ensemble member
        inst_string = '{0}'.format('%04d'%mem)
        if os.path.exists(assim_dir + '/mem'+inst_string+'/') is False: 
            os.makedirs(assim_dir + '/mem'+inst_string+'/')

        os.chdir(assim_dir + '/mem'+inst_string+'/')

        if mem == int(truth_member):
            comd = 'cp ' + case_dir + '/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc truth_data.nc'
        else:
            restarts_in_list.writelines('mem'+inst_string+'/pre_filter_restart.nc \n')
            restarts_out_list.writelines('mem'+inst_string+'/post_filter_restart.nc \n')
            comd = 'cp ' + case_dir + '/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc pre_filter_restart.nc'

        os.system(comd)

        # unpack categories if necessary
        if cat_obs is True:
            if mem != int(truth_member):
                null = unpack_categories('pre_filter_restart.nc')
        
        os.chdir('../')
        mem += 1

    restarts_in_list.close()
    restarts_out_list.close()
    
    symlink_path = 'obs_seq.out'
    obs_file = obs_subdir + '/obs_seq.{0}{1}{2}'.format('%04d' % year, '%02d' % mon, '%02d' % day)

    # Remove existing symlink if necessary
    if os.path.exists(symlink_path) or os.path.islink(symlink_path):
        os.remove(symlink_path)

    # Create the symlink
    os.symlink(obs_file, symlink_path)

    # Set filter namelist settings 
    nml_file = f90nml.read(dart_dir+'/models/cice-scm/scripts/templates/DART_input.nml.template')
    nml_file['filter_nml']['ens_size'] = ensemble_size-1
    nml_file['filter_nml']['num_output_state_members'] = ensemble_size-1
    nml_file['filter_nml']['num_output_obs_members'] = ensemble_size-1
    if date_str == '2011-04-01':
        nml_file['filter_nml']['inf_initial_from_restart'] = [False, False]
        nml_file['filter_nml']['inf_sd_initial_from_restart'] = [False, False]
    else:
        nml_file['filter_nml']['inf_initial_from_restart'] = [True, False]
        nml_file['filter_nml']['inf_sd_initial_from_restart'] = [True, False]
    nml_file['obs_kind_nml']['assimilate_these_obs_types'] = obs_type
    nml_file['algorithm_info_nml']['qceff_table_filename'] = '/glade/work/jrotondo/DART/models/cice-scm/scripts/templates/cice_bounded_qceff_table.csv'
    nml_file['assim_tools_nml']['cutoff'] = 100000.0
    nml_file['model_nml']['model_state_variables'] = state_variables
    nml_file['dart_to_cice_nml']['postprocess'] = postprocessing
    nml_file['obs_impact_tool_nml']['input_filename'] = 'cross_correlations.txt'
    nml_file['obs_impact_tool_nml']['output_filename'] = 'control_impact_runtime.txt'
    nml_file.write('input.nml', force=True)

    # Run the filter
    print('Running the filter for '+date_str+'...')
    comd = './filter > output.filter'
    os.system(comd)
    filt_status = check_completion('output.filter')
    if filt_status < 0:
        print('filter did not complete correctly. Process stopping.')
        sys.exit()

    # Move filter output to case directory
    os.makedirs(case_dir + '/analyses/'+date_str+'/')
    os.makedirs(case_dir + '/forecasts/'+date_str+'/')
    os.makedirs(case_dir + '/output_files/'+date_str+'/')
    shutil.move('output.filter', case_dir + '/output_files/'+date_str+'/')
    # shutil.move('regression_data.txt', case_dir + '/output_files/'+date_str+'/')
    # shutil.move('obs_inc_data.txt', case_dir + '/output_files/'+date_str+'/')
    os.remove('obs_seq.out')

    comd = 'mv input_*.nc ' + case_dir+ '/output_files/'+date_str+'/'
    os.system(comd)
    sleep(2)

    # rename output_priorinf_[mean, sd].nc to input
    os.system('cp output_priorinf_mean.nc input_priorinf_mean.nc')
    os.system('cp output_priorinf_sd.nc input_priorinf_sd.nc')

    comd = 'mv preassim_*.nc ' + case_dir+ '/output_files/'+date_str+'/'
    os.system(comd)
    sleep(2)
    comd = 'mv analysis_*.nc ' + case_dir+ '/output_files/'+date_str+'/'
    os.system(comd)
    sleep(2)
    comd = 'mv output_*.nc ' + case_dir+ '/output_files/'+date_str+'/'
    os.system(comd)
    sleep(2)
    comd = 'mv *_forward_ope_errors* ' + case_dir+ '/output_files/'+date_str+'/'
    os.system(comd)
    shutil.move('obs_seq.final', case_dir+ '/output_files/'+date_str+'/')

    # Do manual post assimilation processes
    mem = 1
    while mem <= ensemble_size:
        inst_string = '{0}'.format('%04d'%mem) 

        # Do everything assimilation related in the assim_dir member directories 
        os.chdir(assim_dir+'/mem'+inst_string+'/')
        if mem == int(truth_member):
            # print('Formatting truth member file...')
            # Save the truth data
            shutil.copy('truth_data.nc', case_dir+'/analyses/'+date_str+'/truth_data_'+inst_string+'.nc' )
            # Put the truth restart back into the case directory
            shutil.move('truth_data.nc', case_dir+'/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc')
        else:
            # Save the restart and state files from the assimilation
            shutil.copy('pre_filter_restart.nc', case_dir+'/analyses/'+date_str+'/pre_filter_restart_'+inst_string+'.nc')
            shutil.copy('post_filter_restart.nc', case_dir+'/analyses/'+date_str+'/post_filter_restart_'+inst_string+'.nc')
            shutil.copy('../input.nml', '.')

            # repack the post_filter_restart.nc file if necessary
            if cat_obs is True:
                null = repack_categories('pre_filter_restart.nc')
                null = repack_categories('post_filter_restart.nc')
            
            #generate a copy of the pre_filter file that will become the postprocessed file 
            shutil.copy('pre_filter_restart.nc', 'postprocessed_restart.nc')

            # Run dart_to_cice postprocessing
            # os.symlink(assim_dir+'/dart_to_cice','dart_to_cice')
            comd = assim_dir + '/dart_to_cice > output.dart_to_cice'
            os.system(comd)
            upst_status = check_completion('output.dart_to_cice')
            if upst_status < 0:
                print('dart_to_cice did not finish correctly. Process stopping.')
                sys.exit()
        
            # Save the post-processed restart file in the analyses folder
            shutil.copy('postprocessed_restart.nc', case_dir+'/analyses/'+date_str+'/postprocessed_restart_'+inst_string+'.nc')
            # Replace the restart file in the case directory 
            shutil.move('postprocessed_restart.nc', case_dir+'/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc')

            # Remove process-specific files 
            os.remove('pre_filter_restart.nc')
            os.remove('post_filter_restart.nc')
            os.remove('output.dart_to_cice')
    
        # Move to the case directory member folder to run the model forward 
        os.chdir(case_dir+'/mem'+inst_string+'/')

        # Edit the namelist file in order to run a forecast with icepack
        icepack_nml = f90nml.read('icepack_in')
        icepack_nml['setup_nml']['ice_ic'] = case_dir + '/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc'
        icepack_nml['setup_nml']['npt'] = 24
        icepack_nml['setup_nml']['runtype_startup'] = False
        icepack_nml.write('icepack_in', force=True)

        # Run Icepack
        comd = './icepack > icepack.out'
        os.system(comd)
        sleep(2)

        if os.path.exists('restart/iced.'+next_day_str+'-00000.nc') is False:
            print('icepack did not create necessary restart file! Process stopping.')
            sys.exit()

        # Clean up 
        shutil.copy('history/icepack.h.{0}{1}{2}'.format('%04d'%year, '%02d'%mon, '%02d'%day)+'.nc', case_dir+'/forecasts/'+date_str+'/icepack.h.'+date_str+'_'+inst_string+'.nc')
        shutil.move('ice_diag.full_ITD', case_dir+'/forecasts/'+date_str+'/ice_diag.full_ITD_'+inst_string)
        os.remove('icepack.out')
        comd = 'rm ice_diag.*'
        os.system(comd)

        # Advance to the next member
        mem += 1 

    print('Done cycling '+case+' for '+date_str+'...')
    # end = perf_counter()

    # print(f'Execution time was {end - start:0.4f} seconds')
    return

###############################################################
## PERFORM CYCLING- DO NOT EDIT BELOW THIS LINE              ##
###############################################################

datelist = pd.date_range(start = first_assim_date, end = end_assim_date).to_pydatetime()
start = perf_counter()

for date in datelist:
    cycle_da(case, source_case, dart_dir, qceff_type, date, truth_member, obs_type, cat_obs=cat_assim)

end = perf_counter()
print(f'Execution time was {(end - start)/60:0.4f} minutes')