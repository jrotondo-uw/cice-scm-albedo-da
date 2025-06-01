###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml
import glob
import os
import shutil
# Removed unused import
from datetime import datetime,timedelta
import sys

import pandas as pd
import xarray as xr
import numpy as np
###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the case from which synthetic observations will be generated
source_case = sys.argv[1]
user = sys.argv[2]

# choose a random ensemble member to use as truth
rand_member = int(sys.argv[3])

# set the relevant directories
# where is this project repo located?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# where is the DART installation on your system?
dart_dir = '/glade/work/'+user+'/DART/'
# where would you like to keep the observations?
obs_dir = project_dir + '/data/processed/synthetic_obs_med_err/'+source_case+'/'+str(rand_member)+'/'
# where is the scratch directory for this user?
scratch_dir = '/glade/derecho/scratch/'+user+'/'

# set dates between which observations are desired 
first_obs_date = datetime(2011,4,1)
last_obs_date = datetime(2011,10,15)

# determine whether observations should be category-based or aggregates
if sys.argv[4] == 'category':
    category = True
else:
    category = False

# set the error type (default or uniform)
error_type   = 'default'


# obs_dir = '/glade/work/'+user+'/DA_obs/TEST/'+case+'/'
if category is True:
    obs_dir= obs_dir + '/itd/'
else:
    obs_dir= obs_dir + '/aggregate/'

###############################################################
## DEFINE HELPER FUNCTIONS HERE                              ##
###############################################################

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    Calculates longitude and latitude coordinates of a point 
    on a sphere with radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R*np.cos(lat_r) * np.cos(lon_r)
    y = R*np.cos(lat_r) * np.sin(lon_r)
    z = R*np.sin(lat_r)

    return x, y, z


def write_blank_obs_seq(lat, lon, year, mon, day, obs_seq_in, ob_types, error = 'default', case_name=None):
    # define status tracker
    status = 0
    data = xr.open_dataset('cice.r.nc')
    data = data.isel(ni=2)

    if error == 'scaled': 
        factors = xr.open_dataset(case_name+'_ensemble_scaling.nc')
        scale_fac = factors.sel(time=year+'-'+mon+'-'+day)
    else:
        scale_fac = 1 #0.01

    # Open file and begin writing to it 
    new_file = open('input.txt', 'w')
    # Enter the number of observations expected 
    num_obs = len(ob_types)
    new_file.writelines(str(num_obs)+'\n')
    # Enter number of copies of data (0 for definition)
    new_file.writelines('0\n')
    # Enter quality control values per field
    new_file.writelines('0\n')

    HL_bounds = [0.0, 0.64, 1.39, 2.47, 4.57]
    HR_bounds = [0.64, 1.39, 2.47, 4.57, 10]
    uniform_errors = [((HR - HL)**2)/12 for HR, HL in zip(HR_bounds, HL_bounds)]

    # Calculate bin widths for albedo error
    interval_widths = np.array(HR_bounds) - np.array(HL_bounds)  # Bin widths
    
    # Assign weights that decrease as albedo increases (e.g., inverse function)
    weights = 1 / (np.array(HR_bounds) + 1)  # Example: More error at low albedo
    weights /= np.sum(weights)  # Normalize so total error sums to 14%

    # Compute weighted errors
    total_error = 0.05  # 14% aggregate error
    albedo_errors = (total_error * weights) * (interval_widths**2) / 12

    fac = 1e-7

    for ob in ob_types:
        # Enter a -1 if there are no more obs
        new_file.writelines('0\n')
        # Enter observation type (string)
        new_file.writelines(ob + '\n')
        # Enter vertical coordinate option (-1, surface)
        new_file.writelines('-1\n')
        # Enter vertical coordinate height
        new_file.writelines('0\n')
        # Enter longitude
        new_file.writelines(str(lon)+'\n')
        # Enter latitude
        new_file.writelines(str(lat)+'\n')
        # Enter time (year month day hour minute second)
        new_file.writelines(year+' '+mon+' '+day+' 00 00 00\n')
        # Error based on ±5% of the albedo radiometric satellite measurement
        error_percentage_albedo = 0.14 # 14% error for upper bound of albedo observation error
        
        if ob == 'SAT_SEAICE_AGREG_CONCENTR':
            A = data.aicen.sum(dim='ncat').values
            error = scale_fac * -0.5*(A**2 - A) 
        elif ob == 'SAT_SEAICE_AGREG_FREEBOARD':
            A = data.aicen.sum(dim='ncat').values
            V = data.vicen.sum(dim='ncat').values
            if A > fac:
                H_bar = V/A
            else:
                H_bar = 0
            x = 0.1*H_bar
            error = scale_fac * 0.1 * x
        elif ob == 'SAT_SEAICE_AGREG_THICKNESS':
            A = data.aicen.sum(dim='ncat').values
            V = data.vicen.sum(dim='ncat').values
            if A > 1e-7:
                H_bar = V/A
            else:
                H_bar = 0
            error = scale_fac * 0.1 * H_bar

        # Aggregate albedo errors
        elif ob == 'SAT_SEAICE_AGREG_ALBEDODIRVIZ':
            alvdr = data.alvdr.values
            error = (scale_fac * error_percentage_albedo * alvdr)

        elif ob == 'SAT_SEAICE_AGREG_ALBEDODIRNIR':
            alidr = data.alidr.values
            error = (scale_fac * error_percentage_albedo * alidr)

        elif ob == 'SAT_SEAICE_AGREG_ALBEDOINDVIZ':
            alvdf = data.alvdf.values
            error = (scale_fac * error_percentage_albedo * alvdf)

        elif ob == 'SAT_SEAICE_AGREG_ALBEDOINDNIR':
            alidf = data.alidf.values
            error = (scale_fac * error_percentage_albedo * alidf)
        ## THIS WILL HAVE TO BECOME MORE FLEXIBLE IF WE EVER THINK ABOUT MORE CATEGORIES
        elif ob in ['SAT_SEAICE_VICE01','SAT_SEAICE_VICE02','SAT_SEAICE_VICE03','SAT_SEAICE_VICE04','SAT_SEAICE_VICE05']:
            idx = int(ob[-1]) - 1
            A = data.aicen.isel(ncat = idx).values
            U = uniform_errors[idx]
            error = A**2 * U
        elif ob in ['SAT_SEAICE_AICE01','SAT_SEAICE_AICE02','SAT_SEAICE_AICE03','SAT_SEAICE_AICE04','SAT_SEAICE_AICE05']:
            idx = int(ob[-1]) - 1
            A = data.aicen.isel(ncat = idx).values
            V = data.vicen.isel(ncat = idx).values
            if A > 1e-7:
                H_bar = V/A
            else:
                H_bar = 0
            U = uniform_errors[idx]
            if H_bar > 0:
                error = (A/H_bar)**2 * U
            else:
                error = fac
        elif ob in ['SAT_SEAICE_ALVDR0001','SAT_SEAICE_ALVDR0002','SAT_SEAICE_ALVDR0003','SAT_SEAICE_ALVDR0004','SAT_SEAICE_ALVDR0005']:
            idx = int(ob[-1]) - 1
            A = data.alvdrn.isel(ncat = idx).values
            U = albedo_errors[idx]
            error = A * total_error * U
        elif ob in ['SAT_SEAICE_ALIDR0001','SAT_SEAICE_ALIDR0002','SAT_SEAICE_ALIDR0003','SAT_SEAICE_ALIDR0004','SAT_SEAICE_ALIDR0005']:
            idx = int(ob[-1]) - 1
            A = data.alidrn.isel(ncat = idx).values
            U = albedo_errors[idx]
            error = A * total_error * U
        elif ob in ['SAT_SEAICE_ALVDF0001','SAT_SEAICE_ALVDF0002','SAT_SEAICE_ALVDF0003','SAT_SEAICE_ALVDF0004','SAT_SEAICE_ALVDF0005']:
            idx = int(ob[-1]) - 1
            A = data.alvdfn.isel(ncat = idx).values
            U = albedo_errors[idx]
            error = A * total_error * U
        elif ob in ['SAT_SEAICE_ALIDF0001','SAT_SEAICE_ALIDF0002','SAT_SEAICE_ALIDF0003','SAT_SEAICE_ALIDF0004','SAT_SEAICE_ALIDF0005']:
            idx = int(ob[-1]) - 1
            A = data.alidfn.isel(ncat = idx).values
            U = albedo_errors[idx]
            error = A * total_error * U
        else:
            error = 0
            status=-1

        if error < fac:
            error = fac
        new_file.writelines(str(error)+'\n')
    # Enter the name of the file to be output 
    new_file.writelines(obs_seq_in)
    new_file.close()

    status = status + 0
    return status

def check_completion(file_string): 

    txt = open(file_string).readlines()
    if ' Finished ... at YYYY MM DD HH MM SS = \n' not in txt:
      print('Process did not finish correctly')
      status = -1
    else:
      os.system('rm '+file_string)
      status = 0

    return status

def unpack_categories(file):

    data_vars = ['aicen','vicen','vsnon', 'alvdrn', 'alidrn', 'alvdfn', 'alidfn']
    ncat = 5
    encd = {'aice01': {'_FillValue': None},
            'aice02': {'_FillValue': None},
            'aice03': {'_FillValue': None},
            'aice04': {'_FillValue': None},
            'aice05': {'_FillValue': None},
            'vice01': {'_FillValue': None},
            'vice02': {'_FillValue': None},
            'vice03': {'_FillValue': None},
            'vice04': {'_FillValue': None},
            'vice05': {'_FillValue': None},
            'vsno01': {'_FillValue': None},
            'vsno02': {'_FillValue': None},
            'vsno03': {'_FillValue': None},
            'vsno04': {'_FillValue': None},
            'vsno05': {'_FillValue': None},
            'alvdr0001': {'_FillValue': None},
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
        ds['aice'+'{:02}'.format(n)]=ds.aicen.isel(ncat=n-1)
        ds['vice'+'{:02}'.format(n)]=ds.vicen.isel(ncat=n-1)
        ds['vsno'+'{:02}'.format(n)]=ds.vsnon.isel(ncat=n-1)
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
        expanded_vars.append('vice'+'{:02}'.format(n))
        expanded_vars.append('aice'+'{:02}'.format(n))
        expanded_vars.append('vsno'+'{:02}'.format(n))
        expanded_vars.append('alvdr000'+'{:01}'.format(n))
        expanded_vars.append('alidr000'+'{:01}'.format(n))
        expanded_vars.append('alvdf000'+'{:01}'.format(n))
        expanded_vars.append('alidf000'+'{:01}'.format(n))

    encd = {'aicen': {'_FillValue': None},
            'vicen': {'_FillValue': None},
            'vsnon': {'_FillValue': None},
            'alvdrn': {'_FillValue': None},
            'alidrn': {'_FillValue': None},
            'alvdfn': {'_FillValue': None},
            'alidfn': {'_FillValue': None}}
            
    ds = xr.open_dataset(file)
    os.rename(file, file+'~')
    aicen=np.zeros((ncat,4))
    vicen=np.zeros((ncat,4))
    vsnon=np.zeros((ncat,4))
    alvdrn=np.zeros((ncat,4))
    alidrn=np.zeros((ncat,4))
    alvdfn=np.zeros((ncat,4))
    alidfn=np.zeros((ncat,4))
    
    for n in range(1,ncat+1):
        aicen[n-1,:]=ds['aice'+'{:02}'.format(n)].values
        vicen[n-1,:]=ds['vice'+'{:02}'.format(n)].values
        vsnon[n-1,:]=ds['vsno'+'{:02}'.format(n)].values
        alvdrn[n-1,:]=ds['alvdr000'+'{:01}'.format(n)].values
        alidrn[n-1,:]=ds['alidr000'+'{:01}'.format(n)].values
        alvdfn[n-1,:]=ds['alvdf000'+'{:01}'.format(n)].values
        alidfn[n-1,:]=ds['alidf000'+'{:01}'.format(n)].values
    
    ds['aicen']=xr.DataArray(data=aicen, dims=['ncat','ni'])
    ds['vicen']=xr.DataArray(data=vicen, dims=['ncat','ni'])
    ds['vsnon']=xr.DataArray(data=vsnon, dims=['ncat','ni'])
    ds['alvdrn']=xr.DataArray(data=alvdrn, dims=['ncat','ni'])
    ds['alidrn']=xr.DataArray(data=alidrn, dims=['ncat','ni'])
    ds['alvdfn']=xr.DataArray(data=alvdfn, dims=['ncat','ni'])
    ds['alidfn']=xr.DataArray(data=alidfn, dims=['ncat','ni'])

    ds=ds.drop(expanded_vars)
    ds.to_netcdf(file, mode= 'w', encoding=encd)

    os.remove(file+'~')

    return len(ds.variables)
###############################################################
## INTERAL PROCESSES                                         ##
###############################################################

if os.path.exists(dart_dir + '/models/cice-scm/work/cice.r.nc') is False:   
    os.symlink(scratch_dir + '/ICEPACK_RUNS/'+source_case+'/mem0001/restart/iced.2011-01-02-00000.nc', dart_dir + '/models/cice-scm/work/cice.r.nc')
else:
    print('dummy cice restart file already exists! Not linking another!')

# set the observation types
if category is True:
    ob_types = ['SAT_SEAICE_ALVDR0001','SAT_SEAICE_ALVDR0002','SAT_SEAICE_ALVDR0003','SAT_SEAICE_ALVDR0004','SAT_SEAICE_ALVDR0005',
                'SAT_SEAICE_ALIDR0001','SAT_SEAICE_ALIDR0002','SAT_SEAICE_ALIDR0003','SAT_SEAICE_ALIDR0004','SAT_SEAICE_ALIDR0005',
                'SAT_SEAICE_ALVDF0001','SAT_SEAICE_ALVDF0002', 'SAT_SEAICE_ALVDF0003', 'SAT_SEAICE_ALVDF0004', 'SAT_SEAICE_ALVDF0005',
                'SAT_SEAICE_ALIDF0001', 'SAT_SEAICE_ALIDF0002', 'SAT_SEAICE_ALIDF0003', 'SAT_SEAICE_ALIDF0004', 'SAT_SEAICE_ALIDF0005']
else:
    ob_types = ['SAT_SEAICE_AGREG_ALBEDODIRVIZ', 'SAT_SEAICE_AGREG_ALBEDODIRNIR', 'SAT_SEAICE_AGREG_ALBEDOINDVIZ', 'SAT_SEAICE_AGREG_ALBEDOINDNIR',
                'SAT_SEAICE_AGREG_THICKNESS', 'SAT_SEAICE_AGREG_CONCENTR', 'SAT_SEAICE_AGREG_FREEBOARD']

# create the obs output directory if it does not already exist
if os.path.exists(obs_dir) is False:
    os.makedirs(obs_dir)

# define an observation file with the latitude and longtiude of the desired observation
location_file = project_dir + '/data/templates/obs_location.nc'
obs_loc_data = xr.open_dataset(location_file)

# define the state variables 
if category is True:
    state_variables = ['aicen', 'QTY_SEAICE_CONCENTR', 'UPDATE', 
                       'vicen', 'QTY_SEAICE_VOLUME', 'UPDATE', 
                       'vsnon', 'QTY_SEAICE_SNOWVOLUME', 'UPDATE',
                       'alvdr0001', 'QTY_SEAICE_ALVDR0001', 'NO_COPY_BACK',
                       'alvdr0002', 'QTY_SEAICE_ALVDR0002', 'NO_COPY_BACK',
                       'alvdr0003', 'QTY_SEAICE_ALVDR0003', 'NO_COPY_BACK',
                       'alvdr0004', 'QTY_SEAICE_ALVDR0004', 'NO_COPY_BACK',
                       'alvdr0005', 'QTY_SEAICE_ALVDR0005', 'NO_COPY_BACK',
                       'alvdr0001', 'QTY_SEAICE_ALIDR0001', 'NO_COPY_BACK',
                       'alidr0002', 'QTY_SEAICE_ALIDR0002', 'NO_COPY_BACK',
                       'alidr0003', 'QTY_SEAICE_ALIDR0003', 'NO_COPY_BACK',
                       'alidr0004', 'QTY_SEAICE_ALIDR0004', 'NO_COPY_BACK',
                       'alidr0005', 'QTY_SEAICE_ALIDR0005', 'NO_COPY_BACK',
                       'alvdf0001', 'QTY_SEAICE_ALVDF0001', 'NO_COPY_BACK',
                       'alvdf0002', 'QTY_SEAICE_ALVDF0002', 'NO_COPY_BACK',
                       'alvdf0003', 'QTY_SEAICE_ALVDF0003', 'NO_COPY_BACK',
                       'alvdf0004', 'QTY_SEAICE_ALVDF0004', 'NO_COPY_BACK',
                       'alvdf0005', 'QTY_SEAICE_ALVDF0005', 'NO_COPY_BACK',
                       'alidf0001', 'QTY_SEAICE_ALIDF0001', 'NO_COPY_BACK',
                       'alidf0002', 'QTY_SEAICE_ALIDF0002', 'NO_COPY_BACK',
                       'alidf0003', 'QTY_SEAICE_ALIDF0003', 'NO_COPY_BACK',
                       'alidf0004', 'QTY_SEAICE_ALIDF0004', 'NO_COPY_BACK',
                       'alidf0005', 'QTY_SEAICE_ALIDF0005', 'NO_COPY_BACK']

else:   
    state_variables = ['alvdr', 'QTY_SEAICE_AGREG_ALBEDODIRVIZ', 'NO_COPY_BACK',
                       'alidr', 'QTY_SEAICE_AGREG_ALBEDODIRNIR', 'NO_COPY_BACK',
                       'alvdf', 'QTY_SEAICE_AGREG_ALBEDOINDVIZ', 'NO_COPY_BACK',
                       'alidf', 'QTY_SEAICE_AGREG_ALBEDOINDNIR', 'NO_COPY_BACK',
                       'aicen', 'QTY_SEAICE_CONCENTR', 'UPDATE',
                       'vicen', 'QTY_SEAICE_VOLUME', 'UPDATE',
                       'vsnon', 'QTY_SEAICE_SNOWVOLUME', 'UPDATE']
    
# set the date range for observations times 
datelist = pd.date_range(start = first_obs_date, end = last_obs_date)
datelist = datelist[~((datelist.month == 2) & (datelist.day == 29))]
datelist = datelist.to_pydatetime()

# cycle through time and generate observations
for t in range(0, datelist.shape[0]):
    count = 0 
    print('Working on ', datelist[t], '...')

    # Parse date information
    today = int(datelist[t].strftime('%Y%m%d')) #curr_date
    yesterday = int((datelist[t] - timedelta(days=1)).strftime('%Y%m%d')) #check_date

    year = datelist[t].strftime('%Y')
    mon =  datelist[t].strftime('%m') 
    day =  datelist[t].strftime('%d')

    year_ystrdy = (datelist[t]-timedelta(days=1)).strftime('%Y')
    mon_ystrdy = (datelist[t]-timedelta(days=1)).strftime('%m')
    day_ystrdy = (datelist[t]-timedelta(days=1)).strftime('%d')

    diff = datelist[t] - datetime(1601,1,1,0)
    nml_file = f90nml.read(dart_dir + '/models/cice-scm/scripts/templates/DART_input_category.nml.template')
    # nml_file = f90nml.read(project_dir + '/data/templates/DART_input.nml.template')
    nml_file['algorithm_info_nml']['qceff_table_filename'] = '/glade/work/jrotondo/DART/models/cice-scm/scripts/templates/cice_bounded_qceff_table.csv'
    nml_file['perfect_model_obs_nml']['init_time_days'] = diff.days
    nml_file['perfect_model_obs_nml']['init_time_seconds'] = diff.seconds
    nml_file['perfect_model_obs_nml']['input_state_files'] = 'cice.r.nc'
    nml_file['obs_kind_nml']['assimilate_these_obs_types'] = ob_types
    nml_file['model_nml']['model_state_variables'] = state_variables
    nml_file['obs_seq_to_netcdf_nml']['obs_sequence_name'] = ''
    nml_file['obs_seq_to_netcdf_nml']['obs_sequence_list'] = 'observations_list.txt'
    nml_file['obs_seq_to_netcdf_nml']['append_to_netcdf'] = True
    nml_file.write('input.nml',force=True)

    # define the truth file for this date (a restart from the source case)
    truth_file = scratch_dir + '/ICEPACK_RUNS/'+source_case+'/mem00{00}/restart/'.format('%02d'%rand_member)+'iced.{0}'.format('%04d'%int(year))+'-{0}'.format('%02d'%int(mon))+'-{0}'.format('%02d'%int(day))+'-00000.nc'
    if os.path.exists('cice.r.nc'):
        os.remove('cice.r.nc')
    os.symlink(truth_file, 'cice.r.nc')
    
    if category is True:
        null = unpack_categories(truth_file)
        null = unpack_categories(dart_dir + '/models/cice-scm/work/cice.r.nc')

    # find corresponding modelled sea ice truth file
    if len(glob.glob(truth_file)) == 0:
        print('Input file for observations not found! No sea ice data for '+str(today)+'.')
        continue
    else:
        print('Extracting observations from:', truth_file)

    # save the only lat, lon (ICEPACK is a single column model)
    lat = obs_loc_data.latitude.values[0]
    lon = obs_loc_data.longitude.values[0]
    # print('Number of locations for requested date:' + str(lat.shape[0]))

    # Begin constructing input file 
    hour_tag = '00000'
    obs_seq_in = 'obs_seq.in_'+str(yesterday)+'_'+hour_tag

    # Construct the input information for create_obs_sequence 
    inpt_status = write_blank_obs_seq(lat, lon, year, mon, day, obs_seq_in, ob_types, error = error_type, case_name=source_case)
    if inpt_status < 0:
        print('Observation types were not recognized! Processing stopping.')
        sys.exit()
    else:
        print('Ready to create observation sequence!')

    # print current directory
    print('Current directory:', os.getcwd())

    # Run create_obs_sequence
    comd = dart_dir + '/models/cice-scm/work/create_obs_sequence < input.txt > output.create_obs_sequence'
    os.system(comd)
    cos_status = check_completion('output.create_obs_sequence')
    if cos_status < 0:
        print('create_obs_sequence did not complete. Processing stopping.')
        sys.exit()
    
    # Run perfect_model_obs
    
    os.symlink(obs_seq_in, 'obs_seq.in')
    if os.path.exists('input_file.nc'):
        os.remove('input_file.nc')
    os.symlink(truth_file, 'input_file.nc')
    comd = dart_dir + '/models/cice-scm/work/perfect_model_obs > output.perfect_model_obs_'+str(yesterday)+'_'+hour_tag
    os.system(comd)
    pmo_status = check_completion('output.perfect_model_obs_'+str(yesterday)+'_'+hour_tag)
    if pmo_status < 0:
        print('perfect_model_obs did not complete. Process stopping.')
        sys.exit()
    
    # Clean up 
    os.rename('obs_seq.out', 'obs_seq.'+str(today))
    if os.path.exists('obs_seq.'+str(today)) is False:
        print('Obs file was not created! Process stopping.')
        sys.exit()
    else:
        if os.path.exists(obs_dir + '/raw_files/') is False:
            os.makedirs(obs_dir + '/raw_files/')
        shutil.move('obs_seq.in_'+str(yesterday)+'_00000', obs_dir+'/raw_files/') 
        os.remove('obs_seq.in')
        # os.remove('input_file.nc')
        os.remove('cice.r.nc')
        # os.remove('obs_seq_flist')
        if category is True:
            null = repack_categories(truth_file)
            null = repack_categories(dart_dir + '/models/cice-scm/work/cice.r.nc')
        
        
        os.remove('input.txt')
        os.remove('forward_op_errors0')
        os.remove('perfect_restart.nc')


# Generate a netcdf file with all the observation sequence files
# Move all generated observation sequence files to obs_dir and clean up
comd = 'mv obs_seq.* '+ obs_dir
os.system(comd)

files = sorted(glob.glob(obs_dir+'/obs_seq.*'))

text_file = open("observations_list.txt", "w")
for file in files:
    n = text_file.write(file)
    n = text_file.write('\n')
text_file.close()

comd = dart_dir + '/models/cice-scm/work/obs_seq_to_netcdf > output.obs_seq_to_netcdf'
os.system(comd)
# Check that the script ran correctly
otncdf_status = check_completion('output.obs_seq_to_netcdf')
if otncdf_status < 0:
    print('obs_seq_to_netcdf did not complete. Process stopping.')
    sys.exit()
else:
    # make a netcdf directory in the obs_dir
    if os.path.exists(obs_dir + '../netcdfs/') is False:
        os.makedirs(obs_dir + '../netcdfs/')
    
    if category is False:
        os.system('mv obs_epoch_001.nc '+obs_dir+'../netcdfs/aggregate_observations.nc')
    else:
        os.system('mv obs_epoch_001.nc '+obs_dir+'../netcdfs/itd_observations.nc')
    
    os.remove('input.nml')
    os.remove('observations_list.txt')

num_files = len(glob.glob(obs_dir + '/obs_seq.*'))
print(str(num_files) + ' days of observations successfully converted. Program finished.')