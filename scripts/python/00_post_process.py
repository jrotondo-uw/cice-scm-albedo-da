###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
import xarray as xr
import glob
import os
import sys
from time import perf_counter

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

case = sys.argv[1]
user = sys.argv[2]

stage = 'restart'

output_path = f'/glade/work/{user}/Projects/cice-scm-da/data/processed/{case}/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

###############################################################
## PERFORM POSTPROCESSING                                    ##
###############################################################

start = perf_counter()

print(f'Working on case: {case}')
print(f'Processing stage: {stage}')
ensemble_dirs = sorted(glob.glob(f'/glade/derecho/scratch/{user}/ICEPACK_RUNS/{case}/mem00??'))

ensemble_size = len(ensemble_dirs)
print(f'Found {ensemble_size} ensemble members.')

DS = []

for i in range(1, ensemble_size + 1):

    restart_path = f'/glade/derecho/scratch/{user}/ICEPACK_RUNS/{case}/mem00{i:02d}/restart/'
    restart_files = sorted(glob.glob(restart_path + 'iced.????-??-??-00000.nc'))

    if not restart_files:
        print(f"[SKIP] No restart files found for member {i} at {restart_path}")
        continue

    try:
        print(f"[INFO] Concatenating {len(restart_files)} files for member {i}...")
        ds_list = [xr.open_dataset(f) for f in restart_files]
        member_ds = xr.concat(ds_list, dim='time')
        member_ds = member_ds.expand_dims('ensemble_member')

        member_output_file = f'{output_path}{stage}_member_00{i:02d}.nc'
        member_ds.to_netcdf(member_output_file)

        DS.append(member_ds)
    except Exception as e:
        print(f"[ERROR] Xarray concat failed for member {i}: {e}")
        continue

# Concatenate across ensemble members
if DS:
    print(f"[INFO] Concatenating {len(DS)} members into ensemble dataset...")
    ens_ds = xr.concat(DS, dim='ensemble_member')
    ens_output_file = f'{output_path}{stage}_ens.nc'
    ens_ds.to_netcdf(ens_output_file)

    # Clean up individual member files
    for f in glob.glob(f'{output_path}{stage}_member_00*.nc'):
        os.remove(f)

    print(f"[SUCCESS] Ensemble file written: {ens_output_file}")
else:
    print("[WARNING] No ensemble data to write.")

end = perf_counter()
print(f'Postprocessing completed. Output path: {output_path}')
print(f'Execution time: {(end - start)/60:0.4f} minutes')