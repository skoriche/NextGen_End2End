############################################################################################
# Author  : Ahmad Jan Khattak
# Contact : ahmad.jan@noaa.gov
# Date    : July 5, 2024
############################################################################################
import os, sys
import pandas as pd
import subprocess
import glob
import yaml
import platform
from src_py import configuration
import json
from pathlib import Path

os_name = platform.system()
workflow_infile   = sys.argv[1]
ngen_cal_basefile = sys.argv[2]

with open(workflow_infile, 'r') as file:
    d = yaml.safe_load(file)

workflow_dir     = d["workflow_dir"]
input_dir        = d["input_dir"]
output_dir       = Path(d["output_dir"])


dformul = d['formulation']
ngen_dir         = dformul["ngen_dir"]
nproc            = int(dformul.get('num_processors_sim', 1))
nproc_adaptive   = int(dformul.get('num_processors_adaptive', True))


dcalib = d['ngen_cal']
ngen_cal_type    = dcalib.get('task_type', None)
#simulation_time  = json.loads(ddformul["simulation_time"])
calibration_time = pd.NaT
validation_time  = pd.NaT

if (ngen_cal_type == 'calibration' or ngen_cal_type == 'calibvalid'):
    calibration_time  = json.loads(dcalib["calibration_time"])
    calib_eval_time  = json.loads(dcalib["calib_eval_time"])

if (ngen_cal_type == 'validation' or ngen_cal_type == 'calibvalid'):
    validation_time  = json.loads(dcalib["validation_time"])
    valid_eval_time  = json.loads(dcalib["valid_eval_time"])

restart_dir = "./"
if (ngen_cal_type == 'restart'):
    restart_dir    = dcalib.get('restart_dir')
    
    if restart_dir is None:
        raise ValueError("ngen_cal_type is restart, however, restart_dir in None. It must be set to a valid directory.")

    if not restart_dir:
        raise FileNotFoundError(f"restart_dir does not exist, provided {restart_dir}.")



def run_ngen_without_calibration():
    
    infile = os.path.join(output_dir, "basins_passed.csv")
    indata = pd.read_csv(infile, dtype=str)

    ngen_exe = os.path.join(ngen_dir, "cmake_build/ngen")

    
    for id, ncats in zip(indata["basin_id"], indata['n_cats']):

        ncats = int(ncats)

        o_dir = output_dir / id
        i_dir = Path(input_dir) / id

        os.chdir(o_dir)
        print ("cwd: ", os.getcwd())
        print ("input_dir: ", i_dir)
        print ("output_dir: ", o_dir)

        

        gpkg_file = Path(glob.glob(str(i_dir / "data" / "*.gpkg"))[0])
        gpkg_name = gpkg_file.stem

        nproc_local = nproc

        file_par = ""
        if (nproc_local > 1):
            nproc_local, file_par = generate_partition_basin_file(ncats, gpkg_file)
        
        print ("Running basin %s on cores %s ********"%(id, nproc_local), flush = True)
        
        realization = glob.glob("json/realization_*.json")

        assert (len(realization) == 1)

        realization = realization[0]
        
        if (nproc_local == 1):
            run_cmd = f'{ngen_exe} {gpkg_file} all {gpkg_file} all {realization}'
        else:
            run_cmd = f'mpirun -np {nproc_local} {ngen_exe} {gpkg_file} all {gpkg_file} all {realization} {file_par}'

        if os_name == "Darwin":
            run_cmd = f'PYTHONEXECUTABLE=$(which python) {run_cmd}'
        
        print (f"Run command: {run_cmd} ", flush = True)
        result = subprocess.call(run_cmd,shell=True)

    
def run_ngen_with_calibration():

    infile = os.path.join(output_dir, "basins_passed.csv")
    indata = pd.read_csv(infile, dtype=str)
   
    for id, ncats in zip(indata["basin_id"], indata['n_cats']):

        ncats = int(ncats)

        o_dir = output_dir / id
        i_dir = Path(input_dir) / id

        os.chdir(o_dir)
        print ("cwd: ", os.getcwd())
        print ("input_dir: ", i_dir)
        print ("output_dir: ", o_dir)


        gpkg_file = Path(glob.glob(str(i_dir / "data" / "*.gpkg"))[0])
        gpkg_name = gpkg_file.stem

        nproc_local = nproc

        

        


        file_par = ""
        if (nproc_local > 1):
            nproc_local, file_par = generate_partition_basin_file(ncats, gpkg_file)
            file_par = os.path.join(o_dir, file_par)

        print ("Running basin %s on cores %s ********"%(id, nproc_local), flush = True)
        

        # Calibration call
        if (ngen_cal_type  == 'calibration' or ngen_cal_type == 'calibvalid'):
            start_time = pd.Timestamp(calibration_time['start_time']).strftime("%Y%m%d%H%M")
        
            #troute_output_file = os.path.join(dir, "outputs/troute", "troute_output_{}.csv".format(start_time))
            #troute_output_file = os.path.join(dir, "outputs/troute", "flowveldepth_{}.csv".format(gpkg_name))
            troute_output_file = os.path.join("./troute_output_{}.nc".format(start_time))
            
            configuration.write_calib_input_files(gpkg_file            = gpkg_file,
                                                  ngen_dir             = ngen_dir,
                                                  output_dir           = o_dir,
                                                  realization_file_par = file_par,
                                                  ngen_cal_basefile    = ngen_cal_basefile,
                                                  troute_output_file   = troute_output_file,
                                                  ngen_cal_type        = 'calibration',
                                                  restart_dir          = restart_dir,
                                                  simulation_time      = calibration_time,
                                                  evaluation_time      = calib_eval_time,
                                                  num_proc             = nproc_local)

            run_command = f"python -m ngen.cal configs/ngen-cal_calib_config.yaml"
            # .yaml file under cat_id/configs
            result = subprocess.call(run_command,shell=True)

        # Validation call
        if (ngen_cal_type  == 'validation' or ngen_cal_type == 'calibvalid'):

            start_time = pd.Timestamp(validation_time['start_time']).strftime("%Y%m%d%H%M")
            troute_output_file = os.path.join("./troute_output_{}.nc".format(start_time))
            
            configuration.write_calib_input_files(gpkg_file            = gpkg_file,
                                                  ngen_dir             = ngen_dir,
                                                  output_dir           = o_dir,
                                                  realization_file_par = file_par,
                                                  ngen_cal_basefile    = ngen_cal_basefile,
                                                  troute_output_file   = troute_output_file,
                                                  ngen_cal_type        = 'validation',
                                                  restart_dir          = restart_dir,
                                                  simulation_time      = validation_time,
                                                  evaluation_time      = valid_eval_time,
                                                  num_proc             = nproc_local)

            run_command = f"python {workflow_dir}/src_py/validation.py configs/ngen-cal_valid_config.yaml"
            # .yaml file under cat_id/configs
            result = subprocess.call(run_command,shell=True)


#####################################################################
def generate_partition_basin_file(ncats, gpkg_file):

    nproc_local = nproc
    json_dir   = "json"

    if (ncats <= nproc_local):
        nproc_local = ncats
    elif(nproc_adaptive):
        nproc_local = min(int(ncats/nproc_local), 20)

    fpar = " "
    
    if (nproc_local > 1):
        fpar = os.path.join(json_dir, f"partition_{nproc_local}.json")
        partition=f"{ngen_dir}/cmake_build/partitionGenerator {gpkg_file} {gpkg_file} {fpar} {nproc_local} \"\" \"\" "
        result = subprocess.call(partition,shell=True)

    return nproc_local, fpar

if __name__ == "__main__":

    if (nproc > 1 and not os.path.exists(f"{ngen_dir}/cmake_build/partitionGenerator")):
        sys.exit("Partitioning geopackage is requested but partitionGenerator does not exit! Quitting...")


    if (not ngen_cal_type in ['calibration', 'validation', 'calibvalid', 'restart']):
        print ("Running NextGen without calibration ...")
        run_ngen_without_calibration()
    else:
        print (f'Running NextGen with {ngen_cal_type}')
        run_ngen_with_calibration()
        
