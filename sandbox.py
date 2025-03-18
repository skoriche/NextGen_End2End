############################################################################################
# Author  : Ahmad Jan Khattak
# Contact : ahmad.jan.khattak@noaa.gov
# Date    : July 16, 2024
############################################################################################

import os, sys
import subprocess
import yaml
import argparse
from pathlib import Path

path = Path(sys.argv[0]).resolve()
workflow_dir = path.parent


from src_py import forcing

from src_py import main
import runner

def Runner(config_workflow, config_calib):
    
    if (args.gpkg):
        print ("Generating geopackages...")
        generate_gpkg = f"Rscript {workflow_dir}/src_r/main.R {config_workflow}"
        status = subprocess.call(generate_gpkg,shell=True)

        if (status):
            sys.exit("Failed during generating geopackge(s) step...")
        else:
            print ("DONE \u2713")

    if (args.forc):
        print ("Generating forcing data...")
        forcings = forcing.ForcingProcessor(config_workflow)
        status = forcings.download_forcing()

        if (status):
            sys.exit("Failed during generating geopackge(s) step...")
        else:
            print ("DONE \u2713")

    if (args.conf):
        print ("Generating config files...")
        main_ = main.Main(config_workflow)
        main_.run()
        #generate_configs = f"python {workflow_dir}/src_py/main.py {config_workflow}"
        #status = subprocess.call(generate_configs,shell=True)

        #if (status):
        #    sys.exit("Failed during generating config files step...")
        #else:
        #    print ("DONE \u2713")
        
    if (args.run):
        print ("Calling Runner...")

        sandbox_runner = runner.SandboxRunner(config_workflow, config_calib)
        sandbox_runner.run()
        #run_command = f"python {workflow_dir}/runner.py {config_workflow} {config_calib}"
        #status = subprocess.call(run_command,shell=True)

        #if (status):
        #    sys.exit("Failed during ngen-cal execution...")
        #else:
        #    print ("DONE \u2713")
    
    print ("**********************************")
    
    

if __name__ == "__main__":
    
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-gpkg", action='store_true', help="generate gpkg files")
        parser.add_argument("-forc", action='store_true', help="generate forcing data")
        parser.add_argument("-conf", action='store_true', help="generate config files")
        parser.add_argument("-run",  action='store_true', help="run nextgen without caliberation")
        parser.add_argument("-i",    dest="workflow_infile",  type=str, required=False,  help="workflow config file")
        parser.add_argument("-j",    dest="calib_infile",     type=str, required=False,  help="caliberation config file")
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if (args.workflow_infile):
        if (os.path.exists(args.workflow_infile)):
            config_workflow = Path(args.workflow_infile).resolve()
        else:
            print ("workflow config file DOES NOT EXIST, provided: ", args.workflow_infile)
            sys.exit(0)
    else:
        config_workflow = f"{workflow_dir}/configs/config_workflow.yaml"

    if (args.calib_infile):
        if (os.path.exists(args.calib_infile)):
            config_calib = Path(args.calib_infile).resolve()
        else:
            print ("caliberation config file DOES NOT EXIST, provided: ", args.calib_infile)
            sys.exit(0)
    else:
        config_calib = f"{workflow_dir}/configs/config_calib.yaml"
    
    if (len(sys.argv) < 2):
        print ("No arguments are provide, printing help()")
        parser.print_help()
        sys.exit(0)
    
    Runner(config_workflow, config_calib)
