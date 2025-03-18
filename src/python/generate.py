############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : October 11, 2023 
############################################################################################

# driver of the script, generates configuration files and realization files by calling
# the respective scripts

import os
import sys
import subprocess
from pathlib import Path
import argparse
import json

from src/python import configuration
from src/python import realization

class Driver:
    def __init__(self, workflow_dir, gpkg_file, forcing_dir, ngen_dir, sim_time, formulation,
                 output_dir, forcing_format, ngen_cal_type, schema):
        self.workflow_dir = workflow_dir
        self.gpkg_file = gpkg_file
        self.forcing_dir = forcing_dir
        self.ngen_dir = ngen_dir
        self.simulation_time = sim_time
        self.formulation_in = formulation.replace(" ", "") # remove space if any
        self.output_dir = output_dir
        self.schema = schema
        self.verbosity = 1
        self.forcing_format = forcing_format
        self.ngen_cal_type = ngen_cal_type
        
        if not os.path.exists(self.gpkg_file):
            sys.exit(f'The gpkg file does not exist: ', self.gpkg_file)

        if not os.path.exists(self.forcing_dir):
            sys.exit(f'The forcing directory does not exist: ', self.forcing_dir)

        self.formulation_supported = [
            "CFE-S",
            "CFE-X",
            "LASAM",
            "NOM,CFE-S",
            "NOM,CFE-X",
            "PET,CFE-S",
            "PET,CFE-X",
            "NON,LASAM",
            "PET,LASAM",
            "NOM,CFE-S,PET",
            "NOM,CFE-X,PET",
            "NOM,CFE-S,SMP,SFT",
            "NOM,CFES-X,SMP,SFT",
            "NOM,LASAM,SMP,SFT",
            "NOM,TOPMODEL",
            "BASELINE,CFE",
            "BASELINE,LAS"
        ]

        if "T-route" in self.formulation_in or "t-route" in self.formulation_in:
            self.formulation_supported_w_troute = [f'{model},T-route' for model in self.formulation_supported]
            
        if self.formulation_in in self.formulation_supported or self.formulation_in in self.formulation_supported_w_troute:
            self.formulation = self.formulation_in
        else:
            str_msg = f"Invalid model option provided or formulation currently not supported: " + self.formulation_in
            sys.exit(str_msg)


        if self.verbosity >= 3:
            print("*******************************************")
            print(self.colors.BLUE)
            print("Running (from driver.py):\n", generate_config_files)
            print("Model option provided: ", args.models_option)
            print("Generating configuration files for model(s) option: ", self.coupled_models_options[args.models_option])
            print(self.colors.ENDC)
            print("*******************************************")

        ConfigGen = configuration.ConfigurationGenerator(workflow_dir = self.workflow_dir,
                                                         gpkg_file = self.gpkg_file,
                                                         forcing_dir = self.forcing_dir,
                                                         output_dir = self.output_dir,
                                                         ngen_dir = self.ngen_dir,
                                                         formulation = self.formulation,
                                                         surface_runoff_scheme = "NASH_CASCADE",
                                                         simulation_time = self.simulation_time,
                                                         verbosity = 1,
                                                         ngen_cal_type = self.ngen_cal_type,
                                                         schema_type = self.schema)

        if self.verbosity >= 3:
            print("Running (from driver.py): \n ", generate_realization_file)
            print(self.colors.ENDC)
            
        if "NOM" in self.formulation:
            print ("caling NOM")
            ConfigGen.write_nom_input_files()
        if "PET" in self.formulation:
            print ("caling PET")
            ConfigGen.write_pet_input_files()
        if "CFE-S" in self.formulation or "CFE-X" in self.formulation:
            print ("caling CFE")
            ConfigGen.write_cfe_input_files()

        if "T-route" in self.formulation_in or "t-route" in self.formulation_in:
            ConfigGen.write_troute_input_files()
        
        result = False
        if result:
            sys.exit("config files could not be generated, check the options provided!")

        if self.verbosity >= 3:
            print("*******************************************")
            print(self.colors.GREEN)
            print("Generating realization file ...")
            
        RealGen = realization.RealizationGenerator(forcing_dir = self.forcing_dir,
                                                   output_dir = self.output_dir,
                                                   ngen_dir = self.ngen_dir,
                                                   formulation = self.formulation,
                                                   surface_runoff_scheme = "NASH_CASCADE",
                                                   simulation_time = self.simulation_time,
                                                   verbosity = 1,
                                                   ngen_cal_type = self.ngen_cal_type,
                                                   forcing_format = self.forcing_format)
        RealGen.write_realization_file()

        if result:
            sys.exit("realization file could not be generated, check the options provided!")
        else:
            print("************* DONE (Baseline realization file successfully generated!) ************** ")
                
    class colors:
        BLUE = '\33[34m'
        BLACK = '\33[30m'
        RED = '\33[31m'
        CYAN = '\033[96m'
        GREEN = '\033[32m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
