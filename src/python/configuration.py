############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : September 28, 2023
############################################################################################

import os
import sys
import argparse
import re
import copy
import glob
import json
import subprocess
import pandas as pd
import geopandas as gpd
import numpy as np
import fiona
import yaml
import platform
from pathlib import Path

try:
    from src.python import schema
except:
    import schema

os_name = platform.system()

class ConfigurationGenerator:
    def __init__(self, workflow_dir, gpkg_file, forcing_dir, output_dir,
                 ngen_dir, formulation, surface_runoff_scheme, simulation_time,
                 verbosity, ngen_cal_type, schema_type = None):
        self.workflow_dir = workflow_dir
        self.gpkg_file = gpkg_file
        self.forcing_dir = forcing_dir
        self.output_dir = output_dir
        self.ngen_dir = ngen_dir
        self.formulation = formulation
        self.surface_runoff_scheme = surface_runoff_scheme
        self.simulation_time = json.loads(simulation_time)
        self.verbosity = verbosity
        #self.json_dir = json_dir
        #self.sim_output_dir = sim_output_dir
        self.ngen_cal_type = ngen_cal_type
        self.schema_type = schema_type


        self.gdf, self.catids = self.read_gpkg_file()
        
        self.soil_class_NWM = self.get_soil_class_NWM()
        
    def get_soil_class_NWM(self):
        nom_soil_file = os.path.join(self.ngen_dir,"extern/noah-owp-modular/noah-owp-modular/parameters/SOILPARM.TBL")
        header = ['index', 'BB', 'DRYSMC', 'F11', 'MAXSMC', 'REFSMC', 'SATPSI', 'SATDK', 'SATDW', 'WLTSMC', 'QTZ', 'BVIC', 'AXAJ', 'BXAJ', 'XXAJ', 'BDVIC', 'BBVIC', 'GDVIC', 'ISLTYP']
        df = pd.read_table(nom_soil_file, delimiter=',', index_col=0, skiprows=3, nrows=19, names=header)
        return df

    def read_gpkg_file(self):
        try:
            gdf_soil = gpd.read_file(self.gpkg_file, layer='divide-attributes')
        except:
            try:
                gdf_soil = gpd.read_file(self.gpkg_file, layer='model_attributes')
            except:
                print("layer 'divide-attributes or model_attributes does not exist!'")
                sys.exit(1)

        gdf_soil.set_index("divide_id", inplace=True)
        gdf_div = gpd.read_file(self.gpkg_file, layer='divides')
        gdf_div = gdf_div.to_crs("EPSG:4326")

        layers = fiona.listlayers(self.gpkg_file)
        flowpath_layer = [layer for layer in layers if 'flowpath' in layer][0]

        if self.verbosity >= 3:
            print("Geopackage layers: ", layers)
            print("\n")

        params = schema.get_schema_model_attributes(gdf_soil)

        gdf_soil['soil_b'] = gdf_soil[params['soil_b']].fillna(16)
        gdf_soil['soil_dksat'] = gdf_soil[params['soil_dksat']].fillna(0.00000338)
        gdf_soil['soil_psisat'] = gdf_soil[params['soil_psisat']].fillna(0.355)
        gdf_soil['soil_smcmax'] = gdf_soil[params['soil_smcmax']].fillna(0.439)
        gdf_soil['soil_smcwlt'] = gdf_soil[params['soil_smcwlt']].fillna(0.066)
        gdf_soil['gw_Zmax'] = gdf_soil[params['gw_Zmax']].fillna(0.01)
        gdf_soil['gw_Coeff'] = gdf_soil[params['gw_Coeff']].fillna(1.8e-05)
        gdf_soil['gw_Expon'] = gdf_soil[params['gw_Expon']].fillna(6.0)
        gdf_soil['soil_slope'] = gdf_soil[params['soil_slope']].fillna(1.0)
        gdf_soil['ISLTYP'] = gdf_soil[params['ISLTYP']].fillna(1).astype(int)
        gdf_soil['IVGTYP'] = gdf_soil[params['IVGTYP']].fillna(1).astype(int)
        gdf_soil['gw_Zmax'] = gdf_soil['gw_Zmax'] / 1000.0
        gdf_soil['gw_Coeff'] = gdf_soil['gw_Coeff'] * 3600 / (7.337700 * 1000 * 1000)
        gdf_soil['elevation_mean'] = gdf_soil[params['elevation_mean']].fillna(4)

        if self.schema_type == 'dangermond':
            gdf_soil['elevation_mean'] = gdf_soil['elevation_mean'] / 100.0

        if 'refkdt' in gdf_soil:
            gdf_soil['soil_refkdt'] = gdf_soil[params['soil_refkdt']].fillna(3.0)
        else:
            gdf_soil['soil_refkdt'] = 3.0

        gdf = gpd.GeoDataFrame(data={'geometry': gdf_div['geometry'].values}, index=gdf_soil.index)
        gdf['soil_b'] = gdf_soil['soil_b'].copy()
        gdf['soil_satdk'] = gdf_soil['soil_dksat'].copy()
        gdf['soil_satpsi'] = gdf_soil['soil_psisat'].copy()
        gdf['soil_slop'] = gdf_soil['soil_slope'].copy()
        gdf['soil_smcmax'] = gdf_soil['soil_smcmax'].copy()
        gdf['soil_wltsmc'] = gdf_soil['soil_smcwlt'].copy()
        gdf['soil_refkdt'] = gdf_soil['soil_refkdt'].copy()
        gdf['max_gw_storage'] = gdf_soil['gw_Zmax'].copy()
        gdf['Cgw'] = gdf_soil['gw_Coeff'].copy()
        gdf['gw_expon'] = gdf_soil['gw_Expon'].copy()
        gdf['ISLTYP'] = gdf_soil['ISLTYP'].copy()
        gdf['IVGTYP'] = gdf_soil['IVGTYP'].copy()
        gdf['elevation_mean'] = gdf_soil['elevation_mean'].copy()

        mask = gdf['soil_b'].gt(0.0)
        min_value = gdf['soil_b'][mask].min()
        mask = gdf['soil_b'].le(0.0)
        gdf.loc[mask, 'soil_b'] = min_value

        mask = gdf['elevation_mean'].le(0.0)
        gdf.loc[mask, 'elevation_mean'] = 1.0

        if "topmodel" in self.formulation:
            gdf['twi'] = gdf_soil[params['twi']]
            gdf['width_dist'] = gdf_soil[params['width_dist']]

        if "CFE-S" in self.formulation or "CFE-X" in self.formulation:
            if self.surface_runoff_scheme == "GIUH" or self.surface_runoff_scheme == 1:
                gdf['giuh'] = gdf_soil[params['giuh']]
            elif self.surface_runoff_scheme == "NASH_CASCADE" or self.surface_runoff_scheme == 2:
                gdf['N_nash_surface'] = gdf_soil[params['N_nash_surface']]
                gdf['K_nash_surface'] = gdf_soil[params['K_nash_surface']]

        if "LASAM" in self.formulation:
            gdf['giuh'] = gdf_soil[params['giuh']]
        
        df_cats = gpd.read_file(self.gpkg_file, layer='divides')
        catids = [int(re.findall('[0-9]+', s)[0]) for s in df_cats['divide_id']]

        return gdf, catids


    def write_nom_input_files(self):
        nom_dir = os.path.join(self.output_dir,"configs/nom")

        self.create_directory(nom_dir)
        start_time = pd.Timestamp(self.simulation_time['start_time']).strftime("%Y%m%d%H%M")
        end_time = pd.Timestamp(self.simulation_time['end_time']).strftime("%Y%m%d%H%M")

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            centroid_x = str(self.gdf['geometry'][cat_name].centroid.x)
            centroid_y = str(self.gdf['geometry'][cat_name].centroid.y)
            soil_type = str(self.gdf.loc[cat_name]['ISLTYP'])
            veg_type = str(self.gdf.loc[cat_name]['IVGTYP'])

            timing = [
                "&timing                                   ! and input/output paths",
                "  dt                 = 3600.0             ! timestep [seconds]",
                f"  startdate          = \"{start_time}\"             ! UTC time start of simulation (YYYYMMDDhhmm)",
                f"  enddate            = \"{end_time}\"             ! UTC time end of simulation (YYYYMMDDhhmm)",
                f"  forcing_filename   = \"{os.path.join(self.forcing_dir, cat_name)}.csv\"         ! file containing forcing data",
                f"  output_filename    = \"output-{cat_name}.csv\"",
                "/\n"
            ]

            params = [
                "&parameters",
                f"  parameter_dir      = \"{os.path.join(nom_dir, 'parameters')}\"  ! location of input parameter files",
                "  general_table      = \"GENPARM.TBL\"                ! general param tables and misc params",
                "  soil_table         = \"SOILPARM.TBL\"               ! soil param table",
                "  noahowp_table      = \"MPTABLE.TBL\"                ! model param tables (includes veg)",
                "  soil_class_name    = \"STAS\"                       ! soil class data source - STAS or STAS-RUC",
                "  veg_class_name     = \"MODIFIED_IGBP_MODIS_NOAH\"   ! vegetation class data source - MODIFIED_IGBP_MODIS_NOAH or USGS",
                "/\n"
            ]

            location = [
                "&location                                         ! for point runs",
                f"  lat              = {centroid_y}                           ! latitude [degrees]  (-90 to 90)",
                f"  lon              = {centroid_x}                           ! longitude [degrees] (-180 to 180)",
                "  terrain_slope    = 0.0                          ! terrain slope [degrees]",
                "  azimuth          = 0.0                          ! terrain azimuth or aspect [degrees clockwise from north]",
                "/ \n"
            ]

            forcing = [
                "&forcing",
                "  ZREF               = 10.0                        ! measurement height for wind speed (m)",
                "  rain_snow_thresh   = 1.0                         ! rain-snow temperature threshold (degrees Celcius)",
                "/ \n"
            ]

            model_opt = [
                "&model_options                                   ! see OptionsType.f90 for details",
                "  precip_phase_option               = 1",
                "  snow_albedo_option                = 1",
                "  dynamic_veg_option                = 4",
                "  runoff_option                     = 3",
                "  drainage_option                   = 8",
                "  frozen_soil_option                = 1",
                "  dynamic_vic_option                = 1",
                "  radiative_transfer_option         = 3",
                "  sfc_drag_coeff_option             = 1",
                "  canopy_stom_resist_option         = 1",
                "  crop_model_option                 = 0",
                "  snowsoil_temp_time_option         = 3",
                "  soil_temp_boundary_option         = 2",
                "  supercooled_water_option          = 1",
                "  stomatal_resistance_option        = 1",
                "  evap_srfc_resistance_option       = 4",
                "  subsurface_option                 = 2",
                "/\n",
            ]

            struct = [
                "&structure",
                f"  isltyp           = {soil_type}               ! soil texture class",
                "  nsoil            = 4               ! number of soil levels",
                "  nsnow            = 3               ! number of snow levels",
                "  nveg             = 27              ! number of vegetation types",
                f"  vegtyp           = {veg_type}               ! vegetation type",
                "  croptype         = 0               ! crop type (0 = no crops; this option is currently inactive)",
                "  sfctyp           = 1               ! land surface type, 1:soil, 2:lake",
                "  soilcolor       = 4               ! soil color code",
                "/\n"
            ]

            init_val = [
                "&initial_values",
                "  dzsnso    =  0.0,  0.0,  0.0,  0.1,  0.3,  0.6,  1.0     ! level thickness [m]",
                "  sice      =  0.0,  0.0,  0.0,  0.0                       ! initial soil ice profile [m3/m3]",
                "  sh2o      =  0.3,  0.3,  0.3,  0.3                       ! initial soil liquid profile [m3/m3]",
                "  zwt       =  -2.0                                        ! initial water table depth below surface [m]",
                "/\n",
            ]

            nom_params = timing + params + location + forcing + model_opt + struct + init_val

            fname_nom = f'nom_config_{cat_name}.input'
            nom_file = os.path.join(nom_dir, fname_nom)
            with open(nom_file, "w") as f:
                f.writelines('\n'.join(nom_params))

        
    def write_cfe_input_files(self,
                              surface_runoff_scheme="NASH_CASCADE",
                              sft_coupled = False):

        cfe_dir = os.path.join(self.output_dir, "configs/cfe")
        self.create_directory(cfe_dir)
        if "CFE-S" in self.formulation:
            surface_water_partitioning_scheme = "Schaake"
        elif "CFE-X" in self.formulation:
            surface_water_partitioning_scheme = "Xinanjiang"
            
        urban_decimal_fraction = 0.0
        delimiter = ","

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            fname = cat_name + '*.txt'

            cfe_params = [
                'forcing_file=BMI',
                f'surface_water_partitioning_scheme={surface_water_partitioning_scheme}',
                f'surface_runoff_scheme={surface_runoff_scheme}',
                'soil_params.depth=2.0[m]',
                f'soil_params.b={self.gdf["soil_b"][cat_name]}[]',
                f'soil_params.satdk={self.gdf["soil_satdk"][cat_name]}[m s-1]',
                f'soil_params.satpsi={self.gdf["soil_satpsi"][cat_name]}[m]',
                f'soil_params.slop={self.gdf["soil_slop"][cat_name]}[m/m]',
                f'soil_params.smcmax={self.gdf["soil_smcmax"][cat_name]}[m/m]',
                f'soil_params.wltsmc={self.gdf["soil_wltsmc"][cat_name]}[m/m]',
                'soil_params.expon=1.0[]',
                'soil_params.expon_secondary=1.0[]',
                f'refkdt={self.gdf["soil_refkdt"][cat_name]}',
                f'max_gw_storage={self.gdf["max_gw_storage"][cat_name]}[m]',
                f'Cgw={self.gdf["Cgw"][cat_name]}[m h-1]',
                f'expon={self.gdf["gw_expon"][cat_name]}[]',
                'gw_storage=0.5[m/m]',
                'alpha_fc=0.33',
                f'soil_storage={self.gdf["soil_smcmax"][cat_name]}[m/m]',
                'K_nash_subsurface=0.03[]',
                'N_nash_subsurface=2',
                'K_lf=0.01[]',
                'nash_storage_subsurface=0.0,0.0',
                'num_timesteps=1',
                'verbosity=0'
            ]

            if self.gdf['soil_b'][cat_name] == 1.0:
                cfe_params[4] = 1.1

            if surface_runoff_scheme == "GIUH" or surface_runoff_scheme == 1:
                giuh_cat = json.loads(self.gdf['giuh'][cat_name])
                giuh_cat = pd.DataFrame(giuh_cat, columns=['v', 'frequency'])
                giuh_ordinates = ",".join(str(x) for x in np.array(giuh_cat["frequency"]))
                cfe_params.append(f'giuh_ordinates={giuh_ordinates}')
            elif surface_runoff_scheme == "NASH_CASCADE" or surface_runoff_scheme == 2:
                cfe_params[2] = 'surface_runoff_scheme=NASH_CASCADE'
                cfe_params.append(f'N_nash_surface={int(self.gdf["N_nash_surface"][cat_name])}[]')
                cfe_params.append(f'K_nash_surface={self.gdf["K_nash_surface"][cat_name]}[h-1]')
                s = [str(0.0),] * int(self.gdf['N_nash_surface'][cat_name])
                s = delimiter.join(s)
                cfe_params.append(f'nash_storage_surface={s}[]')

            if surface_water_partitioning_scheme == 'Xinanjiang':
                cfe_params[1] = 'surface_water_partitioning_scheme=Xinanjiang'
                soil_id = self.gdf['ISLTYP'][cat_name]
                cfe_params.append(f'a_Xinanjiang_inflection_point_parameter={self.soil_class_NWM["AXAJ"][soil_id]}')
                cfe_params.append(f'b_Xinanjiang_shape_parameter={self.soil_class_NWM["BXAJ"][soil_id]}')
                cfe_params.append(f'x_Xinanjiang_shape_parameter={self.soil_class_NWM["XXAJ"][soil_id]}')
                cfe_params.append(f'urban_decimal_fraction={urban_decimal_fraction}')

            if sft_coupled:
                ice_content_threshold = 0.3
                cfe_params.append("sft_coupled=true")
                cfe_params.append(f"ice_content_threshold={ice_content_threshold}")

            fname_cfe = f'cfe_config_{cat_name}.txt'
            cfe_file = os.path.join(cfe_dir, fname_cfe)
            with open(cfe_file, "w") as f:
                f.writelines('\n'.join(cfe_params))

    def write_topmodel_input_files(self):

        topmodel_dir = os.path.join(self.output_dir, "configs/topmodel")
        self.create_directory(topmodel_dir)
        
        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            fname = cat_name + '*.txt'

            topmod = [
                "0",
                f'{cat_name}',
                f"./forcing/{cat_name}.csv",
                f'./{topmodel_dir}/subcat_{cat_name}.dat',
                f'./{topmodel_dir}/params_{cat_name}.dat',
                f'./{topmodel_dir}/topmod_{cat_name}.out',
                f'./{topmodel_dir}/hyd_{cat_name}.out'
            ]

            fname_tm = f'topmod_{cat_name}.run'
            tm_file = os.path.join(topmodel_dir, fname_tm)
            with open(tm_file, "w") as f:
                f.writelines('\n'.join(topmod))

            params = [
                f'Extracted study basin: {cat_name}',
                "0.032  5.0  50.  3600.0  3600.0  0.05  0.0000328  0.002  0  1.0  0.02  0.1"
            ]

            fname_tm = f'params_{cat_name}.dat'
            tm_file = os.path.join(topmodel_dir, fname_tm)
            with open(tm_file, "w") as f:
                f.writelines('\n'.join(params))

            twi_cat = json.loads(self.gdf['twi'][cat_name])
            twi_cat = pd.DataFrame(twi_cat, columns=['v', 'frequency'])
            twi_cat = twi_cat.sort_values(by=['v'], ascending=False)

            width_f = json.loads(self.gdf['width_dist'][cat_name])
            df_width_f = pd.DataFrame(width_f, columns=['v', 'frequency'])
            v_cumm = np.cumsum(df_width_f['frequency'])

            nclasses_twi = len(twi_cat['frequency'].values)
            nclasses_width_function = len(df_width_f['frequency'].values)

            subcat = [
                "1 1 1",
                f'Extracted study basin: {cat_name}',
                f'{nclasses_twi} 1',
                'replace_with_twi',
                f'{nclasses_width_function}',
                'add_width_function',
                '$mapfile.dat'
            ]

            twi_str = ''
            for freq, value in zip(twi_cat['frequency'].values, twi_cat['v'].values):
                twi_str += f"{freq:.6f} {value:.6f}\n"

            subcat[3] = twi_str.strip()

            widthf_str = ''
            for freq, value in zip(v_cumm.values, df_width_f['v'].values):
                widthf_str += f"{freq:.6f} {value:.6f} "

            subcat[5] = widthf_str.strip()

            fname_tm = f'subcat_{cat_name}.dat'
            tm_file = os.path.join(topmodel_dir, fname_tm)
            with open(tm_file, "w") as f:
                f.writelines('\n'.join(subcat))

    def write_sft_input_files(self,
                              surface_water_partitioning_scheme):

        sft_dir = os.path.join(self.output_dir, "configs/sft")
        if not surface_water_partitioning_scheme in ["Schaake", "Xinanjiang"]:
            sys.exit("Runoff scheme should be: Schaake or Xinanjiang")

        ncells = 19
        soil_z = "0.1,0.15,0.18,0.23,0.29,0.36,0.44,0.55,0.69,0.86,1.07,1.34,1.66,2.07,2.58,3.22,4.01,5.0,6.0"
        delimiter = ','
        nsteps_yr = 365 * 24

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            forcing_file = glob.glob(os.path.join(self.forcing_dir, cat_name + '*.csv'))[0]
            df_forcing = pd.read_csv(self.forcing_file, delimiter=',', usecols=['T2D'], nrows=nsteps_yr, index_col=None)
            MAAT = [str(round(df_forcing['T2D'].mean(), 2)),] * ncells
            MAAT = delimiter.join(MAAT)
            soil_id = self.gdf['ISLTYP'][cat_name]

            sft_params = [
                'verbosity=none',
                'soil_moisture_bmi=1',
                'end_time=1.0[d]',
                'dt=1.0[h]',
                f'soil_params.smcmax={self.gdf["soil_smcmax"][cat_name]}[m/m]',
                f'soil_params.b={self.gdf["soil_b"][cat_name]}[]',
                f'soil_params.satpsi={self.gdf["soil_satpsi"][cat_name]}[m]',
                f'soil_params.quartz={self.soil_class_NWM["QTZ"][soil_id]}[]',
                f'ice_fraction_scheme={surface_water_partitioning_scheme}',
                f'soil_z={soil_z}[m]',
                f'soil_temperature={MAAT}[K]'
            ]

            fname_sft = f'sft_config_{cat_name}.txt'
            sft_file = os.path.join(sft_dir, fname_sft)
            with open(sft_file, "w") as f:
                f.writelines('\n'.join(sft_params))

    def write_smp_input_files(self, cfe_coupled, lasam_coupled=False):
        
        smp_dir = os.path.join(self.output_dir, "configs/smp")
        
        soil_z = "0.1,0.15,0.18,0.23,0.29,0.36,0.44,0.55,0.69,0.86,1.07,1.34,1.66,2.07,2.58,3.22,4.01,5.0,6.0"

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            soil_id = self.gdf['ISLTYP'][cat_name]

            smp_params = [
                'verbosity=none',
                f'soil_params.smcmax={self.gdf["soil_smcmax"][cat_name]}[m/m]',
                f'soil_params.b={self.gdf["soil_b"][cat_name]}[]',
                f'soil_params.satpsi={self.gdf["soil_satpsi"][cat_name]}[m]',
                f'soil_z={soil_z}[m]',
                'soil_moisture_fraction_depth=1.0[m]'
            ]

            if cfe_coupled:
                smp_params += ['soil_storage_model=conceptual', 'soil_storage_depth=2.0']
            elif lasam_coupled:
                smp_params += ['soil_storage_model=layered', 'soil_moisture_profile_option=constant',
                               'soil_depth_layers=2.0', 'water_table_depth=10[m]']

            fname_smp = f'smp_config_{cat_name}.txt'
            smp_file = os.path.join(smp_dir, fname_smp)
            with open(smp_file, "w") as f:
                f.writelines('\n'.join(smp_params))

    def write_lasam_input_files(self,
                                sft_coupled=False):
        
        lasam_dir = os.path.join(self.output_dir,"configs/lasam")
        create_directory(lasam_dir)

        lasam_params = os.path.join(self.ngen_dir,"extern/LGAR-C/LGAR-C/data/vG_params_stat_nom_ordered.dat")
        str_sub ="cp -r "+ lasam_params + " %s"%lasam_dir
        out=subprocess.call(str_sub,shell=True)

        sft_calib = "False"
        soil_z = "10.0,15.0,18.0,23.0,29.0,36.0,44.0,55.0,69.0,86.0,107.0,134.0,166.0,207.0,258.0,322.0,401.0,500.0,600.0"
        
        lasam_params_base = [
            'verbosity=none',
            f'soil_params_file={soil_param_file}',
            'layer_thickness=200.0[cm]',
            'initial_psi=2000.0[cm]',
            'timestep=3600[sec]',
            'endtime=1000000000.0[d]',
            'forcing_resolution=3600[sec]',
            'ponded_depth_max=0[cm]',
            'use_closed_form_G=false',
            'layer_soil_type=',
            'wilting_point_psi=15495.0[cm]',
            'field_capacity_psi=340.9[cm]',
            'adaptive_timestep=true',
            'giuh_ordinates='
        ]

        if sft_coupled:
            lasam_params_base.append('sft_coupled=true')
            lasam_params_base.append(f'soil_z={soil_z}[cm]')

        if (sft_coupled and (sft_calib in ["true", "True"])):
            lasam_params_base.append('calib_params=true')
            
        if self.ngen_cal_type in ['calibration', 'validation', 'restart']:
            lasam_params_base.append('calib_params=true')
        
        soil_type_loc = lasam_params_base.index("layer_soil_type=")
            
        giuh_loc_id = lasam_params_base.index("giuh_ordinates=")
            

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            fname = cat_name + '*.txt'

            lasam_params = lasam_params_base.copy()
            lasam_params[soil_type_loc] += str(self.gdf['ISLTYP'][cat_name])

            giuh_cat = json.loads(self.gdf['giuh'][cat_name])
            giuh_cat = pd.DataFrame(giuh_cat, columns=['v', 'frequency'])
            giuh_ordinates = ",".join(str(x) for x in np.array(giuh_cat["frequency"]))

            any_nans = np.any(np.isnan(giuh_cat["frequency"]))
            if any_nans:
                giuh_ordinates = str(1.0)

            lasam_params[giuh_loc_id] += giuh_ordinates

            fname_lasam = f'lasam_config_{cat_name}.txt'
            lasam_file = os.path.join(lasam_dir, fname_lasam)
            with open(lasam_file, "w") as f:
                f.writelines('\n'.join(lasam_params))


    def write_pet_input_files(self):
        pet_dir = os.path.join(self.output_dir,"configs/pet")
        self.create_directory(pet_dir)
                
        df_cats = gpd.read_file(self.gpkg_file, layer='divides')
        df_cats = df_cats.to_crs("EPSG:4326")
        df_cats.set_index("divide_id", inplace=True)
        pet_method = 3

        for catID in self.catids:
            cat_name = 'cat-' + str(catID)
            centroid_x = str(df_cats['geometry'][cat_name].centroid.x)
            centroid_y = str(df_cats['geometry'][cat_name].centroid.y)
            elevation_mean = self.gdf['elevation_mean'][cat_name]

            pet_params = [
                'verbose=0',
                f'pet_method={pet_method}',
                'forcing_file=BMI',
                'run_unit_tests=0',
                'yes_aorc=1',
                'yes_wrf=0',
                'wind_speed_measurement_height_m=10.0',
                'humidity_measurement_height_m=2.0',
                'vegetation_height_m=16.0',
                'zero_plane_displacement_height_m=0.0003',
                'momentum_transfer_roughness_length=0.0',
                'heat_transfer_roughness_length_m=0.0',
                'surface_longwave_emissivity=1.0',
                'surface_shortwave_albedo=0.17',
                'cloud_base_height_known=FALSE',
                'time_step_size_s=3600',
                'num_timesteps=720',
                'shortwave_radiation_provided=1',
                f'latitude_degrees={centroid_y}',
                f'longitude_degrees={centroid_x}',
                f'site_elevation_m={elevation_mean}'
            ]

            fname_pet = f'pet_config_{cat_name}.txt'
            pet_file = os.path.join(pet_dir, fname_pet)
            with open(pet_file, "w") as f:
                f.writelines('\n'.join(pet_params))

    def write_troute_input_files(self):

        troute_basefile = os.path.join(self.workflow_dir, "configs/basefiles/config_troute.yaml")
        troute_dir = os.path.join(self.output_dir,"configs")
        gpkg_name = os.path.basename(self.gpkg_file).split(".")[0]

        if not os.path.exists(troute_basefile):
            sys.exit(f"Sample routing yaml file does not exist, provided is {troute_basefile}")

        with open(troute_basefile, 'r') as file:
            d = yaml.safe_load(file)

        d['network_topology_parameters']['supernetwork_parameters']['geo_file_path'] = self.gpkg_file
        d['network_topology_parameters']['waterbody_parameters']['level_pool']['level_pool_waterbody_parameter_file_path'] = self.gpkg_file
        d['network_topology_parameters']['supernetwork_parameters']['title_string'] = gpkg_name

        dt = 300
        params = self.get_flowpath_attributes(gage_id=self.gpkg_file, full_schema=True)

        columns = {
            'key': params['key'],
            'downstream': params['downstream'],
            'mainstem': params['mainstem'],
            'dx': params['dx'],
            'n': params['n'],
            'ncc': params['ncc'],
            's0': params['s0'],
            'bw': params['bw'],
            'waterbody': params['waterbody'],
            'gages': params['gages'],
            'tw': params['tw'],
            'twcc': params['twcc'],
            'musk': params['musk'],
            'musx': params['musx'],
            'cs': params['cs'],
            'alt': params['alt']
        }

        d['network_topology_parameters']['supernetwork_parameters']['columns'] = columns

        start_time = pd.Timestamp(self.simulation_time['start_time'])
        end_time = pd.Timestamp(self.simulation_time['end_time'])
        diff_time = (end_time - start_time).total_seconds()

        d['compute_parameters']['restart_parameters']['start_datetime'] = start_time.strftime("%Y-%m-%d_%H:%M:%S")

        if self.ngen_cal_type in ['calibration', 'validation', 'calibvalid', 'restart']:
            d['compute_parameters']['forcing_parameters']['qlat_input_folder'] = "./"
        else:
            d['compute_parameters']['forcing_parameters']['qlat_input_folder'] = os.path.join(self.output_dir, "div")

        d['compute_parameters']['forcing_parameters']['qlat_file_pattern_filter'] = "nex-*"
        del d['compute_parameters']['forcing_parameters']['binary_nexus_file_folder']
        d['compute_parameters']['forcing_parameters']['nts'] = int(diff_time / dt)
        d['compute_parameters']['forcing_parameters']['max_loop_size'] = 10000000

        d['compute_parameters']['cpu_pool'] = 10

        if self.ngen_cal_type in ['calibration', 'validation', 'calibvalid', 'restart']:
            stream_output = {
                "stream_output": {
                    "stream_output_directory": "./",
                    'stream_output_time': -1,
                    'stream_output_type': '.nc',
                    'stream_output_internal_frequency': 60
                }
            }
        else:
            stream_output = {
                "stream_output": {
                    'stream_output_directory': os.path.join(self.output_dir, "troute"),
                    'stream_output_time': -1,
                    'stream_output_type': '.nc',
                    'stream_output_internal_frequency': 60
                }
            }

        d['output_parameters'] = stream_output

        with open(os.path.join(troute_dir, "troute_config.yaml"), 'w') as file:
            yaml.dump(d, file, default_flow_style=False, sort_keys=False)


    def get_flowpath_attributes(self,
                                full_schema=False,
                                gage_id=False):

        layers = fiona.listlayers(self.gpkg_file)
        flowpath_layer = [layer for layer in layers if 'flowpath' in layer and not 'flowpaths' in layer][0]
        gdf_fp_attr = gpd.read_file(self.gpkg_file, layer=flowpath_layer)
        params = schema.get_schema_flowpath_attributes(gdf_fp_attr, for_gage_id=gage_id)

        if full_schema:
            return params
        elif gage_id:
            gage_id = params['gages']
            waterbody_id = params['key']
            gdf_fp_cols = gdf_fp_attr[[waterbody_id, gage_id]]
            basin_gage = gdf_fp_cols[gdf_fp_cols[gage_id].notna()]
            basin_gage_id = basin_gage[waterbody_id].tolist()
            return basin_gage_id

    def write_forcing_input_files(self,
                                  forcing_basefile,
                                  forcing_time,
                                  forcing_format):

        if not os.path.exists(forcing_basefile):
            sys.exit(f"Sample forcing yaml file does not exist, provided is {forcing_basefile}")

        with open(forcing_basefile, 'r') as file:
            d = yaml.safe_load(file)

        time_sim = json.loads(forcing_time)
        start_yr = pd.Timestamp(time_sim['start_time']).year
        end_yr   = pd.Timestamp(time_sim['end_time']).year

        if start_yr > end_yr:
            sys.exit(f"end_time ({end_yr})is less than the start_time ({start_yr}")

        if start_yr <= end_yr:
            end_yr = end_yr + 1

        d['gpkg'] = self.gpkg_file
        d["years"] = [start_yr, end_yr]
        d["out_dir"] = os.path.join(os.path.dirname(self.gpkg_file), "forcing")

        out_dir = Path(d['out_dir']) / f'{start_yr}_to_{end_yr}'
        are_identical = out_dir.resolve() == Path(self.forcing_dir).resolve()

        if not are_identical:
            raise RuntimeError(f"Directory mismatch: out_dir={out_dir} is not the same as forcing_dir={self.forcing_dir}.")

        if not os.path.exists(d["out_dir"]):
            os.makedirs("data/forcing")

        if forcing_format == '.csv':
            d['netcdf'] = False

        with open(os.path.join(d["out_dir"], "config_forcing.yaml"), 'w') as file:
            yaml.dump(d, file, default_flow_style=False, sort_keys=False)

        return os.path.join(d["out_dir"], "config_forcing.yaml")

    def create_directory(self, dir_name):
        if os.path.exists(dir_name):
            str_sub = "rm -rf " + dir_name
            out = subprocess.call(str_sub, shell=True)
        os.mkdir(dir_name)

class ConfigurationCalib:
    def __init__(self, gpkg_file, output_dir, ngen_dir, realization_file_par, troute_output_file,
                 ngen_cal_basefile, ngen_cal_type,
                 restart_dir, simulation_time, evaluation_time, num_proc):
        
        self.gpkg_file = gpkg_file
        self.output_dir = output_dir
        self.ngen_dir = ngen_dir
        self.realization_file_par = realization_file_par
        self.simulation_time = simulation_time
        self.evaluation_time = evaluation_time
        #self.verbosity = verbosity
        self.ngen_cal_type = ngen_cal_type
        self.num_proc = num_proc
        self.ngen_cal_basefile = ngen_cal_basefile
        self.troute_output_file = troute_output_file

    def get_flowpath_attributes(self):

        layers = fiona.listlayers(self.gpkg_file)
        flowpath_layer = [layer for layer in layers if 'flowpath' in layer and not 'flowpaths' in layer][0]
        gdf_fp_attr = gpd.read_file(self.gpkg_file, layer=flowpath_layer)
        params = schema.get_schema_flowpath_attributes(gdf_fp_attr, for_gage_id=True)


        gage_id = params['gages']
        waterbody_id = params['key']
        gdf_fp_cols = gdf_fp_attr[[waterbody_id, gage_id]]
        basin_gage = gdf_fp_cols[gdf_fp_cols[gage_id].notna()]
        basin_gage_id = basin_gage[waterbody_id].tolist()
        return basin_gage_id
        
    def write_calib_input_files(self):
        
        conf_dir = os.path.join(self.output_dir, "configs")
        realization = glob.glob(os.path.join(self.output_dir, "json/realization_*.json"))
        print ("OO : ", self.output_dir, realization)
        assert len(realization) == 1

        if not os.path.exists(self.ngen_cal_basefile):
            sys.exit(f"Sample calib yaml file does not exist, provided is {self.ngen_cal_basefile}")

        basin_workflow_dir = os.path.dirname(os.path.dirname(self.ngen_cal_basefile))
        gpkg_name = os.path.basename(self.gpkg_file).split(".")[0]

        with open(self.ngen_cal_basefile, 'r') as file:
            d = yaml.safe_load(file)

        d['general']['workdir']   = self.output_dir.as_posix()
        d['model']['binary']      = os.path.join(self.ngen_dir, "cmake_build/ngen")
        d['model']['realization'] = realization[0]
        d['model']['hydrofabric'] = self.gpkg_file.as_posix()
        d['model']['routing_output'] = self.troute_output_file

        gage_id = self.get_flowpath_attributes()

        if len(gage_id) == 1:
            d['model']['eval_feature'] = gage_id[0]
        else:
            print("more than one rl_gages exist in the geopackage, using max drainage area to filter...")
            div = gpd.read_file(self.gpkg_file, layer='divides')
            df = div[['divide_id', 'tot_drainage_areasqkm']]
            index = df['divide_id'].map(lambda x: 'wb-' + str(x.split("-")[1]))
            df.set_index(index, inplace=True)
            idmax = df['tot_drainage_areasqkm'].idxmax()
            d['model']['eval_feature'] = idmax

        if self.num_proc > 1:
            d['model']['parallel'] = self.num_proc
            d['model']['partitions'] = self.realization_file_par

        if os_name == "Darwin":
            d['model']['binary'] = f'PYTHONEXECUTABLE=$(which python) ' + os.path.join(self.ngen_dir, "cmake_build/ngen")
        else:
            d['model']['binary'] = os.path.join(self.ngen_dir, "cmake_build/ngen")

        if self.ngen_cal_type == 'calibration':
            val_params = {
                'evaluation_start': self.evaluation_time['start_time'],
                'evaluation_stop': self.evaluation_time['end_time']
            }
            d['model']['eval_params'].update(val_params)

        if self.ngen_cal_type == 'validation':
            val_troute_output = {
                'ngen_cal_troute_output': {
                    'validation_routing_output': self.troute_output_file
                }
            }

            try:
                d['model']['plugin_settings'].update(val_troute_output)
            except:
                d['model']['plugin_settings'] = val_troute_output

            val_params = {
                'sim_start': self.simulation_time['start_time'],
                'evaluation_start': self.evaluation_time['start_time'],
                'evaluation_stop': self.evaluation_time['end_time'],
                'objective': "kling_gupta"
            }
            d['model']['val_params'] = val_params

        elif self.ngen_cal_type == 'restart':
            df_par = pd.read_parquet(os.path.join(restart_dir, "calib_param_df_state.parquet"))
            df_params = pd.read_csv(os.path.join(restart_dir, "best_params.txt"), header=None)
            best_itr = str(int(df_params.values[1]))

            best_params_set = df_par[best_itr]
            calib_params = best_params_set.index.to_list()

            for block in d:
                if '_params' in block:
                    for par in d[block]:
                        if par['name'] in calib_params:
                            par['init'] = float(best_params_set[par['name']])

        config_fname = ""
        if self.ngen_cal_type == 'calibration':
            config_fname = "ngen-cal_calib_config.yaml"
        elif self.ngen_cal_type == 'validation':
            config_fname = "ngen-cal_valid_config.yaml"

        with open(os.path.join(conf_dir, config_fname), 'w') as file:
            yaml.dump(d, file, default_flow_style=False, sort_keys=False)
