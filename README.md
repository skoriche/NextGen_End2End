# Next-Generation Framework Sandbox Hub (NextGenSandboxHub)
[NextGen](https://github.com/NOAA-OWP/ngen), Next-Generation Water Resources Modeling Framework, developed by the NOAA's Office of Water Prediction is a standards-based language- and model-agnostic framework, which allows to run a mosaic of surface and subsurface models in a single basin comprised of 10s-100s sub-catchments. 

# Schematic of the NextGenSandboxHub workflow

<div align="center">
<img src="https://github.com/user-attachments/assets/d06b3cf9-6019-4ebd-86f1-e797b4debbae" style="width:800px; height:400px;"/>
</div>

### Running the NextGenSandboxHub workflow
```
python <path_to_NextGenSandboxHub>/main.py option
OPTIONS = [-gpkg -forc -conf -run]
```
Note: These options can be run individually or all together by `path_to/main.py -gpkg -conf -run`. The `-gpkg` is an expensive step, should be run once to get the desired basin geopacakge.


### Configuration steps
Note: The workflow assumes that [ngen](https://github.com/NOAA-OWP/ngen) and other models have already been built, including [t-route](https://github.com/NOAA-OWP/t-route). It is also assumed that both ngen and t-route were built within a Python environment named `.venv_ngen`.
#### Workflow setup
  - `git clone https://github.com/ajkhattak/NextGenSandboxHub && cd NextGenSandboxHub`
  - `git submodule update --init`
  - `source <path_to_venv>/.venv_ngen/bin/activate`
  - `pip install 'extern/ngen-cal/python/ngen_cal[netcdf]'`
  - `pip install -e ./extern/ngen_cal_plugins`
    
#### [Hydrofabric](https://github.com/NOAA-OWP/hydrofabric) installation
Install hydrofabric and related packages (assumes R and Rtools are already installed)
  - Download domain (CONUS or oCONUS) from [lynker-spatial](https://www.lynker-spatial.com/data?path=hydrofabric%2Fv2.2%2F), for instance conus/conus_nextgen.gpkg
  - open configs/config_workflow.yaml [here](configs/config_workflow.yaml) and adjust workflow_dir, input_dir, output_dir, and gpkg_model_params according to your local settings
  - python main.py -gpkg (or open src_r/main.R in RStudio and set infile_config [here](https://github.com/ajkhattak/basin_workflow/blob/nwm-v4-bm/src_r/main.R#L54) and run main.R. This will install the hydrofabric and several other libraries, and if everything goes well, a basin geopackage will be subsetted and stored under `<output_dir>/<basin_id>/data/gage_<basin_id>.gpkg`
    
#### Forcing data
The workflow uses [CIROH_DL_NextGen](https://github.com/ajkhattak/CIROH_DL_NextGen) forcing_prep tool to donwload atmospheric forcing data.
  - `mkdir ~/.venv_forcing`
  - `python -m venv ~/.venv_forcing`
  - `source ~/.venv_forcing/bin/activate`
  - `pip install -r extern/CIROH_DL_NextGen/forcing_prep/requirements.txt`
  - open configs/config_workflow.yaml [here](configs/config_workflow.yaml) and setup `forcings` block
  - `python main.py -forc`

### Generate configuration and realization files
The workflow needs two configuration files, provided [here](https://github.com/ajkhattak/basin_workflow/tree/tnc-dangermond/configs/). Workflow setup and model options and paths need to be adjusted to local settings. Please see the configuration files for further details.



- Option: `-gpkg` downloads geopackage(s) given a gage ID(s), computes TWI, GIUH, and Nash parameters, and append them to the geopackage along with other model parameters (from S3 bucket) as `model-attributes`
- Option: `-conf` generates configuration files for the selected models/basins
- Option: `-run` runs NextGen simulations with and without calibration. The workflow uses [ngen-cal](https://github.com/NOAA-OWP/ngen-cal) for calibration


####################################################################
### Instructions for computing Geomorphological Instantaneous Unit Hydrograph (GIUH) and Topographic Wetness Index (TWI)
This set of scripts offers a workflow for establishing basin-scale simulations from the ground up 
and executing them within the NextGen framework.

- R-based scripts, leveraging the [WhiteBox](https://www.whiteboxgeo.com/manual/wbw-user-manual/book/tool_help.html) tool, a geospatial data analysis software,
  are provided for the computation of GIUH and TWI.
  
- Python-based scripts are provided for generating model(s) (standalone and coupled) configuration files and the
  NextGen realization file.


