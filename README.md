# Data Analysis for the North Slope of Alaska ARM site

This repository is part of the [(AC)3 project](https://wwww.ac3-tr.de), sub-cluster B07. The repository contains scripts to analyze the data from the ARM NSA site in Utqiaǵvik, Alaska.

The data analysis comprise of the cloud remote sensing instrumentation for cloud phase classification. Main instruments considered are:
* KAZR cloud radar,
* MWR microwave radiometer,
* CEIL10m cloud lidar (ceilometer),
* HSRL High spectral resolution Lidar,
* GNDTIR Ground Infrared thermometer,
* INTERPOLATESONDE Radiosonde data interpolated to fix grid,
* RADFLUX shortwave and longwave up- and down-welling radiometers.

The objective is to apply the methodology descrived in [Saavedra Garfias et al. (2023)](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-623/) to the NSA site and characterize cloud properties under the influence of the sea ice conditions observed at the surrounding area of NSA. Sea ice concentration data is obtained by the satellite product ASI by AMSR2 provided by the [https://seaice.uni-bremen.de](University of Bremen).

## Chain of Processing
### 1. Download data
First download the model files needed for Cloudnet, which can be done using the API from [https://cloudnet.fmi.fi](docs.cloudnet.fmi.fi)

### Cloudnet
To convert the ARM data files to be used as input for Cloudnet, user the script: ```Cloudnet4NSA/run_convertor.jl``` by editing the source code and selecting the convination of instruments to use:

Once the NSA input files are converted to Cloudnet complaiant format, run the Cloudnet algorithm using the script: ```Cloudnet4NSA/Process_cloudnet4ARM.jl```. This will use the converted ARM data for lidar, mwr, radar, and model data as indicated in the Dictionary:
```julia
# defining ARM product to be used:
ARMprod = Dict(
    :site => "utqiagvik-nsa", 					#"mosaic",
    :radar => ("KAZR/ARSCL","KAZR/CORGE"),
    :lidar => ("CEIL10m", "HSRL"),
    :mwr => ("MWR/RET","MWR/LOS"),
    :model => "ECMWF",
    :radiosonde => "INTERPOLATEDSONDE",
)
```
In case any of the required instrument data files is missing, the script will skip the Cloudnet calculation for that day and indicate it as a message. Not all ARM sited have model data privided by the Cloudnet servers, for example the ARM site for MOSAiC is not available in the Cloudnet servers, then as alternative it can be use the ARM ```INTERPOLATEDSONDE``` data to recreate Cloudnet input as "model" data with the same parameters.

The default output path for Cloudnet NetCDF files is indicated below in the sub-folder description. The output path as well as the Clounet retrievals can be specified in source code ```Cloudnet4NSA/Process_cloudnet4ARM.jl``` by editing the variable:

```julia
CLNTprod = Dict(
    :categorize => true,   		# for cloud categorization output (must be created for further analysis).
    :classification => false,	# for cloud classification output (must be created for further analysis).
    :lwc => true,				# cloud liquid water content output.
    :iwc => true,				# cloud ice water content output.
    :drizzle => true,			# precipitating water classification  as drizzle output.
    :der => true,				# cloud droplets effective radius output.
    :ier => true,				# cloud ice particle effective radius output.
)
```

### 3. Run atmospheric macro- micro-physical database.
Once Cloudnet categorization output is completed, the relevant atmospheric parameters are computed using the script ```CoupledCloud_Seaice/run_processing2daily_csv.jl```.
In the script the relevant variables are ```winter_jahr::Vector{Int}``` where the range of wintertimes is specified, for example 2023 for the wintertime from November 2023 to April 2024.
From that variable the winter time is assigned as coupled variables for year and month, e.g. ```datum = Date(year, month)``` where  ```year``` ranges from 2023 to 2024 and ```month``` ranges from ```11:4``` representing the range of months from November to April next year.

The output is a daily CSV files with atmospheric parameters in columns at a common temporal resolution of 1 minute.

### 4. Extract Sea Ice data for the region of interest.
Once the Sea Ice data is download for the same period of atmospheric parameters computed in CSV daily files (previous section), the relevant sea ice data is extracted using the script: ```SeaIce/distribution_seaice_winddir.jl``` 

The main parameters to edit in the script are: 
```julia
# Define coordinates for the North Slope Alaska site (or any other place of interest):
nsa_lat = 71.323e0;
nsa_lon = -156.609e0;
# Define the azimuth angles for the sector to consider (avoiding Land):
θₗ₀ = 235e0;
θₗ₁ = 110e0;

# Define the radiuns around the center to calculate Sea Ice (in m):
const R_lim = 50e3;   # radius around NSA location, e.g. 50 km

# det the satellite to use (previously downloaded Sea Ice data):
const SATELLITE = "amsr2";  # "ssmis"; #  (SSMIS is for 2011.11 to 2012.04)

# Define the Sea Ice product to consider, for AMSR2 for example SIC.
PRODUCTS = (:SIC,) # (:DIV, :LF, :SIC)
# Define the wintertime yearss, e.g. 2023:2024 for wintertimes
# from November 2023 to April 2024 and November 2024 to April 2025
winter_jahr = 2023:2024;

```

The script outputs daily Julida Language Data 2 files (JLD2) with a minute time resolution with the information about mean Sea Ice within +/- 3° sector with 50 km radius and azimuth modulated by the maximum water vapour transport direction (optained from the previous CSV daily files). 

### 5. Joining dataset from atmospheric parameters and Sea Ice into a common single CSV file
The Sea Ice and atmospheric paramters data set are merged into a single file using the script ```create_nsa_microphys_db.jl```.

This script needs:
* daily CSV data files after processing for WVT and cloud layers from Cloudnet
* JLD2 sea ice concentration averages at the direction of the WVT

This script will output:
* a yearly CSV data file e.g.  all_nsa_microphys_db_2012-2013.csv
which contains the database from wintertime 2012 Nov to 2013 April with
atmospheric, cloud, and sea ice information merged.

The yearly CSV files can further be merged in a range of years using the script:

## List of files and description:
* ```SeaIce/distribution\_seaice\_winddir.jl``` 
	Computes SIC statistics using AMSR2 ASI h5 daily files. The script needs the input from daily CSV database post-processed Cloudnet files where the azimuth of the maximum WVT profile is stored. The statistics are stored as daily JLD2 files in:
	```SeaIce/data/SIC/2011```

* ```CoupledCloud\_Seaice/run\_processing2daily\_csv.jl``` 
	Process the Cloudnet output files into daily CSV dataset files including atmospheric characterizations like max WVT, direction of WVT, Ri_b, etc. (see paper). Required the Cloudnet categorization, classification, lwc, iwc, der, ier output files.

* ```CoupledCloud\_Seaice/create\_nsa\_microhys\_db.jl```
	Joins the daily CSV cloudnet post-processed files with the SIC daily JLD2 files into a single CSV yearly dataset.

* ```CoupledCloud\_Seaice/join\_years\_microphys\_db.jl``` 
	create a single CSV dataset from multiple years. It needs the yearly CSV files.

## Sub-folders description:

Ancillary data needs to be located in the folders as follow (otherwise full pah need to be specified in Pluto notebook):
* ``` CoupledCloud_Seaice/buffer_data/csv_nsa/```
    * ```all_nsa_microphys_db_2011-2025.csv``` : CSV file with post-processed data base for all years.
* ``` CoupledCloud_Seaice/buffer_data/clima/```
  * ```enso_index.json```: Data file with the ENSO index,
  * ```pdo_index.json``` : Data file with the PDO index,
  * ```norm.daily.ao.cdas.z1000.19500101_current.csv```: Data file with AO index.	 

# Source data location description
The ARM source data is located at the LIM servers in the RemArc working group location:
```/projekt7/remsens/data_new/site-campaign/utquiagvik-nsa```
The data is sorted in the above folder according to ```instrument/{product}/year``` sub-folders, e.g. ```KAZR/ARSCL/2023``` or ```CEIL10m/2023```.

---
(c) Pablo Saavedra Garfias<br>
[pablo.saavedra@uni-leipzig.de](mailto:pablo.saavedra@uni-leipzig.de)<br>
LIM<br>
Faculty of Physics and Geosciences<br>
University of Leipzig<br>
