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

The objective is to apply the methodology descrived in [Saavedra Garfias et al. (2023)](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-623/) to the NSA site and characterize cloud properties under the influence of the sea ice conditions observed at the surrounding area of NSA. Sea ice concentration data is obtained by the satellite product ASI by AMSR2 provided by the [hhtps://seaice.uni-bremen.de](University of Bremen).

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
