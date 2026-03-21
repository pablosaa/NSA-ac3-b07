#!/home/psgarfias/.local/bin/julia
#
#=
Script to join into a daily CSV data files into a singe winterly dataset.
The dataset is comprised of Cloudnet microphysical properties and ARMSR2 sea ice concentration
as a function of the water vapor transport direction.

This script needs:
* daily CSV data files after processing for WVT and cloud layers from Cloudnet
* JLD2 sea ice concentration averages at the direction of the WVT

This script will output:
* a yearly CSV data file e.g.  all_nsa_microphys_db_2012-2013.csv
which contains the database from wintertime 2012 Nov to 2013 April with
atmospheric, cloud, and sea ice information merged.
=#
using StatsPlots
using Statistics
using DataFrames
using CSV
using Dates
using JLD2
using ARMtools
# Statistical analysis for NSA coupled/decoupled micro-physical properties

# Defining the JLD2 files with SIC radius of e.g. 50km, 75km, or 100km
const Rsic = "SIC100km"

# Defining data path
const BASE_PATH = joinpath(homedir(), "LIM/scripts/NSA-ac3-b07")
const DATA_PATH = "/projekt2/ac3data/B07-data/utqiagvik-nsa/" #joinpath(BASE_PATH, "CoupledCloud_Seaice/data")
const LFSIC_PATH = joinpath(DATA_PATH, "SeaIce")
const MIPHY_PATH = joinpath(DATA_PATH, "csv_nsa", Rsic)


# defining the wintertime to process e.g. for year yy:Nov, Dec to yy+1:Jan, Feb, Mar, Apr.
years = (2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)


for jahr ∈ years
    zu_jahr = jahr + 1
    datum = Date(jahr, 11):Day(1):Date(zu_jahr, 4, 30)
    
    DB = DataFrame()

    for heute ∈ datum
        
        yy = year(heute)
        mm = month(heute)
        dd = day(heute)

        # listing all available SeaIce_winddir files 
        lfsic = ARMtools.getFilePattern(LFSIC_PATH, Rsic, yy, mm, dd, fileext=".jld2") #readdir(LFSIC_PATH);

        # llisting all available micro-physical files:
        fn = ARMtools.getFilePattern(DATA_PATH, "csv_nsa", yy, mm, dd, fileext="_I.csv")  #readdir(MIPHY_PATH);


        ## NN = length(list_lfsic)

        ## for lfsic in list_lfsic
        # reading the CSV dataset
        ##		fn = joinpath(MIPHY_PATH, lfsic[8:23]*"_I.csv")
        isnothing(lfsic) && (@warn(" $(heute), sic file $(lfsic) gives nothing. "); continue)
        isnothing(fn) && (@warn("Data file $(fn) gives nothing!"); continue)
               
        cldphys = CSV.read(fn, types=Dict(:coupled=>Bool), DataFrame)

        #reading the jld2 dataset
        seaice = let fn = load_object(lfsic) # joinpath(LFSIC_PATH, lfsic))
            Dict(Symbol(V,db)=>fn[db][V] for db in (:SIC,) for V in (:μ, :σ, :Aμ, :Aσ, :Rμ, :Rσ) ) |> DataFrame
        end
        @assert size(seaice,1)==size(cldphys,1) "the two DB have not the same length $(lfsic)"
        cldphys = hcat(cldphys, seaice) #|> F->filter(row->!isnan(row.lwp), F)
        append!(DB, cldphys, promote=true)
    end

    CSV.write(joinpath(MIPHY_PATH, "yearly", "all_nsa_microphys_db_$(jahr)-$(zu_jahr).csv"), DB)
end
