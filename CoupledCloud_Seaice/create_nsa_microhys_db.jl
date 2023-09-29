#!/home/psgarfias/.local/bin/julia
#
#=
Script to join into a single CSV data file the daily dataset from
Cloudnet microphysical properties and ARMSR2 sea ice concentration.
=#
using StatsPlots
using Statistics
using DataFrames
using CSV
using Dates
using JLD2
using ARMtools
# Statistical analysis for NSA coupled/decoupled micro-physical properties

# Defining data path
const BASE_PATH = joinpath(homedir(), "LIM/scripts/NSA-ac3-b07")
const DATA_PATH = joinpath(BASE_PATH, "CoupledCloud_Seaice/data")
const LFSIC_PATH = joinpath(BASE_PATH, "SeaIce/data")
const MIPHY_PATH = joinpath(DATA_PATH, "csv_nsa")

# defining the wintertime to process e.g. for year yy:Nov, Dec to yy+1:Jan, Feb, Mar, Apr.
years = (2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)

DB = DataFrame()

for jahr ∈ years
    zu_jahr = jahr + 1
    datum = Date(jahr, 11):Day(1):Date(zu_jahr, 4, 30)
    
    for heute ∈ datum
        
        yy = year(heute)
        mm = month(heute)
        dd = day(heute)

        # listing all available SeaIce_winddir files 
        lfsic = ARMtools.getFilePattern(LFSIC_PATH, "SIC", yy, mm, dd, fileext=".jld2") #readdir(LFSIC_PATH);

        # llisting all available micro-physical files:
        fn = ARMtools.getFilePattern(DATA_PATH, "csv_nsa", yy, mm, dd, fileext="_I.csv")  #readdir(MIPHY_PATH);


        ## NN = length(list_lfsic)

        ## for lfsic in list_lfsic
        # reading the CSV dataset
        ##		fn = joinpath(MIPHY_PATH, lfsic[8:23]*"_I.csv")
        isnothing(lfsic) && (@warn(" $(heute), sic file $(lfsic) gives nothing. "); continue)
        isnothing(fn) && (@warn("Data file $(fn) gives nOthing!"); continue)
               
        cldphys = CSV.read(fn, types=Dict(:coupled=>Bool), DataFrame)

        #reading the jld2 dataset
        seaice = let fn = load_object(lfsic) # joinpath(LFSIC_PATH, lfsic))
            Dict(Symbol(V,db)=>fn[db][V] for db in (:SIC,) for V in (:μ, :σ, :Aμ, :Aσ, :Rμ, :Rσ) ) |> DataFrame
        end
        @assert size(seaice,1)==size(cldphys,1) "the two DB have not the same length $(lfsic)"
        cldphys = hcat(cldphys, seaice) #|> F->filter(row->!isnan(row.lwp), F)
        append!(DB, cldphys, promote=true)
    end

    CSV.write(joinpath(DATA_PATH, "all_nsa_microphys_db_$(jahr)-$(zu_jahr).csv"), DB)
end
