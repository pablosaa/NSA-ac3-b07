#!/home/psgarfias/.local/bin/julia
#
#=
Script to join multiple CSV winter_year files into one single file
for statistical analysis.

For instance, when the variable jahre=(2012,2013,2014) the script will
look for yearly CSV files with the ending 2012-2013.csv, 2013-2014.cvs, and
2014-2015.cvs and merged them all into a single CSV file with the ending
like 2012-2014.csv containing the data from all those three winters.
Note: the merging happens in the order of years given by the variable jahre, so don't mess up!
=#
using CSV, DataFrames

PATH_DATA = "/projekt2/ac3data/B07-data/utqiagvik-nsa/csv_nsa/yearly";
# jahre indicates the wintertime period, e.g. 2012 comprises of 2012.11, 2012.12, 2013.1, 2013.2, 2013.3, 2013.4
jahre = (2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)

ff = let dfiles=[]
    tmp=readdir(PATH_DATA, join=true)
    ## filter(x->contains(x,".csv") ,tmp) #for yy ∈ jahre] |> x-> !isempty(x) && push!(dfiles,x)
    foreach(jahre) do yy
        filter(x->contains(x,"$yy-$(yy+1).csv"), tmp) |> x-> !isempty(x) && push!(dfiles, x[1])
    end
    dfiles
end

outff = joinpath(PATH_DATA, "all_nsa_microphys_db_$(jahre[1])-$(jahre[end]+1).csv")
println(outff)

dfs = (DataFrame∘CSV.File).(ff);

df = reduce(vcat, dfs);
CSV.write(outff, df) # joinpath(PATH_DATA, "whole_winters_NSA.csv"), df)
#println(ff)

# end of script
