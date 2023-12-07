#!/home/psgarfias/.local/bin/julia
#
# Script to join multiple CSV winter_year files into one single file
# for statistical analysis.

using CSV, DataFrames

PATH_DATA = "/projekt2/ac3data/B07-data/utqiagvik-nsa/csv_nsa/yearly";
jahre = (2012, 2013) #, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)

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
