#!/home/psgarfias/.local/bin/julia

# Main script to run the conversion of ARM data files to netCDF input CloudNetpy files.
# USAGE:
# > ./run_convertor.jl > status.log 2>&1 & 
# > ./run_convertor.jl jahre="(2012,2013)"
# > ./run_convertor.jl jahre=2025 monaten="(1:4)" tage=1:31

using ARMtools
using Dates
using Printf

using CloudnetTools.ARM

# Declaration of parameter to be used to run
# Possible input parameters are:
# * lidar => "CEIL10m" or "HSRL"
# * radar => "KAZR/ARSCL" or "KAZR/CORGE"
# * radiometer => "MWR/RET" or "MWR/LOS"
# * model => "ECMWF" or "INTERPOLATESONDE"
NSAproduct = Dict(#:radar => ("KAZR/ARSCL", "KAZR/CORGE"),  #",
                  #:ceilometer => "CEIL10m", #
                  #:lidar => "CEIL10m", # "HSRL", #
                  :model => "",
                  #:mwr => ("MWR/RET", "MWR/LOS"), #
                  :radiosonde => "INTERPOLATEDSONDE",
                 );

input_params = Dict(
    :site => "utqiagvik", 
    :campaign => "nsa",
    :products => NSAproduct,
    :data_path => "/projekt7/remsens/data_new/site-campaign/utqiagvik-nsa/",
    :output_path => "/projekt2/ac3data/B07-data/utqiagvik-nsa/CloudNet/input",  #"/home/psgarfias/LIM/data/CloudNet/input"
   );

# Default values to use for range of dates to run the convertor:
jahre = (2023)
monaten = (4)
tage = (5) #:15)

# overwriting the input arguments to julia variables:
# jahre, monaten, tage. If not specified then default values are used.
!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

# Running over all specified products:
for yy ∈ jahre  
    for mm ∈ monaten 
        for dd ∈ tage

            try
                Date(yy, mm, dd);
            catch
                continue
            end
            foreach(input_params[:products]) do (thekey, product)
                listproduct = typeof(product) <: Tuple ? product : (product,)
                for prod in listproduct
                    isempty(prod) && continue

                    local in_params = input_params
                    in_params[:products][thekey] = prod
                    #a = !isempty(prod) && ARMtools.getFilePattern(input_params[:data_path], prod, yy, mm, dd)
                    a = !isempty(prod) && ARM.converter(yy, mm, dd, thekey, in_params);
                    ismissing(a) ? println("$thekey not possible for $(prod) from $(yy).$(mm).$(dd)") : break 
                        
                end
            end
            # --closure over dates loop:
        end
    end
end


# --- end of script.

