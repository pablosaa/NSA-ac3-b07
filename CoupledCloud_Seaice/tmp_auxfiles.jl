# Auxiliary functions
"""
Function to vertically integrate the profiles.
USAGE:
```julia-repl
julia> I₀ₜ = ∫fdh(LWP, H)
julia> I₀ₜ = ∫fdh(LWP, H; h₀=H_bot, hₜ=H_top)
```
WHERE:
* LWP::Matrix 2D variable to integrate along the 2nd dimension (e.g. profile),
* H::Vector  with heights to integrate LWP along the 2nd dimension,
* h₀::Vector the bottom altitude of H to use for integration (optional),
* hₜ::Vector the top altitude of H to use for integration (optional)

RETURN:
* I₀ₜ::Vector with the variable LWP integrated relative to heights H.

(c) 2022 Pablo Saavedra Garfias
_Faculty of Physics and Geosciences_
_Leipzig University_
"""
function ∫fdh(x::AbstractArray, H::AbstractVector; h₀=Real[], hₜ=Real[])
    # getting dimensions:
    nheight, ntime = size(x)

    # checking dimensions:
    nheight != length(H) && @error "Heights muss have same length as matrix 1st dim."
    !isempty(h₀) && ntime != length(h₀) && @error "Optional variable h₀ must have length same as matrix 2nd dim."
    !isempty(hₜ) && ntime != length(hₜ) && @error "Optional variable hₜ must have length same as matrix 2nd dim."
	
    # defining the derivate of H
    δh = Vector{eltype(x)}(undef, length(H)) .= 0
    δh[2:end] = diff(H)
    δh[1] =  δh[2:end] |> minimum

    # integrating variable within limits:
    𝐼₀ₜ = fill(NaN32, ntime)
    for (i, X) ∈ enumerate(eachcol(x))
        # finding the limits for integration:
        i0 = isempty(h₀) ? 1 : !isnan(h₀[i]) ? argmin(abs.(H .- h₀[i])) : nothing
        it = isempty(hₜ) ? nheight : !isnan(hₜ[i]) ? argmin(abs.(H .- hₜ[i])) : nothing

        isnothing(i0) && continue
        isnothing(it) && continue

	lims = i0:it

        Xh = let Xh=X[lims] 
            tmp = ismissing.(Xh) .|| isnan.(Xh)
            all(tmp) && continue
            Xh[tmp] .= 0
            Xh
        end
        #Xh = (collect∘skipmissing)(Xh)
        
	𝐼₀ₜ[i] = Xh'*vec(δh[lims]) 

        
    end
    return 𝐼₀ₜ
end
#----/

# Function to obtain the daily Arctic Oscilation Index and interpolate to any time series:
"""
Function to obtain the daily Arctic Oscilation Index and interpolate to any time series

USAGE:
```julia-repl
julia> using CSV, HTTP, DataFrames
julia> fn = CSV.File(HTTP.get(url).body, header=1);
julia> aoi = aoindex_from_timeseries(fn, ts);
```
WHERE:
* fn::CSV.File the object from the AOI datafile specified by url::String,
* ts::Vector{DateTime} with the time series to interpolate,

RETURN:
* aoi::DataFrame with the files :date=>ts and :aoi the interpolated AO index.

(c) 2023 Pablo Saavedra Garfias
_Faculty of Physics and Geosciences_
_Leipzig University_
"""
function aoindex_from_timeseries(fn, ts::Vector{DateTime})

    aoi_date = DateTime.(fn.year, fn.month, fn.day) .|> Dates.value
    aoi_index= fn.ao_index_cdas
    # preparing variables to interpotale:
    itp = interpolate((aoi_date,), aoi_index, Gridded(Linear()));
    # now interpolating to MOSAiC time series:
    return DataFrame(:date=>ts, :aoi=>itp(Dates.value.(ts)))
end





#----/

