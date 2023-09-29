#!/home/psgarfias/.local/bin/julia

using NCDatasets
using Navigation
using Dates
using Distributions
using JLD2
using GMT
using Printf
using ARMtools
using CSV, DataFrames
using StatsBase

include("aux_math_functions.jl");

SEAICE = include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"));

const PROD_PATH = "/projekt2/ac3data/B07-data/SeaIce";  # joinpath(homedir(), "LIM/data/B07/SeaIce");
const RVNAV_PATH = joinpath(homedir(), "LIM/data/B07/arctic-mosaic");
const LATLON_FILE = joinpath(PROD_PATH, "amsr2", "LongitudeLatitudeGrid-n3125-ChukchiBeaufort.h5");
const DATCSV_PATH = joinpath(homedir(), "LIM/scripts/NSA-ac3-b07/CoupledCloud_Seaice/data");
const DATOUT_PATH = joinpath(homedir(), "LIM/scripts/NSA-ac3-b07/SeaIce/data"); ## old: CoupledCloud_Seaice/data");
const R_lim = 50e3;   # radius around RV polarstern
const MAKEPLOTS = false

# Define coordinates for the North Slope Alaska site:
nsa_lat = 71.323e0;
nsa_lon = -156.609e0;
# Define the angular sector to consider (avoiding Land):
θₗ₀ = 235e0;
θₗ₁ = 110e0;

PRODUCTS = (:SIC,) # (:DIV, :LF, :SIC)

dist_wdir = Dict(data=>Dict() for data ∈ PRODUCTS)

#yy = 2020
#mm = 4
#dd = 15

datum = ((11,2021), (12,2021), (1,2022), (2,2022), (3,2022), (4,2022)) #(11,2021), 
days = (1:31)

!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

for (mm, yy) ∈ datum #, (5,2020)]
#mm = 11; yy=2019;

## # reading 6 hour RV track coordinates to plot:
## RVtrack = if MAKEPLOTS
##         tmp = @sprintf("../RV_polarstern/data/RVpolarstern_track_%04d.jld2", yy)
##         load(tmp, "RV")
##     else
##         nothing
##     end


for dd = days
        heute = try
            Date(yy, mm, dd);
        catch
           continue
        end

for data ∈ PRODUCTS
#let data = :SIC
    # relevant constants specific for data product:
    cbfaktor, cblims, cbbins, cbcolor, cbtrunc, pxsize = if data==:DIV
        1f5, (-1.13, 1.13, 0.15), 0.3, :cork, (-1, 1), 0.07 
    elseif data==:LF
        1, (0, 0.4, 0.025), 0.1, :davos, (0, 1), 0.07
    elseif data==:SIC
        1, (0, 100, 5), 5, :vik, (-1, 0),  0.27
    else
        @error "given data key $(data) not supported"
    end

    ## data == :SIC && (DFsic=CSV.read("data/modis_amsr2_osisaf_statsII.csv", DataFrame))
    
    iplt = 0
    
            # reading RV Polarstern data:
    # @sprintf("../RV_polarstern/data/RVpolarstern_track_%04d.jld2", yy);
    ## RV_PATH = ARMtools.getFilePattern(RVNAV_PATH, "RVpolarstern", yy, mm, dd);
    ## isnothing(RV_PATH) && @error "no RV NAV file found for $(dd).$(mm).$(yy)"
    ## RV = load(RV_PATH, "RV")

    # Reading wind direction of the day
    csv_winddir = let fn = @sprintf("%04d/winddir_%04d%02d%02d_I.csv", yy, yy, mm, dd)
        full_fn = joinpath(DATCSV_PATH, "csv_nsa", fn) 
    end

    !isfile(csv_winddir) && (@warn "$(csv_winddir) not found!"; continue)
    winddir = CSV.read(csv_winddir, DataFrame)
        
    ##    # finding out RV positions on working date:
    ##    iday = findall(==(heute), Date.(RV[:time]))
    ##    isempty(iday) && (@warn "No RV nav for $(heute)";) #continue)
    ##    Nnav = length(iday)

        # loading daily data
    sar = Dict()
    dist_wdir[data] = let tmp=range(cblims[1], stop=cblims[2], step=cblims[3])
        ndat = length(winddir.wvtdir)
        Dict(:bins=>cbfaktor*tmp,
             :freq=>fill(NaN32, length(tmp)-1, ndat),
             :μ => fill(NaN32, ndat),
             :σ => fill(NaN32, ndat),
             :qq=> fill(NaN32, 5, ndat),
             :Afreq=>fill(NaN32, length(tmp)-1, ndat),
             :Aμ => fill(NaN32, ndat),
             :Aσ => fill(NaN32, ndat),
             :Aqq=> fill(NaN32, 5, ndat),
             :Rfreq=>fill(NaN32, length(tmp)-1, ndat),
             :Rμ => fill(NaN32, ndat),
             :Rσ => fill(NaN32, ndat),
             :Rqq=> fill(NaN32, 5, ndat),
                )
    end
    
    dat_coor = if data == :SIC
        lr_filen = ARMtools.getFilePattern(PROD_PATH, "amsr2", yy, mm, dd, fileext=".h5")
        isnothing(lr_filen) && (@warn "No data for $(heute)"; continue)
        SEAICE.read_LatLon_Bremen_product(LATLON_FILE)
    else
            
        lr_filen = SEAICE.load.FilePattern(PROD_PATH, "sentinel-1A_luisa/drift_corrected", yy, mm, dd);
        isnothing(lr_filen) && (@warn "No data for $(heute), $(data) product!"; continue)

        sar_x, sar_y = SEAICE.load.LatLon_AWI_product(lr_filen)
        XX, YY = couple2grid(sar_x, sar_y);
            
        ##sar = SEAICE.load.Data_Divergence_LeadFraction(lr_filen);
        ##XX, YY = couple2grid(sar[:x], sar[:y]);
        lat, lon = xy2latlon(XX, YY, ϕₖ=70.0, λ₀=-45.0);
        Point.(lat, lon) #dat_coor = 
    end
        
    lon, lat = SEAICE.Get_LonLat_From_Point(dat_coor)
        
    println("Working on $(heute) with file $(lr_filen)")

    for i_rv ∈ 1:nrow(winddir) #iday
        
        ## rv_lat, rv_lon = RV[:lat][i_rv], RV[:lon][i_rv]
        Pstern = Point(nsa_lat, nsa_lon)
            
        # finding limits of Box:
        LonLim, LatLim = SEAICE.estimate_box(Pstern, R_lim, δR=10e3)

        idx_box = SEAICE.extract_LonLat_Box(LonLim, LatLim, dat_coor);

        if data==:SIC
            ##fφ = let datum=RV[:time][i_rv]
            ##    df = getSICfixfactor(datum, DFsic)
            ##    try
            ##        df[1, :ratio]
            ##    catch e
            ##        @error "Problem with SIC fix factor in $(datum) and $e"
            ##    end
            ##end
            sar[:SIC] = SEAICE.read_SIC_Bremen_product(lr_filen,
                                                       idx_box)
                                                            # SICPROD="mersic")
        else
            sar = SEAICE.load.Data_Divergence_LeadFraction(lr_filen, idx_box)
        end
        # converting to polar for data within the box
        θ, ρ = SEAICE.LonLat_To_CenteredPolar(Pstern, dat_coor[idx_box]); #!!!! out idx_box
        
        # Tim = findall(RV[:time][i_rv] .≤ winddir[!, :date] .< RV[:time][i_rv+1])
        # Tim = [argmin(abs.(T.-RV[:time][i_rv])) for T ∈ winddir[!,:date]]
        #for iti ∈ Tim
        let iti=i_rv  ## argmin(abs.(winddir[!,:date].-RV[:time][i_rv]))

            # getting wind direction where to extract:
            θᵢ = winddir[iti, :wvtdir]-3; θₑ = winddir[iti, :wvtdir]+3;
            
            idx_wd = findall((θᵢ .≤ θ .≤ θₑ) .& (ρ .≤ 50));  #!!! [idx_box]

            # calculating statistics of box:
            LF = cbfaktor*sar[data]

            μLF, σLF, qqLF = if data==:LF
                stats_𝑁ₗᵤ(filter(≥(0), LF[idx_wd]))
            elseif data==:DIV
                stats_𝑁ₗᵤ(filter(!isnan, LF[idx_wd]), L=nothing, U=nothing)
            elseif data==:SIC
                
                stats_𝑁ₗᵤ(filter(!isnan, LF[idx_wd]), L=0, U=100)
                # make corrections for SIC here!
            else
                @error "data $(data) not supported!"
            end

            # Creating Histogrmas from the wind sector:
            dist_wdir[data][:freq][:, iti] = let Hdat=fit(Histogram,
                                                          LF[idx_wd],
                                                          dist_wdir[data][:bins],
                                                          closed=:right)
                Hdat.weights
            end
            dist_wdir[data][:μ][iti] = μLF
            dist_wdir[data][:σ][iti] = σLF
            dist_wdir[data][:qq][:, iti] = qqLF
            
            # Adding SIC statistics for the whole sector:
            idx_A=findall((θ.≥ θₗ₀ .|| θ.≤ θₗ₁) .& (ρ .≤ 50));
            Aμsic, Aσsic, Aqqsic = stats_𝑁ₗᵤ(filter(!isnan, sar[data][idx_A]), L=0, U=100)
           
            dist_wdir[data][:Aμ][iti] = Aμsic
            dist_wdir[data][:Aσ][iti] = Aσsic
            dist_wdir[data][:Aqq][:, iti] = Aqqsic
            dist_wdir[data][:Afreq][:, iti] = let Hdat=fit(Histogram,
                                                          sar[data][idx_A],
                                                          dist_wdir[data][:bins],
                                                          closed=:right)
                Hdat.weights
            end
            # Adding SIC sector from random direction:
            θᵢ, θₑ = let wvtdir = 359e0rand(1)[1]
                wvtdir-3, wvtdir+3
            end
            idx_rnd = findall((θᵢ .≤ θ .≤ θₑ) .& (ρ .≤ 50));
            Rμsic, Rσsic, Rqqsic = stats_𝑁ₗᵤ(filter(!isnan, sar[data][idx_rnd]), L=0, U=100)
           
            dist_wdir[data][:Rμ][iti] = Rμsic
            dist_wdir[data][:Rσ][iti] = Rσsic
            dist_wdir[data][:Rqq][:, iti] = Rqqsic
            dist_wdir[data][:Rfreq][:, iti] = let Hdat=fit(Histogram,
                                                          sar[data][idx_rnd],
                                                          dist_wdir[data][:bins],
                                                          closed=:right)
                Hdat.weights
            end

            #
            !MAKEPLOTS && continue
            mod(iti,360) != 0 && continue
            # finding the RVtrack point for the given date:
            ## itrack = let DateTrack = Date.(RVtrack[:time])
            ##     findfirst(==(heute), DateTrack)
            ## end

            # creating lines with angle limits to plot
            ##Plon, Plat = SEAICE.create_pair_lines(Point(rv_lat, rv_lon), 50e3, [θᵢ, θₑ])

            # creating circle to plot
            P_circ = SEAICE.Create_Semi_Circle(Pstern, 1f0, 360f0, R_lim=R_lim)
            lon_circ, lat_circ = SEAICE.Get_LonLat_From_Point(P_circ)

    # creating wind lines to plot
            lon_wind, lat_wind = SEAICE.create_pair_lines(Pstern, R_lim, winddir[iti, :wvtdir].+[-3, 3])

            # creating title for plot:
            # NOTE: year≤2018 the version number for AMSR2 is v5. from year 2019 is v5.4
            amsr2ver = yy<2019 ? "v5" : "v5.4"
            plot_name = if data==:SIC
                "/home/psgarfias/quicklooks/NSA/$(data)/$(yy)/"*replace(basename(lr_filen),
                                                                        amsr2ver=>@sprintf("nsa%04dT%s", iplt, Dates.format(winddir[iti,:date], "HH:MM")),
                                                                        ".h5"=>"$(data).png");
                                              
            else
                                              
                "/home/psgarfias/quicklooks/$(data)/"*replace(basename(lr_filen),
                                              "pairs"=>@sprintf("mosaic%04d", iplt),
                                              "deformation_moved.nc"=>"$(data).png");
            end
            
            utc_title_str = @sprintf("%s %3.1f \261 %3.1f on %s", data, μLF, σLF,
                                     Dates.format(winddir[iti, :date], "dd-uuu HH:MM"))
            
            wind_title_str = "@:9:@%14%"*utc_title_str*"@%%@::"

            cbar = if data==:SIC
                GMT.makecpt(cmap=cbcolor, truncate=cbtrunc, range=cblims)
            else
                GMT.makecpt(cmap=cbcolor, truncate=cbtrunc, range=cblims, inverse=:c)
            end
            
            GMT.scatter(lon[idx_box], lat[idx_box], zcolor=LF,
	                proj=(name=:Stereographic, center=[nsa_lon, 90]),  # [mean(LonLim)
                        figsize=8,
                        xaxis=(axes=:Sn, annot=2, ticks="60m", grid=false),
                        yaxis=(axes=:We, annot=0.5, ticks="10m", grid=false),
                        region=(LonLim..., LatLim...), cmap=cbar, markersize=pxsize,
                        marker=:square, title=wind_title_str,
	                par=(FONT_ANNOT_PRIMARY=7, FONT_LABEL=10))
            
            GMT.coast!(proj = (name=:stereographic, center=[nsa_lon, 90], paralles=30),
                       figsize=8, #frame=(axes=:a1g, annot="60m"),
                       res=:high, area=950, land=:lightgreen,
                       shore=:black, show=0)

            GMT.basemap!(region=(LonLim..., LatLim...),
                         proj=(name=:Stereographic, center=[nsa_lon, 90]),
                         figsize=8,
	                 #xaxis=(axes=:Sn, annot=5, ticks="1deg", grid=false),
                         rose = (outside=true, anchor=:TR, width=1., frame=false, offset=(0.5, 0.08), label=true, lc=:lightgreen) #(-0.5, -1.0), label=true)
	                 )
            
            ## GMT.plot!(RVtrack[:lon][1:itrack], RVtrack[:lat][1:itrack], linestyle=:line, lc=:gray, show=0)

            ## GMT.plot!(RV[:lon][1:i_rv], RV[:lat][1:i_rv], linestyle=:line, lc=:blue, show=0)
            GMT.plot!(lon_wind, lat_wind, lc=:orange, linestyle=:line, show=0)
            GMT.colorbar!(pos=(anchor=:CR, length=(4,0.15), offset=(0.2 ,0)), cmap=cbar, frame=(annot="$(cbbins)"), par=(FONT_ANNOT_PRIMARY=8, FONT_LABEL=10))#, xlabel="Lead Fraction", font="9p")
            GMT.plot!(lon_circ, lat_circ, linestyle=:dash, lc=:lightred)
            ##GMT.plot!([Plon...], [Plat...], linestyle=:dash, lc=:lightblue)

            GMT.plot!(nsa_lon, nsa_lat, marker=:star, markersize=0.2, mc=:red, show=0, fmt=:png, savefig=plot_name)

            iplt += 1
        end
    end  # over iday
    
##    end # over daysinmonth
##end # over months and year

end # over variables

output_fn = @sprintf("SIC/%04d/seaice_winddir_%04d%02d%02d.jld2", yy, yy, mm ,dd);

##println(output_tn)
save_object(joinpath(DATOUT_PATH, output_fn), dist_wdir)

end # over days
end # over datum

# end of script
