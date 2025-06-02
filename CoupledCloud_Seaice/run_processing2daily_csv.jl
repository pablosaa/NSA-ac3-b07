#!/home/psgarfias/.local/bin/julia

#= Script to run over all cases for mosaic and create the daily CSV with
 CBH,CTH, CBT, CTT
=#

# general packages:
using Dates
using JLD2
using CSV, DataFrames, HTTP
using ImageFiltering
using Statistics
using Interpolations
using Printf

# my own packages:
using ARMtools
using CloudnetTools
using ATMOStools

# Defining Paths
const CAMPAIGN = "utqiagvik-nsa" # "arctic-mosaic" #
const DATA_PATH = "/projekt2/ac3data/B07-data" #joinpath(homedir(), "LIM/data/B07")
const LPR_PATH = joinpath(DATA_PATH, "LP")
const CLNET_PATH = joinpath(DATA_PATH, CAMPAIGN, "CloudNet/output/V1_68_1")   # "1.9.0"
const CLNET_PRODUCT = "CEIL10m" # "TROPOS/processed/categorize"
const RS_PATH = joinpath("/projekt2/remsens/data_new/site-campaign", CAMPAIGN)
const OUT_CSV = joinpath(DATA_PATH, CAMPAIGN, "csv_nsa") # "/tmp/csv_nsa" #

include("tmp_auxfiles.jl");

# Optional flags:
const ADDAOI = false;

# Reading AO index if flag set to TRUE:
if ADDAOI
    aoi_fn = let infile = "https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.cdas.z1000.19500101_current.csv";
        http_response = HTTP.get(infile);
        CSV.File(http_response.body, header=1);
    end
end


winter_jahr = 2011:2024;
#( (winter_jahr,11), (winter_jahr,12), (,1), (2019,2), (2019,3), (2019,4))
days = (1:31) #,21,22,23,24,25,26,27,28,29,30) #21  #18 #28 #6

!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

datum = [Date(yy, 11)+Month(m) for yy ∈ winter_jahr for m ∈ 0:5]
#datum = [Date(2025,3)] #, Date(2024,3), Date(2024,4)]

for heute in datum
    yy, mm = year(heute), month(heute)
    for dd in days
        try
            Date(yy,mm,dd)
        catch
            continue
        end
        # ******************************************************************
        # READING NEECESSARY DATA FILES:
        # Reading Cloudnet files:
        clnet_filen = ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="categorize.nc");

        class_filen = ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="classification.nc");

        isnothing(clnet_filen) && continue
        isnothing(class_filen) && continue

        clnet = CloudnetTools.readCLNFile(clnet_filen, modelreso=false) # , altfile=class_filen);

        LWC = let nfile=ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="lwc.nc");
            CloudnetTools.readLWCFile(nfile)
        end;

        IWC = let nfile=ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="iwc.nc");
            CloudnetTools.readIWCFile(nfile; inc_rain=false, apply_flag=[0,1,3,6])
        end;


        # *** Reading Radiosonde data from ARM NSA ***
        rs_filen = ARMtools.getFilePattern(RS_PATH, "INTERPOLATEDSONDE", yy, mm, dd)
        !isnothing(rs_filen) ? rs = ARMtools.getSondeData(rs_filen) : (@warn "No Radiosonde $(heute)"; continue)

        # Reading ARM microwave radiometer file:
        mwr = let nfile=ARMtools.getFilePattern(RS_PATH, "MWR/RET", yy, mm, dd)
            tmp = if !isnothing(nfile)
                # MWR/RET surface temperature in K, surface pressure kPa, orig_pwv cm (original PWV)
                ARMtools.getMWRData(nfile, onlyvars=["surface_temp"],addvars=["surface_pres", "orig_pwv"]) 
            else
                # if MWR/RET is not available, trying with RADFLUX data:
                #nfile=ARMtools.getFilePattern(RS_PATH, "MWR/LOS", yy, mm, dd)
                nfile = joinpath(RS_PATH, "RADFLUX")
                tmp = !isnothing(nfile) && ARMtools.read_radflux(nfile,Date(yy,mm,dd), onlyvars=["air_temperature"],addvars=["pressure"])
                Dict(:time=>tmp[:time], :SFT=>tmp[:T_air], :SURFACE_PRES=>tmp[:PRESSURE])
            end
            # MWR/RET is about 2.5 data points per minute, thus
            # interpolating to radiosonde time resolution:
            if !isnothing(nfile)
                Dict(:time=>rs[:time],
                     #:IWV=>CloudnetTools.Interpolate2Cloudnet(rs, tmp[:time], tmp[:ORIG_PWV]),
                     :SFT=>CloudnetTools.Interpolate2Cloudnet(rs, tmp[:time], tmp[:SFT]),
                     :SFPa=>CloudnetTools.Interpolate2Cloudnet(rs, tmp[:time], tmp[:SURFACE_PRES])
                    )
            else
                @warn("No MWR file found $(nfile)")
                nfile
            end
        end;

        #=
        Reading RADLUX LW up and down-welling radiation to estimate surface temperatures from ARM.
        In case RADFLUX data is not available for a given Date, then GND IRT data is used as proxy,
        and when neither RADFLUX for GNDIRT data are available then 2m T from MWR is used as alternative.
        =#
        tir = try
        #    tirpath = joinpath(RS_PATH, "RADFLUX")
        #    radfile = ARMtools.getFilePattern(RS_PATH, "RADFLUX", yy, mm, dd)
        #    tirfile = isnothing(radfile) && ARMtools.getFilePattern(RS_PATH, "GNDIRT", yy, mm, dd)
        #    if !isnothing(radfile)
        #        # RAD FLUX is 1 min resolution:
        #        radflux = ARMtools.read_radflux(tirpath, Date(yy,mm,dd), addvars=["pressure"])
        #        radflux[:IRT] = ATMOStools.T_surf.(radflux[:dw_lw], radflux[:up_lw])  # [K]
        #        radflux[:SFPa] = radflux[:PRESSURE]
        #        radflux
         tirfile = ARMtools.getFilePattern(RS_PATH, "GNDIRT", yy, mm, dd)
         if !isnothing(tirfile)
                # GND IRT is 1 min resolution (but sometimes a day has less than 1440 min).
                tirdat = let tmp = ARMtools.getGNDIRTdata(tirfile)
                    
                    if length(tmp[:time]) ≠ length(rs[:time])
                        tmp = Dict(:time => rs[:time],
                             Dict(k=>CloudnetTools.Interpolate2Cloudnet(rs, tmp[:time], V) for (k,V) ∈ tmp if k≠:time && isa(V, Vector))...)
                       
                    end
                    tmp[:SFPa] = mwr[:SFPa]
                    tmp
                end     
                tirdat
            else
                @warn "On $dd.$mm.$yy using MWR data as surface T proxy!"
                Dict(:time=>mwr[:time], :IRT=> mwr[:SFT], :SFPa=>mwr[:SFPa])
            end
        catch e
            @error("Didn't work with TIR") #i
            println(e)
            nothing
        end
       
        # Calculating Tskin by means of LW down-welling and up-welling radiation:
        radflux = let radfile=ARMtools.getFilePattern(RS_PATH, "RADFLUX", yy, mm, dd)
            if !isnothing(radfile)
                # RAD FLUX is 1 min resolution:
                radpath = joinpath(RS_PATH, "RADFLUX")
                radflux = ARMtools.read_radflux(radpath, Date(yy,mm,dd), addvars=["pressure"])
                radflux[:IRT] = ATMOStools.T_surf.(radflux[:dw_lw], radflux[:up_lw])  # [K]
                radflux[:SFPa] = radflux[:PRESSURE]
                radflux
            else
                radflux = Dict(:time=>rs[:time], :IRT=>fill(NaN32, length(rs[:time])), :SFPa=>fill(NaN32, length(rs[:time])))
            end
            radflux
        end
        ## tir = !isnothing(tir_filen) ? ARMtools.getGNDIRTdata(tir_filen) : Dict(:time=>mwr[:time], :IRT=> mwr[:SFT])

        # *** If TIR is availabe, then add it to radiosonde data as surface level: tir[:IRT] [K]
        # TODO: check what the differences are by using TIR, RADFLUX Tsurf, or T2m ***
        if all(isnan.(radflux[:IRT]))
            ARMtools.attach_Tₛ!(rs, tir[:IRT] .- 273.15, tir[:time]);
        else
            ARMtools.attach_Tₛ!(rs, radflux[:IRT] .- 273.15, radflux[:time]);
        end
        rs[:height][end] < 45 && (rs[:height] .*= 1f3)  # converting km to m 

        # ****************************************************************
        # Calculating atmospheric variables:
        # 1. Richardson number:
        N², Ri = ATMOStools.Ri_N(Float32.(rs[:height]), rs[:U], rs[:V],
                                 rs[:qv], rs[:θ], rs[:T]);
        
        # 2. Planetary boundary layer height:
        PBLH = ATMOStools.estimate_Ri_PBLH(Ri[7:end, :], rs[:height][7:end], ξ_ri=1.0, Hmax=9f3);

        # 3. Vertical gradient of water vapour transport:
        ∇WVT = ATMOStools.calculate_∇WVT(10rs[:Pa], rs[:WSPD], rs[:qv], 1f-3rs[:height]);

        # 4. Maximum ∇WVT height:
        ker2d = ImageFiltering.Kernel.gaussian((2,2));  # smoothing down peaks
        
        WVTH, maxWVT, snrWVT, idxmax,_ = ATMOStools.estimate_WVT_peak_altitude(rs,
                                                               WVT=imfilter(∇WVT, ker2d),
                                                               Hmin=10,
                                                               Hmax=PBLH,
                                                               get_index=true,
                                                              get_maxsnr=true);

        # 5. Cloud base and top heights:
        CBH, CTH, CLB = CloudnetTools.estimate_cloud_layers(clnet, alttime=rs[:time], liquid_base=true);

        # 6. Cloud base and top temperature [°C]
        CBT = CloudnetTools.cloud_temperature(rs, CBH) .+ 273.15; #, clnet=clnet);
        CTT = CloudnetTools.cloud_temperature(rs, CTH) .+ 273.15; #, clnet=clnet);
        CLT = CloudnetTools.cloud_temperature(rs, CLB) .+ 273.15;

        # 7. Finding nearest RS index to Cloudnet
        idxtclnt = [findmin(abs.(T .- rs[:time])) |> x->x[1]>Minute(1) ? -1 : x[2] for T in clnet[:time]] |> x->filter(>(0),x);

        #= finding index of cloud layer which is closer to the max WVTH:
        #  idx is CartesinIndex(data index, layer in)
        =#
        idx = let h=replace(CBH, NaN32=>9f9)
            δH = abs.(h .- WVTH) |> x->argmin(x, dims=2)
        end

        # 8. Calculating layer-base LWP and IWP from Cloudnet:
        # 8.1 Liquid water path integrated over the cloud layeer closest to max_WVT:
        LWP = let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
                X = ∫fdh(1f3LWC[:lwc], LWC[:height],
                                           h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
                𝐼[idxtclnt,i] = X;
            end
            𝐼
        end;
        # 8.2 Ice water path integrated over the cloud layer closest to max_WVT:
        IWP = let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
                X = ∫fdh(1f3IWC[:iwc], IWC[:height],
                                           h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
                𝐼[idxtclnt,i] = X;
            end
            𝐼
        end;
        
        # 8.3 Measured LWP and IWV from the radiometer direcly:
        mwrLWP = let X=fill(NaN32, size(CBH,1))
            X[idxtclnt] = replace(clnet[:LWP], missing=>NaN32)
            vec(X).*1f3  # converting categorize LWP [kg m⁻2] to [g m⁻2]
        end
        #mwrIWV = let X=fill(NaN32, size(CBH,1))
        #    X = replace(mwr[:IWV], missing=>NaN32)
        #    vec(X).*10.0  # Temporal solution to convert cm->kg cm⁻2, future will be used  value from categorize
        #end


        # 9. Calculating layer weighted average effective radius for droplets and ice particles:
        Reff = let fdrop=ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="der.nc")
            fice=ARMtools.getFilePattern(CLNET_PATH, CLNET_PRODUCT, yy, mm, dd, fileext="ier.nc")
            CloudnetTools.readReffFile(fdrop, altfile=fice)
        end;
        
        rₘₑ =  let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
            Nₗ₀ = ∫fdh(Reff[:Nₗ], Reff[:height], h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
            X = ∫fdh(Reff[:dₑᵣ].*Reff[:Nₗ], Reff[:height], h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
            𝐼[idxtclnt,i] = 1f6X./Nₗ₀;
            end
            𝐼
        end;
    
        iₘₑ =  let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
	        δH = (CTH[idxtclnt,i] .- CBH[idxtclnt,i])
	        X = ∫fdh(Reff[:Iₑᵣ], Reff[:height], h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
	        𝐼[idxtclnt,i] = X./δH;
            end
            𝐼
        end;
    
        # 10. Calculating Cloud Temperature Lapse-rate:
        Γₗᵣ = let 𝐼=fill(NaN32, size(CBH))
            LPR = -(rs[:T][2:end,:] .- rs[:T][1:end-1,:])./(rs[:height][2:end] .- rs[:height][1:end-1])
            for i ∈ 1:size(CBH,2)
	        δH = (CTH[:,i] .- CBH[:,i])
	        X = ∫fdh(LPR, rs[:height][1:end-1], h₀=CBH[:,i], hₜ=CTH[:,i])
                𝐼[:,i] = 1f3X./δH;
            end
            𝐼
        end;

        # 11. Calculating the decoupling altitudes:
        decop_hgt = ATMOStools.cloud_decoupling_height(rs[:height], CBH, rs[:θ], θ_thr=0.01);
        topdecop_hgt = ATMOStools.cloud_decoupling_height(rs[:height], CTH, rs[:θ], θ_thr=0.025, topmixlayer=true);

        # 12. creating FLAG for coupling and decoupling:
        
        ϑ_flag = [any(isnan.([D, W, T])) ? missing : (D .≤ W .≤ T) for (D, W, T) in zip(decop_hgt[idx], WVTH, topdecop_hgt[idx])];

        
        # 13. getting wind direction and speed at maximum ∇WVT height:
        wdir, wspd = let WD=fill(NaN32, size(CBH,1)), WS=fill(NaN32, size(CBH,1))
            for (ih, imx) ∈ enumerate(idxmax)
                WD[ih] = rs[:WDIR][imx, ih]
                WS[ih] = rs[:WSPD][imx, ih]
            end
            WD, WS
        end

        # 14. Get interpolated AO index:
        ADDAOI && (aoi = aoindex_from_timeseries(aoi_fn, rs[:time]));

        # 15. storing results:
        wdir_fn=joinpath(OUT_CSV, @sprintf("%04d/winddir_%04d%02d%02d_I.csv", yy, yy, mm, dd))
        # println(names(tir)); println(names(radflux)); println(names(mwr))
        CSV.write(wdir_fn, DataFrame(date=rs[:time],
                                     wvtdir=wdir,
                                     wvtspd=wspd,
                                     wvth=WVTH,
                                     mwvt=maxWVT,
                                     swvt=snrWVT,
                                     pblh=PBLH,
                                     lwp=vec(LWP[idx]),
                                     iwp=vec(IWP[idx]),
                                     Tskin=tir[:IRT], #rs[:T][1,:],
                                     Tsurf=radflux[:IRT],
                                     T2m =mwr[:SFT],
                                     Pa = tir[:SFPa],
                                     cldBT=vec(CBT[idx]),
                                     cldTT=vec(CTT[idx]),
                                     cldLT=vec(CLT[idx]),
                                     cbh=vec(CBH[idx]),
                                     cth=vec(CTH[idx]),
                                     clb=vec(CLB[idx]),
                                     ier=vec(iₘₑ[idx]),
                                     der=vec(rₘₑ[idx]),
                                     lpr=vec(Γₗᵣ[idx]),
                                     decoH=vec(decop_hgt[idx]),
                                     decoT=vec(topdecop_hgt[idx]),
                                     coupled=ϑ_flag,
                                     mwrlwp=mwrLWP)
                                     #mwriwv=mwrIWV) #aoi=aoi.aoi) #
                 )

    end  # over days
end  # over datum
# end of script
