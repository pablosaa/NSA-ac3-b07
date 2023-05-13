#!/home/psgarfias/.local/bin/julia

#= Script to run over all cases for mosaic and create the daily CSV with
 CBH,CTH, CBT, CTT
=#

# general packages:
using Dates
using JLD2
using CSV, DataFrames
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
const CLNET_PATH = joinpath(DATA_PATH, CAMPAIGN, "CloudNet", "output")
const CLNET_PRODUCT = "CEIL10m" # "TROPOS/processed/categorize"
const RS_PATH = joinpath("/projekt2/remsens/data_new/site-campaign", CAMPAIGN)
const OUT_CSV = joinpath("data", "csv_nsa")

include("tmp_auxfiles.jl");

datum =( (2018,11), (2018,12), (2019,1), (2019,2), (2019,3), (2019,4))
days = (1:31) #21  #18 #28 #6

!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

for (yy,mm) in datum
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
            CloudnetTools.readIWCFile(nfile)
        end;


        # Reading ARM microwave radiometer file:
        mwr = let nfile=ARMtools.getFilePattern(RS_PATH, "MWR/RET", yy, mm, dd)
            !isnothing(nfile) && ARMtools.getMWRData(nfile, onlyvars=["time","surface_temp"], addvars=["sonde_times"])
        end;

        # Reading Infrared surface temperatures from ARM
        tir_filen = ARMtools.getFilePattern(RS_PATH, "GNDIRT", yy, mm, dd)
        tir = !isnothing(tir_filen) ? ARMtools.getGNDIRTdata(tir_filen) : Dict(:time=>mwr[:time], :IRT=> mwr[:SFT])


        # Reading Radiosonde data from ARM NSA
        rs_filen = ARMtools.getFilePattern(RS_PATH, "INTERPOLATEDSONDE", yy, mm, dd)
        rs = ARMtools.getSondeData(rs_filen);

        typeof(tir)<:Dict && ARMtools.attach_Tₛ!(rs, tir[:IRT].-273.15, tir[:time]);
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
        LWP = let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
                X = ∫fdh(1f3LWC[:lwc], LWC[:height],
                                           h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
                𝐼[idxtclnt,i] = X;
            end
            𝐼
        end;

        IWP = let 𝐼=fill(NaN32, size(CBH))
            for i ∈ 1:size(CBH,2)
                X = ∫fdh(1f3IWC[:iwc], IWC[:height],
                                           h₀=CBH[idxtclnt,i], hₜ=CTH[idxtclnt,i])
                𝐼[idxtclnt,i] = X;
            end
            𝐼
        end;

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

        # 14. storing results:
        wdir_fn=joinpath(OUT_CSV, @sprintf("%04d/winddir_%04d%02d%02d_I.csv", yy, yy, mm, dd))
        CSV.write(wdir_fn, DataFrame(date=rs[:time],
                                     wvtdir=wdir,
                                     wvtspd=wspd,
                                     wvth=WVTH,
                                     mwvt=maxWVT,
                                     swvt=snrWVT,
                                     pblh=PBLH,
                                     lwp=vec(LWP[idx]),
                                     iwp=vec(IWP[idx]),
                                     Tskin=rs[:T][1,:],
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
                                     coupled=ϑ_flag)
                  )

    end  # over days
end  # over datum
# end of script
