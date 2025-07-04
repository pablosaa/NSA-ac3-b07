#=
 script to create the trend vs SIC statistics for several variables
The data files are pre-calculated and stored at buffer_data under the following structure:
statstype_var_SICxxx-yyy.csv
Where statstype is either ravstat or wavstat,
var can be T2m, lwp, iwp, der, etc.
SICxxx-yyy indicated the range of SIC used e.g. 085-100 is from 85% to 100% and 000-100 is all SIC values.
=#

using Printf
using CSV, DataFrames
using Glob
using Plots
using StatsPlots
using LaTeXStrings

# defining variables to work out for plotting:
Trdvar = :cldTT #:Tsurf #:lwp #:T2m #
Trdstat = :q50 #:μ #

str_lab = Dict(
    :lwp=>(var="LWP", unit="[g m⁻² dec⁻¹]"),
    :iwp=>(var="IWP", unit="[g m⁻² dec⁻¹]"),
    :clb=>(var="CBH", unit="[m dec⁻¹]"),
    :Tsurf=>(var=L"T_\mathrm{skin}", unit="[K dec⁻¹]"),
    :cldTT=>(var=L"CTT", unit="[K dec⁻¹]"),
)

buffer_path = joinpath("buffer_data", String(Trdvar));

add_sic(sic, df, method) = insertcols!(df,1, :sic => sic, :method=>method)
getme_df(df, vstat, vpre, vcop) = filter(d->d.type==vstat && d.pax==vpre && d.coupled==vcop, df)

# Loading the running ave stats:
ravfiles = Dict(i=>f for (i,f) in enumerate(glob("moavstat_$(Trdvar)*.csv", buffer_path)));

ravstat, siclist = let df=DataFrame()

    siclist = Dict()
    coltypes=[Symbol, Float64, Float64, String, Float64, Float64, Float64, Float64, Float64, Symbol, Symbol]


    foreach(ravfiles) do (i,F)
        sicloc = findfirst.(["SIC",".csv"], F) |> x->range(last(x[1])+1, first(x[2])-1) 
        siclim = parse.(Int8, split(F[sicloc], "-"))
        siclist[i]= all(siclim.==100) ? "All" : @sprintf("(%d, %d]", siclim...)
        Fdf = CSV.read(F,types=coltypes, DataFrame)
        Fdf.CI95 = @. (eval∘Meta.parse)(Fdf.CI95)
        add_sic(i, Fdf, :runave)
        df = vcat(df, Fdf)
    end
    df, siclist
end

# ==
# Loading the FFT wavstats:
wavfiles = Dict(i=>f for (i,f) in enumerate(glob("ensostat_$(Trdvar)*.csv", buffer_path)));

wavstat = let df=DataFrame()
    # pax, coupled, type, hat, sig, CI95, pval, snr, r², rmse, chi², cost, θ
    coltypes=[Symbol, Symbol, Symbol, Float64, Float64, String, Float64, Float64, Float64, Float64, Float64, Float64, Float64]

    foreach(wavfiles) do (i,F)
        #sicloc = findfirst.(["SIC",".csv"], F) |> x->range(last(x[1])+1, first(x[2])-1) 
        #siclim = parse.(Int8, split(F[sicloc], "-"))
        
        Fdf = CSV.read(F,types=coltypes, DataFrame)
        Fdf.CI95 = @. (eval∘Meta.parse)(Fdf.CI95)
        add_sic(i, Fdf, :cosfit)
        df = vcat(df, Fdf)
    end
    df
end


# getting the values for SIC from Dictionary keys:
#siclist = [@sprintf("%d<SIC≤%d", parse.(Int8, split(k,"-"))...) for k in keys(ravstat)];
#replace!(siclist,  "(0,100]"=>"All SIC")

mark_Trd = Dict(
    :co=>(run=(:o, 4, [:black :white], stroke(1.5)), wav=(:o, 5, stroke(1.5, :royalblue))),
    :de=>(run=(:^, 4, [:black :white], stroke(1.5)), wav=(:^, 5, stroke(1.5, :tomato))) )

δs = 0.12; # space to separate the symbols for run-ave and enso-fit
enso_col = cgrad(:Spectral_10, 10, categorical=true, rev=true)  # colorscale for period of ENSO
Pₖ = [10.5, 6.3, 5.3, 4.5, 2.9, 2.4];

# creating canvas as background for all 4 subplots:
plt_bg = plot([0],fill(NaN,10), zcolor=(0:0.2:1),
              color=enso_col, colorbar=true, clim=(0,1), #extrema(Pₖ),
              axis=false, grid=false,
              guidefontsize=9,
              ylim = (0,1),
              #ylabel
              #guidefontsize=4,
              xlim = (0,1);
              label= false,
              ann = (1.06, 1.045, text(L"\chi^2", 15)),
              #legend_font_pointsize = 13,
              #legend_background_color = false,
              #legend_foreground_color = false,
              #legend = (1.1,0.95),              
              colorbar_ticks= [0,0.2,0.4,0.6,0.8,1], #["0.0","0.2","0.4","0.6","0.8","1.0"],
              #colorbar_titlefontsize = 10,
              #colorbar_titlefontrotation=-90,
              left_margins=9Plots.mm,
              rigth_margins=-3Plots.mm);

annotate!(-0.11, 0.5, text(L"\frac{\Delta}{\Delta t}"*str_lab[Trdvar].var*"  "*str_lab[Trdvar].unit, 13, :vcenter, rotation=90))

# adding 2x2 frames to be filled with subplots:
plot!([0],[NaN], inset_subplots=[(1, bbox(0.0,0.0,0.42,0.47))], subplot=2, label=false)
plot!([0],[NaN], inset_subplots=[(1, bbox(0.45,0.0,0.42,0.47))], subplot=3, label=false)
plot!([0],[NaN], inset_subplots=[(1, bbox(0.0,0.5,0.42,0.47))], subplot=4, label=false)
plot!([0],[NaN], inset_subplots=[(1, bbox(0.45,0.5,0.42,0.47))], subplot=5, label=false)

# vector to fill with subplots
figtag='a'
plt_Trd = []
sp_idx = 1
for Trdcop in [:co, :de]

    hat_lims = let rng = [filter(d->d.type==Trdstat && d.coupled==Trdcop, ravstat).CI95...,
                          filter(d->d.type==Trdstat && d.coupled==Trdcop, wavstat).CI95...]
        a1 = first.(rng) |> extrema
        a2 = last.(rng) |> extrema
        [a1..., a2...] |> extrema |> e->e.*(1.05, 1.05)
    end

    for Trdpre in [:H, :L]


        #
        global sp_idx += 1 
        # plotting the statistically significant flag:
        df = getme_df(ravstat, Trdstat, Trdpre, Trdcop)
        
        # scatter plot for running average trends:
        ci95mk = ifelse(Trdcop==:co, :o, :^)

        # scatter plot for CI points:
	@df df scatter!(plt_bg[sp_idx], :sic .-δs, :hat,
                        marker=(ci95mk, :grey, 9, stroke(0)),
                        ma=ifelse.(:pval .<0.05, .6, 0),
                        label=ifelse(Trdpre==:H && Trdcop==:co, "CI₉₅",""),
                        title=ifelse(Trdcop==:co, "$(Trdpre)-Pressure",""),
                        framestyle=:box,
                        ann=(0.7, hat_lims[2]*0.9, text("($(figtag))", 10))
                        )
        global figtag += 1
        # adding a horizontal line to mark zero trend:
        hline!(plt_bg[sp_idx], [0], l=:dash, lc=:gray, label=false)

        # scatter plot for running average trends:
        @df df scatter!(plt_bg[sp_idx], :sic .-δs, [:hat :hat],
                        #yerror = :sig, #max.(abs.(:hat)/5, :sig ),
                        yerror=[(:hat .- first.(:CI95), last.(:CI95) .- :hat) 0],
                        #zcolor=fill(NaN, size(:sic)),
                        m=mark_Trd[Trdcop].run,
                        markerstrokecolor=:auto,
                        ylim=hat_lims,
                        label=["" "Eq. 4a"],
                        tickdir=:out, yminorticks=true)

           
        # scatter plot for enso FFT trends:
        df = getme_df(wavstat, Trdstat, Trdpre, Trdcop)
        
        @df df scatter!(plt_bg[sp_idx], :sic .+ δs, :hat,
                        xlabel=ifelse(Trdcop==:de, "SIC", ""),
                        marker=(ci95mk, :grey, 9, stroke(0)),
                        ma=ifelse.(:pval .<0.05, .6, 0),
                        label=false, #ifelse(Trdpre==:H, "95CI",""),
                        ann=(:sic .+ 0.27, :hat.+1.2.*sign.(:sig).*(:sig), text.(round.(inv.(:θ), digits=1), :black, :left, 8, rotation=30) ) )

        vytick = Trdpre == :L && yticks(plt_bg[sp_idx])[1]
        @df  df scatter!(plt_bg[sp_idx], :sic .+ δs , :hat,
                         #yerror= :sig, #max.(abs.(:hat)/5, :sig),
                         yerror=(:hat .- first.(:CI95), last.(:CI95) .- :hat),
                         zcolor=round.(:chi², digits=2),
                         m=mark_Trd[Trdcop].wav,
                         ylim=hat_lims,
                         yticks=ifelse(Trdpre==:H, :auto, (vytick, "")),
                         color=enso_col,
                         clim=(0,1),
                         colorbar=false,
                         label="Eq. 4b",
                         legend=ifelse(Trdpre==:H && Trdcop==:co,
                                       :topright,
                                       ifelse(Trdpre==:L && Trdcop==:de,
                                              :topright,
                                              false)
                                       ),
                         legend_column=1,
                         xtick=(keys(siclist), ifelse(Trdcop==:de, values(siclist), "")),
                         xlims=(0.5,3.7),
                         xrot=0,
                         tickfontsize = 10,
                         top_margin=ifelse(Trdcop==:de, -3, 0)Plots.mm)
            

        push!(plt_Trd, plt_bg[sp_idx])
        
        
    end
end

plot(plt_bg, dpi=600)

#savefig("~/Downloads/quicklooks/nsa/Trends_fun-SIC_tsk.png")

# end of script.
