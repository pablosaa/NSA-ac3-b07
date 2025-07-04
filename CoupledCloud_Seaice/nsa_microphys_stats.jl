### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 0cdaa10d-e72a-4b41-8e41-860317ced277
begin
	using Plots
	using StatsPlots
	using Statistics
	using DataFrames
	using CSV
	using Dates
	using StatsBase
	using LaTeXStrings
	using Distributions
        using LsqFit
        using Printf
	using PrettyTables
	using PlutoUI
        using GLM
        using ATMOStools.CLIMA
end

# ╔═╡ 3d45676a-8105-4e69-a520-605ee9cca2d1
using HypothesisTests

# ╔═╡ 7b9dfaec-9bd3-4608-a932-6c7bf2e6f019
using LinearAlgebra

# ╔═╡ 123e1148-759c-11ed-29d7-e7fc1998111a
html"""<style>
main {
	max-width: 1030px;
	font-size: 13pt;
}
"""

# ╔═╡ f2a22774-dd59-485d-934d-15aa4568d379
md"""
# Statistical analysis for NSA coupled/decoupled micro-physical properties
"""

# ╔═╡ 48c90554-bc23-4e44-a050-4bce892614a2
# Defining data path
begin
    const BASE_PATH = joinpath(pwd(), "buffer_data")
end;

# ╔═╡ d734cb56-6da3-4d2f-8bca-a804bba0ba35
md"""
##### Save Figures? $(@bind SAVEFIG confirm(CheckBox(default=false)))
"""

# ╔═╡ 17af9a85-0159-4c0a-b2f7-2dfb5e15bbbb
md"""
###### Select wintertime years to load: $(@bind wintertime confirm(Select(["2011-2025", "2012-2024", "2012-2022", "2016-2017", "2017-2018", "2012-2013","2018-2019", "2019-2020", "2020-2021"])))
"""

# ╔═╡ 5e79a910-f5b8-46ff-a5cb-d9204ff54cb6
md"""
## Reading the Dataset for the period from $(wintertime)
#### Defining new variables:
"""

# ╔═╡ d56e6ef8-24cc-4afc-a155-5d411c3b9422
DBraw = CSV.read(joinpath(BASE_PATH, "csv_nsa/yearly", "all_nsa_microphys_db_$(wintertime).csv"), header=1, skipto=2, DataFrame);

# ╔═╡ e274212a-ef60-4258-930e-32d6d27f1e36
begin
    # AOi flagging:
    
    DBraw[!, :ΔZ] = let tmp=fill(NaN32, length(DBraw[!, :Pa]))
        paₘ = filter(>(0), DBraw.Pa) |> mean
		taₘ = filter(>(0), DBraw.Tskin) |> mean
        R =  287.05  # [J kg⁻¹ K⁻¹]  Specific gas constant
        g₀ = 9.83  #[m s⁻²]
        ΔZ(T, P, Pₘ) = R/g₀*T*log(Pₘ/P)  # [m]
		
		T₀ = copy(DBraw.Tskin)
		T₀[isnan.(T₀)] .= taₘ    #@. (DBraw.T2m + (DBraw.Tskin+273.15))/2
        tmp = @. ΔZ(T₀, DBraw.Pa, 101.32)
        tmp
    end  #
	# L-H system flagging:
	DBraw[!, :paₓ] = let tmp=fill(0, length(DBraw.ΔZ)) #sign.(DBraw.ΔZ) #
		#paₘ = filter(!isnan, DBraw.Pa) |> mean
        tmp[DBraw.ΔZ .≤ 0] .= -1 
        tmp[DBraw.ΔZ .> 0] .= +1
        tmp
    end
	
   # Cloud thinkness [m]:
   DBraw[!, :δₕ] = (DBraw.cth .- DBraw.clb)
 
end;

# ╔═╡ eb3faf69-cab0-45d7-9484-28ab4963ed9c
md"""
### Creating a new DataSet by filtering the main DataSet according to the following constrains:
* Cloud Top Height < 3.1km
* rejecting missing data of couple/decouple classification,
* Ensure that the detected Cloud Base Height is below Cloud Top Height,
* Extracting data with wind direction from azimuth angles covering land,
* Setting minimum detectable liquid/ice water path to zero,
"""

# ╔═╡ 020695a4-61ee-49f5-8e39-4256ca836183
begin
    # Range of azimuth to filter out corresponding to Land:
    θ₀ = 235e0;
    θ₁ = 110e0;

    DB = let tmp = filter(𝐷-> (𝐷.cth ≤ 3.0f3) && !ismissing(𝐷.coupled), DBraw) #(0< 𝐷.δₕ ≤ 3.5f3)
        # converting type Missing from :coupled to Bool:
		disallowmissing!(tmp, :coupled, error=false)
		
		# Extracting the data with wind direction between the azimuths covering land:
        filter!(𝐷-> 𝐷.wvtdir≥θ₀ || 𝐷.wvtdir≤θ₁, tmp)
		
		# Converting Ice effective radius from m to μm:
        tmp.ier *= 1e6

		# When Liquid and Ice waterpath is below 5 g m⁻² then is considered 0 (retrieval minimum detection):
		#ii = findall((tmp.lwp .< 5) .&& (tmp.iwp .<5))
        #ii = findall(tmp.lwp .> 1000 .|| tmp.lwp .< 5)
        #tmp.lwp[ii] .= NaN32 #0.0
		#ii = findall(tmp.iwp .> 1000 .|| tmp.iwp .< 5)≤ 1f3
        #tmp.iwp[ii] .= NaN32 #0.0
			
		filter!(d-> d.lwp > 5 && d.iwp > 5, tmp) # || 5<d.iwp<1f3, tmp)
		filter!(d-> d.Tskin ≤ 273.9 && d.lwp < 8f2, tmp)
		filter!(d-> !isnan(d.μSIC), tmp)
		tmp
    end;
	# Computing Optical thickness τc:
	DB[:, :τc] = @. 3/2*DB.lwp/DB.der
	DB[:, :τi] = @. 3/2*DB.iwp/DB.ier/0.917
	
	# Categorizing SIC into fixed width bins:
	# 0.0 - 75.0, 0.4 75.0 - 94.0, 0.6 94.0 - 98.0, 0.8: 98.0 - 99.0, 1.0: 99.0 - 100.0
    SIC_bin = (0, 30, 60, 80, 90, 95, 100) #(0.0, 43.941566, 91.19343, 98.2683, 99.8023, 100.0) #(0, 75, 94, 98, 99, 100) #(0:10:100) 

    DB[!, :sicₓ] = let tmp=fill(0f0, length(DB[!, :μSIC]))
        foreach(zip(SIC_bin[1:end-1], SIC_bin[2:end])) do (xb, xt)
            ii = findall((DB.μSIC .> xb) .&& (DB.μSIC .≤ xt))
            tmp[ii] .= mean([xb, xt])
        end
        tmp
    end
	
    # Categorizing data with WINTER fraction, with Nov=0.15, Dec=0.30, Jan=0.45, Feb=0.60, Mar=0.75, Apr=0.90
    jahren = Year.(extrema(DB.date)) |> J->J[1].value:J[2].value  # +1 to include Nov, Dec of last year
    DB[!, :winter] = let tmp = fill(0f0, length(DB.date))
		#yb_mo = [Date(mo>4 ? jj : jj+1,mo) for jj in jahren[1:end-1] for mo in (11,12,1,2,3,4)]
		#yb_fr = Dict(:11=>0, :12=>0.1, :1=>-0.8, :2=>-0.7, :3=>-0.6, :4=>-0.5)
		#foreach(yb_mo) do Djamo
		#	ii = findall((DB.date .≥ firstdayofmonth(Djamo) .&& (DB.date .< lastdayofmonth(Djamo)+Day(1)) ))
		#	tmp[ii] .= year(Djamo) + yb_fr[month(Djamo)]
		#end
        foreach(zip(jahren[1:end-1], jahren[2:end])) do (yb, yt)
            	ii = findall((DB.date .≥ Date(yb,11,1)) .&& (DB.date .≤ Date(yt,5,1) ) )
            tmp[ii] .= yb
        end
        tmp
    end
end

# ╔═╡ 3595bf88-4dbf-4c63-90e5-2feae1b78760
"""
Function to add ```:winter``` to DataFrame.
Wintertime is defined as the period from Nov 1st to Apr 30th of following year, and assigned as the Nov year.
The function requires that the DataFrame has the column ```:date```

```julia> add_winter!(df)```
"""
function add_winter!(df::DataFrame)
	tmp = fill(0f0, length(df.date))
	jahren = Year.(extrema(df.date)) |> J->J[1].value:J[2].value
	foreach(zip(jahren[1:end-1], jahren[2:end])) do (yb, yt)
		ii = findall((df.date .≥ Date(yb,11,1)) .&& (df.date .≤ Date(yt,5,1) ) )
        tmp[ii] .= yb
    end
	idxin = contains.(names(df), "date") |> findfirst
	insertcols!(df, idxin, :winter => tmp) #df[:, :winter] = tmp
	return nothing
end

# ╔═╡ ae59d379-b47d-41fc-8958-4e9e1941b528
# Percentage of RawData with WVT direction comming from the Land sector:
filter(D-> D.wvtdir>θ₁ && D.wvtdir<θ₀ , DBraw) |> D->length(D.date)/length(DBraw.date) # && D.cth≤3f9  && D.mwrlwp>0

# ╔═╡ d1b9d2b9-6b08-4d81-8042-26bddf8f5e54
filter(D-> (D.cth - D.cbh)≤(3f3) , DB).date |> length #&& && D.coupled==(true) && D.lf₀==(false) 

# ╔═╡ 9fe20fcc-8de8-4f4c-aff2-214797c40ba0
Ndat = Dict(:co=>reduce(+, DB.coupled)/length(DB.date), :de=>reduce(+, .!DB.coupled)/length(DB.date), :H=>reduce(+, DB.paₓ .<0)/length(DB.date), :L=>reduce(+, DB.paₓ .>0)/length(DB.date)); pretty_table(HTML, Ndat)

# ╔═╡ 89809ea8-94a4-4158-8f23-b4a2322532b2
md"""
##### Which SIC estimation should be used? : $(@bind sicvar Select([:μSIC, :RμSIC, :AμSIC]))
"""

# ╔═╡ db52e7c0-5a8c-454a-bf1e-df9d377eff25
md"""
### Calculating anciliary varaibles: cloud thickness δₕ, lapse-rate Γ, ice water fraction χᵢ, flag for SIC<95 (LF>0.02)
"""

# ╔═╡ f5a31896-eeb0-4a49-9178-3995fde42a55
# Calculating other variables:
begin
	# SIC threshold to split database:
	sic₀ = 95
        # Cloud lapse-rate [K km⁻¹]:
	DB[!, :Γ] = @. -(DB.cldTT - DB.cldLT)/DB.δₕ/1f-3
	# Cloud phase fraction:
	DB[!, :χᵢ] = DB.iwp./(DB.lwp .+ DB.iwp)
	DB[!, :lf₀] = (DB[!, sicvar] .< sic₀)

end;

# ╔═╡ 838b64cb-5719-4a82-94e0-9bd3f5e18598
md"""
#### General definitions of parameters for plotting:
* farben : two colors for decoupled and coupled data, respectively,
"""

# ╔═╡ e530627f-4fa8-4006-b119-55f9efe72d21
begin
	farben = [:tomato1 :royalblue2]
	hinweistext = ["de. ∀ SIC≥$(sic₀)" "co. ∀ SIC≥$(sic₀)" "de. ∀ SIC<$(sic₀)" "co. ∀ SIC<$(sic₀)"]
	
end

# ╔═╡ 018fc44f-0e6b-4eb2-8b17-4ee61015fb7c
filter(c->!isnan(c.decoH) && c.coupled==(false),DB) |> x-> length(filter(<(5),x.decoH))/length(x.decoH)

# ╔═╡ 16ee6021-e6ec-4b6c-aa91-1da8d8c07411
md"""
### Fitting ``\chi_{ice}`` to empirical function by Coopmann et al. (2018)
`` \chi_{ice}(T; \beta) = \frac{1+tanh(-\beta_1 (T-\beta_2))}{2}``
"""

# ╔═╡ 7cb41328-9a82-4fe4-a65c-b96e424b4e88
md"""
###### Select only SIC<90 ? $(@bind lf002 confirm(CheckBox(true)))
"""

# ╔═╡ 154f4210-a584-4040-9c03-3c9c51c467bd
#begin
#	ΔT = 1
#	Tin = (-50:ΔT:0)
#	NT = length(Tin)
#	CPF = DataFrame(:T=>zeros(2NT), :σT=>zeros(2NT), :μχᵢ=>zeros(2NT), :σχᵢ=>zeros(2NT), :coupled=>zeros(2NT), :N=>zeros(2NT), :Nw=>zeros(2NT), :Ni=>zeros(2NT))
#	for (i, T) in enumerate(Tin)
#		T0 = 273.15 + T-1.2ΔT #T + 273.15
#		T1 = 273.15 + T+1.2ΔT #T + 273.15 + 5
#		
#		## DD=filter(R->T0≤R.cldTT<T1 && R.coupled==(true) && (R.lf₀ || ~lf002) && TG[1]≤R.date≤TG[2], DB)
#		
#		DD0 = filter(R->T0≤R.cldTT<T1 && (R.lf₀ || ~lf002), DB)  ##  && TG[1]≤R.date≤TG[2]
#		
#		# For coupled:
#		DD = filter(R->R.coupled, DD0) 
#		CPF[2i-1, :μχᵢ], CPF[2i-1, :σχᵢ] = Nlu(DD.χᵢ, L=0, U=1)
#		CPF[2i-1, :coupled] = true
#		CPF[2i-1, :T] = T #mean(DD.cldTT)-273.15 #T
#		CPF[2i-1, :σT] = std(DD.cldTT) #T
#		CPF[2i-1, :N] = isempty(DD.χᵢ) || all(isnan.(DD.χᵢ)) ? NaN32 : length(DD.χᵢ)
#		CPF[2i-1, :Nw] = filter(!isnan, DD.lwp) |> length
#		CPF[2i-1, :Ni] = filter(!isnan, DD.iwp) |> length
#
#		# For decoupled:
#		DD=filter(R->!R.coupled, DD0)
#		CPF[2i, :μχᵢ], CPF[2i, :σχᵢ] = Nlu(DD.χᵢ, L=0, U=1)
#		CPF[2i, :coupled] = false
#		CPF[2i, :T] = T #mean(DD.cldTT)-273.15 #T
#		CPF[2i, :σT] = std(DD.cldTT) #T
#		CPF[2i, :N] = isempty(DD.χᵢ) || all(isnan.(DD.χᵢ)) ? NaN32 : length(DD.χᵢ)
#		CPF[2i, :Nw] = filter(!isnan, DD.lwp) |> length |> y->2y
#		CPF[2i, :Ni] = filter(!isnan, DD.iwp) |> length |> y->2y
#		
#	end		
#	
#end

# ╔═╡ 3af8f75f-0c94-4b4a-af8b-14d90585bd25
χice(T, β) = @. 0.5(1+tanh(-β[1]*(T-β[2]))) #1-1/(1+exp(.3(-T-25))) #

# ╔═╡ 263c6114-ee34-4da6-8466-5d7bd50dfd03
#χfit = let tmp = @. !isnan(CPF.μχᵢ) && CPF.coupled==(true)
#	p0 = [0.18, -15]
#	curve_fit(χice, CPF.T[tmp], CPF.μχᵢ[tmp], p0)
#end

# ╔═╡ 3b1bd823-83bb-4630-be82-26b08c90f45a
#begin
#	lfsplt = plot()
#	# Coopman et al. fit:
#	##@df CPF plot!(:T, χice(:T, [0.18, -18]), lw=2, lc=:bone, label="Coopman et al.(2018)")
#	# Fit to coupled data:
#	@df CPF plot!(:T, x->χice(x, χfit.param), lw=2, lc=farben[2], ls=:dash, label="Best fit-model") #[0.09, -25]
#	# Iso-line at -15 C max supersaturation:
#	vline!([-15], ls=:dash, lc=:skyblue, lw=3, label="T = -15 °C")
#	# binned data for decoupled and coupled:
#	@df CPF plot!(:T, :μχᵢ, ribbon=:σχᵢ, lw=0.5, group=:coupled, m=[:^ :o], ms=[5 5], markerstrokewidth=.5, color=farben, fillalpha=.3, xflip=false, xlim=(-51, 1), framestyle=:box, label=["decoupled ± σ" "coupled ± σ"], legend=:bottomleft, xlab="Cloud Top Temperature / °C", ylab="Ice fraction "*L"\chi_{ice}"*" @ Z $zflag", guidefontsize=16,legendfontsize=12, tickfontsize=13, minorticks=true, tickdir=:out)
#	# top histogram: Note, when LF->all ylim=(0, 21.2f3) and legend=:toprigth, when LF>0.02 ylim=(0, 2.2f3) and no legend
#	hist_ytick = (lf002 ? (0:1f4:4.5f4) : (0:5f4:25.1f4))
#	hcpf = @df CPF plot(:T, [:Nw :Ni], group=:coupled, l=[:bar :steppre :bar :steppre], ls=:solid, lw=3, lc=[false farben[1] false farben[2]], bar_width=[0.9 0.4 0.4], fillcolor=[farben farben[2]], fillalpha=[0.9 0.7 0.7], xlim=(-51, 1), yscale=:identity, yticks=hist_ytick, yaxis=(formatter=y->@sprintf("%1.1f",1f-3y)), yguidefontsize=12, ytickfontsize=12, ylim=extrema(:Ni).*(1,1.05), xflip=false, xtickfontcolor=:gray, bottom_margins=-8Plots.mm, ylab="# x 10⁴", tickdir=:out, minorticks=true, ann=(-4, hist_ytick[end-1], text(lf002 ? "(b)" : "(a)", 20)), legend=:topleft, background_color_legend=nothing, foreground_color_legend=nothing, label=["2x  liquid (de)" "2x  Ice (de)" "liquid (co)" "Ice (co)"], top_margins=3Plots.mm) #[false :royalblue2 false :orange] [:royalblue2 :orange :orange], ylim=(0, 1.5f4)
#	chifig = plot(hcpf, lfsplt, layout=@layout([a{0.2h}; b]), size=(700, 550), left_margins=3Plots.mm)
#
#	SAVEFIG ? savefig(chifig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_Xice.png") : chifig
#end

# ╔═╡ bbc95bb1-a296-4409-93be-ecf3c83c1153
md"""
#### Select presure regime: $(@bind zflag confirm(Select([-1=>"H", +1=>"L"], default=+1)))
"""

# ╔═╡ e1f28ccb-e252-4b44-95a3-4ffce7d45455
begin
	cbhfig = @df filter(c->!isnan(c.clb), DB) density(:clb, group=(:lf₀, :coupled), trim=true, color=farben, label=hinweistext, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], lw=[2 2 4 4], xlabel="Cloud liquid base height / m", xlim=(1f2, 6f3), xscale=:log10, xticks=([10,10^2,10^3,10^4], ["10¹","10²","10³","10⁴"]), ylabel="PDF @ Z $zflag", ylim=(0, .005), frame_style=:box, tickdir=:out, minorticks=true, tickfontsize=13, guidefontsize=15, legend=:topright, legendfontsize=15, size=(500,500), right_margin=2Plots.mm , fillalpha=0.4)
	
	SAVEFIG ? savefig(cbhfig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_CBH.png") : cbhfig
end

# ╔═╡ 5c9a28aa-5145-4a70-bc0e-5c299d3ae829
begin
	cdpfig = @df filter(c->c.δₕ>(0), DB) density(:δₕ, group=(:lf₀, :coupled), trim=true, color=farben, label=hinweistext, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], lw=[2 2 4 4], xlabel="Cloud depth / m", xlim=(30, 1.1f4), xscale=:log10, xticks=([10^2,10^3,10^4], ["10²","10³","10⁴"]), ylabel="PDF @ Z $zflag", ylim=(0, .006), frame_style=:box, tickdir=:out, minorticks=true, tickfontsize=13, guidefontsize=15, legend=:topright, legendfontsize=15, size=(500,500), right_margin=2Plots.mm , fillalpha=0.4)
	SAVEFIG ? savefig(cdpfig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_CDP.png") : cdpfig
end

# ╔═╡ 6992ff1a-2c1c-44b0-854d-b260a26909c8
begin
	ctt_plt = @df DB density((:cldTT).-273.15, group=(:lf₀, :coupled), label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], xlim=(-65, 3), xlabel="Cloud top temperature / °C", ylabel="PDF @ Z $zflag", trim=false, lw=[2 2 4 4], frame_style=:box, tickdir=:out,  minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=12, size=(500,500), legend=:topleft)
	SAVEFIG ? savefig(ctt_plt, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_CTT.png") : ctt_plt
end

# ╔═╡ d2accbb3-87d7-484c-82f1-ab6b488b653c
begin
	lprfig = @df filter(c->!isnan(c.lpr), DB) density(:lpr, group=(:lf₀, :coupled), label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], lw=[2 2 4 4], xlim=(-30, 15), xlabel=L"\Gamma_\textrm{cloud}=-\frac{dT}{dh}\,~/~\,\textrm{°C~km^{-1}}", ylabel="PDF @ Z $zflag", ylim=(-0.002, 0.18), trim=true, frame_style=:box, tickdir=:out, minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=15, size=(500,500), legend=:topleft); vline!([6.0],ls=:dash, lw=2, la=0.8, lc=:black, label=L"~\Gamma_\textrm{m}=6~\textrm{K~km^{-1}}")
	
	SAVEFIG ? savefig(lprfig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_LPR.png") : lprfig
end

# ╔═╡ f99cdbfe-2d31-483f-b1a2-664ddffd95f3
begin
	sktfig = @df filter(c->!isnan(c.Tskin), DB) density(:Tskin .- 273.15, group=(:lf₀, :coupled), label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], lw=[2 2 4 4], xlim=(-50, 1.5), xlabel="Surface skin Temperature / °C", ylabel="PDF @ Z $zflag", ylim=(-0.002, 0.15), trim=false, frame_style=:box, tickdir=:out, minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=15, size=(500,500), legend=:topleft) # vline!
	
	SAVEFIG ? savefig(sktfig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_SKT.png") : sktfig
end

# ╔═╡ 1f6657ba-9910-48fc-84a3-399a28d84ef5
begin
	idxco = findall(DB[!, :coupled])
	idxde = findall(.!DB[!, :coupled])
	myedges = range(1, stop=1000.0, length=50) |> collect
	Nlf02co = filter(c->(c.lf₀==(true) && c.coupled==(true)), DB).lwp |> length
	Nlf02de = filter(c->c.lf₀==(true) && c.coupled==(false), DB).lwp |> length
	#marginal
	lwphistco = StatsBase.fit(Histogram, DB[idxco, :lwp], myedges)
	lwphistde = StatsBase.fit(Histogram, DB[idxde, :lwp], myedges)
	iwphistco = StatsBase.fit(Histogram, DB[idxco, :iwp], myedges)
	iwphistde = StatsBase.fit(Histogram, DB[idxde, :iwp], myedges)
		#DB[idxde, :lwp], colorbar_scale=:log10, colorbar=true, colorbar_position=:left)
	#vline!([-5.5])http://localhost:1234/?secret=QNUDmOcxyerror=:σχᵢ, 
	scatter(lwphistco.weights, lwphistde.weights, ms=7, zcolor=(myedges), clim=(1, 1000), xscale=:log10, yscale=:log10, xlim=(.51,1f6),ylim=(.51,1f6), label="Liquid", xlabel="# Coupled cases", ylabel="# Decoupled cases", color=NNcol)
	scatter!(iwphistco.weights, iwphistde.weights, zcolor=(myedges), clim=(1, 800), xscale=:log10, yscale=:log10, xlim=(.51,1f6),ylim=(.51,1f6), ms=6, m=:star6, label="Ice",color=NNcol, colorbar_titlefontsize=13, colorbar_title="\nWater path [g m⁻²]", colorbar_scale=:log10)
	plot!([1, 1f6],[1, 1f6], ls=:dash, lc=:black, lw=2, label="1:1")
	plot!([10 1f2; 1f6 1f6],[1 1; 1f5 1f4], lw=2, ls=:dashdot, lc=[:tomato3 :tomato], label=["10:1" "100:1"], frame_style=:box, tickdir=:out, minorticks=true, tickfontsize=13, guidefontsize=14, legendfontsize=13, legend=:topleft, margin=8Plots.mm, size=(600,500), dpi=600)
	
	#plot(histco); plot!(histde)
end

# ╔═╡ 3473d761-f57f-4d3f-beba-a8a516041da5
begin
	reff_plt = @df filter(c->1<c.der<150, DB) density(:der, group=(:lf₀, :coupled), label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], xlim=(0, 60), xlabel=L"\mathrm{Liquid}~~\overline{r}_{eff}~ /~~\mathrm{\mu~m}", ylabel="PDF @ Z $zflag", trim=false, lw=[2 2 4 4], frame_style=:box, tickdir=:out,  minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=12, size=(500,500), legend=:topright)
	#savefig(reff_plt, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(aoflag)_DEF.png")
end

# ╔═╡ 07d4b3f0-d37e-4b5a-ae8f-46033e57b112
begin
	ieff_plt = @df filter(c->1<c.ier<150, DB) density(:ier, group=(:lf₀, :coupled), label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], xlim=(0, 70), xlabel=L"\mathrm{Ice}~~\overline{r}_{eff}~ /~~\mathrm{\mu~m}", ylabel="PDF @ Z $zflag", trim=false, lw=[2 2 4 4], frame_style=:box, tickdir=:out,  minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=12, size=(500,500), legend=:topleft)
	#savefig(ieff_plt, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI$(zflag)_IEF.png")
end

# ╔═╡ 76c14688-0be1-43a5-932e-352ebac308b4
begin
	aoifig = @df filter(d->!isnan(d.ΔZ), DB) density(:ΔZ, group=(:lf₀, :coupled), bandwidth=3.5, label=hinweistext, lc=farben, linealpha=[.7 .7 1 1], l=[:dash :dash :solid :solid], xlim=(-350, 350), xlabel="ΔZ [m]", ylabel="PDF", trim=true, lw=[2 2 4 4], frame_style=:box, tickdir=:out,  minorticks=true, tickfontsize=13, guidefontsize=15, legendfontsize=12, size=(650,300), legend=:outerright, left_margins=3Plots.mm, bottom_margins=3Plots.mm)
	vline!([15.9], ls=:dash, la=0.5, lc=:black, label="")
	#savefig(aoifig, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_AOI.png")
end

# ╔═╡ 6b11b88f-66e4-4630-a9a0-bc479cda43e1
# Collecting the data into fixed bin sets for LF and SIC:
begin
	## δlf = 0.02
	## lf_bin = range(δlf, step=δlf, stop=.5) |> collect
	## Nbins = length(lf_bin)
	δsic = 10 #mean(diff(sic_bin))
    sic_bin = DB.sicₓ |> unique |> sort  #range(δsic, step=δsic, stop=100) |> collect
    Nbins = length(sic_bin)
	lwp_bin = fill(NaN32, Nbins,2)
	iwp_bin = fill(NaN32, Nbins,2)
	
	#sic_edg = sic_bin .- diff(sic_bin)[1]; push!(sic_edg, 101)
end;

# ╔═╡ ac217227-a672-4f9c-80ff-277c57e473f9
"""
Function to compute average of vector following a given statisitcs methods:
USAGE:
```julia-repl
julia> μ, σ = Nlu(X; stats=:median)
julia> μ, σ = Nlu(X; L=0, U=100, stats=:truncate)
```
WHERE:
* ```X::AbstractVector``` containig the data to average,
* ```stats::Symbol``` (Optional) which method to use: ```:median```, ```:aritmetic```, ```:geometric```, ```:q25```, ```:q75```, ```:truncate``` (default),
* ```L::Real``` and ```U::Real``` are the lower and upper limits for the tats=:truncate methodt (default ```L=0```, ```U=nothing```)
OUTPUT:
* ```μ::Real``` The mean, median, 25th or 75th quantile determined by the option ```stats```,
* ```σ::Real``` The standard deviation, or MAD. for quantiles NaN is returned, for ```:geometric (μ-σₗ, σₕ-μ)``` is returned.
"""
function Nlu(x::AbstractVector; L=0, U=nothing, stats=:truncate)

	x = filter(!isnan, x)
	if isempty(x) && 
            return NaN32, ifelse(stats==:geometric, (0,0), 0)
	end
	
	qq = quantile(x, [.25, .5, .75])
	
	if stats==:aritmetic
		μ = mean(x)
		σ = std(x)
		return μ, σ	
	elseif stats==:median
		
		mad = median(abs.(qq[2] .- x))
		return qq[2], mad  #(qq[1], qq[3]) #
	elseif stats==:q25
		return qq[1], NaN
	elseif stats==:q75
		return qq[3], NaN
	elseif stats==:truncate
		filter!(≥(L), x) 
		𝑁 = fit_mle(Normal, x)
		𝑆 = if isnothing(U)
        	truncated(𝑁, lower=L)
    	elseif L<U
        	truncated(𝑁, lower=L, upper=U)
    	else
        	@error "L must be lower than U, but given $L ≥ $U"
    	end
		return mean(𝑆), std(𝑆)
		
	elseif stats==:geometric
		filter!(>(0), x) 
		X = @. log(x)
		μ = (exp∘mean)(X)
		Xμ = @. (X - μ)^2  # log(x/μ)
		σ = (exp∘std)(X)
		σₗ = μ/σ
		σₕ = μ*σ
		return μ, (μ-σₗ, σₕ-μ)
	else
		@error "optional argument stats=$(stats) not recognized!"
	end
end

# ╔═╡ 49caf4c1-1ca2-483e-bfa7-36b98f3589be
combine(groupby(DB, :winter)) do df
	var=:lwp
	μ, σ = Nlu(df[:, var], stats=:aritmetic)
	u, s = Nlu(df[:, var], stats=:truncate)
	g, h = Nlu(df[:, var], stats=:geometric)
	med, mad = Nlu(df[:, var], stats=:median) #filter(>(0), df.lwp) |> x->(median(x), median(abs.(x .- median(x))))
	(ave = μ, sig = σ, tru = u, tsd = s, med=med, mad=mad, ged=g, gsd=h)
end |> df-> @df df scatter(:winter .+[0 .1 .2 .3], [:ave :tru :med :ged], yerror=[:sig :tsd :mad :gsd])#, xflip=true)
	

# ╔═╡ 8e883748-effd-4898-b6ea-24de538fc5e8
filter(V-> V.coupled==(true), DB) |> F->length(findall(<(10), F.decoH))/length(F.decoH)

# ╔═╡ c5b02891-e561-4667-a07b-a29f52bb66c4
fσₑᵣ(x) = sqrt(mean(x.^2)/length(x))

# ╔═╡ 8846a5e6-ab5a-4156-bc59-b19a2033577f
let edges=[SIC_bin[1].*(1,1)]
	[push!(edges, (SIC_bin[i-1], SIC_bin[i])) for i in 1:Nbins if i>1]
		edges
	end

# ╔═╡ 264627a8-4975-4b75-b4de-136d9861aa0e
begin
    # CLOUD Variables to consider
    #VARIN = [:lwp, :iwp, :der, :ier, :δₕ, :Γ, :χᵢ, :cldTT]

    # auxiliary function to count number of data without NaN:
    Zahlen(N) = filter(!isnan, N) |> length

	# range of acceptable values for variables:
	var_lim = Dict(
			:lwp=>(lim=(0, 1000), stats=:median),
			:iwp=>(lim=(0, 4000), stats=:median),
			:der=>(lim=(0, 150), stats=:median),
			:ier=>(lim=(0, 150), stats=:median),
			:δₕ=>(lim=(0,3.5f3), stats=:median),
			:Γ=>(lim=(-40,40), stats=:median),
			#:χᵢ=>(lim=(-0.1,1.1), stats=:truncate),
			:cldTT=>(lim=(0, 300), stats=:median),
			:clb=>(lim=(0, 3.5f3), stats=:median),
			#:τc=>(lim=(0, 500), stats=:median),
			#:τi=>(lim=(0, 500), stats=:median)
	)

	VARIN = keys(var_lim)
	# auxiliary function to select variable within acceptable limits:
	within(V; lims=(0, Inf)) = ifelse(lims[1]<0, lims[1] ≤ V ≤ lims[2], lims[1] ≤ V < lims[2])

	dat = let tmp = Dict()
		# considering flag of H/L pressure?
		δZ = (-1,1)
		xNbins = length(δZ)*Nbins

		# Initializing variables:
        tmp[:de] = Dict(Symbol(x,y)=>fill(var_lim[y].stats==:geometric && x==:σ ? (0f0, 0f0) : NaN32, xNbins) for x in [:μ, :σ] for y in VARIN) |> DataFrame
		tmp[:co] = Dict(Symbol(x,y)=>fill(var_lim[y].stats==:geometric && x==:σ ? (0f0, 0f0) : NaN32, xNbins) for x in [:μ, :σ] for y in VARIN) |> DataFrame
		
		tmp[:de].NNsic = fill(0, xNbins)
		tmp[:co].NNsic = fill(0, xNbins)
		tmp[:de].coupled = fill(false, xNbins)
		tmp[:co].coupled = fill(true, xNbins)
		
		tmp[:de][!, :sic_bin] = [sic_bin..., sic_bin...]
		tmp[:co][!, :sic_bin] = [sic_bin..., sic_bin...]
		
		tmp[:de][!, :sic_var] = fill(NaN32, xNbins)
		tmp[:co][!, :sic_var] = fill(NaN32, xNbins)

		# pressure flag:
		tmp[:de][!, :Δz] = [fill(δZ[1], Nbins)..., fill(δZ[2], Nbins)...]
		tmp[:co][!, :Δz] = [fill(δZ[1], Nbins)..., fill(δZ[2], Nbins)...]
		
		for z in δZ
			
		for y in VARIN
			
			for i in 2:Nbins

				idx = i + (z==1 ? Nbins : 0)
				isym = Symbol(:μ, y)
				ssym = Symbol(:σ, y)
				# ---
				sic_edgs = (SIC_bin[i-1], SIC_bin[i])
				
				# filling with variables (decoupled):
				
				tmp[:de][idx, isym], tmp[:de][idx, ssym], tmp[:de][idx, :sic_var], tmp[:de][idx, :NNsic] = let
					F = filter(V->within.(V[sicvar], lims=sic_edgs) .&& within.(V[y], lims=var_lim[y].lim) .&& (V.coupled==(false) .&& V.paₓ==(z)), DB)
					Nlu(eval(:($F.$y)); stats=var_lim[y].stats)..., fσₑᵣ(F.σSIC), Zahlen(F.μSIC)
				end
				# filling with variables (coupled):
				tmp[:co][idx, isym], tmp[:co][idx, ssym], tmp[:co][idx, :sic_var], tmp[:co][idx, :NNsic] = let
					F = filter(V->within.(V[sicvar], lims=sic_edgs) .&& within.(V[y], lims=var_lim[y].lim) .&& (V.coupled==(true)) .&& V.paₓ==(z), DB)
					Nlu(eval(:($F.$y)); stats=var_lim[y].stats)..., fσₑᵣ(F.σSIC), Zahlen(F.μSIC)
				end
				#
			end
		end
		end  # over z var
	vcat(tmp[:de], tmp[:co]) |> df->groupby(df, :Δz)
	end
end;

# ╔═╡ bb176db2-fd18-4a59-9021-d2396a030b6d
sic_bin, SIC_bin

# ╔═╡ fa9e471b-3f58-459a-b99c-48cb1087e795
unique(DB.sicₓ)

# ╔═╡ 75609282-c67c-4794-a68b-94d32c004783
md"""
### Fitting data as function of SIC
`` f_s(SIC) = β_1 - β_2\exp[sic^2]``

where ``f_s \in \{\mathrm{lwp},~ \mathrm{iwp},~ r_{eff} \}`` and ``sic \in \{0, \dots, 1\}``
"""

# ╔═╡ 46c4e7cd-7d8a-49b7-abdd-1347692dcfd8
begin
	fᵢ(X, β) = @. β[1]*X^β[2]
    fₛ(X, β; X0=100) = @. β[1] + β[2]*exp((X/X0)^2) #+ β[2]*exp((X/X0)^2) #exp((X/X0)^2) #- β[2]*exp((X/X0)^2 - 1) # 
    fᵧ(X, β; X0=100) = @. β[1] - β[2]*exp((X/X0)^2 -1) #(X/X0)^β[3] #
end

# ╔═╡ acad36c9-c463-47c2-9277-efea563847df
filter(d->!isnan(d.sicₓ), DB) |> df->groupby(df, [:sicₓ, :paₓ, :coupled]) |> gf->combine(gf) do d
	
	Dict([:μclb, :σclb] .=> Nlu(d.clb)) |> DataFrame
end

# ╔═╡ d04fad13-2fe2-4fa2-8528-24de8dafbbad
rfit, R², Chi², Δβ = let ice = :μ,  ## :μLF, 
    tmp = Dict()
	r2 = Dict()
	chi2 = Dict()
	con = Dict()
    foreach(Dict(:co=>true, :de=>false)) do (i, co_status)
		
        r2[i] = Dict()
		chi2[i] = Dict()
        tmp[i] = Dict()
        con[i] = Dict()

		foreach((-1,1)) do varz
			Δz = varz==-1 ? Symbol("H") : Symbol("L")
			r2[i][Δz] = Dict()
			chi2[i][Δz] = Dict()
			tmp[i][Δz] = Dict()
			con[i][Δz] = Dict()

			Y=filter(D->D.coupled==(co_status) && D.Δz==(varz), dat[(Δz=varz,)]) # && D.sic_bin<(100)
			foreach(VARIN) do S
				var = Symbol(ice, S)

				err = replace(String(var), "μ"=>"σ") |> Symbol
				NNice = replace(lowercase(String(ice)), "μ"=>"NN") |> Symbol
				
				## jj = (Y[!, :lf_bin].>0.02) #: (Y[!, :sic_bin].>85) #NNice] > 10) #(Y[!, :lf_bin] > 0.05)
				ii = @. !isnan(Y[!, var]) #&& (Y[!, NNice]>10) && jj
				Yin = Y[ii, var]
				Xin = Y.sic_bin[ii]			
				#𝑆ₑᵣᵣ = @. 1/Y[ii, err]^2 #@. Y[ii, err]/Yin/sqrt(Y[ii, NNice] .-1)
				𝑛 = length(ii)
				ωᵢ, σ² = let
					Y_err = eltype(Y[ii, err]) <: Tuple ? (last.(Y[ii, err]) .- first.(Y[ii,err]))  : Y[ii, err]
					Y_err = @. ifelse(Y_err > 0, Y_err, Y[ii, var]/2)
					σ² = @. Y_err^2
					err_sig = @. inv(σ²)
					err_tot = mean(err_sig)
					AnalyticWeights(err_sig/err_tot), σ²
				end
				#abs(1 - (Y[ii, err]/Yin)^2))  #𝑛*𝑆ₑᵣᵣ/sum(𝑆ₑᵣᵣ) #sqrt(𝑆ₑᵣᵣ/𝑛) #
            	tmp[i][Δz][var], r2[i][Δz][var], chi2[i][Δz][var] = let 𝑓 = fₛ #var==:μΓ ? fᵧ : fₛ
					#println(var, " ", 𝑓)
					Fx = curve_fit(𝑓, Xin, Yin, ωᵢ, [100, 10.0]) #ωᵢ, 
					Yf = 𝑓(Xin, Fx.param)
					n = length(Yf)
					ϵ² = ωᵢ.*(Yin .- Yf).^2 #Fx.resid.^2
					#σ² = Y[ii, err].^2
					p = dof(Fx)
					ŷ = mean(Yin)
					Chi_sq = sum(ϵ²./σ²)/p
					VARᵣₑₛ = sum(ωᵢ.*(Yin .- Yf).^2) #/p
					VARₜₒₜ = sum(ωᵢ.*(Yin .- ŷ).^2) #/(n-1)
					r_sq = cor(Yin, Yf)^2 #cor(Yin, Xin)^2 #1 - VARᵣₑₛ/VARₜₒₜ # 
					Fx, r_sq, Chi_sq
				end
            	 
				try
                	con[i][Δz][var] = confint(tmp[i][Δz][var]) |> X->[[X[1][1], X[2][1]], [X[1][2], X[2][2]]]
				catch e
					println(e)
					println(var, ": ", tmp[i][Δz][var].param, " ", length(ii))
				end
			end
		end  # over variable ice type
    end  # end over co_status
	tmp, r2, chi2, con
end

# ╔═╡ f4f1cf08-4843-46ee-9da2-529b40321c94
let r2=Dict(:coupled=>[], :Z=>[])
for z in (:H, :L)
	
	for c in (:co, :de)
		push!(r2[:coupled], c==:co)
		push!(r2[:Z], ifelse(z==:H, -1, 1))
		for v in keys(rfit[c][z])
			!haskey(r2,v) && (r2[v]=[])
			push!(r2[v], rfit[c][z][v].param[2]) #R²[c][z][v]) # |> atand)
		end
	end
	end
	pretty_table(HTML, r2)
end

# ╔═╡ 848150e8-f3ed-4324-9964-8a6ada5df1ce
@bind ss Select([:co, :de], default=:co)

# ╔═╡ f4e1159d-37ac-47b0-89ff-a9581aa3875f
rfit[ss][:L][:μlwp].param , margin_error(rfit[ss][:L][:μlwp]), fₛ([30,15,5], rfit[ss][:L][:μlwp].param)

# ╔═╡ 045d019e-9817-4bcc-8af8-6882ad43c615
Δβ[:de][:L][:μlwp], rfit[:de][:L][:μlwp].param

# ╔═╡ 47ca6697-a1f3-49ee-a8b8-b2938d1fe1a4
begin
	tmpplt = []
	for mm in (1,2,3,4,11,12)
		kakes = @df filter(d->month(d.date)==(mm) && d.lwp>(0) && d.paₓ==(-1), DB) groupedboxplot(mm<10 ? :winter.+1 : :winter, :lwp, group=:coupled, outliers=false, color=farben, marker=(farben,stroke(0),1), bar_width=0.45, legend=false, ytick=[0,50,100], ylim=(0,170), xlabel=ifelse(mm==4,"year",""),xtick=(jahren, ifelse(mm≠4,"" ,jahren)), xlim=extrema(jahren).+(-0.3,0.3), top_margins=-0Plots.mm)
		push!(tmpplt, kakes)
	end
	plot(tmpplt..., layout=grid(6,1), size=(500,600), dpi=600)
end

# ╔═╡ 628ffc15-49dd-48d9-ac13-3ae2f700dee8
@df filter(d->d.lwp>(0) && d.iwp>(0), DB) scatter(:lwp, :iwp, ms=1, m=:+, xscale=:log10, yscale=:log10)

# ╔═╡ 6913d517-e07c-4962-9c6b-d090b1bb8faa
# Creating the string for winter label like 2023/24 for the wintertime 2023 to 2024:
#strwinter = [@sprintf("%04d/%02d", jj, (jj+1)-2000) for jj in jahren[1:end-1,1]];
strwinter = [@sprintf("%02d/%02d", jj-2000, (jj+1)-2000) for jj in jahren[1:end-1,1]];

# ╔═╡ 49654852-af74-468f-9ddd-cefa4f055907
begin
	var_meta = Dict(
		:lwp=>(unit="g m⁻²", labe="LWP", lege=L"\rm{\overline{LWP}}", lim=(1,150), stats=:geometric),
		:iwp=>(unit="g m⁻²", labe="IWP", lege=L"\rm{\overline{IWP}}", lim=(0.1, 80), stats=:geometric),
		:der=>(unit="μm", labe="Droplet  "*L"r_{eff}", lege=L"\overline{r_{eff}}", lim=(5, 30), stats=:geometric),
		:ier=>(unit="μm", labe="Ice  "*L"r_{eff}", lege=L"\overline{r_{eff}}", lim=(35,55), stats=:geometric),
		:T2m=>(unit="K", labe=L"\rm{T_{2m}}", lege=L"\rm{\overline{T_{2m}}}", lim=(240,274), stats=:aritmetic),
		:Γ=>(unit="K km⁻¹", labe=L"Γ_{\textrm{cloud}}", lege=L"\overline{Γ}_\textrm{cloud}", lim=(-1,12), stats=:aritmetic), #-7
		:δₕ=>(unit="m", labe="Cloud depth ", lege=L"\delta H_\textrm{cloud}", lim=(100, 5.9f3), stats=:aritmetic),
		:τc=>(unit="·", labe="Liquid "*L"\tau", lege=L"\overline{\tau}_\rm{cloud}", lim=(0.0,20), stats=:geometric),
		:τi=>(unit="·", labe="Ice "*L"\tau", lege=L"\overline{\tau}_\rm{cloud}", lim=(0.0,2), stats=:geometric),
		:clb=>(unit="m", labe="Liquid CBH", lege=L"\textrm{CBH}", lim=(10, 2f4), stats=:geometric),
		:cldTT=>(unit="K", labe="Cloud top T", lege=L"\rm{\overline{CTT}}", lim=(245, 273), stats=:aritmetic),
		:Tsurf=>(unit="K", labe="Skin Temp.", lege=L"\textrm{T_{skin}}", lim=(240,277), stats=:aritmetic),
		:μSIC=>(unit="%", labe="SIC", lege=L"\rm{\overline{SIC}}", lim=(10,110), stats=:truncated),
		:σSIC=>(unit="%", labe=L"\rm{\sigma_{SIC}}", lege=L"\rm{\sigma_{SIC}}", lim=(0,30), stats=:geometric),
		:AμSIC=>(unit="%",labe=L"\textrm{SIC}_\bigodot",lege=L"\textrm{SIC}_\bigodot", lim=(0,100), stats=:truncated),
		:AσSIC=>(unit="%",labe=L"\sigma_{\textrm{SIC}_\bigodot}",lege=L"\sigma_{\textrm{SIC}_\bigodot}", lim=(0,30), stats=:geometric),
		:RμSIC=>(unit="%",labe=L"\textrm{SIC_\star}", lege=L"\textrm{\overline{SIC_\star}}", lim=(0,100), stats=:truncated),
		:Pa=>(unit="kPa",labe="P", lege="Pressure", lim=(99,102), stats=:aritmetic),
		:aoi=>(unit="·",labe="AO", lege="AO", lim=(-5,5), stats=:aritmetic),
		:enso=>(unit="·",labe="ENSO", lege="ENSO", lim=(-2,5), stats=:aritmetic),
		:pdo=>(unit="·",labe="PDO", lege="PDO", lim=(-2,5), stats=:aritmetic),
	);
end

# ╔═╡ cd9c2bd6-50d2-4832-b35d-d4c40891d5e8
md"""
#### Subroutine for FFTW funciton:
"""

# ╔═╡ 34a57d0f-b3f3-4ed3-9988-0fe50a7c4764
md"""
SIC range: 	$(@bind siclim Select([(0,100)=>"All", (0,95)=>"(0,95]", (95,100)=>"(95,100]" ]))
Variable: $(@bind varva confirm(Select([:lwp, :iwp, :σSIC, :clb, :δₕ, :Tsurf, :der, :ier, :Γ, :cldTT, :T2m])))
Pressure: $(@bind paxva confirm(Select([:H, :L])))
"""

# ╔═╡ a09bbb85-edcc-4dac-a0b1-115ca73207f7
# [(0,100)=>"All", (0,65)=>"(0,65]", (65,95)=>"(65,95]", (95,98)=>"(95,98]", (98,100)=>"(98,100]"]

# ╔═╡ 6801d499-e203-45e3-b6a1-6e2878779787
"""
Function to simulate the GLM function ```predict(GLM.model)``` for LsqFit.jl package (temporal solution):
```julia-repl
julia> df_pr = predict_curve_fit(X, cfit)
julia> df_pr = predict_curve_fit(X, cfit; lev=0.9, Y, ω, pᵢ)
```

where:
* ```X::Vector``` independent varaible data,
* ```cfit::LsqFit.curve_fit``` the curve fitting output model
Optional arguments:
* ```lev::Real``` level of confident (default 0.95) for 95% CI
* ```Y::Vector``` dependent or model estimated variable, (default empty)
* ```ω::Union{Vector, AnalyticWeights}```, weights for Y (default all ones)
* ```pᵢ::UnitRange{Int}``` indexes of model parameters to consider (default=(1:2))
Output:
* ```DataFrame(prediction, sig_ci)```

Note that prediction is considered as ```y(X, b) = b[1] + b[2]X``` with only parameters ```b=cfit.param[pᵢ]``` are considered. If the model is non-linear, then use with caution and ensure the variable ```X``` is converted to fulfill the linear model.
"""
function predict_curve_fit(X::Vector, cfit; Y::Vector=[], lev::Real=0.95, ω::Union{Vector, AnalyticWeights}=[], pᵢ::UnitRange{Int}=(1:2))::DataFrame

	
	n = length(X)

	# obtaining the weights:
	ω = ifelse(isempty(ω), cfit.wt, ω) |> collect

	# finding if NaN is present in Y data:
	inan = isempty(Y) ?  isnan.(ω) : isnan.(Y)
	
	if length(ω)!=n && !isempty(Y)
		# making ω same length as Y when NaNs present:
		foreach(i->insert!(ω, i, 0), findall(inan))
	end
	
	α = 1- lev
	β = cfit.param
	np = length(β)
	ν = n - length(pᵢ)  # degrees of freedom for 2 parameters
	newx = [ones(n) X]
	# prediction y_hat with length n, same as X:
	y_hat = newx*β[pᵢ]

	# estimating residuals:
	y_res = isempty(Y) ? cfit.resid : ω.*(y_hat .- Y)

	y_pre = ifelse(np==2, 0.0, y_hat .- mean(y_hat))

	Δy_res = (y_res .- y_pre)[.!inan]
	MSE = sum(Δy_res.^2)/ν
	
	# Note: predict from GLM returns DataFrame(prediction=Ym, lower= -cretinterval, upper=2Ym.-cretinterval)
	
	# B = vcov(cfit)[pᵢ, pᵢ] # variance, co-variance matrix of the coefficients.
	#se = diag(newx*B*newx') .|> sqrt
	x_bar = mean(X)
	δx² = @. (X - x_bar)^2
	
	se = MSE.*(1/n .+ δx²/sum(δx²)) .|> sqrt
	ret_ci = se.*quantile(TDist(ν), α/2)
	
	return DataFrame(prediction=y_hat, sig_ci=ret_ci) #lower= ret_ci, upper=ret_ci)
end

# ╔═╡ ecbebaf3-95b9-44e4-bd33-84d278ed1c06
let mdf=filter(d->d.coupled==(false) && !isnan(d.μclb), dat[(Δz=-1,)])
	@df mdf scatter(:sic_bin, :μclb, yerror=:σclb, xflip=true) #Γ δₕ
	wk = mdf.σclb.^2 |> err2 -> inv.(err2)/mean(err2)
	
	fx(x,β) = @. β[1]+β[2]*exp((x/100)^2) #exp((x/100)^β[3]) # β[2]*(x/100)^β[3] #β[1]-β[2]exp((x/100)^2-1 ) #*log((1.01-x/100)) # β[2] *(1-x/100)+β[3]))
	
	fofo = curve_fit(fx, mdf.sic_bin, mdf.μclb, wk, [100, 10.]) #, .100])
	r2 = cor(mdf.μclb, fx(mdf.sic_bin, fofo.param))^2
	pp = fofo.param
	println(pp)
	x = mdf.sic_bin
	#println([ones(length(x)) exp.((x/100))]*fofo.param[1:2])
	fy = fx(x, fofo.param)
	# predict with intervals:
	cf = coef(fofo)
	ci = confidence_interval(fofo, 0.05)    # 5% significance level
	#println(ci)
	#tl, bl = fx(x, [ci[1][1], ci[2][2]]), fx(x, [ci[2][1], ci[1][2]]) #ci[1][1] .+ ci[2][2]*x,   ci[1][2] .+ ci[2][1]*x
	#σp, σm = maximum([tl bl], dims=2) .- 2fy,  2fy .- minimum([tl bl], dims=2)
	ypre = predict_curve_fit(exp.((x/100).^2 ), rfit[:de][:H][:μclb]) #fofo) #, Y=mdf.μΓ)
	yerr_low = @. ifelse(ypre.prediction + 3ypre.sig_ci ≤ 0.0, floor(ypre.prediction, digits=1)/3, 3.0ypre.sig_ci)
	yerr_hig = -3.0ypre.sig_ci
	
	plot!(x,ypre.prediction, ribbon=(yerr_low, yerr_hig), fillalpha=0.3, label="curve_fit", yscale=:log10)
	#plot!(x, fy, label="$(r2)", legend=:bottomright) #; hline!([cf[1]],label="wet")
end

# ╔═╡ aab6bedd-c076-491f-b590-852f179e860b
let
	# general parameters:
	nσ = 1 #0.5
	lf_lim = (-0.02, .55)
#	sic_lim = (7, 101) #(78, 101)
	sic_xin = range(0, 100, length=20) #max.(0.0, sic_bin) # (10:10:100) |> collect #(75:100) (5, (10:10:100)...)
	NNcol = cgrad(:starrynight, 15, categorical=true, alpha=0.5, rev=true, scale=:log10); #:grayyellow
	NNlim =(10, 1.5f3)
	tags = ["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(i)"]
	# **********************************************
	# LWP vs SIC
	a1 = []
	b1 = []
	c1 = []
	Nsicbin = length(SIC_bin)
	proxysic = 1:Nsicbin |> d->[0.15 .+ d..., d...]
	proxybin = range(1, Nsicbin, length=20)
	
	# TO PLOT THE GROUP OF VARIABLES ... CHAGE THE for([...]) LINE:
	# variable options:  #  
	foreach([:μlwp, :μder, :μiwp, :μier]) do var
	#foreach([:μδₕ, :μclb, :μcldTT, :μΓ]) do var
		# For other variables e.g. μder, μier, μcldTT, μlpr, μτc, μτi
		
		vv = replace(String(var), "μ" => "") |> Symbol
		svar = replace(String(var), "μ" => "σ") |> Symbol
		yy0 = diff([var_meta[vv].lim...])[1]
		
		foreach([(-1,:H),(1,:L)]) do (iz, P)
			
			tmpdf = rename(dat[(Δz=iz,)], var=>:var_bin, svar=>:var_std)
			xr2, yr2, ytag = if var==:μΓ && P==:L
				(0.5Nsicbin, 0.15yy0, 0.85yy0)
			elseif var==:μΓ && P==:H
				(0.9Nsicbin, 0.15yy0, 0.85yy0)
			elseif vv==:δₕ || vv==:clb
				(0.8Nsicbin, 0.38yy0, 0.45yy0)
			else
				(0.8Nsicbin, 0.83yy0, 0.85yy0)
			end

			# Scale for the y-axis depending on the variable:
			LogOrLinear = any(var ∈ (:μclb, :μδₕ)) ? :log10 : :identity

			# For coupled:
			𝑓σ = 3.0
			sic_yin = filter(d->d.coupled==true, tmpdf).var_bin
			# calculation prediction +/- CI exp.(sic_xin/100).^2)
			linear_sic_xin = @. exp((sic_xin/100)^2)
			pre_yin = predict_curve_fit(linear_sic_xin, rfit[:co][P][var])
			pre_low = @. ifelse(pre_yin.prediction + 𝑓σ*pre_yin.sig_ci ≤ 0.0 && LogOrLinear==:log10,
					floor(pre_yin.prediction, digits=2),
					-𝑓σ*pre_yin.sig_ci)
			pre_hig = -𝑓σ*pre_yin.sig_ci
		tp = plot(proxybin, pre_yin.prediction,
			ribbon=(pre_low, pre_hig),
			lc=:blue, lw=2, la=.7, fillalpha=0.2, fillcolor=farben[2], label="",
			ann=(xr2, yr2+var_meta[vv].lim[1], text(L"r^2="*@sprintf("\$%3.2f\$\n\$\\chi^2=%3.2f\$", R²[:co][P][var], Chi²[:co][P][var]), 12, color=farben[2])), top_margins=-4Plots.mm)
			
			# For decoupled: exp.((sic_xin/100).^2)
			sic_yin = filter(d->d.coupled==false, tmpdf).var_bin # 𝑓(sic_xin, rfit[:de][P][var].param)
			pre_yin = predict_curve_fit(linear_sic_xin, rfit[:de][P][var])
			pre_low = @. ifelse(pre_yin.prediction + 𝑓σ*pre_yin.sig_ci ≤ 0.0 && LogOrLinear==:log10,
					floor(pre_yin.prediction, digits=2), -𝑓σ*pre_yin.sig_ci)
					pre_hig = -𝑓σ*pre_yin.sig_ci
			
		plot!(proxybin, pre_yin.prediction,
			ribbon=(pre_low, pre_hig),
			lc=:red, lw=2, la=.7, fillalpha=0.2, fillcolor=farben[1], label="",
			ann=(xr2-0.25Nsicbin, yr2+var_meta[vv].lim[1], text(L"r^2="*@sprintf("\$%3.2f\$\n\$\\chi^2=%3.2f\$", R²[:de][P][var], Chi²[:de][P][var]), 12, color=farben[1]), :left), top_margins=ifelse(isempty(c1), 0, -4)Plots.mm, left_margin=-3.5(1+iz)Plots.mm)

		@df tmpdf scatter!(proxysic .- 0.5, :var_bin, yerror=:var_std, group=:coupled, 
			lc=farben, lw=2.0, la=0.5, marker=([:^ :o], 5), mc=farben, msw=1, msc=:grey,
			label=ifelse(length(c1)!=5,"", ["de" "co"]), legend=:topright, legendfontsize=8, legend_background_color=false,
			ann=(6.8, ytag + var_meta[vv].lim[1], popfirst!(tags)),
			xflip=true, xticks=(1:Nsicbin, ""), xlim=(1, Nsicbin+0.1), #sic_xin[1:2:end]
			yscale=LogOrLinear,
			#yformatter=ifelse( any(vv ∈ (:clb, :δₕ)), y->y/1f3, :auto),
		ylabel=ifelse(iz==1,"",var_meta[vv].labe*" [$(var_meta[vv].unit)]"), ylim=var_meta[vv].lim, # max.(0.1*sign(var_meta[vv].lim[1]), var_meta[vv].lim),
			left_margins=ifelse(P==:L, -5, 4)Plots.mm)
		push!(c1, tp)
		end
		
	end
	# retrieving the y-ticks values from H-pressure:
	foreach(enumerate(c1)) do (j, tp)
		iseven(j) && plot!(tp, yticks=(yticks(tp[1])[1], ""))
	end
	plot!.([c1[end-1], c1[end]], xticks=(1:Nsicbin, SIC_bin), xlabel="SIC [%]")  #sic_xin[1:2:end], sic_xin[1:2:end]
	plot!(c1[1], title="High Pressure")
	plot!(c1[2], title="Low Pressure")

	# creating the plots' mosaic with layout 4x2
	cc = plot(c1..., layout=(4,2), tickdir=:out, minorticks=true, guidefontsize=12, tickfontsize=11, size=(800,700), dpi=600, framestyle=:box) #left_margin=4Plots.mm,
	
	SAVEFIG ? savefig(cc, "/home/psgarfias/quicklooks/nsa_yearly/$(wintertime)_HL_MicPhys_SIC.png") : cc
	#savefig(cc, "/home/psgarfias/Downloads/quicklooks/nsa/$(wintertime)_HL_MicPhys_SIC.png")
	#savefig(cc, "/home/psgarfias/Downloads/quicklooks/nsa/$(wintertime)_HL_MacPhys_SIC.png")
end

# ╔═╡ d207a342-af32-4d36-b7d9-41505345ef87
"""
Function to compute t-Test parameters from a givem curve_fit model.
```julia-repl
julia> t_table = t_test4model(X, Y, ϵ_y, model, 𝑓)
julia> t_table = t_test4model(X, Y, ϵ_y, model, 𝑓; pᵢ=(1:4))
julia> t_table = t_test4model(X, Y, ϵ_y, model, 𝑓; fk=2, lev=0.90)
```
where:
* ```X::Vector``` independent variable,
* ```Y::Vector``` dependent variable,
* ```ϵ_y::Vector``` error for the dependent variable,
* ```model::LsqFit.LsqFitResults``` a output of ```curve_fit``` function,
* ```𝑓::Function``` used to obtain the ```model``` best fit.
Optional arguments:
* ```pᵢ::Int``` or ```::Indexes``` the index of the ```curve_fit``` parameter to output (default 2),
* ```fk::Real``` a factor to scate the outputs e.g. ```fk=2``` is ```2σ``` (default 10),
* ```lev::Real``` confidence level to use (default 0.95).
Output:
* ```t_table::NamedTuple{(coef::Real, se::Real, CI::Real, p_value::Real, snr::Real, Chi_sq::Real, R_sq::Real, RMSE::Real}```
or if pᵢ is a range like (1:2),
* ```t_table::NamedTuple{(coef::Vector, se::Vector, CI::Vector, p_value::Vector, snr::Vector, Chi_sq::Real, R_sq::Real, RMSE::Real}```
with the outputs for ```coef```, ```se```, ```CI``` are multiplied by the given optional input factor ```fk```.
"""
function t_test4model(X::Vector, Y::Vector, Yerr::Vector, model::LsqFit.LsqFitResult, 𝑓::Function; fk::Real=10, pᵢ::Union{Int, UnitRange{Int}}=2, lev::Real=0.95)
	μₜ = coef(model)
	σₜ = standard_errors(model)
	ν = dof(model)

	inan = @. !isnan(Y)
	Y = Y[inan]
	X = X[inan]
	Yerr = Yerr[inan]
	
	y̅ = median(Y)
	ŷ = 𝑓(X, model.param)
	
	# RMSE, Chi² and R²:
	ωt = model.wt[inan]
	#ϵ² = model.resid.^2   # include the weights!
	ϵ² = (Y .- ŷ).^2
	n = length(ϵ²)
	IQR = quantile(Y, [.25, .75]) |> diff |> first
	RMSE = √(sum(ϵ²)/n)/IQR

	
	SSE = sum(ϵ²)  #ωt.*
	SST = sum((Y .- y̅).^2) #ωt.*
	R_sq = 1 - SSE/SST

	σ² = eltype(Yerr)<:Tuple ? ((last.(Yerr) .- first.(Yerr))/2).^2 : Yerr.^2
	σ² = @. ifelse(round(σ², digits=3) != 0.000, σ², 1.0)
	
	Chi_sq = sum(ϵ²./σ²)/ν

	
	# T-statistics:
	tT = μₜ./σₜ
	
	CI = confint(model, level=lev)
	p_val = ccdf.(Ref(FDist(1, ν)), abs2.(tT))
	
	return (coef=fk*μₜ[pᵢ], se=fk*σₜ[pᵢ], CI=fk.*CI[pᵢ], p_value=p_val[pᵢ], snr=abs(tT[pᵢ]), Chi_sq=Chi_sq, R_sq=R_sq, RMSE=RMSE )
end

# ╔═╡ 409e45a7-79a9-429b-9889-de843dac49a7
# Preparing strings to be used in plots for legend, labels, etc:
begin
	siclimstr = siclim==true ? "All SIC" : @sprintf("%d < SIC ≤ %d", siclim...) #  ("SIC ∈ (%d,%d]", siclim...)
end

# ╔═╡ 55ea71af-8398-40bf-858d-1a039f20201f
 let fmt = Printf.Format("%s_%s_SIC%03d-%03d.csv")
  	dirfmt = Printf.Format("buffer_data/%s")
  	dirout = Printf.format(dirfmt, varva)
  	!isdir(dirout) && mkdir(dirout)
  	ravpathout = joinpath(dirout, Printf.format(fmt, "moavstat", varva, (siclim==(0,100) ? (100,100) : siclim)...) )
  	#CSV.write(ravpathout, ravstat)
  	wavpathout = joinpath(dirout, Printf.format(fmt, "ensostat", varva, (siclim==(0,100) ? (100,100) : siclim)...) )
  	#CSV.write(wavpathout, wavstat)
 end

# ╔═╡ e868b491-2463-4e0a-a7be-60f896eed14f
"""
Function to return a formated string matrix with either the r² (if f2=true) or the trend ± err, for the given ```df::DataFrame``` with the tabular information of the fit for all statistics and all ENSO/PDO frequencies.
```julia-repl
str = get_trend_str(df, :lwp, :q50, :H; r2=true)
```
"""
function get_trend_str(df::DataFrame, var::Symbol, stat::Symbol, pax::Symbol; vargof=nothing)

	# making String callable 2-var:
	@eval (kk::String)(x, y, U) = Printf.format(Printf.Format(kk), x, y, U)
	
	# making String callable 1-var:
	@eval (kk::String)(x, F) = Printf.format(Printf.Format(kk), x, F)
	
	if isnothing(vargof)
		# Selecting the variable and trend values to show:
		trend_str = select(df, [:hat, :sig, :coupled]=>ByRow((x,y,U)->"\$\\frac{\\Delta}{\\Delta t}\$=%+2.1f\$\\pm\$%2.1f %s"(x, y, "[$(var_meta[var].unit) dec⁻¹]")) => :lege, :type, :coupled, :pax)
		
	elseif vargof==:chi2freq
		trend_str = select(df, [:chi², :θ]=>ByRow((x, P)->"\$\\chi^2\$=%3.2f  @  \$\\nu_k^{-1}\$=%2.1f [yrs]"(x,inv(P)) ) => :lege, :type, :coupled, :pax)
		
	elseif vargof==:r2freq
		trend_str = select(df, [:r², :θ]=>ByRow((x, P)->"\$r^2\$=%3.2f  @  νₖ⁻¹=%2.1f [yrs]"(x,inv(P)) ) => :lege, :type, :coupled, :pax)
		
	elseif vargof==:r2chi2
		trend_str = select(df, [:r², :chi²]=>ByRow((x, P)->"\$r^2\$=%3.2f  &  \$\\chi^2\$=%3.2f"(x,P) ) => :lege, :type, :coupled, :pax)
		
	else
		trend_str = select(df, [:r², :rmse]=>ByRow((x, P)->"\$r^2\$=%3.2f  &  RMSE=%3.2f"(x,P) ) => :lege, :type, :coupled, :pax)
		
	end
		
	filter!(d->d.type==stat && d.pax==pax, trend_str)
	# converting it to format for [de co] legend:
	return reshape(trend_str.lege, 1,2)
end

# ╔═╡ 74c94058-4c91-4a7b-8060-29e4ac6104c2
md"""
## Function to model a time series $y_s(t)$ with trend and periodic component following:
$y_s(t, \beta; \nu_k) = \beta_1 + \beta_2 t + \beta_3\,cos(2\pi\,\nu_k\,t - \beta_4) + \epsilon$
with $\beta$ the parameters for best fit to the data, $\nu_s$ the frequency of the periodic signal. Optionally, it can be chosen $\nu_s$ as a fixed parameter initialized by e.g. the ENSO frequency value of $\nu_s$ = 0.176275, or as a fifth free paramter $\beta_5$ to be fitted.
To optain the time series trend, the coeficient $\beta_2$ in [U/yr] is isolated, with $t$ is in years and $y$ in [U].
"""

# ╔═╡ 6e6c10a4-1fc2-4b82-a619-608100d110f2
𝑦ₛ(t, β; νₛ=nothing) = if !isnothing(νₛ)
		@. β[1] + β[2]t + β[3]cos(2π*νₛ*t - β[4])
	else
		@. β[1] + β[2]t + β[3]cos(2π*β[5]*t - β[4])
	end

# ╔═╡ deef8546-317e-40c9-8992-b721e7e091e9
md"""
### Function to represent the trend of the time series:
$y_t(t, \beta) = \beta_1 + \beta_2 t + \epsilon$
with the parameter $\beta_2$ is the trend of the time series in [U/yr].
"""

# ╔═╡ a06b25e1-6897-4aab-bc5a-32b8bc29ebf7
𝑦ₜ(t, β) = @. β[1] + β[2]t

# ╔═╡ d8b6dba6-aeb8-4d49-bf49-e576338f08e3
"""
Function to fit the time-series data to a cycling function (as model given above).
```julia-repl
julia> enso_fit(Mdf, :lwp, :co, :H; θ=0.176, β₀=[100,2.1,10,0.5])
julia> enso_fit(Mdf, :lwp, :co, :H; fk=10)
```
Where:
* ```Mdf::DataFrame``` is the yearly/winterly data time series classified as :co/:de, :H/:L, for all variables.
* ```var::Symbol``` the variable to use for the fitting, e.g. :lwp, :iwp, etc.
* ```cc::Symbol``` the coupled status, e.g. :co or :de
* ```pa::Symbol``` the pressure system identifier, e.g. :H or :L
* ```θ::Real``` (optional) the frequency to use in the fitting as free parameter, default ```nothing```
* ```β₀::Vector{Real}``` (optional) vector for parameter to initiate the fitting, default ```[50, 1.5, 10, 0.1]```
* ```fk::Real``` (optional) factor to multiply the trend parameter from the fitting, default 1.

In case ```θ``` is not given, then it is no free parameter but rather a fitting parameter.
"""
function enso_fit(mdf::Dict, var::Symbol, cc::Symbol, pa::Symbol; θ = nothing, β₀=[50, 1.5, 10, 0.1], fk=1)
		
	ym(t, β) = 𝑦ₛ(t, β; νₛ=θ)
	
	β₀ = ifelse(!isnothing(θ), β₀, [β₀..., 0.176275])
	
	cc_var = Symbol("Y_", cc)
	cc_err = Symbol("ϵ_", cc) 
	# defining weights:
	xdat, ydat, ωt, σ² = let D = mdf[var][pa][:, cc_err]
		xdat = mdf[var][pa].winter
		ydat = mdf[var][pa][:, cc_var]
		
		err_sig = [(typeof(E)<:Tuple ? ((last(E)+first(E))/2) : E ) for E ∈ D].^2 .|> inv
		
		iout = @. isnan(err_sig) || ~isfinite(err_sig) #findall(x->isnan(x), err_sig)
		tot_sig = all(iout) ? 1.0 : mean(err_sig[.~iout])
		err_sig[iout] .= tot_sig
		W = @. err_sig/tot_sig #abs(1- (err_sig/ydat)) #abs(1- (err_sig/ydat)^2)
		iout = findall(!isnan, ydat)
		
		xdat[iout], ydat[iout], AnalyticWeights(W[iout]), inv.(err_sig[iout])
	end
	
	# Fitting the data to 𝑦ₛ with waights ωt :
	mfit = curve_fit(ym, xdat, ydat, ωt, β₀)

	# extract the fitted trend and its uncertainty and multiply them by factor fk (if given):
	mtest = t_test4model(xdat, ydat, sqrt.(σ²), mfit, ym; fk=fk)
	# combining R_sq and RMSE to a unified parameter to minimize:
	K = let 𝑉=[mtest.Chi_sq, mtest.RMSE^2]
		#, ([(1-R_sq)^2, RMSE^2]) #(sqrt∘sum)([inv(R_sq)^2, RMSE^2])
		#(0.0 ≤ mtest.RMSE ≤1.0) && push!(𝑉, mtest.RMSE^2)
		(-0.2 ≤ mtest.R_sq ≤1.0) && push!(𝑉, 1 - mtest.R_sq)
		(sqrt∘mean)(𝑉)
	end
	
	# Create the DataFrame output:
	outdf = DataFrame(type=var, hat=mtest.coef, sig=mtest.se, CI95=mtest.CI, pval=mtest.p_value, snr=mtest.snr, coupled=cc, pax=pa, r²=mtest.R_sq, rmse=mtest.RMSE, chi²=mtest.Chi_sq , cost=K )
	
	return outdf, mfit
end

# ╔═╡ 72a99203-ff45-477f-8258-b4da6bf96dab
# selection of DataFrame containing only the ENSO/PDO frequency selected:
#### wavstat = combine(groupby(allwavstat, :coupled), df->filter(d->d.θ==θₖ[d.coupled], df))

# ╔═╡ 1e247533-a600-40e1-beeb-f2ade77a4998
# md"""
# #### Finding the maximum r² from all possible fits depending on Θₚ frequencies:
# $(Kmin=filter(d->d.type==:q50 && d.pax==paxva, allwavstat) |> df-> combine(groupby(df, :coupled), :cost=>argmin =>:Kmin) |> df->permutedims(df,1) )
#"""

# ╔═╡ f5988ab9-9fb4-4209-837c-ec6328430307
#md"""
#Coupled period: $(@bind θcₖ Select(Θₚ.νₖ .=> Θₚ.Pₖ; default=Θₚ.νₖ[Kmin.co...]) )
#Decoupled period: $(@bind θdₖ Select(Θₚ.νₖ .=> Θₚ.Pₖ; default=Θₚ.νₖ[Kmin.de...]) )
#"""

# ╔═╡ 6846cf1f-9a49-4731-8552-189b8c095833
#@df filter(d->d.lwp<1f3 && d.iwp<1.5f3, DBraw) histogram(log10.([:lwp :iwp]), yscale=:log10); vline!([log10(5)])
#@df filter(d->d.winter==2016, DB) scatter(:lwp, :iwp, m=:+, ms=1, xscale=:log10, yscale=:log10, minorgrid=true); vline!([5 0.8f3]); hline!([5 3f3]) # filter(d->d.lwp>0 && d.iwp>0, DB)

# ╔═╡ 5d562554-9112-4856-b7a4-5e8ff22ea3d5
"""
Function to obtain the fitting trend using running average of time series.
```julia-repl
julia> rafits = co_de_fit(DB; fitmodel=false, window=10, edges=false)
```
WHERE:
* ```DB::DataFrame``` contains the data to fit and get trends,
* ```fitmodel::Bool``` (default ```false```),
* ```window::Int``` (dafault 5) the size of the average window to use for the averaging,
* ```edges::Bool``` (dafault ```true```) whether or not to include the edges of ```window/2```.
"""
function co_de_fit(newDB; fitmodel=false, window=6, edges=true)

	# Assigning keys for the type of statistics to use for fitting curve: e.g. q50 is median ± mad
	fkey = Dict(:q50=>f->Nlu(f;stats=:median), :q25=>f->Nlu(f;stats=:q25), :q75=>f->Nlu(f;stats=:q75), :μ=>f->Nlu(f;stats=var_meta[varva].stats))

	# Assigning output dictionaries: rfits contains the GLM fitting models, and dffit the DataFrame with values:
	tmpdic = Dict(k=>Dict() for k in keys(fkey))
	rfits = Dict(:co=>Dict(k=>Dict() for k in keys(fkey)),
					:de=>Dict(k=>Dict() for k in keys(fkey)) )
	dffit = Dict(k=>Dict() for k in keys(fkey))

	# ==========
	# Function to get the prediction and confidence interval: [prediction, prediction-lower, upper-prediction]
	function pred_low_hig(X, Y, rfit, W; lev=0.95)
		#garbage = predict(rfit, df, interval = :confidence, level = lev)
		garbage = if typeof(rfit)<:LsqFit.LsqFitResult
			#Ŷ = 𝑦ₜ(df.winter, rfit.param)
			predict_curve_fit(X, rfit; lev=lev, Y=Y, ω=W)
		elseif typeof(rfit)<:StatsModels.TableRegressionModel
			predict(rfit, df, interval = :confidence, level = lev)
		else
			@error("Argument rfit must be ::LsqFit.curve_fit or ::GLM.lm")
		end
		
		#return DataFrames.transform(garbage, names(garbage)=>ByRow((m,lo,hi)->[m-lo,hi-m])=>[:lower, :upper])
		return garbage
	end

	
	# =====
	# Looping over all function of the statistics to use:
	foreach(fkey) do (k, 𝑓) 
		println("worinking on ", k, ":")
		foreach([:H=>-1, :L=>1]) do (px,pz)
		tmdf = groupby(newDB, :winter) |> d->combine(d) do df
			isempty(df) && println(df.winter)
				Y_co = filter(d->d.coupled==(true) && d.paₓ==(pz), df) |> x-> 𝑓(x.Y_var)
				Y_de = filter(d->d.coupled==(false) && d.paₓ==(pz), df) |> x-> 𝑓(x.Y_var)
				(Y_co = Y_co[1], ϵ_co = Y_co[2], Y_de = Y_de[1], ϵ_de = Y_de[2])
		end
		
		# moving average for the coupled/decoupled time series Y_var ~ winter		
		tmdf[!, :S_co] = CLIMA.ave_window(tmdf.Y_co; w=window, edges=edges)
		tmdf[!, :S_de] = CLIMA.ave_window(tmdf.Y_de; w=window, edges=edges)
		
		# Defining the weights for coupled/decoupled time series based on winter's standard deviation:
		scale_it(E,M) = (typeof(E)<:Tuple ? last(E)-first(E) : E)/M |> x->ifelse(isnan(x), 0, x)
		function weight_it(err)
			sig_err = [(typeof(E)<:Tuple ? last(E)-first(E) : E)^2 for E in err] .|> inv
			iout = @. isnan(sig_err) || ~isfinite(sig_err)
			tot_err = all(iout) ? 1.0 : mean(sig_err[.~iout])
			sig_err[iout] .= tot_err
			
			return sig_err./tot_err
		end
		
		tmdf[!, :ω_co] = AnalyticWeights( weight_it(tmdf.ϵ_co)) #(@. abs(1 - scale_it(tmdf.ϵ_co, tmdf.Y_co)^2 ) )
		tmdf[!, :ω_de] = AnalyticWeights( weight_it(tmdf.ϵ_de)) #@. abs(1 - scale_it(tmdf.ϵ_de, tmdf.Y_de)^2 ) )
		
		# Fitting the smoothed time series to yₜ :
		f_co = filter(d->!isnan(d.S_co), tmdf) |> df-> curve_fit(𝑦ₜ, df.winter, df.S_co, df.ω_co, [10, 0.5]) #
		f_de = filter(d->!isnan(d.S_de), tmdf) |> df-> curve_fit(𝑦ₜ, df.winter, df.S_de, df.ω_de, [10, 0.5]) #
		
		#f_co = filter(d->!isnan(d.S_co), tmdf) |> df-> lm(@formula(S_co ~ winter), df; wts=ω_co)
		#f_de = filter(d->!isnan(d.S_de), tmdf) |> df-> lm(@formula(S_de ~ winter), df; wts=ω_de)

		# Assigning rfit
		rfits[:co][k][px] = f_co
		rfits[:de][k][px] = f_de
		
		# For coupled
		tmdf = hcat(tmdf, pred_low_hig(tmdf.winter, tmdf.Y_co, f_co, tmdf.ω_co) |> x->rename(x, [:lin_co, :err_co]) )
		# For decoupled
		tmdf = hcat(tmdf, pred_low_hig(tmdf.winter, tmdf.Y_de, f_de, tmdf.ω_de) |> x->rename(x, [:lin_de, :err_de]) )
			
		dffit[k][px] = tmdf
		end
	end
	
	return dffit, rfits #tmdf, coeff_co, coeff_de, Dict(:co=>f_co, :de=>f_de)
	
end

# ╔═╡ 76063b20-0848-4fd3-9e53-d76149d3d115
mdf, rfits = let var=varva
	low_threshold = 0.0 # var==:Tskin ? 1.5*var_meta[var].lim[1] : 0.0
    newDB = filter(d->d[var]>(low_threshold) && (siclim[1] < d.μSIC ≤ siclim[2]), DB)
    rename!(newDB, Dict(var=>:Y_var))
	# computing the data fitting using running window average method:
	co_de_fit(newDB; fitmodel=true, window=10, edges=true)
end;

# ╔═╡ fb63d204-e954-4813-bdc5-b02b67c080c5
let df=DataFrame(x=mdf[:q50][paxva].winter, y=mdf[:q50][paxva].S_de, Y=mdf[:q50][paxva].Y_de, e=mdf[:q50][paxva].ϵ_de)
	#(x=(1:12),y= 20.1 .+ 0.86*(1:12) .+0.8randn(12), e=rand(12))
	filter!(d->!isnan(d.y), df)
	
	lev = 0.95
	mfit = lm(@formula(y ~ x), df)
	#println("type of lm fit: ", typeof(mfit)<:StatsModels.TableRegressionModel)
	R_cholU = [1.6248076809271923 0.0; 3.059411708155671 25.495097567963924]
	n = length(df.x)
	newx = [ones(n) df.x]
	
	chol = cholesky!(mfit.model.pp)
	#print("pp=", chol isa CholeskyPivoted, " ")
	ip = invperm(chol.p)
    chol.U[ip, ip] |> println

	p = coef(mfit)
	lmint = confint(mfit, level=lev)
	#println("lm CI:", lmint)
	

	garbage = predict(mfit, df, level=lev, interval=:confidence) #DataFrame(prediction=
	insertcols!(garbage, 1, :x=>df.x)

	garbage[!, :lm] = newx*p #f(df.x, p)

	#println("pp=",cholesky!(mfit.model.pp).U)
	#println(mfit)
	
	residvar = ones(size(newx,2)) * deviance(mfit)/dof_residual(mfit)
	retvariance = (newx/R_cholU).^2 * residvar
	retinterval = quantile(TDist(dof_residual(mfit)), (1. - lev)/2)*sqrt.(retvariance)
	#println("retinterval", retinterval, "  residvar", residvar)
	garbage[!, :lm_lo] = garbage.lm .- retinterval
	garbage[!, :lm_up] = garbage.lm .+ retinterval
	
	# For curve_fit
	f(x,p) = @. p[1] + p[2]*x
	
	wwe = ifelse.(df.e .> 0.0, inv.(df.e.^2), 0.0) |> w->w./mean(w)
	#println("x=", df.x, " y=",df.y, " we=",wwe)
	cfit = curve_fit(f, df.x, df.y, wwe, [1,0.5]) # inv.(df.e.^2).*df.Y, 
	#println("type of curve_fit: ", typeof(cfit)<:LsqFit.LsqFitResult)
	cint = confidence_interval(cfit, 0.95)
	cresvar = ones(size(newx,2)) * sum(cfit.resid.^2)/dof(cfit)
	cretvariance = (newx/R_cholU).^2 * cresvar
	cretinterval = quantile(TDist(dof(cfit)), (1. - lev)/2)*sqrt.(cretvariance)
	
	#print("  curve_fit residvar=", cresvar)
	
	m = f(df.x, cfit.param)
	
	cstd = stderror(cfit)
	#println("coef=",cfit.param, "err=", cstd, "cint=", cint)
	
	garbage[!, :cfit] = m
	garbage[!, :cfit_lo] = m .+ cretinterval #f(df.x, cfit.param .- cstd)
	garbage[!, :cfit_up] = m .- cretinterval #f(df.x, cfit.param .+ cstd)
	
	qq = predict_curve_fit(df.x, cfit)
	#println(qq)
	@df df scatter(:x, :Y, label="data", mc=:gray) # yerror=:e,
	@df garbage plot!(df.x, :prediction, ribbon=(:prediction .- :lower, :upper .- :prediction), fillalpha=0.2, lc=:blue, fillcolor=:blue, label=round(10p[2], digits=3))
	#@df garbage plot!(df.x, :lm, ribbon=(:lm_lo, :lm_up), fillalpha=0.2, lc=:green, fillcolor=:green)
	#:cfit],   (:cfit_lo, :cfit_up)], fillalpha=0.2, linecolor=[:blue :green :red], fillcolor=[:blue :green :red])
	@df qq plot!(df.x, m, ribbon=:sig_ci, label="predic_") # [m.+:lower m.-:upper])
	#@df mdf[:q50][:H] plot!(:winter, :lin_de, ribbon=:err_de, label=round(10rfits[:de][:q50][:H].param[2], digits=3))
end

# ╔═╡ 73f5e706-d9c4-4d3e-b9bd-cd24c1b29dba
# here loop over rfits[c][k][pa] where c is (:co, :de), k is (:q50, :μ, ...), and pa is (:H, :L) then Vpa is the LsqFit model to use:
ravstat= [DataFrame((:type , :hat, :sig, :CI95, :pval, :snr, :chi², :r², :rmse, :coupled, :pax) .=> [k, values(t_test4model(mdf[k][pa].winter, mdf[k][pa][:, Symbol(:Y_, c)], mdf[k][pa][:, Symbol(:ϵ_, c)], Vpa, 𝑦ₜ))... , c, pa]) for (c,F) in rfits for (k,v) in F for (pa,Vpa) in v] |> D->reduce(vcat, D, cols=:union)

# ╔═╡ 6d4a626d-69a3-4538-8d6e-30e22a3e51ed
begin
	lwpts = let var=:T2m
		fk=1; fb=0 #-273.15;
		
		newDB = rename(filter(d->d[var]>var_meta[var].lim[1] && d.paₓ==1, DB), Dict(var=>:Y_var))  # && (d.aoₓ!=(0))
		newDB[!, :Y_var] .*=fk
		newDB[!, :Y_var] .+=fb
		tmpplot = @df newDB groupedboxplot(:winter, :Y_var, group=:coupled, whisker_range=0.5, bar_width=0.5, outliers=false, notch=true, la=0.7, fillcolor=farben, fillalpha=0.7, label="", xtick=(jahren, strwinter), xrot=33, xlabel="Wintertime", ylabel=var_meta[var].labe*" $(var_meta[var].unit)", size=(650,550), left_margin=5Plots.mm, bottom_margin=7Plots.mm, tickdir=:out, yminorticks=true, ytickfontsize=14, xtickfontsize=10, legendfontsize=12, guidefontsize=15, 
		ylim = var_meta[var].lim,
		);
		# Calculating mean of LWP co & de:
		#mdf, rfits = co_de_fit(newDB)
		coeff_co, coeff_de  = let subdf = filter(d->d.type==:μ && d.pax==:L, ravstat)
			dfco = filter(d->d.coupled==:co, subdf)
			dfde = filter(d->d.coupled==:de, subdf)
			(dfco.hat[1], dfco.CI95[1]...), (dfde.hat[1], dfde.CI95[1]...)
		end
		#@df mdf[:olr] scatter!(tmpplot, :winter.+[-0.1 0.1], [:Y_de :Y_co], mc=farben, label=false, m=:square) 
		@df mdf[:μ][:L] plot!(tmpplot, :winter.+[-0.1 0 0.1 0], [:Y_de :lin_de :Y_co :lin_co], seriestype=[:scatter :line],
		ribbon=[:err_de :err_co], fillcolor=[false farben[1] false farben[2]], fillalpha=0.3,
		m=[:^ :none :o :none], lc=[false farben[1] false farben[2]], lw=[0 2], ms=4, markerstrokewidth=1, mc=reshape(repeat(farben,2),1,4), 
		legend=:topleft,
		label=["(de)   "*var_meta[var].lege @sprintf("%+3.1f %s/yr @CI₉₅ [%+3.1f %3.1f]", coeff_de[1], var_meta[var].unit, coeff_de[2], coeff_de[3]) "(co)   "*var_meta[var].lege @sprintf("%+3.1f %s/yr @CI₉₅ [%+3.1f %3.1f]", coeff_co[1], var_meta[var].unit, coeff_co[2], coeff_co[3])], dpi=400)
		tmpplot
	end
end

# ╔═╡ e533eb73-dfaf-4928-9be6-72340edca990
@df DataFrame((:co, :de).=>[rfits[c][:q50][paxva].resid for c in (:co, :de)]) density([:co, :de], bandwidth=1)

# ╔═╡ 9e498f61-a297-4529-a7a1-1692e6cff222
CLIMA.ave_window(sind.(0:15:270); w=10, edges=false)

# ╔═╡ 2cc166c0-464b-41a4-be91-44d936f1eedf
@df filter(d->!isnan(d.lwp) && d.Tskin>(-100), DB) groupedboxplot(:winter, :Tskin, group=:paₓ, bar_width=0.4, outliers=false, fillcolor=farben, label=["de" "co"], xtick=jahren, ylabel="LWP [g m⁻²]", size=(850,400), left_margin=3Plots.mm, title="Pressure level L")

# ╔═╡ 4a5f18c6-3c6a-40f5-ad18-d1c46048ed3d
# ICE WATER PATH

# ╔═╡ c1c5b1ac-fd5c-4734-aa9a-cfc18b640eea
@df filter(d->!isnan(d.iwp), DB) groupedboxplot(:winter, :iwp, group=:coupled, bar_width=0.4, outliers=false, fillcolor=farben, label=["de" "co"], xtick=jahren, ylabel="IWP [g m⁻²]", size=(850,400), left_margin=3Plots.mm, title="AO index -")

# ╔═╡ 9a86e6d3-9e4c-4d14-b53f-81aba7e0e362
@df filter(d->!isnan(d.σSIC), DB) groupedboxplot(:winter, :σSIC, group=:paₓ, bar_width=0.4, outliers=false, fillcolor=farben, xtick=jahren, ylabel="IWP [g m⁻²]", size=(850,400), left_margin=3Plots.mm, title="Pressure level L") #; hline!([-2 2])

# ╔═╡ da33cd4a-7eec-4437-b31e-eabd91f5a8ae
md"""
## Time series Sea Ice concentration around NSA 50 km radius
"""

# ╔═╡ 4c56ee27-9a31-47e7-9959-05c342b0c4d2
begin
	let var=:AμSIC
		BW = 7
		pltsic = @df filter(d->d.AμSIC≥(0), DB) violin((:winter).-.015, :AμSIC, side=:left, label="\$\\oslash\$ 50km", xtick=(jahren, strwinter), xrot=35, bandwidth=BW, trim=false, lw=0, color=:gray)
		
		@df filter(d->d.μSIC≥(0), DB) violin!((:winter).+.015, :μSIC, side=:right, bandwidth=BW, trim=false, lw=0, color=:navy, label="\$ f(\\mathrm{WVT})\$", xtick=(jahren, strwinter), xrot=30, xlabel="Wintertime", ylabel="<SIC>  [$(var_meta[var].unit)]", xtickfontsize=9, ylim=(0, 100), ytickfontsize=14, minorticks=true, tickdir=:out, box=true, legend=:right, guidefontsize=14, size=(850,400), dpi=600, left_margin=3Plots.mm, bottom_margins=5.5Plots.mm)
	pltsic
	end
end

# ╔═╡ 2a4020e5-e353-4ca6-8b18-057c495209a6
tt = heatmap(repeat(vec([NNlim...]),1,2), color=NNcol, clim=NNlim, colorbar_scale=:log10, #colorbar_title=" No. occurence", 
colorbar_titlefontsize=15, tickfontsize=13, colorbar_tickswidth=4, xlim=(-4,-1), framestyle=:none, top_margins=1.5Plots.mm)

# ╔═╡ 93fbf31a-333b-4ce1-b4ec-e601146be8e2
md"""
## Time series for meteorological and climatological variables:
"""

# ╔═╡ 1d4fc3cc-d3c3-4b98-8b01-fd99a986a766
"""
Function to add ```:we``` to DataFrame including the fraction of the week to wintertime.
The value of ```:we``` will be then ```:winter+week(yy)/Δw``` where ```yy``` is the year and ```Δw``` id the period of winter time. By default the period is from Nov 1st to Apr 30th.

```julia-repl
julia> add_week_winter!(df)
julia> add_week_winter!(df; T0=Date(2012,11), Tf=Date(2013,04,30))
```
with ```T0``` and ```Tf``` indicating the starting and final Date of the dataset to be considered.
"""
function add_week_winter!(df::DataFrame; T0=Date(2012,11), Tf=Date(2013,04,30))
	w0 = week(T0)
	wf = week(Tf)
	Δw = w0 - wf
	w = week.(df.date)
	"winter" ∉ names(df) && add_winter!(df)
	insertcols!(df, 2, :we => @. df[:, :winter] + round((ifelse(w<w0, 51. + w, w) - 44.)/Δw, digits=3) )
	#df[:, :we] = 
	return nothing
end

# ╔═╡ 2f35259b-b297-4e2f-a535-7d48d0d01873
"""
Function to convert a DataFrame and make the data weekly median values, with added running average.

```julia-repl
julia> WeDF = make_it_weekly(df::DataFrame, vars::Vector{Symbol}; rename_it=[:idx=>:enso], nanrow=true)
julia> WeDF = make_it_weekly(df::DataFrame, vars::Vector{Symbol}; nanrow=true, W=6)
```
Output is the ```df::DataFrame``` but reduced to weekly values ```WeDF::DataFrame```.
Optional input arguments:
```nanrow::Bool``` if true then add a NaN row at the final of the continues winter period (default true).
```W::Numeric``` the running average window in number of data points (default w=4)
```rename_it::Vector{Tuple(Symbol, Symbol)}``` to rename the variable names e.g. [:idx=>:enso, :aoi=>:AOi]
"""
function make_it_weekly(dt::DataFrame, vars::Vector{Symbol}; rename_it=[], nanrow=true, W=4)
	Nvars = length(vars)
	tmp = combine(groupby(filter(d->d.winter≠0, dt), [:we, :winter]),
		vars .=> x->Nlu(x, stats=:median)[1], renamecols=false)
	
	# running average with windows w=4 over variables:
	transform!(groupby(tmp, :winter), vars .=> x->CLIMA.ave_window(Vector(x); w=W))

	# if rename_it has the type [:idx=>:aoi] then it will rename the dataframe column :idx to :aoi
	if !isempty(rename_it)
		vars_tmp = [Symbol(k,"_function")=>Symbol(v,"_function") for (k,v) in rename_it]
		rename_it = vcat(rename_it, vars_tmp)
		rename!(tmp, rename_it...)
	end
	
	# adding a NaN row when weeks has a jump, i.e. from Apr 30th to Nov 1st:
	if nanrow
		inan = findall(>(0), diff(tmp.winter))
		foreach(i->insert!(tmp, i+1, (tmp.winter[i]+0.99, tmp.winter[2], fill(NaN, 2Nvars)...)) , inan)
	end
	return tmp
end

# ╔═╡ ac5548e8-093b-496e-b742-eb12f0588a6c
# Loading database for PDO index: pdo_index.json or  ersst.v5.pdo.dat
pdo = let dt_path = joinpath(pwd(), "buffer_data", "clima", "pdo_index.json")
    dt=CLIMA.load_climate_index(dt_path; Tlim=(DateTime(1990,1), DateTime(2025,03)));
	add_week_winter!(dt)
	dt
end;

# ╔═╡ ac0eafdc-8b05-480a-ac82-4b7a9bdf4530
# Loading database for ENSO index: enso_index.json or enso_index.json
enso = let dt_path = joinpath(pwd(), "buffer_data", "clima", "enso_index.json")
    dt=CLIMA.load_climate_index(dt_path; Tlim=(DateTime(1990,1), DateTime(2025,03)));
	add_week_winter!(dt)
	dt
end;

# ╔═╡ 4cf8a75a-9838-4049-90ee-92028f9cf545
# Loading database for Arctic Oscilation index: 
# 
aoi = let dt_path = joinpath(pwd(), "buffer_data", "clima", "norm.daily.ao.cdas.z1000.19500101_current.csv")
        #"https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.cdas.z1000.19500101_current.csv" 
    df=CLIMA.load_climate_index(dt_path; Tlim=(DateTime(1990,1), DateTime(2025,03)));
	add_week_winter!(df)
	df
end;

# ╔═╡ 3766c2a6-4540-485a-b98d-07622c1f159b
# Merging datasets for time-series plot:
begin
	vars = [:T2m, :Pa, :AμSIC, :AσSIC, :paₓ]
	Nvars = length(vars)
	
	dfts = let dt = DBraw[:, [:date, vars...]] 
		add_week_winter!(dt)
		dt[:, :T2m] .-= 273.15

		# weekly average all quantities in 'vars':
		tmpL = make_it_weekly(dt, vars)

		# Adding the AOi climate index:
		tmpR = make_it_weekly(aoi, [:idx]; rename_it=[:idx=>:aoi])
		tmpL = outerjoin(tmpL, tmpR, on=[:we, :winter]) |> d->sort(d, :we)
		
		# Adding the ENSO climate variables index:
		tmpR = make_it_weekly(enso, [:idx]; rename_it=[:idx=>:enso], W=8)
		tmpL = outerjoin(tmpL, tmpR, on=[:we, :winter]) |> d->sort(d, :we)

		# Adding the PDO climate variables index:
		tmpR = make_it_weekly(pdo, [:idx]; rename_it=[:idx=>:pdo], W=8)
		outerjoin(tmpL, tmpR, on=[:we, :winter]) |> d->sort(d, :we)
	end	
end;

# ╔═╡ 2385ec55-60b2-44df-8c6f-a4de7f6da443
let df=dfts
	gr_vars = [:T2m, :Pa, :AμSIC, :AσSIC, :enso, :pdo] #
	
	nn = length(gr_vars) #-3 Nvars-1  # -1 because paₓ is not plotted.
	xwinter = (2011:2024)
	HLcolor = cgrad(:roma, 2, categorical=true)
	met = []
	label_letter = 'a'
	for vv in gr_vars
	 	
	 	vv == :paₓ && continue
	 	
		kakes = select(df, [:we, vv, :paₓ, Symbol(vv, "_function")] .=> [:we, :var, :pa, :rav])
		
		yylims = let tt = extrema(filter(!isnan, skipmissing(kakes.var)))
			dtt = diff([tt...,]).*(-0.03, +0.03)
			dtt .+= tt
		end
				
		# adding to the plot vertical lines to indicate the initial of wintertime:
		tmp = vline([xwinter...], lc=:gray, ls=:dash, label="")
	 	any(vv ∈ (:enso, :pdo, :aoi)) && hline!(tmp, [0], lc=:gray, ls=:dash, label="")

		# plotting the scattered data symbol colored by H- and L-pressure:
		@df kakes scatter!(:we, :var, markersize=3, markerstrokewidth=0, markeralpha=0.4,
			zcolor=replace(:pa, missing=>NaN), color=HLcolor, clim=(-1,1), label="", 
			ylabel=var_meta[vv].labe*" [$(var_meta[vv].unit)]",
			yguidefontsize=11, ytickfontsize=11, colorbar=false,
			#ann=(2013, maximum(skipmissing(:var)), "($(label_letter))"),
		)
		label_letter +=1  # adding the label letter for the next plot
		
		# adding the running average of variables on top as time series:
		@df dropmissing(kakes, :rav) plot!(tmp, :we, :rav,
			lw=1, lc=:black,
			xticks=(xwinter, ""), xlim=(2010.9, 2025.), ylim=yylims, 
			tickdir=:out, label="", bottom_margins=-4Plots.mm)
		
		# collecting each plot for variable into met[]:
	 	push!(met, tmp)
		
	 	# adding BoxPlots to the right of the time series for L and H-pressure:	 && !ismissing(d.pa) 
		tmp = @df filter(d->!isnan(d.var) && d.pa≠(0), dropmissing(kakes)) groupedboxplot(:pa, :var,
			group=:pa, outliers=false, notch=true, wister_width=:half, bar_width=3, 
			label="", fillcolor=[HLcolor[1] HLcolor[2]], fillalpha=0.4, 
			xticks=([-1.5,1.5],["H","L"]), xlim=(-4,4), 
			ylim=yylims, ymirror=true, ytickfontsize=11, yguideposition=:right, 
			tickdir=:out, left_margins=-7Plots.mm, bottom_margins=-4Plots.mm)
		
	 	push!(met, tmp)
	 end
	plot!(met[2nn-1], xticks=(xwinter, strwinter), xrot=30, xlabel="Wintertime [+2000 year]", bottom_margins=+1.5Plots.mm)
	
	plot(met..., layout=grid(nn,2, widths=(.85,.15)), size=(800,700), dpi=400)
	#savefig("/home/psgarfias/Downloads/quicklooks/nsa/TimeSeries_$(wintertime)_clima.png")
end

# ╔═╡ 06311e0c-129d-4485-8599-70e89c4806e0
let tmp=groupby(DB, :winter)
	plx = []
	
	Dsic = combine(tmp) do df
		#p=density(filter(!isnan, df.AμSIC) , trim=true, bandwidth=3) #AμSIC
		#nx = [p.series_list[1].plotattributes[:x]..., 100]
		siclin = range(1,99,length=30)
		s = StatsBase.fit(Histogram, filter(!isnan, df.μSIC), siclin) #|> StatsBase.normalize
		
		h = StatsBase.fit(Histogram, filter(!isnan, df.AμSIC), siclin) #
		Nh = maximum(h.weights)
		Ns = maximum(s.weights)
		#h = StatsBase.normalize(h)
		fNN(d) = [d..., NaN32] |> x->x/maximum(filter(!isnan,x))
		
		( x = h.edges[1], #[p.series_list[1].plotattributes[:x]..., 100],
		 y = fNN(h.weights), #[p.series_list[1].plotattributes[:y]..., NaN32], |> d->d/maximum(d)
		s = fNN(s.weights),
		 Nh = fill(Nh, length(h.edges[1])),
		Ns = fill(Ns, length(s.edges[1]))
		)
		
	end
	#print(ceil.(unique(Dsic.N)/1f3))
	foreach([isodd, iseven]) do ff
		
		p1 = @df filter(d->ff(d.winter), Dsic) plot(:x, :winter .- [:y :s], seriestype=:step, lw=[1 2], fillalpha=[0.2 0.4], fillrange=:winter, label="", color=[:red :gray], yminorgrid=true, yaxis=:flip,
		ygridalpha=0.7, ytickdir=:out, ytick=(unique(:winter), ""), ylabel=ff==isodd ? "Normalized distribution" : "",
		#ytick=(unique(:winter), ["    $(w)" for w in strwinter if ff(parse(Int, w[1:4]))]), 
		ann=(2, unique(:winter).-1.18, text.(["(20$(w))" for w in strwinter if ff(parse(Int, w[1:2]))], 8, :left)),
		ylim=extrema(:winter).-[1.2, 0], yrot=90,
		xlim=(0,100), xtickdir=:out, xminorgrid=true, xgridlinewidth=2, xgridlinestyle=:dash, xtickfontsize=12, xguidefontsize=14, xlabel="SIC [%]")
		
		let df = filter(d->ff(d.winter), Dsic)
			dh = unique(df.Nh)
			ds = unique(df.Ns)
			
			foreach(enumerate(zip(dh,ds))) do (i,tt)
				annotate!(p1, -6, unique(df.winter)[i].-[0, 0.6, 1.2], text.(["0", "0.5", "1"], 7, :left))
				annotate!(p1, 30, unique(df.winter)[i].-1.1, text(tt[1], 10, :left, :red))
				annotate!(p1, 50, unique(df.winter)[i].-1.1, text(tt[2], 10, :left, :gray))
			end
		end
		push!(plx, p1)
	end
	
	plot(plx..., layout=(1,2), size=(700,600), right_margin=2Plots.mm)
	savefig("/home/psgarfias/Downloads/quicklooks/nsa/$(wintertime)_SIC_2hist.png")
end

# ╔═╡ 67552873-5276-4955-9dff-358eaf9c390a
md"""
### Calculating the FFT for climate variables and NSA observables:
* First estimate the FFT frequencies for ENSO using all timeseries,
* then plot only Wintertime data along with best fit,
* plot frequencies for ENSO and NSA variables:
"""

# ╔═╡ 2d686363-8f9a-48dc-bfba-9535e7b73d70
# combine(groupby(aoi, :winter), :idx=>mean, renamecols=false) |> df->filter(d->d.winter≠0, df) |> df-> 
enso_fft = let df=enso
    yfft = CLIMA.FourierFrequencies(df.date, df.idx; P=Year)
    yfft
end;

# ╔═╡ 230869e6-0f06-4757-aac3-31dbf883f334
md"""
### Defining the free parameter Θₚ to consider as frequency of signal oscilation: e.g. ENSO periods Pₖ:
$(Θₚ = filter(k->k.YdBₖ≥18 && 2< k.Pₖ < 15, enso_fft) |> df->(νₖ=df.νₖ, Pₖ = round.(inv.(df.νₖ), digits=1)); )
"""

# ╔═╡ 82ccaaa0-f5d7-4c1f-93bc-37da45b2b048
# Calculating Fitting parameters for every Climatological significant cycle (e.g. ENSO)
# ====
allwavstat, efits = let tmp=DataFrame[]
	efit = Dict()
#	#tmp = [enso_fit(mdf, var, cc, pax; fk=10) for cc in (:de, :co) for var in (:q50, :μ, :q25, :q75) for pax in (:H, :L)]
	 for cc in (:de, :co)
		 vardic = Dict() #cc=>
		 for var in (:q50, :μ, :q25, :q75)
			 !haskey(vardic, var) && (vardic[var]=Dict(:H=>[], :L=>[], :θ=>[]))
			 for pax in (:H, :L)
				 for (k, νₖ) ∈ enumerate(Θₚ.νₖ)
					
				 	pax==:H && push!(vardic[var][:θ], Θₚ.νₖ[k] )
					 df, tmpdic = enso_fit(mdf, var, cc, pax; θ=νₖ, fk=10) #,  θ=θₖ[cc]
	 				push!(vardic[var][pax], tmpdic)
					 df[:, :θ] .= Θₚ.νₖ[k]
				 	push!(tmp, df)
					
				 end
#		 		#vardic[var]=Dict(pax=>ffit)
	 		end
	 	end
		efit[cc]=vardic
		
	 end
	reduce(vcat, tmp, cols=:union), efit
end;

# ╔═╡ 4cb72574-9b2a-4482-b12e-642eb137c0d1
Lmin = filter(d->d.type==:q50, allwavstat) |> df->combine(groupby(df, [:pax, :coupled]), :cost =>argmin => :Kmin) |> df->transform(df, :Kmin => ByRow(i->Θₚ.νₖ[i]) => :θmin)

# ╔═╡ e13860c4-cd7a-4d37-9dc4-cb2c01f177ad
inv.(Lmin.θmin)

# ╔═╡ cbeeec75-64bb-4507-aa33-df9b1b90cad1
# selection of DataFrame containing only the ENSO/PDO frequency selected:
wavstat = combine(groupby(allwavstat, [:pax, :coupled])) do Gdf 
	pp, cc = (Gdf.pax[1], Gdf.coupled[1])
	
	Imin = filter(d->d.pax==pp && d.coupled==cc, Lmin).θmin[1]
	filter(d->d.θ==Imin, Gdf)
end;

# ╔═╡ a1941c4d-28f9-4f25-a3ad-350324003a33
# Plotting the trend for different statistical variables, coupling and pressure systems:
begin
	ttplt=plot()
	let xvarstr=[:q25, :μ, :q50, :q75]
		xtikstr=["¹    ²\n"*L"Q_1", "¹    ²\n "*L"\overline{\mu}_g", "¹    ²\n "*L"\mu_{1/2}", "¹    ²\n "*L"Q_3"]
		ragdf = groupby(ravstat, :pax)
		wagdf = groupby(wavstat, :pax)
		
		invfarben = circshift(farben, (0,1))
		for (i,df) in enumerate([ragdf[(pax=paxva,)], wagdf[(pax=paxva,)]])
		
		df = transform(df, :type =>ByRow(x-> findall(==(x), xvarstr)[1]-1/3+i/5) =>:xstr)
			# plotting the statistically significant flag:
			@df df scatter!(ttplt, :xstr, :hat, group=:coupled, m=[:o :^], mc=:gray, ma=ifelse.(:pval .<0.05, .5,0), ms=10, label=ifelse(i==1 && siclim==true, ["p<0.05" ""],"") )

			# plotting the trend for both methods: 1) running average, 2) ENSO fit: (-first.(:CI95)+:hat, last.(:CI95)-:hat)
			@df df scatter!(ttplt, :xstr, :hat, yerror=(:sig), group=:coupled, m=[:o :^], ms=5, mc=invfarben,lw=2, lc=invfarben, xtick=((1:4), xtikstr),
			xlabel=ifelse(siclim[1]==0, "Statistic used for trend estimation",""), yminorticks=true, #ylim=ifelse(paxva==:L,(0,15),(-10,10)), 
			ylabel=ifelse(paxva==:H, L"\frac{\Delta}{\Delta t}"*var_meta[varva].labe*" [$(var_meta[varva].unit) decade⁻¹]", ""), 
			legendtitle=siclimstr, legendforegroundcolor=false, legendbackgroundcolor=false,
				guidefontsize=13, tickfontsize=13, label=ifelse(i==1 && siclim==true, ["co ± σ" "de ± σ"],""), tickdir=:out, framestyle=:box, bottom_margin=5Plots.mm, top_margin=1Plots.mm, size=(450,300))
			
		end
	end	
	ttplt
	#[findall(ravstat) for (ityp, typ) in enumerate(ravstat.type)]
end

# ╔═╡ 59a25151-b64c-44f5-b3a0-0ad7d5697fc8
let Ldeco = filter(d->d.pax==paxva, Lmin) |> df->permutedims(select(df,[:coupled, :θmin]), 1)
md"""
Coupled period: $(@bind θcₖ Select(Θₚ.νₖ .=> Θₚ.Pₖ; default=Ldeco.co[1]))
___ Decoupled period: $(@bind θdₖ Select(Θₚ.νₖ .=> Θₚ.Pₖ; default=Ldeco.de[1]))
"""
end

# ╔═╡ e4ef2843-07dc-4f7f-ade6-acef5612cb52
θₖ = Dict(:co=>θcₖ, :de=>θdₖ); println(θₖ)

# ╔═╡ d2c0d96e-ead0-4ba2-adc0-585cf161bb64
begin
	# plotting scattering data for coupled and decoupled:
	@df mdf[:q50][paxva] scatter(:winter .+[0.1 -0.1], [:Y_de :Y_co], yerror=[:ϵ_de :ϵ_co], mc=farben, lw=2, la=0.4, lc=farben,  marker=([:^ :o], stroke(0.01), 6), label=ifelse(siclim==(20,100), ["decoupled:  "*L"\mu_{1/2}\pm \sigma_{\textrm{mad}}" "coupled:   "*L"\mu_{1/2}\pm \sigma_{\textrm{mad}}"], ""),
	legend_column=2, legendfontsize=10, legendforegroundcolor=false, legendbackgroundcolor=false, legend_position=:top, 
	gridlinewidth=.5, frame=:box, tickdir=:out, yminorticks=true, ytickfontsize=12, yguidefontsize=14,
	xtickfontsize=11, xguidefontsize=11, title=siclimstr*" & "*String(paxva)*"-pressure"
	)

	# Estimaint the limits for the Y-Axis based on the data and given lims in 'var_meta[].lim' variable:
	yye_lims = extrema([var_meta[varva].lim...]) #, extrema(mdf[:q50][paxva].Y_co)..., extrema(mdf[:q50][paxva].Y_de)...])
	yye_lims = (yye_lims[1], yye_lims[2]*1.00)
	# finding the index of frequency to show and getting curve fit parameters:
		
	efit_co, efit_de = let Ico = argmin(abs.(efits[:co][:q50][:θ] .- θₖ[:co]))
		Ide = argmin(abs.(efits[:de][:q50][:θ] .- θₖ[:de]))
		efits[:co][:q50][paxva][Ico], efits[:de][:q50][paxva][Ide]
	end
	βc = efit_co.param; #..., θₖ]
	βd = efit_de.param; #..., θₖ]
	# Calculating uncertainty region:
	
	# plotting the periodic fitted curve:
	plot!(2011:0.1:2024, [w->max.(1, 𝑦ₛ(w, βd; νₛ=θₖ[:de])), w->𝑦ₛ(w, βc; νₛ=θₖ[:co])],
		lw=2, la=0.7, lc=farben, ls=:dash,
		label=get_trend_str(wavstat, varva, :q50, paxva; vargof=:chi2freq))
	
	# plotting the trend line (2nd coefficient from fitted curve):
	@df let w=mdf[:q50][paxva].winter
		
		n = length(w)
		newx = w * [0 1]; newx[:,1].=1;
		Ŷde = 𝑦ₜ(w, βd[1:2])
		Ŷco = 𝑦ₜ(w, βc[1:2])

		Y_de = mdf[:q50][paxva].Y_de
		
		ret_de = predict_curve_fit(w, efit_de, Y=Y_de) |> df->rename(df, [:lin_de, :err_de]) #, :up_de])
		transform!(ret_de, [:lin_de, :err_de] => ByRow( (p,s)-> (ifelse.(p+s≤0, 0.95floor(p, digits=1), -s), -s) ) =>:sig_de)
		
		Y_co = mdf[:q50][paxva].Y_co
		ret_co = predict_curve_fit(w, efit_co, Y=Y_co) |> df->rename(df, [:lin_co, :err_co]) #, :up_co])
		transform!(ret_co, [:lin_co, :err_co] => ByRow( (p,s)-> (ifelse.(p+s≤0, floor(0.9p, digits=1), -s), -s) ) =>:sig_co)
		
		ttcat = hcat(DataFrame(winter=w), ret_de, ret_co)
		#println(ttcat)
		ttcat
	end	plot!(:winter, [:lin_de :lin_co], ribbon=[(first.(:sig_de), last.(:sig_de)) (first.(:sig_co), last.(:sig_co))], lc=farben, lw=2, la=0.9, ls=:solid, fillalpha=0.2, fillcolor=farben, label=get_trend_str(wavstat, varva, :q50, paxva),
	xtick=(:winter,strwinter), xrot=30, xlabel = "Wintertime [+2000 year]", #yscale=:log10, 
	ylim=ifelse(varva==:lwp, (1, 500), yye_lims), #ifelse(varva != :Γ, yye_lims, (0, 10)),
	#ylim=(100, 4500), yscale=:log10, # for Γ_cloud (0, 10) #
	ylabel = var_meta[varva].labe*" [$(var_meta[varva].unit)]", bottom_margins=+2Plots.mm, legend=(0.1, 0.96)) #
	
end

# ╔═╡ 8e7dea67-d49a-442d-a76e-6bbd7172ebcc
begin
	stat_labe = Dict(:q50=>"median ", :q25=>"Q1 ", :q75=>"Q3 ", :μ=>"mean ");
	# Plotting time series with original data points:
	@df mdf[:q50][paxva] scatter(:winter .+[.05 -.05], [:Y_de :Y_co], yerror=[:ϵ_de :ϵ_co], mc=farben, lw=2, la=0.6, lc=farben,  marker=([:^ :o], 6), mscolor=:grey, msw=2, label=ifelse(siclim==(10,100), ["de: "*L"\mu_{1/2}\pm \sigma_{m}" "co: "*L"\mu_{1/2}\pm \sigma_{m}"], ""), legend_column=2, legendfontsize=10, 
		legendforegroundcolor=false, legendbackgroundcolor=false,
	gridlinewidth=.5, frame=:box, tickdir=:out, yminorticks=true, ytickfontsize=12, yguidefontsize=14,
	xtickfontsize=11, xguidefontsize=11,
	title=siclimstr*" & "*String(paxva)*"-pressure")

	# Plotting time series with smoothed data points:
	#@df mdf[:q50][paxva] scatter!(:winter .+[.05 -.05], [:S_de :S_co], marker=:x, ms=5, mc=farben, label=false)
	
	# Plotting fitted trend lines:
	
	@df  let df = mdf[:q50][paxva]
		
		transform!(df, [:lin_co, :err_co] => ByRow( (p,s)-> (ifelse.(p+s≤0, floor(0.9p, digits=1), -s), -s) ) =>:sig_co)
		transform!(df, [:lin_de, :err_de] => ByRow( (p,s)-> (ifelse.(p+s≤0, floor(0.9p, digits=1), -s), -s) ) =>:sig_de)
		#ttcat = hcat(DataFrame(winter=w), ret_de, ret_co)
		
		df
	end plot!(:winter, [:lin_de :lin_co], ribbon=[(first.(:sig_de), last.(:sig_de)) (first.(:sig_co), last.(:sig_co))], lw=2, la=0.9, fillalpha=0.2, lc=farben, fillcolor=farben, label=get_trend_str(ravstat, varva, :q50, paxva; vargof=:r2chi2).*"\n".*get_trend_str(ravstat, varva, :q50, paxva),
	xtick=(:winter, strwinter), xrot=30, xlabel = "Wintertime [+2000 year]",
	yscale=ifelse(any(varva ∈ (:δₕ, :clb, :lwp)), :log10, :identity), ylabel = var_meta[varva].labe*" [$(var_meta[varva].unit)]", ylim=ifelse(varva==:lwp, (5, 500), yye_lims),
	#ylim=ifelse(varva!=:Γ, var_meta[varva].lim.*(1,1.0), (0, 10)), # for Γ_cloud (0, 10) #var_meta[varva].lim.*(1,1.0)
	legend=(0.09, 0.9), bottom_margins=+2Plots.mm, left_margins=+2Plots.mm)
	
end

# ╔═╡ 72d1dc70-cc6e-46b1-a9e5-1140a24cf6fc
begin
	kakesplot=[]
	foreach(enumerate([:de, :co])) do (i, cc)
		Vy = Symbol(:Y_, cc)
		Verr = Symbol(:ϵ_, cc)
		Vrng = !isempty(kakesplot) && yticks(kakesplot[end])[1][1] #range(yye_lims[1], step=5, stop= round(yye_lims[2]/5)*5) #..., length=7) .|> round
		Anny = diff([yye_lims...,])[1]*0.03 + yye_lims[1]
		
		dat = rename(mdf[:q50][paxva], Vy=>:Yvar, Verr=>:Yerr)
		tmp = @df dat scatter(:winter, :Yvar, yerror=:Yerr, mc=farben[i], lc=farben[i], lw=2, la=.7, #lc=farben[i],
			marker=(ifelse(i==1, :^, :o), stroke(1.5), 6),
			label=ifelse(siclim==(10,100), ifelse(i==1, "decoupled:  "*L"\mu_{1/2}\pm \sigma_{\textrm{mad}}", "coupled:   "*L"\mu_{1/2}\pm \sigma_{\textrm{mad}}"), ""),
			title=ifelse(cc==:de, siclimstr*" &", String(paxva)*"-pressure"), title_position=ifelse(cc==:de, :right, :left),
			legend_column=1, legendfontsize=11, legendforegroundcolor=false, legendbackgroundcolor=false, legend_position=:topleft,
			gridlinewidth=.5, frame=:box, tickdir=:out, 
			yminorticks=true, ytickfontsize=12, yguidefontsize=14, yticks=ifelse(i==1,:auto, (Vrng,"")),
			ylims=ifelse(varva != :Γ, yye_lims, (0, 10)), # for Γ_cloud (0, 10) # yye_lims, 
			xtickfontsize=9, xguidefontsize=11, left_margins=ifelse(i==1, +5, -3)Plots.mm	)
			# 
		# retrieving fitting parameters:
		β = ifelse(i==1, efit_de.param, efit_co.param);
		efit_ln = ifelse(i==1, efit_de, efit_co)
		
		# plotting the periodic fitted curve:
		plot!(tmp, 2011:0.1:2024, w->𝑦ₛ(w, β; νₛ=θₖ[cc]), lw=2, la=0.7, lc=farben[i], ls=:dash, label=get_trend_str(wavstat, varva, :q50, paxva; vargof=:chi2freq)[i])

		@df let w=mdf[:q50][paxva].winter
			Y = ifelse(i==1, mdf[:q50][paxva].Y_de, mdf[:q50][paxva].Y_co)
			n = length(w)
			newx = w * [0 1]; newx[:,1].=1;
			Ŷ = 𝑦ₜ(w, β[1:2])
			
			ret_ln = predict_curve_fit(w, efit_ln, Y=Y) |> df->rename(df, [:lin_y, :err_y])
			transform!(ret_ln, [:lin_y, :err_y] => ByRow( (p,s)-> (ifelse.(p+s≤0, 0.95floor(p, digits=1), -s), -s) ) =>:sig_de)
			
			hcat(DataFrame(winter=w), ret_ln)
		end plot!(tmp, :winter, :lin_y, ribbon=(first.(:sig_de), last.(:sig_de)), lc=farben[i], lw=2, la=0.9, ls=:solid, fillalpha=0.2, fillcolor=farben[i], label=get_trend_str(wavstat, varva, :q50, paxva)[i], xtick=(:winter,strwinter),
	xrot=30, xlabel = "Wintertime [+2000 year]",
		yscale=ifelse(any(varva ∈ (:δₕ, :clb, :iw3p, :lwp)), :log10, :identity), ylabel = ifelse(i==1, var_meta[varva].labe*" [$(var_meta[varva].unit)]", ""), ylim=ifelse(varva==:lwp, (5, 500), yye_lims), 
		legend=(0.1, 0.95) ) # || varva==:δₕ ylims=(5, 4500), 
		
		push!(kakesplot, tmp)
		# adding the plot to the mosaic final plot:
	end
	plot(kakesplot..., layout=(1,2), size=(900,400), dpi=600, bottom_margins=6Plots.mm)
end

# ╔═╡ 53243445-6771-45bf-9b7e-4943cf8cef20
let dfcode = filter(d->d.type==:q50 && d.pax==paxva, allwavstat)
	invfarben=reverse(farben) # reshape(repeat(farben[j],3), 1,6)
	Y_vals = (0:0.25:1.2)
	scores_plt = []
	foreach(enumerate([:de, :co])) do (j, cc)
		tmp = @df filter(d->d.coupled==cc,dfcode) plot(:θ, [1 .- :r² :rmse :chi²], m=[:o :square :x], mc=farben[j], ms=[3 4 6], stroke=8, la=0.5, l=:dash, lc=farben[j], 
		xlabel=ifelse(j==2, "Signal period νₖ⁻¹ [years]", ""), xticks=(:θ, ifelse(cc==:de, "", round.(inv.(:θ), digits=1))), xflip=true, 
		yminorticks=true, 
		label=["1 - r²" "nRMSE" "χ²" "" "" ""], tickdir=:out, legend=:outerleft, legend_title=cc, legendforegroundcolor=false, tickfontsize=11, guidefontsize=12,
		title=ifelse(cc==:de, var_meta[varva].labe*": "*siclimstr*" & "*String(paxva)*"-pressure", ""),
		bottom_margin=ifelse(cc==:de, -4Plots.mm, 0Plots.mm)
		)
		vline!(tmp, [θₖ[cc]], lc=farben[j], lw=5, la=0.4, yticks=Y_vals, ylim=(0, 1.25), label=false, frame=:box)
		push!(scores_plt, tmp)
	end
	plot(scores_plt..., layout=(2,1), left_margins=-10Plots.mm)
end

# ╔═╡ fc318b5f-baa4-4683-a31e-c869543043e8
pdo_fft = let df=pdo 
    yfft = CLIMA.FourierFrequencies(df.date, df.idx; P=Year)
    yfft
end;

# ╔═╡ 74657565-b504-4360-9073-3d67e7d52099
aoi_fft = let df = filter(d->d.date>Date(1992,1), aoi)
	yfft = CLIMA.FourierFrequencies(df.date, df.idx; P=Year)
        yfft	
end;

# ╔═╡ 24d9a027-28f9-4eec-9be4-d6aba07ce202
let allvars = (:aoi, :pdo, :enso)
	winterplts=[]
	winterts = []
	for clivar in allvars
	tmpdf=select(dfts,[:winter, clivar] .=> [:winter, :var])
	filter!(d->!ismissing(d.var) && !isnan(d.var), tmpdf)
	Wtoms = 26 # (DateTime(2001,4,30)-DateTime(2000,11,1))/Millisecond(Dates.toms(Week(1)))
	fwedate(x) = modf(x) |> v->DateTime(last(v), 11,1) + Week(round(Wtoms*first(v)))
	tmpdf = combine(groupby(tmpdf, :winter), :var => (x->Nlu(x; stats=:aritmetic)[1]), renamecols=false)
	transform!(tmpdf, :winter => ByRow(fwedate) => :datum)
	
	
	yfft = CLIMA.FourierFrequencies(tmpdf.datum, Float32.(tmpdf.var); P=Year) #bar(:νₖ, :Yfft, xticks=:Pₖ, xflip=true)
	tsplt = @df tmpdf plot(:datum, :var, m=:o, xlim=(DateTime(2012,10), DateTime(2024,5)))
		push!(winterts, tsplt)
	psplt = @df filter(d->d.Pₖ>(1), yfft) bar(:νₖ, :Yₖ, bar_width=.001, fa=[1 0.5], fillcolor=[:grey :orange], xrot=45,
	 	xticks=(:νₖ, @. Printf.format(Printf.Format("%3.1f"), :Pₖ)), xflip=true, xlim=(1/12, 1/2), label=String(clivar), legend=:left)
		push!(winterplts, psplt)
	end
	plot(winterplts..., layout=(3,1))
	
end

# ╔═╡ 3d57920f-b46e-4cf9-8c4d-0b1624a166dd
let vars=Dict(:ENSO=>enso_fft, :AO=>aoi_fft, :PDO=>pdo_fft)
	label_letter = 'a'
	tmp = []
	Plims = (1.3, 45)
	foreach(vars) do (k,vv)
		sub_df = filter(d->Plims[1]<d.Pₖ<Plims[2], vv)
		xString_Ticks = if k==:ENSO
			(sub_df.νₖ, @. Printf.format(Printf.Format("%3.1f"), sub_df.Pₖ))
		else
			(sub_df.νₖ, "")
		end
		pltj = @df sub_df plot(:νₖ, :Yₖ, seriestype=:stem, m=:+, lw=0.15(:Pₖ[1]), fillcolor=:grey,
			xrot=45, xlim=(1/Plims[2], 1/Plims[1]), xticks=xString_Ticks, xflip=true, xtickdir=:out, xguidefontsize=12, xtickfontsize=10,
			ylim=(0, 335), ytickdir=:out, yguidefontsize=12,
			ann=(1/1.32, 290, text("($(label_letter)) "*String(k), halign=:left) ),
			legend=false, top_margins=ifelse(label_letter=='a',0,-5)Plots.mm)
		
		k==:ENSO && scatter!(pltj, Θₚ.νₖ, [9], ms=5, m=:^, mc=:red, label="")
		push!(tmp, pltj)
		label_letter +=1
	end
	plot(tmp..., layout=grid(3,1), xlabel=["" "" "Signal period (νₖ⁻¹) [year]"], ylabel=["" "Power spectrum per Δνₖ" ""], bottom_margin=+3Plots.mm)
	#savefig("/home/psgarfias/Downloads/quicklooks/nsa/ClimaIndex_FFT.png")
end

# ╔═╡ c35b719f-742b-4c85-b83d-71e6c98575cf
let vvstats = (
		T2m =Dict(:arg=>Dict(:stats=>:aritmetic,), :S=>1),
		Γ = Dict(:arg=>Dict(:stats=>:aritmetic,), :S=>2),
		δₕ = Dict(:arg=>Dict(:stats=>:median, :L=>0), :S=>3),
		lwp = Dict(:arg=>Dict(:stats=>:median,), :S=>4),
		iwp = Dict(:arg=>Dict(:stats=>:median,), :S=>5),
		μSIC = Dict(:arg=>Dict(:stats=>:aritmetic, :L=>0, :U=>100), :S=>6)
	)
	vars = keys(vvstats) #(:T2m, :Γ, :δₕ, :lwp, :iwp, :μSIC)
	
	Wtoms = 26 # (DateTime(2001,4,30)-DateTime(2000,11,1))/Millisecond(Dates.toms(Week(1)))
	fwedate(x) = modf(x) |> v->DateTime(last(v), 11,1) + Week(round(Wtoms*first(v)))
	# coupled==ifelse(cc==:co, true, false), DB) #d.we≠0 &&
	tmp = Dict(cc=>let df=filter(d-> d.paₓ==ifelse(cc==:H, -1, 1), DB)
				
		#Dict(vv=>combine(groupby(DB, :winter), vv => (x->Nlu(x; vvstats[vv]...)[1]) => :var) for vv in vars)
		Dict(vv=>begin
			tmpdf = combine(groupby(df, :winter), vv => (x->Nlu(x; vvstats[vv][:arg]...)[1]) => :var)
			filter(d->!isnan(d.var), transform(tmpdf, :winter => ByRow(fwedate) => :datum))
		end for vv in vars )
		end
	for cc in (:H, :L) ) #(:de, :co) )
	
	#ytp = Dict(cc=>Dict(vv=>CLIMA.FourierFrequencies(dd.datum, dd.var; P=Year) for (vv, dd) in tmp[cc]) for cc in (:de, :co) )

	
	ytp = Dict(vv=>let ttcat = DataFrame()
		foreach(tmp) do (cc,dd)
			
			ytp0 = CLIMA.FourierFrequencies(dd[vv].datum, dd[vv].var; P=Year)
			ytp0[:, :Δz] .= ifelse(cc==:H, -1, 1)
			ttcat = vcat(ttcat, ytp0)
		end
		ttcat
	end for vv in vars )
		
	HLcolor = cgrad(:roma, 2, categorical=true)
	tmplt = []
	labcha = 'a'
	foreach(vars) do vv
	  	kakes = plot()
	  	
	  	@df filter(d->d.νₖ≤(0.6), ytp[vv]) bar!(kakes, :νₖ .+ [-0.01ones(7)...,0.01ones(7)...], :Yₖ, bar_width=[0.012 0.012], group=:Δz,
		lc=[HLcolor[1] HLcolor[2]], color=[HLcolor[1] HLcolor[2]], fillalpha=0.4,
		xticks=(:νₖ, @. Printf.format(Printf.Format("%3.1f"), :Pₖ)), xflip=true, xrot=45, xlim=(0.01, 0.55),
	  	ann=(0.35, maximum(:Yₖ), "($(labcha)) "*var_meta[vv].labe), label=ifelse(vv==:T2m, ["H" "L"], "") )
	  	
	  	push!(tmplt, kakes)
		labcha +=1
	  end
	 plot(tmplt..., layout=grid(2,3), xlabel=["" "" "" "" "Signal period (νₖ⁻¹) [years]" ""], ylabel=["|FFT|² power stectrum [dB Δνₖ]" "" "" "|FFT|² power spectrum [dB Δνₖ]"  "" ""], guidefontvalign=[:bottom :top])
end
	#|> dd->plot(dd.winter, dd.var, marker=:o); plot!(2012:0.5:2024, T->cos(2π*enso_fft.νₖ[3]*T+1.9)) #

# ╔═╡ 1cb9ec34-8474-4af3-8877-87b396262ccf
begin
	paₘ = 101.32 #filter(!isnan, DBraw.Pa) |> mean
	taₘ = filter(!isnan, DBraw.T2m) |> mean
    R =  287.05  # [J kg⁻¹ K⁻¹]  Specific gas constant
    g₀ = 9.83  #[m s⁻²]
    ΔZ(Tₐ, P) = R/g₀*(Tₐ)*log(paₘ/P)  # [m] +taₘ/2
	Φ = [ΔZ(x, y) for x in (220:2:280), y in (95:110)]
	heatmap((220:2:280), (95:110), Φ', color=:RdBu, clim=(-300,300)); vline!([taₘ]); hline!([paₘ])
end

# ╔═╡ 3b2571ec-9e79-42c6-bcc8-c8e395506e17
begin
	tmp=plot(); ccol=cgrad(:jet, 6)
	[plot!(tmp, (88:117), ΔZ.(t, (88:117)), lc=ccol[k], label="$(t)") for (k,t) in enumerate(225:5:275)]
	vline!(tmp,[paₘ], label=false, lc=:black); hline!([-50 50], l=:dash, lc=:black, label=false, xminorticks=true)
end

# ╔═╡ 98df9486-a69b-4280-ba6c-2b12106b989d
groupby(DBraw, :coupled) |> gf->combine(gf) do df
	pa=filter(!isnan, df.Pa)
	dz=filter(!isnan, df.ΔZ) 
	(Pa=quantile(pa, (.25, .5, .75)), Dz=quantile(dz, (.25, .5, .75)), N=length(pa))
end

# ╔═╡ ad29355c-59a0-494e-96fa-805853b6166d
filter(!isnan, DBraw.Pa) |> d->quantile(d, (.25, .5, .75))

# ╔═╡ 45c00bda-bc4c-43fd-9976-897cad5fe77e
@df filter(d->d.ΔZ>(-9999) && isfinite(d.ΔZ), DBraw) density(:ΔZ, group=:coupled, lc=farben, trim=true, xscale=:identity, label=["D"  "C"], xlim=(-300,300)) # ; vline!([-38 15.9 70.9 ])ΔZ xlim=(235,278), 

# ╔═╡ dca14b2e-bb63-4f1d-97f5-ac221c40b62e
let kakes = @df filter(d->!isnan(d.μSIC), DB) fit(Histogram, :μSIC, collect(0:10:100))
	sic_cdf = cumsum(kakes.weights)/sum(kakes.weights)
	sic_xx = kakes.edges[1] #[1:end-1]
	cdf_plt = plot(sic_xx[1:end-1], sic_cdf)
	per = (0:0.2:1) #[.25, .5, .75, 1] #[0.16, 0.25, 0.5, 1]
	idxii = [1]
	siclims = quantile(DB.μSIC, per)
	println(siclims)
	for i ∈ 1:length(siclims)-1
		println(siclims[i], "-", siclims[i+1], ": ", findall(siclims[i] .≤ DB.μSIC .≤ siclims[i+1]) |> length)
		vline!(cdf_plt, [siclims[i]], label="$(siclims[i]) $(per[i])", l=:dash, lc=cgrad(:jet, categorical=true, 6)[i])
		hline!(cdf_plt, [per[i]], label=false, l=:dash, lc=cgrad(:jet, categorical=true, 6)[i], yticks=(0:0.2:1))
	end
	# foreach(enumerate(per)) do (j, pp)
	# 	ii = argmin(abs.(pp .- sic_cdf))
	# 	push!(idxii, ii)
	# 	println(pp, " $(sic_xx[idxii[end-1]]) - $(sic_xx[idxii[end]]) ", findall(sic_xx[idxii[end-1]] .< DB.μSIC .≤sic_xx[idxii[end]]) |> length )
	# 	vline!(cdf_plt, [sic_xx[ii]], label="$(pp) at $(sic_xx[ii])")
		
	# end
	cdf_plt
end

# ╔═╡ 0ba23046-f3b8-4c26-83d3-e06744c50587
let df25 = filter(d-> Date(2024,11,1) < d.date < Date(2024,12,30), DB)
@df df25 plot(:date, [:lwp 100*(:T2m .- 273).+500])
hline!([500], label=false)
hline!([Nlu(df25.lwp, stats=:median)[1] ], label="median", ylim=(0, 600))
end

# ╔═╡ 1238f26a-6a0d-42d9-a40d-256c671d5398
groupby(DB, :coupled) |> gdf->combine(gdf) do df
	
	ntot = filter(d->!isnan(d.clb), df) |> d->length(d.clb)
	nlls = filter(d->!isnan(d.clb) && d.clb≤100, df) |> d->length(d.clb)
	perlls = nlls/ntot
	println(df.coupled[1], perlls)
end

# ╔═╡ 8dede2cb-914a-4d46-aafe-f6a47bc2ed4d
findall(!isnan, DB.Tskin) |> length

# ╔═╡ e4578098-49d1-4ccc-a5f3-e6e9319b591d
groupby(DB, [:winter, :paₓ]) |> gb->combine(gb, :coupled=>sum) |> tt->plot(tt.winter, tt.coupled_sum, group=tt.paₓ)

# ╔═╡ 26fd5073-f21d-4bf8-ba24-b00b6074dac1
groupby(DB, :winter) |> gb->combine(gb, :coupled=>sum) |> tt->plot(tt.winter, tt.coupled_sum)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ATMOStools = "f4bd9c7d-eeda-4f89-9af3-b8c6f828feac"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
GLM = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
ATMOStools = "~0.1.0"
CSV = "~0.10.14"
DataFrames = "~1.7.0"
Distributions = "~0.25.108"
GLM = "~1.9.0"
HypothesisTests = "~0.11.2"
LsqFit = "~0.15.0"
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
StatsBase = "~0.34.3"
StatsPlots = "~0.15.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "f823fb2db1be1a4024caabc5590afb57a39d860f"

[[deps.ATMOStools]]
deps = ["CSV", "DataFrames", "Dates", "FFTW", "HTTP", "JSON3", "Printf", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "a496b2b9766c9e3f5dd97be26477f298829aa569"
repo-rev = "main"
repo-url = "git@github.com:pablosaa/ATMOStools.jl.git"
uuid = "f4bd9c7d-eeda-4f89-9af3-b8c6f828feac"
version = "0.1.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "0b4190661e8a4e51a842070e7dd4fae440ddb7f4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.118"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "b10bdafd1647f57ace3885143936749d61638c3b"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "273bd1cd30768a2fddfa3fd63bbc746ed7249e5f"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Rmath", "Roots", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "6c3ce99fdbaf680aa6716f4b919c19e902d67c9c"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "1d322381ef7b087548321d3f878cb4c9bd8f8f9b"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.1"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "40acc20cfb253cf061c1a2a2ea28de85235eeee1"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "3cebfc94a0754cc329ebc3bab1e6c89621e791ad"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.20"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "650a022b2ce86c7dcfbdecf00f78afeeb20e5655"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "24be21541580495368c35a6ccef1454e7b5015be"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.11"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "3a7c7e5c3f015415637f5debdf8a674aa2c979c4"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.1"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "9022bcaa2fc1d484f1326eaa4db8db543ca8c66d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.4"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "3b1dcbf62e469a67f6733ae493401e53d92ff543"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─123e1148-759c-11ed-29d7-e7fc1998111a
# ╟─f2a22774-dd59-485d-934d-15aa4568d379
# ╠═0cdaa10d-e72a-4b41-8e41-860317ced277
# ╠═48c90554-bc23-4e44-a050-4bce892614a2
# ╟─d734cb56-6da3-4d2f-8bca-a804bba0ba35
# ╠═17af9a85-0159-4c0a-b2f7-2dfb5e15bbbb
# ╟─5e79a910-f5b8-46ff-a5cb-d9204ff54cb6
# ╠═d56e6ef8-24cc-4afc-a155-5d411c3b9422
# ╠═e274212a-ef60-4258-930e-32d6d27f1e36
# ╟─eb3faf69-cab0-45d7-9484-28ab4963ed9c
# ╠═020695a4-61ee-49f5-8e39-4256ca836183
# ╟─3595bf88-4dbf-4c63-90e5-2feae1b78760
# ╠═ae59d379-b47d-41fc-8958-4e9e1941b528
# ╠═d1b9d2b9-6b08-4d81-8042-26bddf8f5e54
# ╠═9fe20fcc-8de8-4f4c-aff2-214797c40ba0
# ╟─89809ea8-94a4-4158-8f23-b4a2322532b2
# ╠═db52e7c0-5a8c-454a-bf1e-df9d377eff25
# ╠═f5a31896-eeb0-4a49-9178-3995fde42a55
# ╟─838b64cb-5719-4a82-94e0-9bd3f5e18598
# ╠═e530627f-4fa8-4006-b119-55f9efe72d21
# ╠═018fc44f-0e6b-4eb2-8b17-4ee61015fb7c
# ╟─16ee6021-e6ec-4b6c-aa91-1da8d8c07411
# ╟─7cb41328-9a82-4fe4-a65c-b96e424b4e88
# ╠═154f4210-a584-4040-9c03-3c9c51c467bd
# ╠═3af8f75f-0c94-4b4a-af8b-14d90585bd25
# ╠═263c6114-ee34-4da6-8466-5d7bd50dfd03
# ╠═3b1bd823-83bb-4630-be82-26b08c90f45a
# ╠═bbc95bb1-a296-4409-93be-ecf3c83c1153
# ╠═e1f28ccb-e252-4b44-95a3-4ffce7d45455
# ╠═5c9a28aa-5145-4a70-bc0e-5c299d3ae829
# ╠═6992ff1a-2c1c-44b0-854d-b260a26909c8
# ╠═d2accbb3-87d7-484c-82f1-ab6b488b653c
# ╠═f99cdbfe-2d31-483f-b1a2-664ddffd95f3
# ╠═1f6657ba-9910-48fc-84a3-399a28d84ef5
# ╠═3473d761-f57f-4d3f-beba-a8a516041da5
# ╠═07d4b3f0-d37e-4b5a-ae8f-46033e57b112
# ╠═76c14688-0be1-43a5-932e-352ebac308b4
# ╠═6b11b88f-66e4-4630-a9a0-bc479cda43e1
# ╟─ac217227-a672-4f9c-80ff-277c57e473f9
# ╠═49caf4c1-1ca2-483e-bfa7-36b98f3589be
# ╠═8e883748-effd-4898-b6ea-24de538fc5e8
# ╠═c5b02891-e561-4667-a07b-a29f52bb66c4
# ╠═8846a5e6-ab5a-4156-bc59-b19a2033577f
# ╠═264627a8-4975-4b75-b4de-136d9861aa0e
# ╠═bb176db2-fd18-4a59-9021-d2396a030b6d
# ╠═fa9e471b-3f58-459a-b99c-48cb1087e795
# ╟─75609282-c67c-4794-a68b-94d32c004783
# ╠═46c4e7cd-7d8a-49b7-abdd-1347692dcfd8
# ╠═acad36c9-c463-47c2-9277-efea563847df
# ╠═d04fad13-2fe2-4fa2-8528-24de8dafbbad
# ╠═f4f1cf08-4843-46ee-9da2-529b40321c94
# ╠═f4e1159d-37ac-47b0-89ff-a9581aa3875f
# ╠═848150e8-f3ed-4324-9964-8a6ada5df1ce
# ╠═045d019e-9817-4bcc-8af8-6882ad43c615
# ╠═47ca6697-a1f3-49ee-a8b8-b2938d1fe1a4
# ╠═628ffc15-49dd-48d9-ac13-3ae2f700dee8
# ╠═ecbebaf3-95b9-44e4-bd33-84d278ed1c06
# ╠═aab6bedd-c076-491f-b590-852f179e860b
# ╠═6913d517-e07c-4962-9c6b-d090b1bb8faa
# ╠═49654852-af74-468f-9ddd-cefa4f055907
# ╠═6d4a626d-69a3-4538-8d6e-30e22a3e51ed
# ╠═fb63d204-e954-4813-bdc5-b02b67c080c5
# ╟─cd9c2bd6-50d2-4832-b35d-d4c40891d5e8
# ╟─34a57d0f-b3f3-4ed3-9988-0fe50a7c4764
# ╠═a09bbb85-edcc-4dac-a0b1-115ca73207f7
# ╠═76063b20-0848-4fd3-9e53-d76149d3d115
# ╟─6801d499-e203-45e3-b6a1-6e2878779787
# ╟─d207a342-af32-4d36-b7d9-41505345ef87
# ╠═73f5e706-d9c4-4d3e-b9bd-cd24c1b29dba
# ╠═409e45a7-79a9-429b-9889-de843dac49a7
# ╠═8e7dea67-d49a-442d-a76e-6bbd7172ebcc
# ╠═55ea71af-8398-40bf-858d-1a039f20201f
# ╠═a1941c4d-28f9-4f25-a3ad-350324003a33
# ╟─e868b491-2463-4e0a-a7be-60f896eed14f
# ╟─74c94058-4c91-4a7b-8060-29e4ac6104c2
# ╟─6e6c10a4-1fc2-4b82-a619-608100d110f2
# ╟─deef8546-317e-40c9-8992-b721e7e091e9
# ╟─a06b25e1-6897-4aab-bc5a-32b8bc29ebf7
# ╠═3d45676a-8105-4e69-a520-605ee9cca2d1
# ╠═7b9dfaec-9bd3-4608-a932-6c7bf2e6f019
# ╠═d8b6dba6-aeb8-4d49-bf49-e576338f08e3
# ╟─230869e6-0f06-4757-aac3-31dbf883f334
# ╠═82ccaaa0-f5d7-4c1f-93bc-37da45b2b048
# ╠═72a99203-ff45-477f-8258-b4da6bf96dab
# ╠═1e247533-a600-40e1-beeb-f2ade77a4998
# ╠═4cb72574-9b2a-4482-b12e-642eb137c0d1
# ╠═e13860c4-cd7a-4d37-9dc4-cb2c01f177ad
# ╠═cbeeec75-64bb-4507-aa33-df9b1b90cad1
# ╠═f5988ab9-9fb4-4209-837c-ec6328430307
# ╟─59a25151-b64c-44f5-b3a0-0ad7d5697fc8
# ╠═e4ef2843-07dc-4f7f-ade6-acef5612cb52
# ╠═6846cf1f-9a49-4731-8552-189b8c095833
# ╠═d2c0d96e-ead0-4ba2-adc0-585cf161bb64
# ╠═72d1dc70-cc6e-46b1-a9e5-1140a24cf6fc
# ╠═53243445-6771-45bf-9b7e-4943cf8cef20
# ╠═e533eb73-dfaf-4928-9be6-72340edca990
# ╠═5d562554-9112-4856-b7a4-5e8ff22ea3d5
# ╠═9e498f61-a297-4529-a7a1-1692e6cff222
# ╠═2cc166c0-464b-41a4-be91-44d936f1eedf
# ╠═4a5f18c6-3c6a-40f5-ad18-d1c46048ed3d
# ╠═c1c5b1ac-fd5c-4734-aa9a-cfc18b640eea
# ╠═9a86e6d3-9e4c-4d14-b53f-81aba7e0e362
# ╟─da33cd4a-7eec-4437-b31e-eabd91f5a8ae
# ╠═4c56ee27-9a31-47e7-9959-05c342b0c4d2
# ╠═2a4020e5-e353-4ca6-8b18-057c495209a6
# ╟─93fbf31a-333b-4ce1-b4ec-e601146be8e2
# ╟─1d4fc3cc-d3c3-4b98-8b01-fd99a986a766
# ╟─2f35259b-b297-4e2f-a535-7d48d0d01873
# ╠═ac5548e8-093b-496e-b742-eb12f0588a6c
# ╠═ac0eafdc-8b05-480a-ac82-4b7a9bdf4530
# ╠═4cf8a75a-9838-4049-90ee-92028f9cf545
# ╠═3766c2a6-4540-485a-b98d-07622c1f159b
# ╠═2385ec55-60b2-44df-8c6f-a4de7f6da443
# ╠═06311e0c-129d-4485-8599-70e89c4806e0
# ╟─67552873-5276-4955-9dff-358eaf9c390a
# ╠═2d686363-8f9a-48dc-bfba-9535e7b73d70
# ╠═fc318b5f-baa4-4683-a31e-c869543043e8
# ╠═74657565-b504-4360-9073-3d67e7d52099
# ╠═24d9a027-28f9-4eec-9be4-d6aba07ce202
# ╠═3d57920f-b46e-4cf9-8c4d-0b1624a166dd
# ╠═c35b719f-742b-4c85-b83d-71e6c98575cf
# ╠═1cb9ec34-8474-4af3-8877-87b396262ccf
# ╠═3b2571ec-9e79-42c6-bcc8-c8e395506e17
# ╠═98df9486-a69b-4280-ba6c-2b12106b989d
# ╠═ad29355c-59a0-494e-96fa-805853b6166d
# ╠═45c00bda-bc4c-43fd-9976-897cad5fe77e
# ╠═dca14b2e-bb63-4f1d-97f5-ac221c40b62e
# ╠═0ba23046-f3b8-4c26-83d3-e06744c50587
# ╠═1238f26a-6a0d-42d9-a40d-256c671d5398
# ╠═8dede2cb-914a-4d46-aafe-f6a47bc2ed4d
# ╠═e4578098-49d1-4ccc-a5f3-e6e9319b591d
# ╠═26fd5073-f21d-4bf8-ba24-b00b6074dac1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
