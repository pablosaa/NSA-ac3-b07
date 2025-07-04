#= script to create the Figure 1 for NSA paper:
Frame1: Beafurt-Chucki seas are with 50 km around NSA in the center,
Frame2: zoom-up to the north slope of Alaska,
Frame3: 55km square with Sea Ice Concentration values for given date.
=#

using NCDatasets
using CSV, DataFrames
using Navigation
using Dates
using Printf
using GMT 
using ARMtools
using Statistics
import ATMOStools as ATM
import SEAICEtools as SEAICE

#const SEAICE=include("/home/psgarfias/LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl")

# Defining parameters and Paths:
const CAMPAIGN = "utqiagvik-nsa";
# variables definition for SeaIce:
const DATA_PATH = joinpath(homedir(), "LIM/data"); #data") #
const SENSOR = "amsr2/HDF";
const PROC_PATH = joinpath(DATA_PATH, CAMPAIGN, "SeaIce");
const LATLON_FILE = joinpath(PROC_PATH, SENSOR, "LongitudeLatitudeGrid-n3125-ChukchiBeaufort.h5");

# defining Geo-location variables:
#nsa_lat = 71.323e0
#nsa_lon = -156.609e0
nsa_coor = Point(71.323e0, -156.609e0)

R_lim = 50f0

θ₀ = -133f0
θ₁ = 113f0
LonLim = (-159.5, -153.5) #(-158.8, -154.2)
LatLim = (70.7, 72.25) #(70.95, 72)

dd = 19 #23
mm = 01 #03
yy = 2019

# Reading the DB:
nsa = let dat = CSV.read("/home/psgarfias/LIM/scripts/NSA-ac3-b07/CoupledCloud_Seaice/buffer_data/csv_nsa/yearly/all_nsa_microphys_db_2011-2025.csv", DataFrame);
    filter!(d-> Date(d.date)==(Date(yy,mm,dd)), dat)
end
# getting the wind directions:
wvtdir = select(nsa, [:date, :wvtdir, :mwvt]) |> df -> filter(d->Minute(d.date).value==0, df);
wvtdirstr = Dates.format.(wvtdir.date, DateFormat("HH:MM"))

windlines = let dd1=SEAICE.create_pair_lines(nsa_coor, 1.15e3R_lim, wvtdir.wvtdir)
    dd2 = SEAICE.create_pair_lines(nsa_coor,
                                   0.9e3(R_lim.-25wvtdir.mwvt./maximum(wvtdir.mwvt)), wvtdir.wvtdir)
    txt=SEAICE.create_pair_lines(nsa_coor, 1.25e3R_lim, wvtdir.wvtdir)
    first(dd1)[2:2:end] = first(dd2)[1:2:end]
    last(dd1)[2:2:end] = last(dd2)[1:2:end]
    θ = mod.(270.0 .- wvtdir.wvtdir, 360.0)
    (lon=first(dd1), lat=last(dd1), dir=θ, len=wvtdir.mwvt./maximum(wvtdir.mwvt),
     txx=first(txt), txy=last(txt))
end

# Defining and reading SIC data:
sic_coor = SEAICE.read_LatLon_Bremen_product(LATLON_FILE);

idx_box = SEAICE.extract_LonLat_Box(LonLim, LatLim, sic_coor);

θ_all, ρ_all = SEAICE.LonLat_To_CenteredPolar(nsa_coor, sic_coor);

idx_sector = SEAICE.Get_Sector_Indexes(θ₀, θ₁, θ_all, ρ_all, R_lim=50f0);

sic_filen = ARMtools.getFilePattern(PROC_PATH, SENSOR, yy, mm, dd, fileext=".h5");

# Reading corresponding SIC data for the selected sector with idx indexes:
SIC = SEAICE.read_SIC_Bremen_product(sic_filen, idx_box);

# Creating canvas for Geo-Location:
#GMT.scatter(ρ_all[idx_box][iirad], SIC[iirad], region=(0,52, 0, 100), marker=:circle, show=0, fmt=:png, figsize=(6, 3))

lat = map(p->p.ϕ, sic_coor[idx_box]);
lon = map(p->p.λ, sic_coor[idx_box]);
# creating a semi-circle for visualization
P_circ = SEAICE.Create_Semi_Circle(nsa_coor, θ₀+360, θ₁, R_lim=51e3)
lon_circ, lat_circ = SEAICE.Get_LonLat_From_Point(P_circ)

# defining the colorbar scheme for SIC:
sic_bar = GMT.makecpt(cmap=:vik, truncate=(-.8, 0), range=[70, 80, 85, 90, 92, 94, 96, 98,100], continues=true);

GMT.scatter(lon, lat, color=sic_bar, zcolor=SIC, markersize=0.18, marker=:square,
	    region = (LonLim..., LatLim...), frame=:none, figsize=9,
	    proj = (name=:stereographic, center = [-156.5, 90], paralles=30),
	    xlabel="lon", ylabel="@:2:lat",
	    title="", show=0)

#GMT.coast(region = (-200.0, -100.0, 55.0, 85.0), show=0)
GMT.coast!(region = (LonLim..., LatLim...), #(-180.0, -130.0, 68.0, 75.0),
	   proj =(name=:stereographic, center = [-156.5, 90], paralles=30),
	   figsize=9,
           xaxis = (annot="60m", ticks="60m", grid="60m"),
           yaxis = (annot="30m", ticks="30m", grid="30m"),
           offset = (4, 8),
           map_scale="jTL+c40+w50+f+o0.5/0.5+lkm", 
	   res=:high, area=950, land=:lightbrown, background=0, #figscale="1:1900000",
	   rose = (inside=true,
                   anchor=:TR,
                   width=1.5,
                   color=:green,
                   offset=(-0.8, -0.03),
                   label=true),
           shore=:thinnest,
           inset=(

               coast, R=[-183, -120, 59, 77],
               frame = :a1g,
               center = [-156.5, 90],
               proj=:merc,
               shore=:thinnest,
               land=:lightbrown,
               water=:azure,
               rect=(0.35, :red),
               res=:middle),

           show=0)

GMT.text!(text="NSA", font=(8,"Helvetica", :red), region_justify=:C, offset=(-0.17, 2.4), noclip=true)
GMT.text!(text="ALASKA", font=(11,"Helvetica-Bold", "white=thin"), region_justify=:C, offset=(-0.25, 0.1), noclip=true)

# adding SIC colorbar:
GMT.colorbar!(pos=(anchor=:CR, length=(4, 0.2),offset=(0.35, -1.0), triangles=:b), color=sic_bar, frame=(ylabel="SIC [%]  ",), show=0)
#GMT.plot!(windlines.lon, windlines.lat, ls=:dash, zcolor=wvtdir.mwvt, color=wvt_bar, lw=0.5, alpha=10, show=0)

wvt_bar = GMT.makecpt(cmap=:jet, truncate=(0, 1), range=(0, 20, 2), continues=true, nobg=true)
GMT.colorbar!(pos=(anchor=:TC, length=(5, 0.15),offset=(-1.0, 0.4)), color=wvt_bar, frame=(ylabel="@~\321@~@-z@-WVT [kg m@+-2@+ s@+-1@+]",), show=0)
foreach(1:length(wvtdir.mwvt)) do i
    j = 2(i-1)+1
    GMT.plot!(windlines.lon[j:j+1], windlines.lat[j:j+1], 
              arrow=(len=0.09, start=:tail, shape=0.0005, angle=45),
              level=(data=wvtdir.mwvt[i], nofill=true), lw=0.5,
              color=wvt_bar, show=0)
end

# Including clock hours to the WVT lines:
idxin = [(1:2:7)..., (11:4:33)..., 45]
GMT.text!(mat2ds([windlines.txx[idxin] windlines.txy[idxin]], wvtdirstr[(idxin.÷2 .+1)]), font=(6.5, "Courier"), rotate=30, show=0)

GMT.plot!(lon_circ, lat_circ, lc=:lightred, lw=0.9, ls=:dash, alpha=20, show=0)
GMT.plot!(nsa_coor.λ, nsa_coor.ϕ, symbol=(symb=:+, size=.3), mc=:red, size=(600,500), show=1, savefig="NSA_SIC_50km.png")

# move SIC[%] a bit higher.
# Hour of lines non-linearly separated and rotated.
# add km to scale.
# change color of NSA to red.
# change the limits of main frame around NSA.

# end of script.
