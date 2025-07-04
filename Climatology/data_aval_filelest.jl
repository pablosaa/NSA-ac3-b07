# script to correlate list of files for synergy:

import ARMtools as ARM
using Dates
using DataStructures, DataFrames
using Plots

const SITE="utqiagvik-nsa"
const BASE_PATH=joinpath(homedir(),"LIM/remsens");

# Define list of directories to test:
PATH_LIST = Dict(
    1=>joinpath(homedir(), "LIM/data/B07", SITE,"CloudNet/input/"),
    2=>joinpath(homedir(), "LIM/data/B07", SITE,"CloudNet/input/"),
    3=>joinpath(homedir(), "LIM/data/B07", SITE,"CloudNet/input/"),
    4=>"",
    5=>"",
    6=>joinpath(homedir(), "LIM/data/B07", SITE,"CloudNet/input/"),
    7=>joinpath(homedir(), "LIM/data/B07") )

ALTE_LIST = Dict(
    1=>"KAZR/CORGE",
    2=>"MWR/LOS",
    3=>"",
    4=>"",
    5=>"GNDIRT",
    6=>"",
    7=>"SeaIce/ssmis")

DATA_LIST = OrderedDict(
    1=>"KAZR/ARSCL",
    2=>"MWR/RET",
    3=>"CEIL10m",
    4=>"INTERPOLATEDSONDE",
    5=>"RADFLUX",
    6=>"ECMWF",
    7=>"SeaIce/amsr2")

# Cloudnet path: joinpath(homedir(), "LIM/data/B07", SITE,""),
# Cloudnet alt path:, "CloudNet/1.9.0/output/CEIL10m"
# Cloudnet in DATA_LIST: "CloudNet/1.8.0/output/CEIL10m",

TimeSpan=(Date(2011,11):Date(2025,4,30));

# Defining global variables:
corrlist = DataFrame()
Nlist = length(DATA_LIST)

# creating short names:
#reshape([split(last(k),"/")[1] for k in DATA_LIST], 1, Nlist)
LABEL_LIST = reshape([replace(v, "/"=>": ")*ifelse(isempty(ALTE_LIST[k]), "", "/"*last(split(ALTE_LIST[k],"/"))) for (k,v) in DATA_LIST], 1, Nlist)

#corrlist = Dict(k=>0 for k in DATA_LIST) |> DataFrame
for td ∈ TimeSpan
    
    yy, mm, dd = year(td), month(td), day(td)
    5 ≤ mm ≤ 10 && continue
    
    tmp = joinpath(BASE_PATH, SITE)
    rowdf = OrderedDict{String, Any}("date"=>td)
    foreach(DATA_LIST) do (di, dl)
        path_file = isempty(PATH_LIST[di]) ? tmp : PATH_LIST[di]
        nfiles, faktor = let alte = ARM.getFilePattern(path_file, dl, yy, mm, dd)
            fak = 0
            if isnothing(alte) && !isempty(ALTE_LIST[di])
                @info "On $(td) - $(dl) alternative used!"
                alte = ARM.getFilePattern(path_file, ALTE_LIST[di], yy, mm, dd)
                fak=0.5
            end
            alte, fak
        end
        rowdf[LABEL_LIST[di]] = isnothing(nfiles) ? NaN16 : (faktor+di)
    end
    global corrlist = vcat(corrlist, DataFrame(rowdf))
end

# Analyzing part:
# Adding column with winter value:
jahren = Year.(extrema(corrlist.date)) |> J->J[1].value:J[2].value

let tmp = fill(0, length(corrlist.date))
    foreach(zip(jahren[1:end-1], jahren[2:end])) do (yb, yt)
        ii = findall((corrlist.date .≥ Date(yb,11,1)) .&& (corrlist.date .≤ Date(yt,5,1) ) )
        tmp[ii] .= yb
    end
    insertcols!(corrlist, 2, :winter=>tmp)
end
wintertime = unique(corrlist.winter)
Nwinters = length(wintertime)


# Plotting availability of data for selected time span:
colinstru = palette(:tab20, Nlist)

pltaval=[]
foreach(1:2) do j
    ninstru = ceil(Int, Nlist/2)
    ji = ninstru*(j-1) + 1
    je = min(j*ninstru, Nlist)
    nsub = je-ji+1
    labinstru = let tmp=[]

        foreach(ji:je) do id
            push!(tmp, "  $(id) $(LABEL_LIST[id])")
        end
        reshape(tmp, 1, nsub)
    end
    
    push!(pltaval,
          plot([1], fill(NaN16, 1, nsub),
               color=colinstru[ji:je]', label=labinstru, legend=:bottom,
               legend_font_pointsize=6, foreground_color_legend=nothing,
               axis=false, grid=false,
               lw=10, bottom_margins=-9Plots.mm, top_margins=2Plots.mm))

end

for i in (1:Nwinters)

    DB = filter(d->d.winter==wintertime[i], corrlist)
    tm_tick = DB[1, :date]:Week(1):DB[end, :date];
    str_tick = i≥(Nwinters-1) ? Dates.format.(tm_tick, "u/dd") : ""
    lineval, alphaval = let tmp=Matrix(DB[!, 3:end]) .|> modf
        aval = getindex.(tmp, 1) |> V->replace(x->isnan(x) ? 0 : -3x, V)
        aval .+= 2.5
        lval = getindex.(tmp, 2)
        lval, aval
    end
    
    tmp = plot(DB[!,:date], lineval,
               lw=alphaval, color=colinstru[:]', legend=false,
               xticks=(tm_tick, str_tick), xrot=60, xminorticks=true, xtickdir=:out,
               xlim=extrema(DB.date), xtickfontsize=5,
               ylabel=wintertime[i], ytickfontsize=4, ytickdir=:out,
               yguidefontsize=7, ylim=(0, Nlist+1), yticks=(1:Nlist),
               top_margins=i>2 ? -6Plots.mm : -1Plots.mm,
               left_margins=-3Plots.mm,
              )
    push!(pltaval, tmp)
end

# Final composite graph:
plot(pltaval..., layout=grid(8,2), dpi=600)

# end of script
