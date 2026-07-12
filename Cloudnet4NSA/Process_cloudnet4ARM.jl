#!/home/psgarfias/.local/bin/julia


# Script to process ARM data with CloudNetpy
using CloudnetTools.ACTRIS
using CloudnetTools.Vis
using Printf

site_campaign = "utqiagvik-nsa"; # "arctic-mosaic";  # 
const CLNTVER="V1_68_1"

CLOUDNET_PATH = "/home/psgarfias/LIM/CloudNet"
DATA_PATH = "/projekt2/ac3data/B07-data/$(site_campaign)/CloudNet/input"; #"/home/psgarfias/LIM/remsens/CloudNet/input/";
CATE_PATH = "/projekt2/ac3data/B07-data/$(site_campaign)/CloudNet/output/$(CLNTVER)"; # "/home/psgarfias/LIM/data/CloudNet/output/";
PLOT_PATH = "/projekt2/ac3data/B07-data/$(site_campaign)/CloudNet/plot/$(CLNTVER)"; #"/home/psgarfias/LIM/data/CloudNet/plots/";

# defining range of processing:
years = (2025) #:2024)
months = (11,12)
days = (1:31)

!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

println("** Calculating for the periods of years=$years, months=$months and days=$days")
##__PAC__PIPE__##
# **** WICHTIG ******
# definitions for Cloudnet output products:
PROD_TYPE = "CEIL10m"; #or "HSRL"; # 

# defining ARM product to be used:
ARMprod = Dict(
    :site => "utqiagvik-nsa", #"mosaic",
    :radar => ("KAZR/ARSCL","KAZR/CORGE"),
    :lidar => PROD_TYPE,
    :mwr => ("MWR/RET","MWR/LOS"),
    :model => "ECMWF", # "ICON-IGLO-108-119",
    :radiosonde => "INTERPOLATEDSONDE",
)

# alternatice model name: e.g. const model_alt = "mosaic_icon-iglo-108-119";
const NAME_ALT_MODEL = "arm-nsa_"*lowercase(ARMprod[:model]);

CLNTprod = Dict(
    :categorize => true,
    :classification => false,
    :lwc => true,
    :iwc => true,
    :drizzle => true,
    :der => true,
    :ier => true,
)

PLTNTprod = Dict(
    :categorize => false,
    :classification => false,
    :lwc => false,
    :iwc => false,
    :drizzle => false,
    :der => false,
    :ier => false,
)


##
for yy in years
    for mm in months
        for dd in days
                                
            datestr = @sprintf("%04d/%04d%02d%02d_%s_", yy, yy, mm, dd, ARMprod[:site])

            # Creating Dictionary for input files (model change to NAME_ALT_MODEL if given):
            data_files = Dict(k => let tmp=""
                                  for fn in (typeof(V)<:Tuple ? V : (V,))
                                      tmp = joinpath(DATA_PATH, fn, datestr*String(k)*".nc")
                                      k == :model && (tmp=replace(tmp, ARMprod[:site]*"_model"=>NAME_ALT_MODEL))
                                      k == :radiosonde && (tmp=replace(tmp, "radiosonde"=>"ecmwf") )
                                      !isfile(tmp) ? (tmp=nothing) : break
                                  end
                                  tmp
                              end for (k, V) in ARMprod if k!=:site)
            
            # Cheking availability for all files
            # map(isfile, values(data_files)) |> all
            aval_flag = filter(d-> first(d) != :radiosonde && isnothing(last(d)), data_files)
            
            !isempty(aval_flag) && (println("Skipping $dd.$mm.$yy due to lack of $(keys(aval_flag))"); continue)

            # file names for CloudNet output:
            OUTPUT_PATH = joinpath(CATE_PATH, PROD_TYPE)

            # *******************************************
            # File names for retrived products:
            clnet_file = Dict(k => joinpath(OUTPUT_PATH, datestr*String(k)*".nc") for (k,flag) in CLNTprod if flag)

            # Checking if output files' path exist, otherwise create it:
            foreach([dirname(clnet_file[k]) for (k, B) in CLNTprod if B]) do tmp_dir
                !isdir(tmp_dir) && (mkpath(tmp_dir); println("$(tmp_dir) created!"))
            end

            
           ### Starting the Py Block:
	    if CLNTprod[:categorize]
                cateuuid = ACTRIS.categorize_it(data_files, clnet_file[:categorize])
                if isnothing(cateuuid) && isfile(data_files[:radiosonde])
                    @info("Trying with Radiosonde as :model")
                    data_files[:model] = data_files[:radiosonde]
                    cateuuid = ACTRIS.categorize_it(data_files, clnet_file[:categorize])
                    isnothing(cateuuid) && @warn("Unsuccess :model neither ECMWF or RADIOSONDE")
                end
            end

            isfile(clnet_file[:categorize]) && ACTRIS.generate_products(clnet_file, CLNTprod);


            #=
            #ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python")
            pushfirst!(PyVector(pyimport("sys")."path"), CLOUDNET_PATH);
            pyimport("sys").executable

            =#
            # File names for plot of products:
            ~PLTNTprod[:categorize] && continue

            VIS_PATH = joinpath(PLOT_PATH, PROD_TYPE)
            #pltnet_file = Dict(k => joinpath(VIS_PATH, datestr*String(k)*".png") for (k,flag) in PLTNTprod if flag)
            pltnet_file = replace(clnet_file[:categorize], "categorize.nc"=>"classify-data.png")

            isfile(clnet_file[:categorize]) && Vis.show_measurements(clnet_file[:categorize], savefig=pltnet_file, showclassific=true);
#
##
##            catcnet = pyimport("cloudnetpy.categorize")
##            global procnet = pyimport("cloudnetpy.products")
##	   
##		for (K, V) in CLNTprod
##			ex = Meta.parse("uuid=procnet.generate_$K(clnet_file[:categorize], clnet_file[$(K)])")
##			try
##				V && eval(ex)
##			catch e
##				println("\e[1m\e[38;2;230;30;30;249m","* Error trying to run $K from $dd.$mm.$yy");
##				println(e)
##			end
##		end 
##            #try
            #    uuid = catcnet.generate_categorize(data_files, clnet_file[:categorize])
            #    uuid = procnet.generate_classification(clnet_file[:categorize],
            #                                           clnet_file[:classification])
            #    uuid = procnet.generate_iwc(clnet_file[:categorize], clnet_file[:iwc])
            #    uuid = procnet.generate_lwc(clnet_file[:categorize], clnet_file[:lwc])
            #    uuid = procnet.generate_drizzle(clnet_file[:categorize],
            #                                    clnet_file[:drizzle])
		#uuid = procnet.generate_der(clnet_file[:categorize], clnet_file[:der])
		#uuid = procnet.generate_ier(clnet_file[:categorize], clnet_file[:ier])

            #catch e
            #end
            # Plotting:
            #
            #try
            #	pltid = CloudnetTools.Vis.show_classific(clnet_file[:categorize],
            #                                         SITENAME="NSA Utqiagvik",
            #                                         savefig=pltnet_file[:categorize]);
            #catch e
		# println("\e[1m\e[38;2;230;30;30;249m","* Error trying to visualize  data from $dd.$mm.$yy");
            #end
            # closure over date loops
        end
    end
end


# **** WICHTIG!! ******************            
# **** Template to user categorize directly without the storage of temporal nc files:
# utils = pyimport("cloudnetpy.utils")
# catcnet = pyimport("cloudnetpy.categorize")
# radar = catcnet.radar.Radar(radar_template_file)
# radar.data["Z"].data = Ze  or set!(radar.data, "Z", Ze)
# radar.data["v"].data = v, etc.
# lidar = catcnet.lidar.Lidar(lidar_template_file)
# 
# time, height = utils.time_grid(), radar.height
# wl=utils.get_wl_band(radar.radar_frequency)
#
# **************************************


#from cloudnetpy.plotting import generate_figure

#py"""
#input_files = {
#    'radar': $radar_file,
#    'lidar': $lidar_file,
#    'model': $model_file,
#    'mwr': $mwr_file
#}
#
#from cloudnetpy.categorize import generate_categorize
#uuid = generate_categorize(input_files, $categorize_file)
##
#from cloudnetpy.products import generate_classification
#uuid = generate_classification($categorize_file, $classific_file)
##
#from cloudnetpy.products import generate_iwc
#generate_iwc($categorize_file, $iwc_file)
##
#from cloudnetpy.products import generate_lwc
#generate_lwc($categorize_file, $lwc_file)
##
#from cloudnetpy.products import generate_drizzle
#generate_drizzle($categorize_file, $drizzle_file)
#
#### Visualization:
#
#from cloudnetpy.plotting import generate_figure
#generate_figure($categorize_file, ['Z','beta','lwp'], show=False, max_y=7, image_name=$zbetalwp_plotf)
#
#generate_figure($classific_file, ['target_classification', 'detection_status'], show=False, max_y=7, image_name=$classific_plotf)
#
#generate_figure($iwc_file, ['iwc', 'iwc_error', 'iwc_retrieval_status'], show=False, max_y=5, image_name=$iwc_plotf)
#
#generate_figure($lwc_file, ['lwc', 'lwc_error', 'lwc_retrieval_status'], show=False, max_y=5, image_name=$lwc_plotf)
#
#generate_figure($drizzle_file, ['Do', 'mu', 'S'], max_y=3, show=False, image_name=$drizzle_plotf)
#
#"""
#
## end of Script
