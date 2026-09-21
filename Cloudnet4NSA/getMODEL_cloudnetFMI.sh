#!/usr/bin/csh
# ------------------------------------------------------
# Script to retrieve ECMWF model files from the API
# cloudnet.fmi.fi API site for ARM NSA 
# ------------------------------------------------------
#
# (c) 2020 Pablo Saavedra Garfias
# pablo.saavedra@uni-leipzig.de
# (AC)3 Project, B07
# LIM
# University of Leipzig
# ------------------------------------------------------
#

#set WEB_SITE='http://devcloudnet.fmi.fi/cnet'
#https://cloudnet.fmi.fi/api/download/product/4edf4172-3633-46bc-952b-b263684a91da/20121226_arm-nsa_ecmwf.nc
#https://cloudnet.fmi.fi/api/download/product/1ce9fc63-eafd-4721-9fa9-26eafc604e3c/20121225_arm-nsa_ecmwf.nc
## .... 
set WEB_SITE='https://cloudnet.fmi.fi/api/'
set SITE='arm-nsa'
set TYPE='calibrated'
set LOG='ecmwfstatus.log'
set MODEL='ecmwf' # 'icon-iglo-108-119'
set OUT_PATH=/projekt2/ac3data/B07-data/utqiagvik-nsa/CloudNet/input/ECMWF
#/home/psgarfias/LIM/data/arctic-mosaic/ICONG
#set OUT_PATH=/home/psgarfias/LIM/data/utqiagvik-nsa/ECMWF
set TEMP_PATH=/tmp

echo '================' >& ${TEMP_PATH}/${LOG}
/bin/date >>& ${TEMP_PATH}/${LOG}
echo '================' >>& ${TEMP_PATH}/${LOG}

set years=(2026)
set months=(1 2 3 4)

setenv LC_TIME en_EN.utf8
foreach year ($years)
    foreach month ($months)
	set day=1
	while ( $day <= 31 )
	    set strday=`date --date="$year-$month-$day" +%d`
	    set strmonth=`date --date="$year-$month-$day" +%m`

	    if ("${strday}" == "") then
		break
	    endif
	    #> curl "https://cloudnet.fmi.fi/api/model-files?site=arm-nsa&date=2012-12-10&model=ecmwf" | jq '.[]["downloadUrl"]' | xargs -n1 curl -O
	    
	    set filename=${year}${strmonth}${strday}_${SITE}_${MODEL}.nc

	    echo 'curl "'${WEB_SITE}'model-files?site='${SITE}'&date='${year}'-'${strmonth}'-'${strday}'&model='${MODEL}'" | jq '\''.[]["downloadUrl"]'\'' | xargs -n1 curl -O' >& fmi.in #
	set outfile=${OUT_PATH}/${year}

	    echo 'mv *.nc '${outfile} >>& fmi.in
	    	    # calling wget to retrieve the file
	    set cmd=(sh fmi.in)
	    ${cmd}
	
	    @ day++
	end
	
	
	#echo ${fullpath}
	
    end
end

/bin/date >>& ${TEMP_PATH}/${LOG}
echo '================' >>& ${TEMP_PATH}/${LOG}

exit 0

#set strmm='xxx' #`/usr/bin/date --date="$year-$month-$day" +%b | awk '{print tolower($0)}'`
#set strmm=`date --date="$year-$month-$day" +%b 2>/dev/null | awk '{print tolower($0)}'`
#wget ${fullpath} -a ${TEMP_PATH}/seaicestatus.log &&  /home/psgarfias/LIM/scripts/SeaIce/h4toh5 ${filename}

#echo ${webaddress}'?getfile' >& nclink.in
#set outfile=${OUT_PATH}${ncfile}

#set cmd=(/usr/bin/wget --user=superuser --password=15k-Nimbus --spider -i nclink.in -O ${outfile} -a ${TEMP_PATH}/ceilimstatus.log) 

#${cmd}

# End of script
