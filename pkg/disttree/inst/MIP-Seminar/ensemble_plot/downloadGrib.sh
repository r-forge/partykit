# -------------------------------------------------------------------
# - NAME:        downloadGrib.sh
# - AUTHOR:      Reto Stauffer
# - DATE:        2018-03-13
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2018-03-13, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2018-03-13 19:33 on marvin
# -------------------------------------------------------------------


# Downloading, filtering, and subsetting some data.
ftpbase="http://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod"
ftppath="gefs.%Y%m%d/%H/pgrb2/gep<member>.t%Hz.pgrb2f<step>"
outtemplate="RAW/GFS_%Y%m%d%H_<step>_<member>.grb2"
nctemplate="NetCDF/GFS_%Y%m%d%H_<step>_<member>.nc"
if [ ! -d RAW ] ; then mkdir RAW ; fi
if [ ! -d NetCDF ] ; then mkdir NetCDF ; fi

member=1
step=36  # 36h-forecast
step=240 # 156h-forecast
date="20180307"

# Subset specification
lonW=-40
lonE=60
latS=30
latN=60

# No filtering the grib data using grib_filter
# This here is the filter file with the rules.
filterfile=".filter"
cat > $filterfile <<EOF
# Only surface temperature and precipitation
switch(shortName) {
   case "tp":
      print "found [shortName]";
      write;
   case "2t":
      print "found [shortName]";
      write;
   default:
}
EOF

while [ $member -le 20 ] ; do
    # Verbose output
    printf " * Processing member %d\n" ${member}

    # String (two digit member)
    memstr=`printf "%02d" ${member}`
    stepstr=`printf "%02d" ${step}`

    # Create URL
    tmp=`date --date "${date}" "+${ftppath}" | \
         sed "s/<member>/"${memstr}"/g" | \
         sed "s/<step>/"${step}"/g"`

    # Output file
    outfile=`date --date "${date}" "+${outtemplate}" | \
         sed "s/<member>/"${memstr}"/g" | \
         sed "s/<step>/"${step}"/g"`

    # Downloading grib file if required
    if [ ! -f $outfile ] ; then
        url=`printf "%s/%s" "${ftpbase}" "${tmp}"`
        printf "   Downloading %s\n" ${url}
        printf "   Save data as %s\n" ${outfile}

        wget ${url} -O ${outfile}
    fi

    # Output file
    ncfile=`date --date "${date}" "+${nctemplate}" | \
         sed "s/<member>/"${memstr}"/g" | \
         sed "s/<step>/"${step}"/g"`

    # Create NetCDF file if required
    if [ ! -f ${ncfile} ] ; then

        # Create netcdf files
        sudo grib_filter ${filterfile} ${outfile} -o tmp.grb2
        sudo wgrib2 tmp.grb2 -small_grib ${lonW}:${lonE} ${latS}:${latN} small.grb2 && \
            sudo grib_to_netcdf small.grb2 -o ${ncfile} && sudo rm small.grb2 tmp.grb2 && \
            sudo chown retos:retos ${ncfile}

    fi

    # Increase member number
    let member=$member+1
done













