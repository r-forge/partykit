# -------------------------------------------------------------------
# - NAME:        Spatial.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2018-03-13
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2018-03-13, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2018-03-14 12:46 on marvin
# -------------------------------------------------------------------


    library( "raster" )
    library( "ncdf4" )


    step <- 36
    step <- 240

    ncfiles <- list.files("NetCDF",sprintf("^GFS_[0-9]{10}_%02d_[0-9]{2}.nc$",step))
    if ( length(ncfiles)==0 ) stop("No netcdf files found")
    ncfiles <- sprintf("NetCDF/%s",sort(ncfiles))


    # Helper function (matrix + lon/lat vectors -> raster)
    data2raster <- function( x, lons, lats, flip = TRUE ) {
        delta <- min(diff(lons)) / 2
        x <- t(x); if ( flip ) x <- x[nrow(x):1,]
        r <- raster( x, xmn = min(lons)-delta, xmx = max(lons)+delta,
                        ymn = min(lats)-delta, ymx = max(lats)+delta )
        r
    }

    res <- list(t2m=list(),tp=list())
    for ( file in ncfiles ) {
        # Open NetCDF file
        nc <- nc_open( file )
        # Extract member from file name
        mem <- as.integer(gsub("\\.nc","",regmatches(file,regexpr("[0-9]{2}.nc",file))))
        # Loading longitude/latitude arrays
        lats <- ncvar_get(nc,"latitude")
        lons <- ncvar_get(nc,"longitude")
        # Extracting data
        tmp_t2m <- ncvar_get( nc, "t2m" ) - 273.15
        tmp_tp  <- ncvar_get( nc, "tp" )
        # Convert to raster and store
        res$t2m[[sprintf("member_%d",mem)]] <- data2raster( tmp_t2m, lons, lats )
        res$tp[[sprintf("member_%d",mem)]]  <- data2raster( tmp_tp,  lons, lats )
        # Close NetCDF connection
        nc_close( nc )
    }


    t2m <- brick( res$t2m )
    tp  <- brick( res$tp  )

    t2m <- disaggregate( t2m, fact=5, method="bilinear" )
    tp  <- disaggregate( tp,  fact=5, method="bilinear" )

    zlim_t2m <- range(cellStats(t2m,range))
    zlim_tp <- range(cellStats(tp,range))


    library( "colorspace" )
    plot( t2m, zlim=zlim_t2m, col=diverge_hcl(51) )
    plot( tp, zlim=zlim_tp,   col=rev(sequential_hcl(51)) )



    # Orography
    nc <- nc_open( "NetCDF/orog.nc" )

    lats <- ncvar_get(nc,"latitude")
    lons <- ncvar_get(nc,"longitude")

    # Pick orography
    orog    <- data2raster( ncvar_get( nc, "orog" ), lons, lats, FALSE )
    # Convert to raster and store

    library( "sp" )
    library( "maps" )

    innsbruck <- SpatialPoints( data.frame( lon = 11.38, lat = 47.3 ) )

    plot( crop(orog,extent( 0, 30, 40, 55 )), col = colorspace::terrain_hcl(51),
         main = expression(paste("GFS/GEFS Model Topography 1",degree,"x1",degree)))
    map( add = TRUE )

    # Adding Innsbruck
    innsbruck_height <- extract( orog, innsbruck, method = "bilinear" )
    points( innsbruck, pch = 19, cex = 3 )
    text( innsbruck, sprintf("%d",round(innsbruck_height)), adj=c(-.5,.5), cex = 2 )


    nc_close( nc )






