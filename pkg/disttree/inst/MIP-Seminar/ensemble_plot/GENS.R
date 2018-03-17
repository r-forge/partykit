# -------------------------------------------------------------------
# - NAME:        GENS.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2018-03-13
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2018-03-13, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2018-03-13 17:43 on marvin
# -------------------------------------------------------------------


    library( "zoo" )

    Sys.setenv("TZ"="UTC")

# -------------------------------------------------------------------
# GFS ensemble forecast
# -------------------------------------------------------------------

    GFS   <- read.table( "GENS_00_innsbruck-flughafen.txt", header = TRUE )
    init  <- min(as.POSIXct(GFS$timestamp-GFS$step*3600,origin="1970-01-01"))
    title <- sprintf("%s\n%s",
                "Global Forecast System (GFS) Ensemble Forecast",
                sprintf("Forecast initialized %s",strftime(init,"%Y-%m-%d %H:%M UTC")))

    getZoo <- function( x, variable ) {
        x <- subset( x, varname == variable )
        if ( nrow(x) == 0 ) stop("Ups, variable seems not to exist at all!")
        # Else create zoo
        x <- zoo( subset(x,select=-c(varname,timestamp,step)),
                  as.POSIXct(x$timestamp,origin="1970-01-01") )
        x
    }
    GFS_t2m  <- getZoo( GFS, "tmp2m" )
    GFS_rain <- getZoo( GFS, "apcpsfc" )

    # Two POSIXct vectors for the axis
    main  <- seq(min(as.POSIXct(as.Date(index(GFS_t2m)))),max(index(GFS_t2m)),by=86400)
    minor <- seq(min(index(GFS_t2m)),max(index(GFS_t2m)),by=3*3600)

    par(mar=c(0.2,0,0.2,0),oma=c(5,5,5,5))
    layout( matrix(1:2,ncol=1) )
    plot( GFS_t2m, screen=1, xaxs="i", xaxt="n" )
    abline( v = main, lty = 2 )
    mtext( side = 3, line = 2, cex = 1.2, font = 2, title )
    plot( GFS_rain, screen=1, xaxs="i", xaxt="n", yaxs="i", ylim=range(GFS_rain)*c(0,1.05) )
    abline( v = main, lty = 2 )

    axis( side = 1, at = minor, strftime(minor,"%H") )
    axis( side = 1, at = main + 42300,  strftime(main,"%b %d"),
          line = 2, col.ticks=NA, col.axis="gray30", col = NA )



# -------------------------------------------------------------------
# ECMWF Example Forecast
# -------------------------------------------------------------------





