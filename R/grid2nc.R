#     grid2nc.R Export climate4R grids to netCDF
#
#     Copyright (C) 2021 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Grid 2 netCDF export
#' @description Export a loadeR grid to NetCDF
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncatt_put ncvar_put nc_close
#' @param data A climate4R grid data object (\url{http://www.meteo.unican.es/climate4R})
#' @param NetCDFOutFile Name of the file created by the function. (default to \code{"out.nc4"})
#' @param missval Missing value codification (default to \code{1e20})
#' @param globalAttributes Optional. A list of global attributes included in the NetCDF file. Same format as \code{varAttributes}.
#' @param varAttributes Optional. List of attributes to be included in the variable written in the NetCDF file. Default to \code{NULL}.
#' It has the format \code{list("name_of_attribute1" = list(attribute1, prec of attribute1), "name_of_attribute2" = list(attribute1, prec of attribute1)} etc.
#' @param prec Precision to write the attribute. If not specified, the written precision is given by the variable's
#' corresponding attribute. This can be overridden by specifying the following argument values:
#'  \code{"short"}, \code{"float"}, \code{"double"}, or \code{"text"}.
#' @param compression Integer value between 1 (less compression) and 9 (highest compression), enabling 
#' compression for the variable as it is written to the file. Note that turning compression on forces the 
#' output file to netcdf v4 format (the default), not compatible with older software 
#' that can only handle version 3 files. 
#' @param shuffle Logical. Turns on (\code{TRUE}) or off (\code{FALSE}, the default) the shuffle filter.
#'  According to netcdf documentation, turning the shuffle filter on can improve compression for integer variables.
#'   Note that turning the shuffle filter on forces the output file to netcdf v4 format (the default), 
#'   not compatible with older software that can only handle version 3 files. 
#' @param verbose Optional. If set to \code{TRUE}, switches \code{\link[ncdf4]{nc_create}} to verbose mode
#' @return A NetCDF-4 file with the variable and attributes defined in the inputs.
#' @references
#' \itemize{
#' \item David Pierce \email{dpierce@@ucsd.edu}, Interface to Unidata netCDF (version 4 or earlier) format data files, http://dwpierce.com/software
#' }
#' @author S. Herrera, W. Franssen and J. Bedia
#' @export
#' @examples \dontrun{
#' library(loadeR.2nc)
#' data(tx) # A climate4R grid
#' # Name of output file:
#' fileName <- "/tmp/tasmax_WFDEI_JJA_W2001_2010.nc"
#' # Including a global attribute:
#' globalAttributeList <- list("title" = "Daily mean maximum temperature data from the WFDEI Database over Southern Europe used as an example of netCDF export from climate4R",
#'                             "institution" = "SantanderMetGroup, <http://www.meteo.unican.es>", 
#'                             "source" = "climate4R <https://github.com/SantanderMetGroup/climate4R>")
#' # Including variable attributes. The last one not really correct, but for illustration purpose:
#' varAttributeList <- list("name" = list("tasmax", prec = "text"),
#'                          "description" = list("Daily maximum near surface air temperature", prec = "text"),
#'                          "level" = list(2L, prec = "short"))
#' # Create file:
#' grid2nc(data = tx,
#'         NetCDFOutFile = fileName,
#'         missval = 1e20,
#'         prec = "float",
#'         globalAttributes = globalAttributeList,
#'         varAttributes = varAttributeList)
#' }


grid2nc <- function(data,
                    NetCDFOutFile = "out.nc4",
                    missval = 1e20,
                    globalAttributes = NULL,
                    varAttributes = NULL,
                    prec = NULL,
                    compression = 4,
                    shuffle = FALSE,
                    verbose = FALSE) {
      
      prec <- match.arg(prec, choices =  c("float", "short", "double", "text"))
      
      # Global attribute defs
      if (is.null(varAttributes[["name"]])) {
            tmpStdName <- data$Variable$varName
      } else {
            message("[", Sys.time(), "] The original \'name\' attribute was overriden by the value specified in the variable attr list")
            tmpStdName <- varAttributes[["name"]][[1]]
      }
      if (is.null(varAttributes[["units"]])) {
            tmpUnits <- attributes(data$Variable)$"units"   
      } else {
            message("[", Sys.time(), "] The original \'units\' attribute was overriden by the value specified in the variable attr list")
            tmpUnits <- varAttributes[["units"]][[1]]
      }
      if (is.null(varAttributes[["long_name"]])) {
            tmpLongName <- attr(data$Variable, "longname")   
      } else {
            message("[", Sys.time(), "] The original \'longname\' attribute was overriden by the value specified in the variable attr list")
            tmpLongName <- varAttributes[["long_name"]][[1]]
      }
      if (is.null(tmpLongName)) tmpLongName <- tmpStdName ## Assign shortname if longname is missing
      
      time.index <- grep("^time$", attr(data$Data, "dimensions"))
      lon.index <- grep("^lon$", attr(data$Data, "dimensions"))
      lat.index <- grep("^lat$", attr(data$Data, "dimensions"))
      member.index <- grep("^member$", attr(data$Data, "dimensions"))
      datesList <- as.POSIXct(data$Dates$start, tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
      if (all(is.na(datesList))) datesList <- as.POSIXct(data$Dates$start, tz = "GMT")
      times <- (as.double(datesList) - as.double(datesList[1])) / 86400
      dimtime <- ncdim_def("time", paste("days since", data$Dates$start[1]), times, unlim = FALSE, calendar = "gregorian", create_dimvar = TRUE)
      dimlon  <- ncdim_def("lon", units = "degrees_east", data$xyCoords$x, longname = "longitude", create_dimvar = TRUE)
      dimlat  <- ncdim_def("lat", units = "degrees_north", data$xyCoords$y, longname = "latitude", create_dimvar = TRUE)
      if (length(member.index) > 0) {
            dimens  <- ncdim_def("member", units = "member", 0:(dim(data$Data)[member.index] - 1), longname = "realization", create_dimvar = TRUE)
            perOrdered <- c(lon.index, lat.index, member.index, time.index)
            dimOrdered <- list(dimlon, dimlat, dimens, dimtime)
      } else {
            perOrdered <- c(lon.index,lat.index,time.index)
            dimOrdered <- list(dimlon,dimlat,dimtime)
      }
      dataOrdered <- aperm(data$Data, perOrdered)
      var <- ncvar_def(name = tmpStdName,
                       units = tmpUnits,
                       dim = dimOrdered,
                       missval = missval,
                       longname = tmpLongName,
                       prec = prec,
                       compression = compression,
                       shuffle = shuffle)
      ncnew <- nc_create(NetCDFOutFile, var, verbose = verbose)
      on.exit(nc_close(ncnew))
      ncatt_put(ncnew, "time", "standard_name","time")
      ncatt_put(ncnew, "time", "axis","T")
      ncatt_put(ncnew, "time", "_CoordinateAxisType","Time")
      ncatt_put(ncnew, "time", "_ChunkSize",1)
      ncatt_put(ncnew, "lon", "standard_name","longitude")
      ncatt_put(ncnew, "lon", "_CoordinateAxisType","Lon")
      ncatt_put(ncnew, "lat", "standard_name","latitude")
      ncatt_put(ncnew, "lat", "_CoordinateAxisType","Lat")
      if (length(member.index) > 0) {
            ncatt_put(ncnew, "member", "standard_name","realization")
            ncatt_put(ncnew, "member", "_CoordinateAxisType","Ensemble")
            ncatt_put(ncnew, "member", "ref","http://www.uncertml.org/samples/realisation")
      }
      ncatt_put(ncnew, tmpStdName, "missing_value", missval)
      if (!is.null(varAttributes)) {
            sapply(1:length(varAttributes), function(x) {
                  prec1 <- if (is.null(varAttributes[[x]][["prec"]])) {
                        NA
                  } else {
                        varAttributes[[x]][["prec"]]
                  }
                  ncatt_put(nc = ncnew, 
                            varid = var$name,
                            attname =  names(varAttributes)[x],
                            attval = as.character(varAttributes[[x]]),
                            prec = prec1)
            })
      }
      if (is.null(varAttributes[["description"]])) ncatt_put(ncnew, var$name, "description", attributes(data$Variable)$"description")
      if (is.null(varAttributes[["long_name"]])) ncatt_put(ncnew, var$name, "long_name", attributes(data$Variable)$"longname")
      z <- attributes(data$Variable$level)
      if (!is.null(z)) ncatt_put(ncnew, var$name, "level", z)
      if (!is.null(globalAttributes)) {      
            sapply(1:length(globalAttributes), function(x) ncatt_put(ncnew, 0, names(globalAttributes)[x], as.character(globalAttributes[[x]])))
      }
      # Bias-corrected products
      if (length(attr(data$Data, "correction")) > 0) {
            ncatt_put(ncnew, 0, "product", "Bias-Correction")
            ncatt_put(ncnew, 0, "bc_method", attr(data$Data, "correction"))
      }
      # Downscaled products
      if (length(attr(data$Data, "downscaling:method")) > 0) {
            ncatt_put(ncnew, 0, "product", "Downscaling")
            ncatt_put(ncnew, 0, "downscaling_method", attr(data$Data, "downscaling:method"))
      }
      if (length(attr(data, "dataset")) > 0) {
            ncatt_put(ncnew, 0, "dataset", attr(data, "dataset"))
      }
      if (length(attr(data, "source")) > 0) {
            ncatt_put(ncnew, 0, "source", attr(data, "source"))
      }
      ncatt_put(ncnew, 0, "Origin", "NetCDF file created by loadeR.2nc <https://github.com/SantanderMetGroup/loadeR.2nc>")
      ncatt_put(ncnew, 0, "Conventions", "CF-1.4")
      ncvar_put(ncnew, var, dataOrdered)
      message("[", Sys.time(), "] NetCDF file written in: ", NetCDFOutFile)
}



