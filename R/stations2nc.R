#     Copyright (C) 2016 Santander Meteorology Group (http://www.meteo.unican.es)
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

#' @title Station data 2 netCDF export
#' @description Export a loadeR grid to NetCDF
#' @importFrom ncdf4 ncdim_def 
#' @importFrom ncdf4 ncvar_def 
#' @importFrom ncdf4 nc_create 
#' @importFrom ncdf4 ncatt_put 
#' @importFrom ncdf4 ncvar_put
#' @importFrom ncdf4 nc_close
#' @param data A grid data object coming from \code{\link[loadeR]{loadGridData}} or \code{\link[transformeR]{interpGrid}} 
#'  or the function \code{\link[loadeR.ECOMS]{loadECOMS}} of package \pkg{loadeR.ECOMS}. 
#' @param NetCDFOutFile Name of the file created by the function. (default to "out.nc4")
#' @param missval Missing value codification (default to 1e20)
#' @param globalAttributes Optional. A list of global attributes included in the NetCDF file. Same format as \code{varAttributes}.
#' @param varAttributes Optional. List of attributes to be included in the variable written in the NetCDF file. (NULL, default).
#' It has the format \code{list("name_of_attribute1" = attribute1, "name_of_attribute2" = attribute2} etc.
#' @param prec Precision to write the attribute. If not specified, the written precision is the same as the variable whose 
#' attribute this is. This can be overridden by specifying this argument with a value
#'  of \code{"short"}, \code{"float"}, \code{"double"}, or \code{"text"}.
#' @param compression If set to an integer between 1 (less compression) and 9 (highest compression), enables 
#' compression for the variable as it is written to the file. Turning compression on forces the 
#' output file to netcdf v4 format (the default), not compatible with older software 
#' that can only handle version 3 files. 
#' @param shuffle Logical. Turns on (\code{TRUE}) or off (\code{FALSE}, the default) the shuffle filter.
#'  According to netcdf documentation, turning the shuffle filter on can improve compression for integer variables.
#'   Turning the shuffle filter on forces the output file to netcdf v4 format (the default), 
#'   not compatible with older software that can only handle version 3 files. 
#' @param standardName Standard name of the variable.
#' @param verbose Optional. If set to \code{TRUE}, switches \code{\link[ncdf4]{nc_create}} to verbose mode
#' @return A NetCDF-4 file with the variable and attributes defined in the inputs.
#' @references
#' \itemize{
#' \item David Pierce \email{dpierce@@ucsd.edu}, Interface to Unidata netCDF (version 4 or earlier) format data files, http://dwpierce.com/software
#' }
#' @author S. Herrera, W. Franssen, J. Bedia and M. Iturbide
#' @export
#' @examples \dontrun{
#' require(transformeR)
#' data("VALUE_Iberia_tas")
#' # Name of output file:
#' fileName <- "tas_VALUE_DJF_Iberia.nc4"
#' # Including a global attribute:
#' globalAttributeList <- list("institution" = "Santander MetGroup, http://www.meteo.unican.es/")
#' # Including two variable attributes:
#' varAttributeList <- list(var_attr1 = "one_attribute", var_attr2 = "another_attribute")
#' # Create file:
#' stations2nc(data = VALUE_Iberia_tas,
#'         NetCDFOutFile = fileName,
#'         missval = 1e20,
#'         prec = "float",
#'         globalAttributes = globalAttributeList,
#'         varAttributes = varAttributeList)
#' }

stations2nc <- function(data,
                    NetCDFOutFile = "out.nc4",
                    missval = 1e20,
                    globalAttributes = NULL,
                    varAttributes = NULL,
                    prec = "float",
                    compression = 4,
                    shuffle = FALSE,
                    standardName = NULL,
                    verbose = FALSE) {
	tmpStdName <- data$Variable$varName
	tmpUnits <- attributes(data$Variable)$"units"
	if (is.null(tmpUnits)) {tmpUnits <- ""}
	time.index <- grep("^time$", attr(data$Data, "dimensions"))
	loc.index <- grep("^loc$", attr(data$Data, "dimensions"))
	member.index <- grep("^member$", attr(data$Data, "dimensions"))
	datesList <- as.POSIXct(data$Dates$start, tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
	if (all(is.na(datesList))) datesList <- as.POSIXct(data$Dates$start, tz = "GMT")
	times <- (as.double(datesList) - as.double(datesList[1])) / 86400
	dimtime <- ncdim_def("time", paste("days since", data$Dates$start[1]), times, unlim = FALSE, calendar = "gregorian", create_dimvar = TRUE)
	dimSta  <- ncdim_def("station", "", 1:length(data$Metadata$station_id), create_dimvar = FALSE)
	dimnchar <- ncdim_def("id_strlen", "", 1:max(nchar(data$Metadata$station_id)), create_dimvar=FALSE )
	nName <- suppressWarnings(tryCatch({ncdim_def("name_strlen", "", 1:max(nchar(data$Metadata$name)), create_dimvar=FALSE )}, error = function(e){NA}))
	nProj <- ncdim_def("projection_strlen", "", 1:1, create_dimvar=FALSE )
	if (length(member.index) > 0) {
		dimens  <- ncdim_def("member", units = "member", 0:(dim(data$Data)[member.index] - 1), longname = "realization", create_dimvar = TRUE)
		perOrdered <- c(loc.index, member.index, time.index)
		dimOrdered <- list(dimSta, dimens, dimtime)
	} else {
		perOrdered <- c(loc.index,time.index)
		dimOrdered <- list(dimSta,dimtime)
	}
	dataOrdered <- aperm(data$Data, perOrdered)
	varLon <- ncvar_def("lon", "degrees_east" , dim = list(dimSta), longname = "station longitude", prec = "double")
	varLat <- ncvar_def("lat", "degrees_north", dim = list(dimSta), longname = "station latitude",  prec = "double")
	varProj <- ncvar_def("projection", "", dim = list(nProj), prec = "char")
	varStation <- ncvar_def("station_id", "", dim = list(dimnchar, dimSta), longname = "station identifier", prec= "char")
	var <- ncvar_def(data$Variable$varName, units = tmpUnits, dim = dimOrdered, missval, compression = compression, shuffle = shuffle, prec = prec) 
	varHeight <- NULL
	varName <- NULL
	if(!is.null(data[["Metadata"]][["altitude"]])) varHeight <- ncvar_def("alt", "m", dim = list(dimSta), missval=NA, longname = "height", prec= "double")
	if(!is.null(data[["Metadata"]][["name"]])) varName <- ncvar_def("station_name", "", dim = list(nName, dimSta), longname = "station name", prec= "char")
	vars <- list(varLat, varLon, varStation, varName, varHeight, varProj, var) ##  ... , force_v4=FALSE
	vars <- vars[unlist(lapply(vars, function(e) !is.null(e)))]
	ncnew <- nc_create(NetCDFOutFile, vars, verbose = verbose)
	ncatt_put(ncnew, data$Variable$varName, "standard_name", standardName)
	ncatt_put(ncnew, "time", "standard_name","time")
	ncatt_put(ncnew, "time", "axis","T")
	ncatt_put(ncnew, "time", "_CoordinateAxisType","Time")
	ncatt_put(ncnew, "time", "_ChunkSize",1)
	ncatt_put(ncnew, "station_id", "cf_role","timeseries_id")
	ncatt_put(ncnew, "projection", "EPSG_code","EPSG:4326")
	ncatt_put(ncnew, "lon", "standard_name","longitude")
	ncatt_put(ncnew, "lat", "standard_name","latitude")
	if(!is.null(data[["Metadata"]][["altitude"]])) ncatt_put(ncnew, "alt", "standard_name","altitude")
	if(!is.null(data[["Metadata"]][["altitude"]])) ncatt_put(ncnew, "alt", "missing_value", NA, prec = "double")
	if (length(member.index) > 0) {
		ncatt_put(ncnew, "member", "standard_name","realization")
		ncatt_put(ncnew, "member", "_CoordinateAxisType","Ensemble")
		ncatt_put(ncnew, "member", "ref","http://www.uncertml.org/samples/realisation")
	}
	ncatt_put(ncnew, data$Variable$varName, "missing_value", missval, prec = prec)
	if (!is.null(varAttributes)) {
		sapply(1:length(varAttributes), function(x) ncatt_put(ncnew, var$name, names(varAttributes)[x], as.character(varAttributes[[x]])))
	}
	ncatt_put(ncnew, var$name, "description", attr(data$Variable, "description"))
	ncatt_put(ncnew, var$name, "longname", attr(data$Variable, "longname"))
	ncatt_put(ncnew, var$name, "coordinates", "lat lon")
	z <- attributes(data$Variable$level)
	if (!is.null(z)) ncatt_put(ncnew, var$name, "level", z)
	if (!is.null(globalAttributes)) {      
		sapply(1:length(globalAttributes), function(x) ncatt_put(ncnew, var$name, names(globalAttributes)[x], as.character(globalAttributes[[x]])))
	}
	# Bias-corrected products
	if (length(attr(data$Variable, "correction")) > 0) {
		ncatt_put(ncnew, 0, "product", "Bias-Correction")
		ncatt_put(ncnew, 0, "bc_method", attr(data$Variable, "correction"))
	}
	# Downscaled products
	if (length(attr(data$Variable, "downscaling:method")) > 0) {
		ncatt_put(ncnew, 0, "product", "Downscaling")
		ncatt_put(ncnew, 0, "downscaling_method", attr(data$Variable, "downscaling:method"))
	}
	if (length(attr(data$Variable, "dataset")) > 0) {
	  ncatt_put(ncnew, 0, "dataset", attr(data$Variable, "dataset"))
	}
	if (length(attr(data, "source")) > 0) {
		ncatt_put(ncnew, 0, "source", attr(data, "source"))
	}
	ncatt_put(ncnew, 0, "Origin", "NetCDF file created by loadeR.2nc: https://github.com/SantanderMetGroup/loadeR.2nc")
	ncatt_put(ncnew, 0, "Conventions", "CF-1.4")
	ncatt_put(ncnew, 0, "coordinates", "projection")
	ncatt_put(ncnew, 0, "featureType", "timeSeries")
	ncvar_put(ncnew, var, dataOrdered)
	ncvar_put(ncnew, varStation, data$Metadata$station_id)
	if(!is.null(data[["Metadata"]][["name"]])) ncvar_put(ncnew, varName, data$Metadata$name)
	if(!is.null(data[["Metadata"]][["altitude"]])) ncvar_put(ncnew, varHeight, data$Metadata$altitude)
	ncvar_put(ncnew, varLon, data$xyCoords$x)
	ncvar_put(ncnew, varLat, data$xyCoords$y)
	nc_close(ncnew)
	message("[", Sys.time(), "] NetCDF file written in: ", NetCDFOutFile)
}

#end
