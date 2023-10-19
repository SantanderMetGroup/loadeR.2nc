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

#' @title climate4R grid 2 netCDF export
#' @description Export a climate4R grid to NetCDF
#' @param data A climate4R grid data object (\url{http://www.meteo.unican.es/climate4R})
#' @param NetCDFOutFile Name of the file created by the function. (default to \code{"out.nc4"})
#' @param missval Missing value codification (default to \code{1e20})
#' @param globalAttributes Optional. A list of global attributes included in the NetCDF file. Same format as \code{varAttributes}.
#' @param varAttributes Optional. List of attributes to be included in the variable written in the NetCDF file. Default to \code{NULL}.
#' It has the format \code{list("name_of_attribute1" = attribute1, "name_of_attribute2" = attribute2} etc.
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
#' @param gridNorthPole Vector with the longitude and latitude. Default is c("39.25","-162.0")
#' @param coordBounds Default is NULL.
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncatt_put ncvar_put nc_close
#' @import transformeR
#' @details 
#' 
#' \strong{netCDF attributes}
#' Some attributes are automatically inferred by the function from the climate4R input grid metadata.
#' These can be replaced by user-defined attributes via the global and/or variable attribute lists.
#' It is highly advised that an authoritative attribute convention is followed, such as CF and/or ESIP:
#' \itemize{
#' \item \url{https://cfconventions.org/}
#' \item \url{https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#creator_url}
#' }
#' This conventions are indicated by default in the global attributes of the output. 
#' 
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
#' fileName <- "tasmax_WFDEI_JJA_W2001_2010.nc4"
#' # Including a global attribute:
#' globalAttributeList <- list("institution" = "SantanderMetGroup, http://www.meteo.unican.es/")
#' # Including two variable attributes:
#' varAttributeList <- list(var_attr1 = "one_attribute", var_attr2 = "another_attribute")
#' # Create file:
#' grid2nc(data = tx,
#'         NetCDFOutFile = fileName,
#'         missval = 1e20,
#'         prec = "float",
#'         globalAttributes = globalAttributeList,
#'         varAttributes = varAttributeList)
#' # Inspection fo the file: Requires nectcdf-bin installed in your system:
#' system("ncdump -h tasmax_WFDEI_JJA_W2001_2010.nc4")
#' }

grid2nc <- function(data,
                    NetCDFOutFile = "out.nc4",
                    missval = 1e20,
                    globalAttributes = NULL,
                    varAttributes = NULL,
                    prec = "float",
                    compression = 4,
                    shuffle = FALSE,
                    verbose = FALSE,
                    gridNorthPole = c("39.25","-162.0"),
                    coordBounds = NULL) {
   
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
      message("[", Sys.time(), "] The original \'long_name\' attribute was overriden by the value specified in the variable attr list")
      tmpLongName <- varAttributes[["long_name"]][[1]]
   }
   if (is.null(tmpLongName)) tmpLongName <- tmpStdName ## Assign shortname if longname is missing
   
   time.index <- grep("^time$", attr(data$Data, "dimensions"))
   lon.index <- grep("^lon$", attr(data$Data, "dimensions"))
   lat.index <- grep("^lat$", attr(data$Data, "dimensions"))
   member.index <- grep("^member$", attr(data$Data, "dimensions"))
   var.index <- grep("^var$", attr(data$Data, "dimensions"))
   startList <- transformeR::getRefDates(data, "start")
   datesList <- as.POSIXct(startList, tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
   if (all(is.na(datesList))) datesList <- as.POSIXct(startList, tz = "GMT")
   if ((length(grep(":00:00 GMT",startList)) == 0) & ((length(grep("GMT",startList)) != 0))) {
      startList <- gsub(startList, pattern = "00:00:00 GMT", replacement = "GMT")
   } 
   daysincestring <- if (grepl("\\d{2}-\\d{2} GMT", startList[1])) { # A not very elegant arrangement for netcdf-java API compliance
       gsub(pattern = " GMT", replacement = "T00:00:00", startList[1])        
   } else {
       startList[1]
   }
   times <- (as.double(datesList) - as.double(datesList[1])) / 86400
   dimtime <- ncdim_def("time", paste("days since", daysincestring), times, unlim = FALSE,
                        calendar = "gregorian", create_dimvar = TRUE)
   if (!is.null(attr(data$xyCoords, "projection")) & attr(data$xyCoords, "projection") == "RotatedPole") {
      dimlon  <- ncdim_def("rlon", units = "degrees", data$xyCoords$x,
                           longname = "longitude in rotated pole grid", create_dimvar = TRUE)
      dimlat  <- ncdim_def("rlat", units = "degrees", data$xyCoords$y,
                           longname = "latitude in rotated pole grid", create_dimvar = TRUE)
   } else if (!is.null(attr(data$xyCoords, "projection")) & (length(grep("LambertConformal",
                                                                         attr(data$xyCoords, "projection"),
                                                                         ignore.case = TRUE)) > 0)) {
      dimlon  <- ncdim_def("x", units = "m", data$xyCoords$x,
                           longname = "x coordinate of projection", create_dimvar = TRUE)
      dimlat  <- ncdim_def("y", units = "m", data$xyCoords$y,
                           longname = "y coordinate of projection", create_dimvar = TRUE)
   } else {
      dimlon  <- ncdim_def("lon", units = "degrees_east", data$xyCoords$x,
                           longname = "longitude", create_dimvar = TRUE)
      dimlat  <- ncdim_def("lat", units = "degrees_north", data$xyCoords$y,
                           longname = "latitude", create_dimvar = TRUE)
   }
   if (length(member.index) > 0) {
      dimens  <- ncdim_def("member", units = "", 1:(dim(data$Data)[member.index]), create_dimvar = FALSE)
      dimnchar <- ncdim_def("nchar", "", 1:max(nchar(data$Members)), create_dimvar = FALSE)
      perOrdered <- c(lon.index, lat.index, member.index, time.index)
      dimOrdered <- list(dimlon, dimlat, dimens, dimtime)
   } else {
      perOrdered <- c(lon.index, lat.index,time.index)
      dimOrdered <- list(dimlon, dimlat,dimtime)
   }
   if (transformeR::isMultigrid(data)) {
      if (var.index < min(perOrdered)) perOrdered <- perOrdered - 1
   }
   if (!is.null(coordBounds)) {
      dimBounds  <- ncdim_def("vertices", units = "", c(1:4), create_dimvar = FALSE)
   }
   if (transformeR::isMultigrid(data)) {
      dataOrdered <- lapply(1:length(tmpStdName), function(v){
         data.var <- transformeR::subsetGrid(data, var = tmpStdName[v])
         data.var <- aperm(data.var$Data, perOrdered)
      })
      var <- lapply(1:length(tmpStdName), function(v) {
         var.var <- ncvar_def(tmpStdName[v], units = tmpUnits[v], dim = dimOrdered,
                              missval = missval, longname = tmpStdName[v],
                              compression = compression,
                              shuffle = shuffle, prec = prec)
      })
   } else {
      dataOrdered <- aperm(data$Data, perOrdered)
      var <- ncvar_def(data$Variable$varName, units = tmpUnits, dim = dimOrdered,
                       missval = missval, longname = tmpStdName,
                       compression = compression, shuffle = shuffle, prec = prec)
   }
   if (!is.null(attr(data$xyCoords, "projection")) & attr(data$xyCoords, "projection") == "RotatedPole"){
      varProj <- ncvar_def("rotated_pole", units = "", dim = list(), prec = "char")
      varLon <- ncvar_def("lon", units = "degrees_east", dim = list(dimlon, dimlat),
                          longname = "longitude", prec = "double")
      varLat <- ncvar_def("lat", units = "degrees_north", dim = list(dimlon, dimlat),
                          longname = "latitude", prec = "double")
      if (!is.null(coordBounds)) {
         varLonBounds <- ncvar_def("lon_vertices", units = "degrees_east", dim = list(dimBounds, dimlon, dimlat),
                                   longname = "longitude", prec = "double")
         varLatBounds <- ncvar_def("lat_vertices", units = "degrees_north", dim = list(dimBounds, dimlon ,dimlat),
                                   longname = "latitude", prec = "double")
      }
      if (!transformeR::isMultigrid(data)) {
         var <- list(var, varLon, varLat, varProj)
      } else {
         var[[length(var) + 1]] <- varLon
         var[[length(var) + 1]] <- varLat
         var[[length(var) + 1]] <- varProj
      }
      if (length(member.index) > 0) {
         if (is.character(data$Members)) {
            varMem <- ncvar_def("member", units = "", dim = list(dimnchar,dimens), prec = "char")
         } else {
            varMem <- ncvar_def("member", units = "1", dim = list(dimens), prec = "int")
         }
         if (!is.null(coordBounds)) {
            var[[length(var) + 1]] <- varLonBounds
            var[[length(var) + 1]] <- varLatBounds
         }
         var[[length(var) + 1]] <- varMem
      } else {
         if (!is.null(coordBounds)){
            var[[length(var) + 1]] <- varLonBounds
            var[[length(var) + 1]] <- varLatBounds
         }
      }
   } else if (!is.null(attr(data$xyCoords, "projection")) & (length(grep("LambertConformal",
                                                                         attr(data$xyCoords, "projection"),
                                                                         ignore.case = TRUE)) > 0)){
      varProj <- ncvar_def("LambertConformal", units = "", dim = list(), prec = "char")
      varLon <- ncvar_def("lon", units = "degrees_east", dim = list(dimlon,dimlat),
                          longname = "longitude", prec = "double")
      varLat <- ncvar_def("lat", units = "degrees_north", dim = list(dimlon,dimlat),
                          longname = "latitude", prec = "double")
      if (!is.null(coordBounds)){
         varLonBounds <- ncvar_def("lon_vertices", units = "degrees_east",
                                   dim = list(dimBounds,dimlon,dimlat), longname = "longitude", prec = "double")
         varLatBounds <- ncvar_def("lat_vertices", units = "degrees_north",
                                   dim = list(dimBounds,dimlon,dimlat), longname = "latitude", prec = "double")
      }
      if (!transformeR::isMultigrid(data)) {
         var <- list(var, varLon, varLat, varProj)
      } else {
         var[[length(var) + 1]] <- varLon
         var[[length(var) + 1]] <- varLat
         var[[length(var) + 1]] <- varProj
      }
      if (length(member.index) > 0) {
         if (is.character(data$Members)) {
            varMem <- ncvar_def("member", units = "",
                                dim = list(dimnchar, dimens), prec = "char")
         } else {
            varMem <- ncvar_def("member", units = "1", dim = list(dimens), prec = "int")
         }
         if (!is.null(coordBounds)) {
            var[[length(var) + 1]] <- varLonBounds
            var[[length(var) + 1]] <- varLatBounds
         }
         var[[length(var) + 1]] <- varMem
      } else {
         if (!is.null(coordBounds)){
            var[[length(var) + 1]] <- varLonBounds
            var[[length(var) + 1]] <- varLatBounds
         }
      }
   } else {
      if (length(member.index) > 0) {
         if (is.character(data$Members)){
            varMem <- ncvar_def("member", units = "",
                                dim = list(dimnchar, dimens), prec = "char")
         } else {
            varMem <- ncvar_def("member", units = "1",
                                dim = list(dimens), prec = "int")
         }
         if (transformeR::isMultigrid(data)) {
            var[[length(var) + 1]] <- varMem
         } else {
            var <- list(var, varMem)
         }
      }
   }
   ncnew <- nc_create(NetCDFOutFile, var, verbose = verbose)
   ncatt_put(ncnew, "time", "standard_name","time")
   ncatt_put(ncnew, "time", "axis","T")
   ncatt_put(ncnew, "time", "_CoordinateAxisType","Time")
   ncatt_put(ncnew, "time", "_ChunkSize",1)
   if (!is.null(attr(data$xyCoords, "projection")) & attr(data$xyCoords,
                                                          "projection") == "RotatedPole") {
      ncatt_put(ncnew, "rlon", "standard_name","grid_longitude")
      ncatt_put(ncnew, "rlon", "axis","X")
      ncatt_put(ncnew, "rlon", "_CoordinateAxisType","GeoX")
      ncatt_put(ncnew, "rlat", "standard_name","grid_latitude")
      ncatt_put(ncnew, "rlat", "axis","Y")
      ncatt_put(ncnew, "rlat", "_CoordinateAxisType","GeoY")
      ncatt_put(ncnew, varProj$name, "grid_mapping_name","rotated_latitude_longitude")
      ncatt_put(ncnew, varProj$name, "grid_north_pole_latitude", gridNorthPole[1])
      ncatt_put(ncnew, varProj$name, "grid_north_pole_longitude", gridNorthPole[2])
      ncatt_put(ncnew, varLon$name, "standard_name","longitude")
      ncatt_put(ncnew, varLon$name, "_CoordinateAxisType","Lon")
      ncatt_put(ncnew, varLat$name, "standard_name","latitude")
      ncatt_put(ncnew, varLat$name, "_CoordinateAxisType","Lat")
   } else if (!is.null(attr(data$xyCoords, "projection")) & (length(grep("LambertConformal",
                                                                         attr(data$xyCoords, "projection"),
                                                                         ignore.case = TRUE)) > 0)) {
      ncatt_put(ncnew, "x", "standard_name","projection_x_coordinate")
      ncatt_put(ncnew, "x", "axis","X")
      ncatt_put(ncnew, "y", "standard_name","projection_y_coordinate")
      ncatt_put(ncnew, "y", "axis","Y")
      # ncatt_put(ncnew, "rlon", "_CoordinateAxisType","GeoX")
      # ncatt_put(ncnew, "rlat", "_CoordinateAxisType","GeoY")
      ncatt_put(ncnew, varProj$name, "grid_mapping_name","lambert_conformal_conic")
      ncatt_put(ncnew, varProj$name, "latitude_of_projection_origin",
                as.numeric(attr(data$xyCoords, "latitude_of_projection_origin")))
      ncatt_put(ncnew, varProj$name, "longitude_of_central_meridian",
                as.numeric(attr(data$xyCoords, "longitude_of_central_meridian")))
      if (length(grep("par2", attr(data$xyCoords, "projection"), ignore.case = TRUE)) > 0 ) {
         ch <- strsplit(attr(data$xyCoords, "projection"), ",")
         par1 <- strsplit(ch[[1]][grep("par1", ch[[1]])], "=")
         par2 <- strsplit(ch[[1]][grep("par2", ch[[1]])], "=")
         ncatt_put(ncnew, varProj$name, "standard_parallel", c(as.numeric(par1[[1]][2]), as.numeric(par2[[1]][2])))
      } else {
         ncatt_put(ncnew, varProj$name, "standard_parallel", as.numeric(attr(data$xyCoords, "standard_parallel")))
      }
      if (length(grep("falseEasting", attr(data$xyCoords, "projection"), ignore.case = TRUE)) > 0 ){
         ch <- strsplit(attr(data$xyCoords, "projection"), ",")
         par1 <- strsplit(ch[[1]][grep("falseEasting",ch[[1]])],"=")
         ncatt_put(ncnew, varProj$name, "false_easting",as.numeric(par1[[1]][2]))
      }
      if (length(grep("falseNorthing", attr(data$xyCoords, "projection"), ignore.case = TRUE)) > 0 ){
         ch <- strsplit(attr(data$xyCoords, "projection"), ",")
         par1 <- strsplit(gsub("}", " ", ch[[1]][grep("falseNorthing",ch[[1]])]), "=")
         ncatt_put(ncnew, varProj$name, "false_northing",as.numeric(par1[[1]][2]))
      }
      ncatt_put(ncnew, varLon$name, "standard_name","longitude")
      ncatt_put(ncnew, varLon$name, "_CoordinateAxisType","Lon")
      ncatt_put(ncnew, varLat$name, "standard_name","latitude")
      ncatt_put(ncnew, varLat$name, "_CoordinateAxisType","Lat")
   } else {
      ncatt_put(ncnew, "lon", "standard_name","longitude")
      ncatt_put(ncnew, "lon", "_CoordinateAxisType","Lon")
      ncatt_put(ncnew, "lat", "standard_name","latitude")
      ncatt_put(ncnew, "lat", "_CoordinateAxisType","Lat")
   }
   if (length(member.index) > 0) {
      ncatt_put(ncnew, "member", "standard_name","realization")
      ncatt_put(ncnew, "member", "_CoordinateAxisType","Ensemble")
      ncatt_put(ncnew, "member", "ref","http://www.uncertml.org/samples/realisation")
      Nvar <- 1
      if (transformeR::isMultigrid(data)) {
         Nvar <- length(tmpStdName)
      }
      for (v in c(1:Nvar)) {
         ncatt_put(ncnew, var[[v]]$name, "missing_value", missval, prec = prec)
         if (!is.null(varAttributes)) {
            sapply(1:length(varAttributes), function(x) {
               ncatt_put(ncnew, var[[v]]$name,
                         names(varAttributes)[x],
                         as.character(varAttributes[[x]]))
            })
         }
         if (transformeR::isMultigrid(data)){
            ncatt_put(ncnew, var[[v]]$name, "description", attributes(data$Variable)$"description"[[v]])
            ncatt_put(ncnew, var[[v]]$name, "longname", attributes(data$Variable)$"longname"[[v]])
         } else {
            ncatt_put(ncnew, var[[v]]$name, "description", attributes(data$Variable)$"description")
            ncatt_put(ncnew, var[[v]]$name, "longname", attributes(data$Variable)$"longname")
         }
         if (!is.null(attr(data$xyCoords, "projection")) & attr(data$xyCoords,
                                                                "projection") == "RotatedPole") {
            ncatt_put(ncnew, var[[v]]$name, "grid_mapping", "rotated_pole")
         }
      }
   } else {
      if ((!is.null(attr(data$xyCoords,
                         "projection"))) & ((attr(data$xyCoords,
                                                  "projection") == "RotatedPole") | (length(grep("LambertConformal",
                                                                                                 attr(data$xyCoords,
                                                                                                      "projection"),
                                                                                                 ignore.case = TRUE)) > 0))) {
         Nvar <- 1
         if (transformeR::isMultigrid(data)){
            Nvar <- length(tmpStdName)
         }
         for (v in c(1:Nvar)){
            ncatt_put(ncnew, var[[v]]$name, "missing_value", missval, prec = prec)
            if (!is.null(varAttributes)) {
               sapply(1:length(varAttributes), function(x) ncatt_put(ncnew, var[[v]]$name,
                                                                     names(varAttributes)[x],
                                                                     as.character(varAttributes[[x]])))
            }
            if (transformeR::isMultigrid(data)){
               ncatt_put(ncnew, var[[v]]$name, "description", attributes(data$Variable)$"description"[[v]])
               ncatt_put(ncnew, var[[v]]$name, "longname", attributes(data$Variable)$"longname"[[v]])
            } else {
               ncatt_put(ncnew, var[[v]]$name, "description", attributes(data$Variable)$"description")
               ncatt_put(ncnew, var[[v]]$name, "longname", attributes(data$Variable)$"longname")
            }
            if (!is.null(attr(data$xyCoords, "projection")) & attr(data$xyCoords, "projection") == "RotatedPole"){
               ncatt_put(ncnew, var[[v]]$name, "grid_mapping", "rotated_pole")
            } else if (!is.null(attr(data$xyCoords,
                                    "projection")) & (length(grep("LambertConformal",
                                                                  attr(data$xyCoords, "projection"),
                                                                  ignore.case = TRUE)) > 0)){
               ncatt_put(ncnew, var[[v]]$name, "grid_mapping", "lambert_conformal_conic")
            }
         }
      } else {
         if (transformeR::isMultigrid(data)) {
            Nvar <- length(tmpStdName)
            for (v in c(1:Nvar)) {
               ncatt_put(ncnew, var[[v]]$name, "missing_value", missval, prec = prec)
               if (!is.null(varAttributes)) {
                  sapply(1:length(varAttributes), function(x) {
                     ncatt_put(ncnew, var[[v]]$name, names(varAttributes)[x], as.character(varAttributes[[x]]))
                  })
               }
               ncatt_put(ncnew, var[[v]]$name, "description", attributes(data$Variable)$"description"[[v]])
               ncatt_put(ncnew, var[[v]]$name, "long_name", attributes(data$Variable)$"longname"[[v]])
            }
         } else {
            ncatt_put(ncnew, var$name, "missing_value", missval, prec = prec)
            if (!is.null(varAttributes)) {
               sapply(1:length(varAttributes), function(x) {
                  ncatt_put(ncnew, var$name, names(varAttributes)[x], as.character(varAttributes[[x]]))
               })
            }
            ncatt_put(ncnew, var$name, "description", attributes(data$Variable)$"description")
            ncatt_put(ncnew, var$name, "long_name", attributes(data$Variable)$"longname")
         }
      }
   }
   ## ncatt_put(ncnew, data$Variable$varName, "missing_value", missval)
   z <- attributes(data$Variable$level)
   if (!is.null(z)) ncatt_put(ncnew, var$name, "level", z)
   if (is.null(globalAttributes)) globalAttributes <- attributes(data)[-1]
   if (!is.null(globalAttributes)) {      
      indGlobal <- which(names(globalAttributes) != "_NCProperties") ## _NCProperties is a Global Attribute imposed by the writting software,
      if (length(indGlobal) > 0){
         sapply(indGlobal, function(x) {
            ncatt_put(ncnew, 0, names(globalAttributes)[x], as.character(globalAttributes[[x]]))}
         )
      }
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
   ncatt_put(ncnew, 0, "acknowledgement", "NetCDF file created by loadeR.2nc, a package of the 'climate4R' open-source framework <https://github.com/SantanderMetGroup/climate4R>")
   if (is.null(globalAttributes[["Conventions"]])) {
      ncatt_put(ncnew, 0, "Conventions", "CF-1.8, ACDD-1.3")   
   }
   
   if (((!is.null(attr(data$xyCoords,
                       "projection"))) & ((attr(data$xyCoords,
                                                "projection") == "RotatedPole") | (length(grep("LambertConformal",
                                                                                               attr(data$xyCoords,
                                                                                                    "projection"),
                                                                                               ignore.case = TRUE)) > 0)))  | is.character(data$Members)){
      if (transformeR::isMultigrid(data)) {
         for (v in c(1:length(tmpStdName))) {
            ncvar_put(ncnew, var[[v]], dataOrdered[[v]])
         }
      } else {
         v <- 1
         ncvar_put(ncnew, var[[v]], dataOrdered)
      }
      if ((!is.null(attr(data$xyCoords,
                         "projection"))) & ((attr(data$xyCoords,
                                                  "projection") == "RotatedPole") | (length(grep("LambertConformal",
                                                                                                 attr(data$xyCoords,
                                                                                                      "projection"),
                                                                                                 ignore.case = TRUE)) > 0))) {
         ncvar_put(ncnew, var[[length(tmpStdName)+1]], t(data$xyCoords$lon))
         ncvar_put(ncnew, var[[length(tmpStdName)+2]], t(data$xyCoords$lat))
         if (!is.null(coordBounds)) {
            ncvar_put(ncnew, var[[length(tmpStdName)+4]], aperm(coordBounds$lon, c(3,2,1)))
            ncvar_put(ncnew, var[[length(tmpStdName)+5]], aperm(coordBounds$lat, c(3,2,1)))
         }
      }
      if (is.character(data$Members)){
         ncvar_put(ncnew, var[[length(var)]], data$Members)
      }
   } else if (transformeR::isMultigrid(data)) {
      for (v in c(1:length(tmpStdName))) {
         ncvar_put(ncnew, var[[v]], dataOrdered[[v]])
      }
   } else {
      ncvar_put(ncnew, var, dataOrdered)
   }
   nc_close(ncnew)
   message("[", Sys.time(), "] NetCDF file written in: ", NetCDFOutFile)
}
