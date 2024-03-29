% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid2nc.R
\name{grid2nc}
\alias{grid2nc}
\title{climate4R grid 2 netCDF export}
\usage{
grid2nc(
  data,
  NetCDFOutFile = "out.nc4",
  missval = 1e+20,
  globalAttributes = NULL,
  varAttributes = NULL,
  prec = "float",
  compression = 4,
  shuffle = FALSE,
  verbose = FALSE,
  gridNorthPole = c("39.25", "-162.0"),
  coordBounds = NULL
)
}
\arguments{
\item{data}{A climate4R grid data object (\url{http://www.meteo.unican.es/climate4R})}

\item{NetCDFOutFile}{Name of the file created by the function. (default to \code{"out.nc4"})}

\item{missval}{Missing value codification (default to \code{1e20})}

\item{globalAttributes}{Optional. A list of global attributes included in the NetCDF file. Same format as \code{varAttributes}.}

\item{varAttributes}{Optional. List of attributes to be included in the variable written in the NetCDF file. Default to \code{NULL}.
It has the format \code{list("name_of_attribute1" = attribute1, "name_of_attribute2" = attribute2} etc.}

\item{prec}{Precision to write the attribute. If not specified, the written precision is given by the variable's
corresponding attribute. This can be overridden by specifying the following argument values:
 \code{"short"}, \code{"float"}, \code{"double"}, or \code{"text"}.}

\item{compression}{Integer value between 1 (less compression) and 9 (highest compression), enabling 
compression for the variable as it is written to the file. Note that turning compression on forces the 
output file to netcdf v4 format (the default), not compatible with older software 
that can only handle version 3 files.}

\item{shuffle}{Logical. Turns on (\code{TRUE}) or off (\code{FALSE}, the default) the shuffle filter.
According to netcdf documentation, turning the shuffle filter on can improve compression for integer variables.
 Note that turning the shuffle filter on forces the output file to netcdf v4 format (the default), 
 not compatible with older software that can only handle version 3 files.}

\item{verbose}{Optional. If set to \code{TRUE}, switches \code{\link[ncdf4]{nc_create}} to verbose mode}

\item{gridNorthPole}{Vector with the longitude and latitude. Default is c("39.25","-162.0")}

\item{coordBounds}{Default is NULL.}
}
\value{
A NetCDF-4 file with the variable and attributes defined in the inputs.
}
\description{
Export a climate4R grid to NetCDF
}
\details{
\strong{netCDF attributes}
Some attributes are automatically inferred by the function from the climate4R input grid metadata.
These can be replaced by user-defined attributes via the global and/or variable attribute lists.
It is highly advised that an authoritative attribute convention is followed, such as CF and/or ESIP:
\itemize{
\item \url{https://cfconventions.org/}
\item \url{https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#creator_url}
}
This conventions are indicated by default in the global attributes of the output.
}
\examples{
\dontrun{
library(loadeR.2nc)
data(tx) # A climate4R grid
# Name of output file:
fileName <- "tasmax_WFDEI_JJA_W2001_2010.nc4"
# Including a global attribute:
globalAttributeList <- list("institution" = "SantanderMetGroup, http://www.meteo.unican.es/")
# Including two variable attributes:
varAttributeList <- list(var_attr1 = "one_attribute", var_attr2 = "another_attribute")
# Create file:
grid2nc(data = tx,
        NetCDFOutFile = fileName,
        missval = 1e20,
        prec = "float",
        globalAttributes = globalAttributeList,
        varAttributes = varAttributeList)
# Inspection fo the file: Requires nectcdf-bin installed in your system:
system("ncdump -h tasmax_WFDEI_JJA_W2001_2010.nc4")
}
}
\references{
\itemize{
\item David Pierce \email{dpierce@ucsd.edu}, Interface to Unidata netCDF (version 4 or earlier) format data files, http://dwpierce.com/software
}
}
\author{
S. Herrera, W. Franssen and J. Bedia
}
