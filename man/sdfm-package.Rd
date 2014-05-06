\name{sdfm-package}
\alias{sdfm-package}
\alias{sdfm}
\docType{package}
\title{
Processing of SDF and SMI chemistry files
}
\description{
Processing of SDF and SMI chemistry files
}
\details{
\tabular{ll}{
Package: \tab sdfm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-16\cr
License: \tab GPL3\cr
}
Processing of SDF and SMI chemistry files
}
\author{
Pau Carrió
}
\references{

}
\keyword{ package }
\seealso{
}
\examples{
 
# data <- readSdf( "mySDF.sdf" )
# view fields in the sdf
# a special field STRUCTURE stores the structure information
# colnames(data)
# create new fields with
# data$activity <- ( -log(data$FIELD1) )
# save sdf
# writeSdf( data, "myNewSDF.sdf")

}
