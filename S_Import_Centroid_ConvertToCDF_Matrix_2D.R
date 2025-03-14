#Load packages
library("mzR")
library("data.table")
library("parallel")
library(RNetCDF)
library(stringr)

#Set paths
mzML_dir <- "F:/PhD/Matrix/IMS/RPLCxHILIC_mzML/RPLCxHILIC_Wheat/SecondAttempt"
cdf_dir <- "F:/PhD/Matrix/IMS/RPLCxHILIC_CDF/Matrix/RPLCxHILIC/SecondAttempt"
          
setwd(mzML_dir)
fileList <- dir(pattern = ".mzML")

#Get functions
source("C:/Users/mht541/Documents/2D_LC_IMS_R/2DIMS/Import_centroid_cdf/functions.R")  # Assuming functions.R contains the centroiding_parallel and MakeVector functions

#specify number of cores
max_cores <- 10
SampleVec <- (length(fileList) -1)
# Loop over fileList
for (k in SampleVec) { #
  
  #Set working directory
  setwd(mzML_dir)
  
  # Load data
  sample_file <- fileList[k]
  centroided_data <- process_centroiding(sample_file, max_cores)
  
  # Read meta data
  metadata <- read_metadata(sample_file, centroided_data$idx)
  
  # Save to CDF
  save_to_cdf(sample_file, cdf_dir, metadata, centroided_data)
  
  #Remove these variables from the environment 
  rm(list = intersect(ls(), c("centroided_data", "metadata")))
  gc()
  }

## #Set paths
mzML_dir <- "F:/PhD/Matrix/IMS/RPLCxHILIC_mzML/RPLCxHILIC_Wheat/SecondAttempt/MS3"
cdf_dir <- "F:/PhD/Matrix/IMS/RPLCxHILIC_CDF/Matrix/RPLCxHILIC/SecondAttempt/MS3"

setwd(mzML_dir)
fileList <- dir(pattern = ".mzML")

# Loop over fileList
for (k in SampleVec) { #
  
  #Set working directory
  setwd(mzML_dir)
  
  # Load data
  sample_file <- fileList[k]
  centroided_data <- process_centroiding(sample_file, max_cores)
  
  # Read meta data
  metadata <- read_metadata(sample_file, centroided_data$idx)
  
  # Save to CDF
  save_to_cdf(sample_file, cdf_dir, metadata, centroided_data)
  
  #Remove these variables from the environment 
  rm(list = intersect(ls(), c("centroided_data", "metadata")))
}

