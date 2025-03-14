# Functions

##################################            CENTROIDING          ###################
############# MakeVector2 USED
MakeVectorChunked <- function(pks, chunk_size = 10000) {
  # Initialize an empty data frame to store the result
  result <- data.frame()
  
  # Total number of elements in `pks`
  total <- length(pks)
  # Process `pks` in chunks
  for (start_idx in seq(1, total, by = chunk_size)) {
    # Determine the end index for the current chunk
    end_idx <- min(start_idx + chunk_size - 1, total)
    # Extract the current chunk
    chunk <- pks[start_idx:end_idx]
    # Convert each element of the chunk to a data frame and bind it to the result
    #result[start_idx:end_idx] <- rbind(result, rbindlist(lapply(chunk, as.data.frame)))
    
    result <- rbind(result, rbindlist(lapply(chunk, as.data.frame)))
    # Clear the processed chunk to free memory
    pks[start_idx:end_idx] <- 0
    #gc()  # Run garbage collection to reclaim memory
    
    # Print progress
    cat(sprintf("Processed %.1f%% of %d chunks\n", (end_idx / total) * 100, ceiling(total / chunk_size)))
    
  }
  
  return(result)
}


####### Read spectra_cluster_apply ## USED
read_spectra_clusterapply <- function(sample_file, max_cores = 4) {
  # Get the number of spectra
  mz <- openMSfile(sample_file, backend = "pwiz")
  num_spectra <- length(mz)
  close(mz)
  
  # Set up parallel processing
  num_cores <- min(detectCores()-1, max_cores)
  chunk_size <- ceiling(num_spectra / num_cores)
  chunks <- split(seq_len(num_spectra), ceiling(seq_len(num_spectra) / chunk_size))
  
  # Create a cluster
  cl <- makeCluster(num_cores)
  
  # Load necessary libraries in each worker
  clusterEvalQ(cl, {
    library(mzR)
    library(MALDIquant)
  })
  
  # Define the process_chunk function without `%>%`
  process_chunk <- function(chunk, sample_file) {
    mz <- openMSfile(sample_file, backend = "pwiz")
    spectra_list <- lapply(chunk, function(scan_idx,halfWindowSize = 1) {
      spectrum_data <- peaks(mz, scan_idx)  # Retrieve spectrum data
      
      if (nrow(spectrum_data) >= 2 * halfWindowSize + 1) {
        # Create spectrum object and detect peaks if the spectrum has non-zero intensity
        spectrum <- createMassSpectrum(mass = spectrum_data[, 1], intensity = spectrum_data[, 2])
        detected_peaks <- detectPeaks(spectrum, halfWindowSize = 1, SNR = 0)
        return(matrix(c(detected_peaks@mass, detected_peaks@intensity), ncol = 2, byrow = FALSE))
      } else {
        return(matrix(,nrow = 0, ncol = 2, byrow = FALSE))
      } 
    })
    close(mz)
    return(spectra_list)
  }
  
  # Export variables and the function to cluster workers
  clusterExport(cl, varlist = c("sample_file", "process_chunk"), envir = environment())
  
  # Use clusterApply to process chunks
  results <- clusterApply(cl, chunks, function(chunk) {
    process_chunk(chunk, sample_file)
  })
  
  # Stop the cluster
  stopCluster(cl)
  print("Make into one vector")
  # Combine results into a single list
  all_spectra <- vector("list", num_spectra)
  for (i in seq_along(chunks)) {
    all_spectra[chunks[[i]]] <- results[[i]]
  }
  
  return(all_spectra)
}



########################          Pipeline ###################################

# Centroiding process
process_centroiding <- function(sample_file, max_cores) {
  cat('Centroiding file:  ', sample_file,'\n')
  
  # Perform centroiding
  pks_cent <- read_spectra_clusterapply(sample_file, max_cores)
  
  # Count peaks in each spectrum
  pks_count <- unlist(lapply(pks_cent, nrow))
  
  # Filter non-empty spectras
  idx <- which(pks_count > 0)
  pks_count <- pks_count[idx]
  pks_cent <- pks_cent[idx]  # Keep only non-empty spectra
  
  # Convert to vector with chunking
  pks_cent <- MakeVectorChunked(pks_cent)
  
  # Return both values as a list
  return(list(pks_cent = pks_cent, pks_count = pks_count,idx = idx))
}

#Read meta data
read_metadata <- function(sample_file, idx) {
  cat('Reading meta data of file:  ', sample_file,'\n')
  
  # Open MS file
  mz <- openMSfile(sample_file, backend = "pwiz")
  
  # Read header for specified scans
  hdr <- header(mz, scans = c(1:max(idx)))
  
  # Extract metadata fields
  retentionTime <- hdr$retentionTime
  acquisitionNum <- hdr$acquisitionNum[idx]
  msLevel <- hdr$msLevel
  numberRows <- nrow(hdr)
  
  # Check for ion mobility drift time
  if (any(!is.na(hdr$ionMobilityDriftTime))) {
    ionMobilityDriftTime <- hdr$ionMobilityDriftTime
  } else {
    ionMobilityDriftTime <- NULL
  }
  
  # Return metadata as a list
  return(list(
    retentionTime = retentionTime,
    acquisitionNum = acquisitionNum,
    ionMobilityDriftTime = ionMobilityDriftTime,
    msLevel = msLevel,
    numberRows = numberRows
  ))
}

## Saving to CDF
save_to_cdf <- function(sample_file, cdf_dir, metadata, centroid_data, chunk_size = 10000) {
  
  cat('Saving file:  ', sample_file, '  to CDF\n')
  
  # Check if the directory exists and if we have write permissions
  if (!dir.exists(cdf_dir)) {
    stop("The directory does not exist: ", cdf_dir)
  }
  
  if (file.access(cdf_dir, mode = 2) != 0) {
    stop("Permission denied: Cannot write to the directory")
  }
  
  # Set working directory
  setwd(cdf_dir)
  
  # Define output CDF file name with absolute path
  cdf_filename <- file.path(cdf_dir, paste0(substring(sample_file, 1, nchar(sample_file) - 5), ".cdf"))
  
  # Remove existing file if necessary
  if (file.exists(cdf_filename)) {
    cat("Removing existing CDF file: ", cdf_filename, "\n")
    file.remove(cdf_filename)
  }
  
  # Create NetCDF file
  nc <- create.nc(cdf_filename, format = "netcdf4")
  
  # Extract metadata
  numberRows <- metadata$numberRows
  retentionTime <- metadata$retentionTime
  acquisitionNum <- metadata$acquisitionNum
  msLevel <- metadata$msLevel
  numberOfNNZScans <- length(acquisitionNum) 
  
  # Extract centroiding results
  pks_cent <- centroid_data$pks_cent
  pks_count <- centroid_data$pks_count
  
  # Define dimensions
  dim.def.nc(nc, "retentionTime", numberRows)
  dim.def.nc(nc, "intensity", sum(pks_count))
  dim.def.nc(nc, "NonZeroScans", numberOfNNZScans)
  
  # Create variables
  var.def.nc(nc, "scan_acquisition_time", "NC_DOUBLE", "retentionTime")
  var.def.nc(nc, "scan_acquisition_number", "NC_INT", "NonZeroScans")
  var.def.nc(nc, "point_count", "NC_DOUBLE", "NonZeroScans")
  var.def.nc(nc, "mslevel", "NC_INT", "retentionTime")
  var.def.nc(nc, "mass_values", "NC_DOUBLE", "intensity")
  var.def.nc(nc, "intensity_values", "NC_INT", "intensity")
  
  # Optional: Define ion mobility drift time if available
  if (exists("ionMobilityDriftTime", where = metadata)){
    var.def.nc(nc, "drift_acquisition_time", "NC_DOUBLE", "retentionTime")
  }
  
  # Fill values
  var.put.nc(nc, "scan_acquisition_time", retentionTime)
  var.put.nc(nc, "scan_acquisition_number", acquisitionNum)
  var.put.nc(nc, "point_count", pks_count)
  var.put.nc(nc, "mslevel", msLevel)
  
  # Optional: Store ion mobility drift time
  if (exists("ionMobilityDriftTime", where = metadata)){
    var.put.nc(nc, "drift_acquisition_time", as.double(metadata$ionMobilityDriftTime))
  }
  
  # Chunked writing of mass and intensity values
  num_chunks <- ceiling(length(pks_cent$V1) / chunk_size)
  
  for (i in 1:num_chunks) {
    start <- (i - 1) * chunk_size + 1
    count <- min(chunk_size, length(pks_cent$V1) - start + 1)
    
    var.put.nc(nc, "mass_values", pks_cent$V1[start:(start + count - 1)], start = start)
    var.put.nc(nc, "intensity_values", pks_cent$V2[start:(start + count - 1)], start = start)
  }
  
  # Close the NetCDF file
  close.nc(nc)
  
  cat("CDF file saved:", cdf_filename, '\n')
}

