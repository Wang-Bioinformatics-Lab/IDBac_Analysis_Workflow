# Load the MALDIquant package
library(MALDIquant)
library(MALDIquantForeign)

process_mzML_file <- function(input_file, output_file, contains_ms2) {
  # Check if MS1+MS2 code path
  if (contains_ms2) {
    # Simply, reoutput
    cat("Input spectra contain MS2 data. Skipping baseline correction and peak detection.\n")
    file.copy(input_file, output_file)
    return()
  }

  # Read the mzML file
  spectra <- importMzMl(input_file)

  # Check if any intensities are greater than zero
  any_nonzero_intensity <- any(sapply(spectra, function(x) any(x@intensity > 0)))

  if (!any_nonzero_intensity) {
    # Print warning message
    cat("Warning: All intensities in the spectra are zero. Skipping baseline correction and peak detection.\n")

    # Copy the input file to the output file (to avoid errors in the next step)
    file.copy(input_file, output_file)
    
    return()
  }

  cat("Performing baseline correction and peak detection on MS1 spectra.\n")

  spectra <- MALDIquant::smoothIntensity(spectra,
                                         method = 'SavitzkyGolay',
                                         halfWindowSize = 20L)
                                         
  # Perform baseline subtraction using the SNIP algorithm
  spectra_baseline_corrected <- removeBaseline(spectra, method="TopHat")

  # Perform peak detection
  peaks <- detectPeaks(spectra_baseline_corrected,
                      method='SuperSmoother',
                      halfWindowSize=20L,
                      SNR=1)

  # Export the processed spectra as an mzML file
  exportMzMl(peaks, output_file)
}

# Get the input and output file paths from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
contains_ms2 <- args[3]

contains_ms2 <- tolower(args[3]) %in% c("yes", "y", "true")

process_mzML_file(input_file, output_file, contains_ms2)