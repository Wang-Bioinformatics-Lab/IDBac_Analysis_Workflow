# Load the MALDIquant package
library(MALDIquant)
library(MALDIquantForeign)

process_mzML_file <- function(input_file, output_file) {
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

  # Perform baseline subtraction using the SNIP algorithm
  spectra_baseline_corrected <- removeBaseline(spectra, method="SNIP", iterations=100)

  # Perform peak detection
  peaks <- detectPeaks(spectra_baseline_corrected)

  # Export the processed spectra as an mzML file
  exportMzMl(peaks, output_file)
}

# Get the input and output file paths from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

process_mzML_file(input_file, output_file)
