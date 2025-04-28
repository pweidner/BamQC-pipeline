	    library(optparse)
	    library(AneuFinder)
	    library(dplyr)
	    library(BSgenome.Hsapiens.UCSC.hg38)

# Define command-line options using optparse
	    option_list <- list(
				  make_option(c("-i", "--inputfolder"), type="character",
					                    help="Input folder containing BAM/other files for Aneufinder", metavar="folder"),
				  make_option(c("-o", "--outputfolder"), type="character",
					                    help="Output folder to store all Aneufinder results", metavar="folder"),
				  make_option(c("--nCPU"), type="integer", default=2,
					                    help="Number of CPU threads to use [default %default]", metavar="n")
				  )
	    opt <- parse_args(OptionParser(option_list=option_list))

	    # Create the persistent output folder if it does not exist.
	    if (!dir.exists(opt$outputfolder)) {
		      dir.create(opt$outputfolder, recursive=TRUE)
	    }

	    cat("Running Aneufinder on input folder:\n", opt$inputfolder, "\n")
	    cat("Results will be saved in:\n", opt$outputfolder, "\n")

	    # Load the hg38 BSgenome
	    hg38 <- BSgenome.Hsapiens.UCSC.hg38

	    # Run Aneufinder with the desired parameters.
	    result <- tryCatch({
		        Aneufinder(
				         inputfolder = opt$inputfolder,
					       outputfolder = opt$outputfolder,
					       numCPU = opt$nCPU,
					             method = "edivisive",
					             correction.method = "GC",
						           GC.bsgenome = hg38,
						           refine.breakpoints = FALSE
							       )
	    }, error = function(e) {
		        cat("Error during Aneufinder run:", conditionMessage(e), "\n")
	        stop(e)
	    })

	    cat("Aneufinder run completed successfully.\n")

	    # Create a marker file to indicate successful completion.
	    marker_file <- file.path(opt$outputfolder, "done.txt")
	    write("Aneufinder run completed successfully.", file = marker_file)
	    cat("Marker file written to", marker_file, "\n")
