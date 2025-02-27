################################################################################
# This script was executed while setting up the R project for the analysis
# 
# It should be visited time-to-time to upgrade the package library and/or to
# install new packages.
# Once new packages are installed, a snapshot of latest package library is must
# 
# renv::init() will download and install the latest version of package from CRAN
# if you need a specific version of a package, you will need to install it
# manually yourself for the first time.
# For later instances, the renv.lock will have enough metadata to download and
# install the package on it's own.
# 
# 
# Reference: # https://rstudio.github.io/renv/articles/renv.html
# 
# 
# 
# Updated: 12/09/2024
# 
# The script can be executed line by line to setup Renv for the RStudio project
# Please read the comments carefully to understand the proceedings
# 
# 
# Written & implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################


# Clearing the R session
rm(list = ls())


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 | length(args) > 1) {
  stop("Please provide specified number of command-line argument.", call. = FALSE)
}


## R version check
r_ver <- paste0(unlist(unname(R.version$major)),
                ".",
                unlist(unname(R.version$minor)))

# checking the R version used in the session
if(r_ver != "4.1.2") {
  message("The R version (4.1.2) used in the original analysis do not match with R version (", r_ver, ") you are using.
  Please use the correct R version to fully replicate the environment for reproducibility.")
}


#### FOR COLLABORATOR ####
# If you want to use the renv functionality to reproduce the analysis

# below command creates a replica of R session environment used during the 
# original analysis
# This should install all the required packages and their dependencies from the 
# sources mentioned in the renv.lock file
Sys.setenv("GITHUB_PAT"= args)
message("Restoring R environment")
# run to restore packages to version recorded in lock file provided by
# the developer. This command takes significant amount of time to execute
renv::restore()

# write the current .libPaths() for renv to a .txt file to use later
writeLines(.libPaths(), "renv_library_paths.txt")


# #### FOR DEVELOPER ####
# # DO NOT RUN/UNCOMMENT IF YOU ARE A COLLABORATOR
# 
# # Getting started: to convert the project to use renv call below function
# # The command below adds files to your directory where the R project resides
# 
# renv::init()
# 
# # project library, ‘renv/library’ (contains packages currently used by the project)
# # lockfile, ‘renv.lock’ (metadata of each package)
# # project R profile, ‘.Rprofile’ (configure  R session to use the project library)
# 
# 
# # we want to capture all the packages installed in the project library
# renv::settings$snapshot.type("implicit")
# 
# # saving the metadata file (renv.lock)
# renv::snapshot()
