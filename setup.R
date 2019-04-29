# Clear all variables
rm(list = ls())



# Set repository so that required packages can be downloaded
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)



# Install path
p <- "~/R/x86_64-redhat-linux-gnu-library/3.5/"
if(!dir.exists(p)){
  dir.create(path = p, recursive = TRUE)
}



# BiocManager package version
bioc_version <- "3.8"



# Gather required packages
packages <- c("ape", 
              "bigmemory", "biganalytics", "BiocManager", 
              "curl", "compiler", "car", 
              "data.table", "devtools", "DataCombine", "DBI", 
              "EMMREML", 
              "genetics", "gplots", "grid", 
              "httr", 
              "lme4", "LDheatmap", 
              "rvest", 
              "scatterplot3d", "sparklyr", 
              "plyr", "dplyr", "tidyverse", 
              "yaml")

# Check packages and install them if needed
invisible(lapply(packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "https://cran.cnr.berkeley.edu/", lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
bioc_packages <- c("zlibbioc", "snpStats", "multtest")

# Check packages and install them if needed
invisible(lapply(bioc_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, version = bioc_version, lib.loc = p, lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
github_packages <- c("ggbiplot")

# Check packages and install them if needed
invisible(lapply(github_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install_github("vqv/ggbiplot", lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# Print this after all packages are successfully installed
loaded_packages <- sessionInfo()
print(names(loaded_packages$otherPkgs))



# Create these folders
foldernames <- c("raw_data", "BLUP_BLUE", "reference_files", "output")
for (i in 1:length(foldernames)) {
  temp <- file.path("..", foldernames[i])
  if(!dir.exists(temp)){
    dir.create(temp)
    if(dir.exists(temp)){
      print(paste(foldernames[i], " folder has been created!!!", sep = ""))
    } else{
      print(paste(foldernames[i], " folder cannot be created!!!", sep = ""))
    }
  } else{
    print(paste(foldernames[i], " folder exists!!!", sep = ""))
  }
}



# Create Demo.yaml file 
filename <- "Demo.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}

# Create required_data.yaml file 
filename <- "required_data.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}

# Create Arabidopsis1001.yaml file 
filename <- "Arabidopsis1001.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}



# Copy this file into raw directory
filename <- "raw_add_A.R"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", "raw_data", filename))
if (file.exists(file.path("..", "raw_data", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}

