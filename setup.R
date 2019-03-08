# Clear all variables
rm(list = ls())



# Gather required packages
packages <- c("bigmemory", "biganalytics", "BiocManager", 
             "curl", "compiler", "car", 
             "data.table", "devtools", "DataCombine", 
             "httr", 
             "lme4", 
             "rvest", 
             "sparklyr", 
             "tidyverse", 
             "yaml")

# Check packages and install them if needed
invisible(lapply(packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "https://cran.cnr.berkeley.edu/")
    library(x, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
bioc_packages <- c("multtest")

# Check packages and install them if needed
invisible(lapply(bioc_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x)
    library(x, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
github_packages <- c("ggbiplot")

# Check packages and install them if needed
invisible(lapply(github_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install_github("vqv/ggbiplot")
    library(x, character.only = TRUE)
  }
}))



# Print this after all packages are successfully installed
loaded_packages <- sessionInfo()
print(names(loaded_packages$otherPkgs))



foldernames <- c("raw_data", "reference_files", "output")
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

filename <- "required_data.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}



filename <- "raw_add_A.R"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", "raw_data", filename))
if (file.exists(file.path("..", "raw_data", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}