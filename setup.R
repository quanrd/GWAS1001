# Clear all variables
rm(list = ls())



# Gather required packages
packages <- c("bigmemory", "biganalytics", "BiocManager", 
             "curl", "compiler", "car", 
             "data.table", "devtools", 
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



if(!dir.exists(file.path("..", "raw_data"))){
  dir.create(file.path("..", "raw_data"))
  if(dir.exists(file.path("..", "raw_data"))){
    print("raw_data folder has been created!!!")
  } else{
    print("raw_data folder cannot be created!!!")
  }
} else{
  print("raw_data folder exists!!!")
}

if(!dir.exists(file.path("..", "reference_files"))){
  dir.create(file.path("..", "reference_files"))
  if(dir.exists(file.path("..", "reference_files"))){
    print("reference_files folder has been created!!!")
  } else{
    print("reference_files folder cannot be created!!!")
  }
} else{
  print("reference_files folder exists!!!")
}

if(!dir.exists(file.path("..", "output"))){
  dir.create(file.path("..", "output"))
  if(dir.exists(file.path("..", "output"))){
    print("output folder has been created!!!")
  } else{
    print("output folder cannot be created!!!")
  }
} else{
  print("output folder exists!!!")
}



dat <- readLines(con = "required_data.YAML")
writeLines(dat, con = file.path("..", "required_data.YAML"))
if (file.exists(file.path("..", "required_data.YAML"))) {
  print("required_data.YAML has been created!!!")
} else{
  print("required_data.YAML cannot be created!!!")
}



dat <- readLines(con = "raw_add_A.R")
writeLines(dat, con = file.path("..", "raw_data", "raw_add_A.R"))
if (file.exists(file.path("..", "raw_data", "raw_add_A.R"))) {
  print("raw_add_A.R has been created!!!")
} else{
  print("raw_add_A.R cannot be created!!!")
}