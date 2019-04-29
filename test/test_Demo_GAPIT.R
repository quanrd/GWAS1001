## Clean the workspace
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



# Older version of GAPIT
# source("https://raw.githubusercontent.com/yenon118/GWAS_Cli_App/master/emma.txt")
# source("https://raw.githubusercontent.com/yenon118/GWAS_Cli_App/master/gapit_functions_20160415.txt")

# Import GAPIT Library
# source("http://www.zzlab.net/GAPIT/GAPIT.library.R")

# Import EMMA
source("http://www.zzlab.net/GAPIT/emma.txt")

# Import FarmCPU
source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")

# Import GAPIT
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")


dat <- read.table(file = file.path("/home/ycth8/data/GWAS1001/BLUP_BLUE/Demo/mdp_traits.txt"), 
                header = TRUE,
                stringsAsFactors = FALSE, 
                check.names = FALSE)

if(nrow(dat) >= 10){
  print(dat[1:10, ])
} else{
  print(dat)
}

if(dir.exists("/home/ycth8/data/GWAS1001/reference_files/Demo/")){
  print("dir exists")
}else{
  print("dir does not exist")
}

start_column <- 2

for (i in start_column:ncol(dat)){
  Results <- GAPIT(
    Y = dat[,c(1,i)],
    file.Ext.G = "hmp.txt",
    file.G = "mdp_genotype_chr",
    file.from = 1,
    file.to = 10,
    file.path = "/home/ycth8/data/GWAS1001/reference_files/Demo/"
  )
}