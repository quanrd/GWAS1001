## Clean the workspace
rm(list=ls())


## Import library
library(sparklyr)
library(tidyverse)
library(yaml)
library(data.table)
library(lme4)
library(car)


# Source R files
source("func_remove_duplicates.R")
source("func_outlier_removal.R")
source("func_boxcox_transformation.R")
source("func_generate_BLUP.R")


## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)


## Import configuration file and read in raw data and reference files
yaml_filename <- file.path("../required_data.yaml")
if(file.exists(yaml_filename)){
  yaml_dat <- try(read_yaml(yaml_filename))
  
  if(exists("yaml_dat")){
    
    ## Import raw data
    if(!is.null(yaml_dat$raw_data)){
      if(file.exists(yaml_dat$raw_data)){
        if(endsWith(yaml_dat$raw_data, ".csv")){
          raw_data <- try(read.csv(yaml_dat$raw_data, header = TRUE, stringsAsFactors = FALSE))
        } else if(endsWith(yaml_dat$raw_data, ".txt")){
          raw_data <- try(read.table(yaml_dat$raw_data, header = TRUE, stringsAsFactors = FALSE))
        }
      } else{
        print("The raw data file does not exists.")
      }
    } else{
      print("The raw_data parameter is NULL.")
    }
    
    ## Create output folder
    if(!is.null(yaml_dat$output)){
      if(!dir.exists(yaml_dat$output)){
        dir.create(path = yaml_dat$output, showWarnings = TRUE, recursive = TRUE)
      } else{
        print("The output folder exists.")
      }
    } else{
      print("The output parameter is NULL.")
    }

    ## Import reference files
    # Genotype file
    if(!is.null(yaml_dat$genotype_file)){
      if(file.exists(yaml_dat$genotype_file)){
        if(endsWith(yaml_dat$genotype_file, ".csv")){
          genotype_file <- try(read.csv(yaml_dat$genotype_file, header = TRUE, stringsAsFactors = FALSE))
        } else if(endsWith(yaml_dat$genotype_file, ".txt")){
          genotype_file <- try(read.table(yaml_dat$genotype_file, header = TRUE, stringsAsFactors = FALSE))
        }
      } else{
        print("The Genotype file does not exists.")
      }
    } else{
      print("The genotype_file parameter is NULL.")
    }
    
    # SNPs file
    if(!is.null(yaml_dat$SNPs_file)){
      if(file.exists(yaml_dat$SNPs_file)){
        if(endsWith(yaml_dat$SNPs_file, ".csv")){
          SNPs_file <- try(read.csv(yaml_dat$SNPs_file, header = TRUE, stringsAsFactors = FALSE))
        } else if(endsWith(yaml_dat$SNPs_file, ".txt")){
          SNPs_file <- try(read.table(yaml_dat$SNPs_file, header = TRUE, stringsAsFactors = FALSE))
        }
      } else{
        print("The SNPs file does not exists.")
      }
    } else{
      print("The SNPs_file parameter is NULL.")
    }

    # Ecotype file
    if(!is.null(yaml_dat$ecotype_file)){
      if(file.exists(yaml_dat$ecotype_file)){
        if(endsWith(yaml_dat$ecotype_file, ".csv")){
          ecotype_file <- try(read.csv(yaml_dat$ecotype_file, header = TRUE, stringsAsFactors = FALSE))
        } else if(endsWith(yaml_dat$ecotype_file, ".txt")){
          ecotype_file <- try(read.table(yaml_dat$ecotype_file, header = TRUE, stringsAsFactors = FALSE))
        }
      } else{
        print("The ecotype file does not exists.")
      }
    } else{
      print("The ecotype_file parameter is NULL.")
    }

    # Protein file
    if(!is.null(yaml_dat$protein_file)){
      if(file.exists(yaml_dat$protein_file)){
        if(endsWith(yaml_dat$protein_file, ".csv")){
          protein_file <- try(read.csv(yaml_dat$protein_file, header = TRUE, stringsAsFactors = FALSE))
        } else if(endsWith(yaml_dat$protein_file, ".txt")){
          protein_file <- try(read.table(yaml_dat$protein_file, header = TRUE, stringsAsFactors = FALSE))
        }
      } else{
        print("The protein file does not exists.")
      }
    } else{
      print("The protein_file parameter is NULL.")
    }
    
  } else{
    print("The yaml data was not read into the work space!!!")
  }
  
} else{
  print("The yaml file does not exists!!!")
}

## Main
if(identical(args, character(0))){
  print("No arguments found!!!")
} else{
  
  # removeDuplicates
  if(all("-removeDuplicates" %in% args)){
    index <- match("-removeDuplicates", args)
    print(paste(index, ": removeDuplicates", sep = ""))
    
    if(!is.null(yaml_dat$raw_data) & exists("raw_data") & 
       !is.null(yaml_dat$output) & dir.exists(yaml_dat$output) & 
       !is.null(yaml_dat$by_column)){
      
      # Using customized function to remove duplicates
      raw_data <- remove_duplicates(dat = raw_data, by_column = yaml_dat$by_column)
      
      write.csv(x = raw_data, file = paste(yaml_dat$output, "raw_data.csv", sep = ""), row.names = TRUE, na = "")
    }
  }
  
  # outlierRemoval
  if(all("-outlierRemoval" %in% args)){
    index <- match("-outlierRemoval", args)
    print(paste(index, ": outlierRemoval", sep = ""))
    
    if(!is.null(yaml_dat$raw_data) & exists("raw_data") & 
       !is.null(yaml_dat$output)  & dir.exists(yaml_dat$output) & 
       !is.null(yaml_dat$by_column) & !is.null(yaml_dat$start_column)){
      
      # Using customized function to remove outliers
      results <- outlier_removal(dat = raw_data, by_column = yaml_dat$by_column, start_column = yaml_dat$start_column)
      
      if(is.list(results)){
        
        raw_data <- results$Outlier_removed_data
        
        capture.output(results$Outliers_residuals, file = paste(yaml_dat$output, "Outliers_residuals.txt", sep = ""))
        write.csv(x = results$Outlier_data, file = paste(yaml_dat$output, "Outlier_data.csv", sep = ""), row.names = FALSE, na = "")
        write.csv(x = results$Outlier_removed_data, file = paste(yaml_dat$output, "Outlier_removed_data.csv", sep = ""), row.names = FALSE, na = "")
      } else{
        
        raw_data <- results
        
        write.csv(x = results, file = paste(yaml_dat$output, "Outlier_removed_data.csv", sep = ""), row.names = FALSE, na = "")
        print("No outlier has been found!!!")
      }
      
    }
  }
  
  # boxcoxTransformation
  if(all("-boxcoxTransformation" %in% args)){
    index <- match("-boxcoxTransformation", args)
    print(paste(index, ": boxcoxTransformation", sep = ""))
    
    if(!is.null(yaml_dat$raw_data) & exists("raw_data") & 
       !is.null(yaml_dat$output)  & dir.exists(yaml_dat$output) & 
       !is.null(yaml_dat$by_column) & !is.null(yaml_dat$start_column)){
      
      # Using customized function to perform box-cox transformation
      results <- boxcox_transformation(dat = raw_data, by_column = yaml_dat$by_column, start_column = yaml_dat$start_column)
      
      if(is.list(results)){
        
        raw_data <- results$Boxcox_transformed_data
        
        capture.output(results$Lambda_values, file = paste(yaml_dat$output, "Lambda_values.txt", sep = ""))
        write.csv(x = results$Boxcox_transformed_data, file = paste(yaml_dat$output, "Boxcox_transformed_data.csv", sep = ""), row.names = FALSE, na = "")
      } else{
        
        raw_data <- results
        
        print("No lambda found!!! Data returned without transformed!!!")
      }
      
    }
  }
  
  # generateBLUP
  if(all("-generateBLUP" %in% args)){
    index <- match("-generateBLUP", args)
    print(paste(index, ": generateBLUP", sep = ""))
    
    if(!is.null(yaml_dat$raw_data) & exists("raw_data") & 
       !is.null(yaml_dat$output)  & dir.exists(yaml_dat$output) & 
       !is.null(yaml_dat$by_column) & !is.null(yaml_dat$start_column)){
      
      # Using customized function to generate BLUP
      results <- generate_BLUP(dat = raw_data, by_column = yaml_dat$by_column, start_column = yaml_dat$start_column)
      
      if(is.list(results)){
        BLUP <- results$BLUP
        
        write.csv(x = results$BLUP, file = paste(yaml_dat$output, "BLUP.csv", sep = ""), row.names = FALSE, na = "")
      } else{
        print("No BLUP generated!!!")
      }
      
      
    }
  }
  
  # FarmCPU
  if(all("-farmCPU" %in% args)){
    index <- match("-farmCPU", args)
    print(paste(index, ": farmCPU", sep = ""))
  }

}

