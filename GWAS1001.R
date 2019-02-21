## Clean the workspace
rm(list = ls())


## Import library
library(sparklyr)
library(tidyverse)
library(yaml)
library(data.table)
library(lme4)
library(car)
library(bigmemory)
library(biganalytics)
require(compiler)

#Import GAPIT
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#source("gapit_functions.txt")

#Import FarmCPU
source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")
#source("FarmCPU_functions.txt")

# Source R files
source("func_remove_duplicates.R")
source("func_outlier_removal.R")
source("func_boxcox_transformation.R")
source("func_generate_BLUP.R")
source("func_farming_with_FarmCPU.R")


## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the YAML file path is in the args
if (!identical(args, character(0)) &
    length(args) > 0 & file.exists(file.path(args[1]))) {
  print("YAML file exists!!!")
  
  if (length(args) > 1) {
    print(args[2:length(args)])
    
    ## Import configuration file and read in raw data and reference files
    yaml_dat <- try(read_yaml(file.path(args[1])))
    
    #######################################################################
    ## Read in all file and data that are specified in the YAML file
    #######################################################################
    
    if (exists("yaml_dat")) {
      
      ## Import raw data
      if (!is.null(yaml_dat$raw_data)) {
        if (file.exists(yaml_dat$raw_data)) {
          if (endsWith(yaml_dat$raw_data, ".csv")) {
            raw_data <-
              try(read.csv(yaml_dat$raw_data,
                           header = TRUE,
                           stringsAsFactors = FALSE))
          } else if (endsWith(yaml_dat$raw_data, ".txt")) {
            raw_data <-
              try(read.table(yaml_dat$raw_data,
                             header = TRUE,
                             stringsAsFactors = FALSE))
          }
        } else{
          print("The raw data file does not exists.")
        }
      } else{
        print("The raw_data parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$by_column)) {
        by_column <- yaml_dat$by_column
      } else{
        print("The by_column parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$start_column)) {
        start_column <- yaml_dat$start_column
      } else{
        print("The start_column parameter is NULL.")
      }
      
      ## Import BULP data
      if (!is.null(yaml_dat$BLUP)) {
        if (file.exists(yaml_dat$BLUP)) {
          if (endsWith(yaml_dat$BLUP, ".csv")) {
            BLUP <-
              try(read.csv(yaml_dat$BLUP,
                           header = TRUE,
                           stringsAsFactors = FALSE))
          } else if (endsWith(yaml_dat$BLUP, ".txt")) {
            BLUP <-
              try(read.table(yaml_dat$BLUP,
                             header = TRUE,
                             stringsAsFactors = FALSE))
          }
        } else{
          print("The raw data file does not exists.")
        }
      } else{
        print("The BLUP parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$BLUP_by_column)) {
        BLUP_by_column <- yaml_dat$BLUP_by_column
      } else{
        print("The BLUP_by_column parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$BLUP_start_column)) {
        BLUP_start_column <- yaml_dat$BLUP_start_column
      } else{
        print("The BLUP_start_column parameter is NULL.")
      }
      
      ## Import reference files
      # Genotype file
      if (!is.null(yaml_dat$genotype_file)) {
        if (file.exists(yaml_dat$genotype_file)) {
          if (endsWith(yaml_dat$genotype_file, ".csv")) {
            genotype_file <-
              try(read.csv(
                yaml_dat$genotype_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          } else if (endsWith(yaml_dat$genotype_file, ".txt")) {
            genotype_file <-
              try(read.table(
                yaml_dat$genotype_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          }
        } else{
          print("The Genotype file does not exists.")
        }
      } else{
        print("The genotype_file parameter is NULL.")
      }
      
      # SNPs file
      if (!is.null(yaml_dat$SNPs_file)) {
        if (file.exists(yaml_dat$SNPs_file)) {
          if (endsWith(yaml_dat$SNPs_file, ".csv")) {
            SNPs_file <-
              try(read.csv(yaml_dat$SNPs_file,
                           header = TRUE,
                           stringsAsFactors = FALSE))
          } else if (endsWith(yaml_dat$SNPs_file, ".txt")) {
            SNPs_file <-
              try(read.table(yaml_dat$SNPs_file,
                             header = TRUE,
                             stringsAsFactors = FALSE))
          }
        } else{
          print("The SNPs file does not exists.")
        }
      } else{
        print("The SNPs_file parameter is NULL.")
      }
      
      # Ecotype file
      if (!is.null(yaml_dat$ecotype_file)) {
        if (file.exists(yaml_dat$ecotype_file)) {
          if (endsWith(yaml_dat$ecotype_file, ".csv")) {
            ecotype_file <-
              try(read.csv(
                yaml_dat$ecotype_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          } else if (endsWith(yaml_dat$ecotype_file, ".txt")) {
            ecotype_file <-
              try(read.table(
                yaml_dat$ecotype_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          }
        } else{
          print("The ecotype file does not exists.")
        }
      } else{
        print("The ecotype_file parameter is NULL.")
      }
      
      # Protein file
      if (!is.null(yaml_dat$protein_file)) {
        if (file.exists(yaml_dat$protein_file)) {
          if (endsWith(yaml_dat$protein_file, ".csv")) {
            protein_file <-
              try(read.csv(
                yaml_dat$protein_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          } else if (endsWith(yaml_dat$protein_file, ".txt")) {
            protein_file <-
              try(read.table(
                yaml_dat$protein_file,
                header = TRUE,
                stringsAsFactors = FALSE
              ))
          }
        } else{
          print("The protein file does not exists.")
        }
      } else{
        print("The protein_file parameter is NULL.")
      }
      
      ## Create output folder
      if (!is.null(yaml_dat$output)) {
        output <- file.path(yaml_dat$output)
        if (!dir.exists(output)) {
          dir.create(
            path = output,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The output folder exists.")
        }
      } else{
        print("The output parameter is NULL.")
      }
      
    } else{
      print("The yaml data was not read into the work space!!!")
    }
    
    #######################################################################
    ## Check args and run action based on args
    #######################################################################
    
    # removeDuplicates
    if (all("-removeDuplicates" %in% args)) {
      index <- match("-removeDuplicates", args)
      print(paste(index, ": removeDuplicates", sep = ""))
      
      if (exists("raw_data") &
          exists("by_column") &
          exists("start_column") & dir.exists(output)) {
        folder_path <- file.path(output, "removeDuplicates")
        
        if (!dir.exists(folder_path)) {
          dir.create(
            path = folder_path,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The removeDuplicates folder exists.")
        }
        
        # Using customized function to remove duplicates
        raw_data <-
          remove_duplicates(dat = raw_data, by_column = by_column)
        
        write.csv(
          x = raw_data,
          file = file.path(folder_path, "raw_data.csv"),
          row.names = TRUE,
          na = ""
        )
      }
    }
    
    # outlierRemoval
    if (all("-outlierRemoval" %in% args)) {
      index <- match("-outlierRemoval", args)
      print(paste(index, ": outlierRemoval", sep = ""))
      
      if (exists("raw_data") &
          exists("by_column") &
          exists("start_column") & dir.exists(output)) {
        folder_path <- file.path(output, "outlierRemoval")
        
        if (!dir.exists(folder_path)) {
          dir.create(
            path = folder_path,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The outlierRemoval folder exists.")
        }
        
        # Using customized function to remove outliers
        results <-
          outlier_removal(dat = raw_data,
                          by_column = by_column,
                          start_column = start_column)
        
        if (is.list(results)) {
          raw_data <- results$Outlier_removed_data
          
          capture.output(
            results$Outliers_residuals,
            file = file.path(folder_path, "Outliers_residuals.txt")
          )
          write.csv(
            x = results$Outlier_data,
            file = file.path(folder_path, "Outlier_data.csv"),
            row.names = FALSE,
            na = ""
          )
          write.csv(
            x = results$Outlier_removed_data,
            file = file.path(folder_path, "Outlier_removed_data.csv"),
            row.names = FALSE,
            na = ""
          )
        } else{
          raw_data <- results
          
          write.csv(
            x = results,
            file = file.path(folder_path, "Outlier_removed_data.csv"),
            row.names = FALSE,
            na = ""
          )
          print("No outlier has been found!!!")
        }
        
      }
    }
    
    # boxcoxTransformation
    if (all("-boxcoxTransformation" %in% args)) {
      index <- match("-boxcoxTransformation", args)
      print(paste(index, ": boxcoxTransformation", sep = ""))
      
      if (exists("raw_data") &
          exists("by_column") &
          exists("start_column") & dir.exists(output)) {
        folder_path <- file.path(output, "boxcoxTransformation")
        
        if (!dir.exists(folder_path)) {
          dir.create(
            path = folder_path,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The boxcoxTransformation folder exists.")
        }
        
        # Using customized function to perform box-cox transformation
        results <-
          boxcox_transformation(dat = raw_data,
                                by_column = by_column,
                                start_column = start_column)
        
        if (is.list(results)) {
          raw_data <- results$Boxcox_transformed_data
          
          write.csv(
            x = results$Lambda_values,
            file = file.path(folder_path, "Lambda_values.csv"),
            row.names = FALSE,
            na = ""
          )
          write.csv(
            x = results$Boxcox_transformed_data,
            file = file.path(folder_path, "Boxcox_transformed_data.csv"),
            row.names = FALSE,
            na = ""
          )
        } else{
          raw_data <- results
          
          print("No lambda found!!! Data returned without transformed!!!")
        }
        
      }
    }
    
    # generateBLUP
    if (all("-generateBLUP" %in% args)) {
      index <- match("-generateBLUP", args)
      print(paste(index, ": generateBLUP", sep = ""))
      
      if (exists("raw_data") &
          exists("by_column") &
          exists("start_column") & dir.exists(output)) {
        folder_path <- file.path(output, "generateBLUP")
        
        if (!dir.exists(folder_path)) {
          dir.create(
            path = folder_path,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The generateBLUP folder exists.")
        }
        
        # Using customized function to generate BLUP
        results <-
          generate_BLUP(dat = raw_data,
                        by_column = by_column,
                        start_column = start_column)
        
        if (is.list(results)) {
          BLUP <- results$BLUP
          BLUP_by_column <- 1
          BLUP_start_column <- 2
          
          write.csv(
            x = results$BLUP,
            file = file.path(folder_path, "BLUP.csv"),
            row.names = FALSE,
            na = ""
          )
        } else{
          print("No BLUP generated!!!")
        }
        
      }
    }
    
    # FarmCPU
    if (all("-farmCPU" %in% args)) {
      index <- match("-farmCPU", args)
      print(paste(index, ": farmCPU", sep = ""))
      
      if (exists("BLUP") &
          exists("BLUP_by_column") &
          exists("BLUP_start_column") &
          exists("genotype_file") & exists("SNPs_file") &
          exists("ecotype_file") &
          exists("protein_file") & dir.exists(output)) {
        folder_path <- file.path(output, "FarmCPU")
        
        if (!dir.exists(folder_path)) {
          dir.create(
            path = folder_path,
            showWarnings = TRUE,
            recursive = TRUE
          )
        } else{
          print("The FarmCPU folder exists.")
        }
        
        # Using customized function to generate BLUP
        results <-
          farming_with_FarmCPU(
            dat = BLUP,
            by_column = BLUP_by_column,
            start_column = BLUP_start_column,
            output_path = folder_path,
            genotype_file = genotype_file,
            SNPs_file = SNPs_file,
            ecotype_file = ecotype_file,
            protein_file = protein_file
          )
        
      }
    }
    
    # GAPIT
    if (all("-GAPIT" %in% args)) {
      index <- match("-GAPIT", args)
      print(paste(index, ": GAPIT", sep = ""))
    }
    
  } else{
    print("No action is required to perform!!!")
    
    # Print help file here
  }
  
} else{
  print("YAML file does not exists!!!")
}
