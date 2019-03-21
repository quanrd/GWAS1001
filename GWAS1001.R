## Clean the workspace
rm(list = ls())


# Set repository so that required packages can be downloaded
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


## Import library
library(ape)
library(bigmemory)
library(biganalytics)
library(car)
library(data.table)
library(DataCombine)
library(DBI)
library(EMMREML)
library(genetics)
library(gplots)
library(lme4)
# library(LDheatmap)
library(multtest)
library(scatterplot3d)
library(sparklyr)
library(tidyverse)
library(yaml)

require(compiler)


#Import GAPIT Library
source("http://www.zzlab.net/GAPIT/GAPIT.library.R")

#Import GAPIT
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#source("gapit_functions.txt")

#Import FarmCPU
source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")
#source("FarmCPU_functions.txt")

# Source R files
source("func_read_file.R")
source("func_remove_duplicates.R")
source("func_outlier_removal.R")
source("func_boxcox_transformation.R")
source("func_generate_BLUP.R")
source("func_farming_with_FarmCPU.R")


## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

cat(rep("\n", 2));print("-------------------- GWAS1001 Start --------------------")

# Check if the YAML file path is in the args
if (!identical(args, character(0)) & length(args) > 0 & file.exists(file.path(args[1]))) {

  print("YAML file exists!!!")
  
  # If there are more than one argument
  if (length(args) > 1) {

    print(args[2:length(args)])
    
    ## Import configuration file and read in raw data and reference files
    yaml_dat <- tryCatch({ 
                  read_yaml(file.path(args[1])) 
                }, error = function(e) { 
                  print("The yaml file is invalid!!!")
                  quit(status = -1)
                })
    
    #######################################################################
    ## Read in all file and data that are specified in the YAML file
    #######################################################################
    
    if (exists("yaml_dat")) {
      
      ## Import raw data
      raw_data <- read_file(file_path = yaml_dat$raw_data)
      if (is.null(raw_data)) {
        print("The raw_data parameter is NULL.")
      } else{
        print("raw_data has been loaded into memory.")
      }
      
      if (!is.null(yaml_dat$by_column)) {
        by_column <- yaml_dat$by_column
        print(paste("by_column: ", by_column, sep = ""))
      } else{
        print("The by_column parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$start_column)) {
        start_column <- yaml_dat$start_column
        print(paste("start_column: ", start_column, sep = ""))
      } else{
        print("The start_column parameter is NULL.")
      }
      
      ## Import BULP data
      BLUP <- read_file(file_path = yaml_dat$BLUP)
      if (is.null(BLUP)) {
        print("The BLUP parameter is NULL.")
      } else{
        print("BLUP has been loaded into memory.")
      }
      
      if (!is.null(yaml_dat$BLUP_by_column)) {
        BLUP_by_column <- yaml_dat$BLUP_by_column
        print(paste("BLUP_by_column: ", BLUP_by_column, sep = ""))
      } else{
        print("The BLUP_by_column parameter is NULL.")
      }
      
      if (!is.null(yaml_dat$BLUP_start_column)) {
        BLUP_start_column <- yaml_dat$BLUP_start_column
        print(paste("BLUP_start_column: ", BLUP_start_column, sep = ""))
      } else{
        print("The BLUP_start_column parameter is NULL.")
      }
      
      ## Import GAPIT reference files
      # GAPIT hapmap file
      GAPIT_hapmap <- read_file(file_path = yaml_dat$GAPIT_hapmap)
      if (is.null(GAPIT_hapmap)) {
        print("The GAPIT_hapmap parameter is NULL.")
      } else{
        print("GAPIT_hapmap has been loaded into memory.")
      }

      # GAPIT model
      if (!is.null(yaml_dat$GAPIT_model)) {
        GAPIT_model <- yaml_dat$GAPIT_model
        print(paste("GAPIT_model: ", GAPIT_model, sep = ""))
      } else{
        print("The GAPIT_model parameter is NULL.")
      }

      # GAPIT p.Value FDR threshold
      if (!is.null(yaml_dat$GAPIT_p_value_fdr_threshold)) {
        GAPIT_p_value_fdr_threshold <- as.numeric(yaml_dat$GAPIT_p_value_fdr_threshold)
        print(paste("GAPIT_p_value_fdr_threshold: ", GAPIT_p_value_fdr_threshold, sep = ""))
      } else{
        print("The GAPIT_p_value_fdr_threshold parameter is NULL.")
      }
      
      # GAPIT LD_number
      if (!is.null(yaml_dat$GAPIT_LD_number)) {
        GAPIT_LD_number <- as.numeric(yaml_dat$GAPIT_LD_number)
        print(paste("GAPIT_LD_number: ", GAPIT_LD_number, sep = ""))
      } else{
        print("The GAPIT_LD_number parameter is NULL.")
      }

      ## Import FarmCPU reference files
      # FarmCPU genotype file
      FarmCPU_genotype <- read_file(file_path = yaml_dat$FarmCPU_genotype)
      if (is.null(FarmCPU_genotype)) {
        print("The FarmCPU_genotype parameter is NULL.")
      } else{
        print("FarmCPU_genotype has been loaded into memory.")
      }
      
      # FarmCPU SNPs file
      FarmCPU_SNPs <- read_file(file_path = yaml_dat$FarmCPU_SNPs)
      if (is.null(FarmCPU_SNPs)) {
        print("The FarmCPU_SNPs parameter is NULL.")
      } else{
        print("FarmCPU_SNPs has been loaded into memory.")
      }
      
      # FarmCPU ecotype file
      FarmCPU_ecotype <- read_file(file_path = yaml_dat$FarmCPU_ecotype)
      if (is.null(FarmCPU_ecotype)) {
        print("The FarmCPU_ecotype parameter is NULL.")
      } else{
        print("FarmCPU_ecotype has been loaded into memory.")
      }
      
      # FarmCPU filter threshold p_value_threshold and p_value_fdr_threshold
      if (!is.null(yaml_dat$FarmCPU_p_value_threshold)) {
        FarmCPU_p_value_threshold <- as.numeric(yaml_dat$FarmCPU_p_value_threshold)
        print(paste("FarmCPU_p_value_threshold: ", FarmCPU_p_value_threshold, sep = ""))
      } else{
        print("The FarmCPU_p_value_threshold parameter is NULL.")
      }
      if (!is.null(yaml_dat$FarmCPU_p_value_fdr_threshold)) {
        FarmCPU_p_value_fdr_threshold <- as.numeric(yaml_dat$FarmCPU_p_value_fdr_threshold)
        print(paste("FarmCPU_p_value_fdr_threshold: ", FarmCPU_p_value_fdr_threshold, sep = ""))
      } else{
        print("The FarmCPU_p_value_fdr_threshold parameter is NULL.")
      }
      
      # FarmCPU LD_number
      if (!is.null(yaml_dat$FarmCPU_LD_number)) {
        FarmCPU_LD_number <- as.numeric(yaml_dat$FarmCPU_LD_number)
        print(paste("FarmCPU_LD_number: ", FarmCPU_LD_number, sep = ""))
      } else{
        print("The FarmCPU_LD_number parameter is NULL.")
      }

      ## Haploview
      # Hapmap_numeric file
      hapmap_numeric <- read_file(file_path = yaml_dat$hapmap_numeric)
      if (is.null(hapmap_numeric)) {
        print("The hapmap_numeric parameter is NULL.")
      } else{
        print("hapmap_numeric has been loaded into memory.")
      }
      
      ## Match Gene Start and Gene Stop
      # gff file
      gff <- read_file(file_path = yaml_dat$gff)
      if (is.null(gff)) {
        print("The gff parameter is NULL.")
      } else{
        print("gff has been loaded into memory.")
      }
      
      ## Create output folder
      if (!is.null(yaml_dat$output)) {
        if (!dir.exists(yaml_dat$output)) {
          dir.create(path = yaml_dat$output, showWarnings = TRUE, recursive = TRUE)
          if (dir.exists(yaml_dat$output)) {
            output <- file.path(yaml_dat$output)
          }
        } else{
          output <- file.path(yaml_dat$output)
          print("The output folder exists.")
        }
      } else{
        print("The output parameter is NULL.")
        quit(status = -1)
      }
      
    } else{
      print("The yaml data was not read into the work space!!!")
      quit(status = -1)
    }
    
    #######################################################################
    ## Check args and run action based on args
    #######################################################################
    
    # removeDuplicates
    if (all("-removeDuplicates" %in% args)) {
      index <- match("-removeDuplicates", args)
      print(paste(index, ": removeDuplicates", sep = ""))
      
      if (exists("raw_data") & exists("by_column") & exists("start_column") & dir.exists(output)) {
        
        folder_path <- file.path(output, "removeDuplicates")
        
        if (!dir.exists(folder_path)) {
          dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
        } else{
          print("The removeDuplicates folder exists.")
        }
        
        # Using customized function to remove duplicates
        raw_data <- remove_duplicates(dat = raw_data, by_column = by_column)
        
        write.csv(x = raw_data, file = file.path(folder_path, "raw_data.csv"), row.names = TRUE, na = "")
      }
    }
    
    # outlierRemoval
    if (all("-outlierRemoval" %in% args)) {
      index <- match("-outlierRemoval", args)
      print(paste(index, ": outlierRemoval", sep = ""))
      
      if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {
        
        folder_path <- file.path(output, "outlierRemoval")
        
        if (!dir.exists(folder_path)) {
          dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
        } else{
          print("The outlierRemoval folder exists.")
        }
        
        # Using customized function to remove outliers
        results <- outlier_removal(dat = raw_data, by_column = by_column, start_column = start_column)
        
        if (is.list(results)) {
          raw_data <- results$Outlier_removed_data
          
          capture.output( results$Outliers_residuals, file = file.path(folder_path, "Outliers_residuals.txt"))
          write.csv(x = results$Outlier_data, file = file.path(folder_path, "Outlier_data.csv"), row.names = FALSE, na = "" )
          write.csv( x = results$Outlier_removed_data, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
        } else{
          raw_data <- results
          
          write.csv(x = results, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
          print("No outlier has been found!!!")
        }
        
      }
    }
    
    # boxcoxTransformation
    if (all("-boxcoxTransformation" %in% args)) {
      index <- match("-boxcoxTransformation", args)
      print(paste(index, ": boxcoxTransformation", sep = ""))
      
      if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {
        
        folder_path <- file.path(output, "boxcoxTransformation")
        
        if (!dir.exists(folder_path)) {
          dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
        } else{
          print("The boxcoxTransformation folder exists.")
        }
        
        # Using customized function to perform box-cox transformation
        results <- boxcox_transformation(dat = raw_data, by_column = by_column, start_column = start_column)
        
        if (is.list(results)) {
          raw_data <- results$Boxcox_transformed_data
          
          write.csv(x = results$Lambda_values, file = file.path(folder_path, "Lambda_values.csv"), row.names = FALSE, na = "")
          write.csv(x = results$Boxcox_transformed_data, file = file.path(folder_path, "Boxcox_transformed_data.csv"), row.names = FALSE, na = "")
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
      
      if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {
        
        folder_path <- file.path(output, "generateBLUP")
        
        if (!dir.exists(folder_path)) {
          dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
        } else{
          print("The generateBLUP folder exists.")
        }
        
        # Using customized function to generate BLUP
        results <-
          generate_BLUP(dat = raw_data, by_column = by_column, start_column = start_column)
        
        if (is.list(results)) {
          BLUP <- results$BLUP
          BLUP_by_column <- 1
          BLUP_start_column <- 2
          
          write.csv( x = results$BLUP, file = file.path(folder_path, "BLUP.csv"), row.names = FALSE, na = "" )
        } else{
          print("No BLUP generated!!!")
        }
        
      }
    }
    
    # GAPIT
    if (all("-gapit" %in% args)) {
      index <- match("-gapit", args)
      print(paste(index, ": gapit", sep = ""))

      if (exists("BLUP") & !is.null(BLUP) & exists("BLUP_by_column") & exists("BLUP_start_column") & 
          exists("GAPIT_hapmap") & !is.null(GAPIT_hapmap) & 
          exists("hapmap_numeric") & !is.null(hapmap_numeric) & exists("gff") & !is.null(gff) & 
          exists("GAPIT_LD_number") & GAPIT_LD_number >= 0 & dir.exists(output)) {
            
        folder_path <- file.path(output, "GAPIT")
        
        if (!dir.exists(folder_path)) {
          dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
        } else{
          print("The GAPIT folder exists.")
        }

      }
    }

    # FarmCPU
    if (all("-farmCPU" %in% args)) {
      index <- match("-farmCPU", args)
      print(paste(index, ": farmCPU", sep = ""))
      
      if (exists("BLUP") & !is.null(BLUP) & exists("BLUP_by_column") & exists("BLUP_start_column") & 
          exists("FarmCPU_genotype") & !is.null(FarmCPU_genotype) & exists("FarmCPU_SNPs") & !is.null(FarmCPU_SNPs) & 
          exists("hapmap_numeric") & !is.null(hapmap_numeric) & exists("gff") & !is.null(gff) & 
          exists("FarmCPU_LD_number") & FarmCPU_LD_number >= 0 & dir.exists(output)) {
        
        # Check BLUP file and genotype file
        if (sum(match(BLUP[,1], FarmCPU_genotype[,1]), na.rm = TRUE) == 0) {
          print(paste0("Elements in ", colnames(BLUP)[1], " column of BLUP file does not able to match with any element in ", 
                       colnames(FarmCPU_genotype)[1], " column of genotype file."))
          quit(status = -1)
        }
        
        folder_path <- file.path(output, "FarmCPU")
        
        if (!dir.exists(folder_path)) {
          dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
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
            p_value_threshold = FarmCPU_p_value_threshold, 
            p_value_fdr_threshold = FarmCPU_p_value_fdr_threshold,
            ld_number = FarmCPU_LD_number, 
            genotype = FarmCPU_genotype,
            SNPs = FarmCPU_SNPs,
            hapmap_numeric = hapmap_numeric,
            gff = gff
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

cat(rep("\n", 2));print("-------------------- GWAS1001 Stop --------------------")
