
## farming with FarmCPU
farming_with_FarmCPU <- function(dat, by_column = 1, start_column = 2, output_path, 
                                 p_value_threshold = NA, p_value_fdr_threshold = NA, ld_number = 0, 
                                 genotype, SNPs, hapmap, gff) {
  
  #######################################################################
  ## Create folders to store outputs
  #######################################################################
  
  auto_save_path <- file.path(getwd(), output_path, "farmCPU_auto_output")
  farmCPU_manhattan_plot_save_path <- file.path(getwd(), output_path, "farmCPU_Manhattan_Plot")
  farmCPU_qq_plot_save_path <- file.path(getwd(), output_path, "farmCPU_QQ_Plot")
  farmCPU_significant_save_path <- file.path(getwd(), output_path, "farmCPU_significant")
  ped_and_info_save_path <- file.path(getwd(), output_path, "Haploview_PEDandINFO")
  ld_data_save_path <- file.path(getwd(), output_path, "Haploview_LD_data")
  ld_plot_save_path <- file.path(getwd(), output_path, "Haploview_LD_plot")
  haplotypes_gabriel_blocks_save_path <- file.path(getwd(), output_path, "Haploview_Haplotypes_gabriel_blocks")
  gff_save_path <- file.path(getwd(), output_path, "GFF")
  
  temp <- c(auto_save_path, farmCPU_manhattan_plot_save_path, farmCPU_qq_plot_save_path, farmCPU_significant_save_path, 
            ped_and_info_save_path, ld_data_save_path, ld_plot_save_path, haplotypes_gabriel_blocks_save_path, gff_save_path)
  
  for (i in 1:length(temp)) {
    if (!dir.exists(temp[i])){
      try(dir.create(temp[i]))
    }
  }
  
  #######################################################################
  ## Set the working directory
  #######################################################################
  
  # Get current working directory
  current_working_directory <- getwd()
  
  setwd(auto_save_path)
  
  #######################################################################
  ## farming with FarmCPU and other operations
  #######################################################################
  
  # farming with FarmCPU
  gwas_result_list <- list()
  for (i in start_column:ncol(dat)){
    farmCPU_result <- FarmCPU(
      Y = dat[,c(1,i)],
      GD = genotype,
      GM = SNPs,
      MAF.calculate = TRUE,
      maf.threshold = 0.05,
      method.bin="optimum",
      threshold.output=1
    )
    
    if(file.exists(file.path(paste("FarmCPU.", colnames(dat)[i], ".Manhattan.Plot.Genomewise.pdf", sep = "")))){
      system(paste("convert", file.path(paste("FarmCPU.", colnames(dat)[i], ".Manhattan.Plot.Genomewise.pdf", sep = "")), 
                   file.path(farmCPU_manhattan_plot_save_path, paste("FarmCPU.", colnames(dat)[i], ".Manhattan.Plot.Genomewise.png", sep = "")), sep = " "))
    }
    
    if(file.exists(file.path(paste("FarmCPU.", colnames(dat)[i], ".QQ-Plot.pdf", sep = "")))){
      system(paste("convert", file.path(paste("FarmCPU.", colnames(dat)[i], ".QQ-Plot.pdf", sep = "")),
                   file.path(farmCPU_qq_plot_save_path, paste("FarmCPU.", colnames(dat)[i], ".QQ-Plot.png", sep = "")), sep = " "))
    }
    
    gwas_result <- as.data.frame(farmCPU_result$GWAS)
    gwas_result <- gwas_result[!is.na(gwas_result$P.value),]
    gwas_result$P.value.fdr <- p.adjust(gwas_result$P.value, method = "fdr")
    
    if(!is.na(p_value_threshold)){
      gwas_result <- gwas_result[gwas_result$P.value <= p_value_threshold,]
    } else if(!is.na(p_value_fdr_threshold)){
      gwas_result <- gwas_result[gwas_result$P.value.fdr <= p_value_fdr_threshold,]
    }
    
    if(nrow(gwas_result) == 0){
      gwas_result[1,1] <- NA
    }
    
    gwas_result$Trait <- colnames(dat)[i]
    gwas_result$Method <- "FarmCPU"
    
    gwas_result_list[[colnames(dat)[i]]] <- gwas_result
  }
  
  #######################################################################
  ## Reset the working directory
  #######################################################################
  
  # Set back the working directory
  setwd(current_working_directory)
  
  #######################################################################
  ## Run Haploview here
  #######################################################################
  
  for (i in 1:length(gwas_result_list)) {
    temp_gwas_result <- gwas_result_list[[i]]
    
    temp_gwas_result$LD_number <- ld_number
    temp_gwas_result$LD_start <- temp_gwas_result$Position-ld_number
    temp_gwas_result$LD_end <- temp_gwas_result$Position+ld_number
    temp_gwas_result$Haploblock_number <- NA
    temp_gwas_result$Haploblock_start <- NA
    temp_gwas_result$Haploblock_stop <- NA

    # For each row in each table
    j <- 1
    while(j <= nrow(temp_gwas_result)){

      # Clear temp_hapmap because we might use it later
      temp_hapmap <- data.frame()

      if(!is.na(temp_gwas_result$LD_start[j]) & !is.na(temp_gwas_result$LD_end[j])){

        # Change colnames of Hapmap
        colnames(hapmap)[1] <- "Chromosome"
        colnames(hapmap)[2] <- "Positions"

        # Create a temporary hapmap from chromosome, LD start, and LD stop
        temp_hapmap <- hapmap[(hapmap$Positions >= temp_gwas_result$LD_start[j] & 
                              hapmap$Positions <= temp_gwas_result$LD_end[j] & 
                              hapmap$Chromosome == temp_gwas_result$Chromosome[j]), ]

        # If column names are integer, the system will append "X" in front of each column name
        # This step is to remove the "X" from all integer column names
        if(all(startsWith(colnames(temp_hapmap)[3:ncol(temp_hapmap)], "X"))){
          colnames(temp_hapmap)[3:ncol(temp_hapmap)] <- gsub("[[:alpha:]]", "", colnames(temp_hapmap)[3:ncol(temp_hapmap)])
        }

        if(nrow(temp_hapmap) > 0){
          # Put positions to info file
          info <- data.frame(temp_hapmap$Positions, temp_hapmap$Positions)

          # Remove Chromosome and Positions columns
          temp_hapmap <- temp_hapmap[, c(-1, -2)]

          # Prepare ped file
          ibd <- paste("IBD",1:ncol(temp_hapmap), sep = "")
          ped <- rbind(temp_hapmap[0,], ibd, 0, 0, 7, 1, temp_hapmap[c(1:nrow(temp_hapmap)),])
          ped <- t(ped)

          filename <- paste0(temp_gwas_result$Trait[j], "_", temp_gwas_result$SNP_ID[j], "_", 
                            temp_gwas_result$Chromosome[j], "_", temp_gwas_result$LD_start[j], "-", 
                            temp_gwas_result$LD_end[j])

          # Name ped and info file
          ped_file_name <- paste0(filename, ".ped")
          info_file_name <- paste0(filename, ".info")

          # Save info and ped file
          write.table(info, file.path(ped_and_info_save_path, info_file_name), sep = '\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
          write.table(ped, file.path(ped_and_info_save_path, ped_file_name), sep = '\t', row.names = TRUE, col.names = FALSE, quote=FALSE)
        
          # Prepare the Haploview command
          ld_data_command <- paste("java -jar ", file.path(current_working_directory, "Haploview.jar"), 
                                   " -n -out ", file.path(ld_data_save_path, filename), 
                                   " -pedfile ", file.path(ped_and_info_save_path, ped_file_name), 
                                   " -info ", file.path(ped_and_info_save_path, info_file_name), 
                                   " -skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05", 
                                   sep = "")
          
          # Run Haploview
          system(ld_data_command)

          # Move all the plots into their specific folders
          if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD.PNG", sep = "")))){
            system(paste("mv", file.path(ld_data_save_path, paste(filename, ".LD.PNG", sep = "")), ld_plot_save_path, sep = " "))
          }
          if(file.exists(file.path(ld_data_save_path, paste(filename, ".GABRIELblocks", sep = "")))){
            system(paste("mv", file.path(ld_data_save_path, paste(filename, ".GABRIELblocks", sep = "")), haplotypes_gabriel_blocks_save_path, sep = " "))
          }

          # Read in LD_data and gabriel_block_string
          if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD", sep = "")))){
            LD_data <- try(read.table(file.path(ld_data_save_path, paste(filename, ".LD", sep = "")), check.names = FALSE, header = TRUE))
          }
          if(file.exists(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = "")))){
            gabriel_block_string <- try(readLines(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = ""))))
          }

          if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD", sep = ""))) & 
            file.exists(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = ""))) & 
            !identical(gabriel_block_string, character(0))){

            # Get Haploblock start and stop
            ld <- sort(unique(c(LD_data[,1], LD_data[,2])))

            # Parse gabriel blocks data
            gbb <- list()
            for (k in 1:length(gabriel_block_string)) {
              if(grepl("MARKERS: ", gabriel_block_string[k], ignore.case = TRUE)){
                gbb <- append(gbb, strsplit(x = gsub(".*MARKERS: ", "" , gabriel_block_string[k]), split = " "))
              }
            }

            # Get all the start and stop of markers from ld and gbb
            for(m in 1:length(gbb)){
              for (n in 1:length(gbb[[m]])) {
                gbb[[m]][n] <- as.numeric(ld[as.integer(gbb[[m]][n])])
              }
            }

            # Write the number of haploblocks to the corresponding row of Haploblock_number column
            temp_gwas_result$Haploblock_number[j] <- as.numeric(length(gbb))

            # # Put all the markers start and stop to the gwas results
            # for(m in 1:length(gbb)){
            #   if(m == 1){
            #     temp_gwas_result$Haploblock_start[j] <- as.numeric(gbb[[m]][1])
            #     temp_gwas_result$Haploblock_stop[j] <- as.numeric(gbb[[m]][length(gbb[[m]])])
            #   } else if(m > 1){
            #     temp_gwas_result <- InsertRow(temp_gwas_result, NewRow = temp_gwas_result[j,], RowNum = j+1)
            #     j <- j + 1
            #     temp_gwas_result$Haploblock_start[j] <- as.numeric(gbb[[m]][1])
            #     temp_gwas_result$Haploblock_stop[j] <- as.numeric(gbb[[m]][length(gbb[[m]])])
            #   }
            #   # Remove any row that contains all NA
            #   temp_gwas_result <- temp_gwas_result[rowSums(is.na(temp_gwas_result)) != ncol(temp_gwas_result),]
            # }

            # Write the Haploblock_start and Haploblock_stop that enclose the position
            for(m in 1:length(gbb)){
              if(temp_gwas_result$Position[j] >= as.numeric(gbb[[m]][1]) & temp_gwas_result$Position[j] <= as.numeric(gbb[[m]][length(gbb[[m]])]) ){
                temp_gwas_result$Haploblock_start[j] <- as.numeric(gbb[[m]][1])
                temp_gwas_result$Haploblock_stop[j] <- as.numeric(gbb[[m]][length(gbb[[m]])])
              }
            }

          }
        }
      }

      # Update counter for the while loop
      j <- j + 1
    }

    # Remove any row that contains all NA
    temp_gwas_result <- temp_gwas_result[rowSums(is.na(temp_gwas_result)) != ncol(temp_gwas_result),]

    gwas_result_list[[i]] <- temp_gwas_result
  }

  #######################################################################
  ## Save all GWAS Results
  #######################################################################
  
  for (i in 1:length(gwas_result_list)) {
    gwas_result_filename <- paste("FarmCPU.", names(gwas_result_list)[i], ".GWAS.Results.csv", sep = "")
    write.csv(gwas_result_list[[i]], file.path(farmCPU_significant_save_path, gwas_result_filename), row.names = FALSE)
  }
  
  #######################################################################
  ## Find gene based on Haploblock start and stop or LD start and stop
  #######################################################################

  for (i in 1:length(gwas_result_list)) {
    temp_gwas_result <- gwas_result_list[[i]]

    temp_gwas_result$Gene_name <- NA
    temp_gwas_result$Gene_start <- NA
    temp_gwas_result$Gene_stop <- NA
    temp_gwas_result$Gene_description <- NA

    # For each row in each table
    j <- 1
    while(j <= nrow(temp_gwas_result)){

      temp_gff <- data.frame()

      if(nrow(gff) > 1 & !is.na(temp_gwas_result$Haploblock_start[j]) & !is.na(temp_gwas_result$Haploblock_stop[j])){
        # Match with chromosome
        temp_gff <- gff[as.character(gff[,1]) %in% as.character(temp_gwas_result$Chromosome[j]), ]

        # Match gff with Haploblock start and stop
        if(nrow(temp_gff) > 0){
          temp_gff <- temp_gff[(
            temp_gff[,4] < temp_gwas_result$Haploblock_start[j] & 
            temp_gff[,4] < temp_gwas_result$Haploblock_stop[j] & 
            temp_gff[,5] > temp_gwas_result$Haploblock_start[j]
          ) | (
            temp_gff[,4] < temp_gwas_result$Haploblock_stop[j] & 
            temp_gff[,5] > temp_gwas_result$Haploblock_start[j] & 
            temp_gff[,5] > temp_gwas_result$Haploblock_stop[j]
          ) | (
            temp_gff[,4] > temp_gwas_result$Haploblock_start[j] & 
            temp_gff[,5] < temp_gwas_result$Haploblock_stop[j]
          ) | (
            temp_gff[,4] < temp_gwas_result$Haploblock_start[j] & 
            temp_gff[,4] < temp_gwas_result$Haploblock_stop[j] & 
            temp_gff[,5] > temp_gwas_result$Haploblock_start[j] & 
            temp_gff[,5] > temp_gwas_result$Haploblock_stop[j]
          ),]
        }

      } else if(nrow(gff) > 1 & !is.na(temp_gwas_result$LD_start[j]) & !is.na(temp_gwas_result$LD_end[j])){
        # Match with chromosome
        temp_gff <- gff[as.character(gff[,1]) %in% as.character(temp_gwas_result$Chromosome[j]), ]

        # Match gff with LD start and stop when Haploblock start and stop are NA
        if(nrow(temp_gff) > 0){
          temp_gff <- temp_gff[(
            temp_gff[,4] < temp_gwas_result$LD_start[j] & 
            temp_gff[,4] < temp_gwas_result$LD_end[j] & 
            temp_gff[,5] > temp_gwas_result$LD_start[j]
          ) | (
            temp_gff[,4] < temp_gwas_result$LD_end[j] & 
            temp_gff[,5] > temp_gwas_result$LD_start[j] & 
            temp_gff[,5] > temp_gwas_result$LD_end[j]
          ) | (
            temp_gff[,4] > temp_gwas_result$LD_start[j] & 
            temp_gff[,5] < temp_gwas_result$LD_end[j]
          ) | (
            temp_gff[,4] < temp_gwas_result$LD_start[j] & 
            temp_gff[,4] < temp_gwas_result$LD_end[j] & 
            temp_gff[,5] > temp_gwas_result$LD_start[j] & 
            temp_gff[,5] > temp_gwas_result$LD_end[j]
          ),]
        }

      }

      # If the results after matching have at least 1 row, write all the results to temp_gwas_result
      if(nrow(temp_gff) > 0){
        for(m in 1:nrow(temp_gff)){

          if(m == 1){
            temp_gwas_result$Gene_name[j] <- temp_gff[m,9]
            temp_gwas_result$Gene_start[j] <- temp_gff[m,4]
            temp_gwas_result$Gene_stop[j] <- temp_gff[m,5]
            temp_gwas_result$Gene_description[j] <- temp_gff[m,11]
          } else if(m > 1){
            temp_gwas_result <- InsertRow(temp_gwas_result, NewRow = temp_gwas_result[j,], RowNum = j+1)
            j <- j + 1
            temp_gwas_result$Gene_name[j] <- temp_gff[m,9]
            temp_gwas_result$Gene_start[j] <- temp_gff[m,4]
            temp_gwas_result$Gene_stop[j] <- temp_gff[m,5]
            temp_gwas_result$Gene_description[j] <- temp_gff[m,11]
          }

          # Remove any row that contains all NA
          temp_gwas_result <- temp_gwas_result[rowSums(is.na(temp_gwas_result)) != ncol(temp_gwas_result),]
        }
      }

      # Update the counter
      j <- j + 1
    }

    # Remove any row that contains all NA
    temp_gwas_result <- temp_gwas_result[rowSums(is.na(temp_gwas_result)) != ncol(temp_gwas_result),]

    gwas_result_list[[i]] <- temp_gwas_result
  }

  #######################################################################
  ## Save all GWAS Results
  #######################################################################
  
  for (i in 1:length(gwas_result_list)) {
    gwas_result_filename <- paste("FarmCPU.", names(gwas_result_list)[i], ".GWAS.Results.csv", sep = "")
    write.csv(gwas_result_list[[i]], file.path(farmCPU_significant_save_path, gwas_result_filename), row.names = FALSE)
  }
  
  return(0)
}