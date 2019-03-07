
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

    for(j in 1:nrow(temp_gwas_result)){

      if(!is.na(temp_gwas_result$LD_start[j]) & !is.na(temp_gwas_result$LD_end[j])){

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

          ibd <- paste("IBD",1:ncol(temp_hapmap), sep = "")

          ped <- rbind(temp_hapmap[0,], ibd, 0, 0, 7, 1, temp_hapmap[c(1:nrow(temp_hapmap)),])
          ped <- t(ped)

          filename <- paste0(temp_gwas_result$Trait[j], "_", temp_gwas_result$SNP_ID[j], "_", 
                            temp_gwas_result$Chromosome[j], "_", temp_gwas_result$LD_start[j], "-", 
                            temp_gwas_result$LD_end[j])
          
          ped_file_name <- paste0(filename, ".ped")
          info_file_name <- paste0(filename, ".info")

          write.table(info, file.path(ped_and_info_save_path, ped_file_name), sep = '\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
          write.table(ped, file.path(ped_and_info_save_path, info_file_name), sep = '\t', row.names = TRUE, col.names = FALSE, quote=FALSE)
        }

      }
    }

    temp_gwas_result_row_number <- nrow(temp_gwas_result)

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