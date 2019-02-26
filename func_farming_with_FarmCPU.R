
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