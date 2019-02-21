
## farming with FarmCPU
farming_with_FarmCPU <- function(dat, by_column = 1, start_column = 2, output_path, genotype_file, SNPs_file, ecotype_file, protein_file) {
  
  #######################################################################
  ## Create folders to store outputs
  #######################################################################
  
  auto_save_path <- file.path(output_path, "farmCPU_auto_output")
  if (!dir.exists(auto_save_path)){
    try(dir.create(auto_save_path))
  }
  
  farmCPU_manhattan_plot_save_path <- file.path(output_path, "farmCPU_Manhattan_Plot")
  if (!dir.exists(farmCPU_manhattan_plot_save_path)){
    try(dir.create(farmCPU_manhattan_plot_save_path))
  }
  
  farmCPU_qq_plot_save_path <- file.path(output_path, "farmCPU_QQ_Plot")
  if (!dir.exists(farmCPU_qq_plot_save_path)){
    try(dir.create(farmCPU_qq_plot_save_path))
  }
  
  farmCPU_significant_save_path <- file.path(output_path, "farmCPU_significant")
  if (!dir.exists(farmCPU_significant_save_path)){
    try(dir.create(farmCPU_significant_save_path))
  }
  
  ped_and_info_save_path <- file.path(output_path, "Haploview_PEDandINFO")
  if (!dir.exists(ped_and_info_save_path)){
    try(dir.create(ped_and_info_save_path))
  }
  
  ld_data_save_path <- file.path(output_path, "Haploview_LD_data")
  if (!dir.exists(ld_data_save_path)){
    try(dir.create(ld_data_save_path))
  }
  
  ld_plot_save_path <- file.path(output_path, "Haploview_LD_plot")
  if (!dir.exists(ld_plot_save_path)){
    try(dir.create(ld_plot_save_path))
  }
  
  haplotypes_gabriel_blocks_save_path <- file.path(output_path, "Haploview_Haplotypes_gabriel_blocks")
  if (!dir.exists(haplotypes_gabriel_blocks_save_path)){
    try(dir.create(haplotypes_gabriel_blocks_save_path))
  }
  
  gff_save_path <- file.path(output_path, "GFF")
  if (!dir.exists(gff_save_path)){
    try(dir.create(gff_save_path))
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
  
  # farming with FarmCPU here
  
  
  #######################################################################
  ## Reset the working directory
  #######################################################################
  
  # Set back the working directory
  setwd(current_working_directory)
  
  return(0)
}