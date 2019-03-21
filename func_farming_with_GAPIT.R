## farming with GAPIT
farming_with_GAPIT <- function(dat, by_column = 1, start_column = 2, output_path, 
                                model, p_value_fdr_threshold = NA, ld_number = 0, 
                                hapmap_numeric, gff) {

    #######################################################################
    ## Create folders to store outputs
    #######################################################################

    auto_save_path <- file.path(getwd(), output_path, "GAPIT_auto_output")
    farmCPU_manhattan_plot_save_path <- file.path(getwd(), output_path, "GAPIT_Manhattan_Plot")
    farmCPU_qq_plot_save_path <- file.path(getwd(), output_path, "GAPIT_QQ_Plot")
    farmCPU_significant_save_path <- file.path(getwd(), output_path, "GAPIT_significant")
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
    ## farming with GAPIT and other operations
    #######################################################################
    # gapit_result <- GAPIT(
    #     Y = myY[,c(1,2)],
    #     KI = NULL,
    #     CV = NULL,
    #     G = NULL,
    #     GD = myGD,
    #     GM = myGM,
    #     model = "GLM",
    #     PCA.total = 5,
    #     file.output = FALSE
    # )

    #######################################################################
    ## Reset the working directory
    #######################################################################

    # Set back the working directory
    setwd(current_working_directory)

    #######################################################################
    ## Run Haploview here
    #######################################################################

    # Some code here

    #######################################################################
    ## Save all GWAS Results
    #######################################################################

    # Some code here

    #######################################################################
    ## Find gene based on Haploblock start and stop or LD start and stop
    #######################################################################

    # Some code here

    #######################################################################
    ## Save all GWAS Results
    #######################################################################

    # Some code here
    
    print("In GAPIT")

    return(0)
}