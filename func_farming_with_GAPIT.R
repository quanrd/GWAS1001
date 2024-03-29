## farming with GAPIT
farming_with_GAPIT <- function(dat, 
                                by_column = 1, 
                                start_column = 2, 
                                output_path, 
                                p_value_threshold = NA,
                                p_value_fdr_threshold = NA,
                                ld_number = 0,
                                KI = NULL,
                                CV = NULL,
                                G = NULL,
                                GD = NULL,
                                GM = NULL,
                                file.Ext.G = NULL,
                                file.Ext.GD = NULL,
                                file.Ext.GM = NULL,
                                file.G = NULL,
                                file.GD = NULL,
                                file.GM = NULL,
                                file.path = NULL,
                                file.from = 0,
                                file.to = 0,
                                model = NULL,
                                SNP.MAF = 0,
                                PCA.total = 0,
                                Model.selection = FALSE,
                                SNP.test = FALSE,
                                file.output = FALSE) {

    #######################################################################
    ## Create folders to store outputs
    #######################################################################

    auto_save_path <- file.path(output_path, "GAPIT_auto_output")
    GAPIT_manhattan_plot_save_path <- file.path(output_path, "GAPIT_Manhattan_Plot")
    GAPIT_qq_plot_save_path <- file.path(output_path, "GAPIT_QQ_Plot")
    GAPIT_significant_save_path <- file.path(output_path, "GAPIT_significant")
    ped_and_info_save_path <- file.path(output_path, "Haploview_PEDandINFO")
    ld_data_save_path <- file.path(output_path, "Haploview_LD_data")
    ld_plot_save_path <- file.path(output_path, "Haploview_LD_plot")
    haplotypes_gabriel_blocks_save_path <- file.path(output_path, "Haploview_Haplotypes_gabriel_blocks")
    gff_save_path <- file.path(output_path, "GFF")

    temp <- c(auto_save_path, GAPIT_manhattan_plot_save_path, GAPIT_qq_plot_save_path, GAPIT_significant_save_path, 
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

    # farming with GAPIT
    combined_gwas_result = foreach(i = start_column:ncol(dat), .combine = rbind) %dopar% {
        gapit_result <- GAPIT(
            Y = dat[,c(1,i)],
            KI = KI,
            CV = CV,
            G = G,
            GD = GD,
            GM = GM,
            file.Ext.G = file.Ext.G,
            file.Ext.GD = file.Ext.GD,
            file.Ext.GM = file.Ext.GM,
            file.G = file.G,
            file.GD = file.GD,
            file.GM = file.GM,
            file.path = file.path,
            file.from = file.from,
            file.to = file.to,
            model = model,
            SNP.MAF = SNP.MAF,
            PCA.total = PCA.total,
            group.to = nrow(dat),
            group.from = nrow(dat),
            Model.selection = Model.selection,
            SNP.test = SNP.test,
            file.output = file.output
        )

        if(file.exists(file.path(paste("GAPIT.", model,".", colnames(dat)[i], ".Manhattan.Plot.Genomewise.pdf", sep = "")))){
            system(paste("convert", file.path(paste("GAPIT.", model,".", colnames(dat)[i], ".Manhattan.Plot.Genomewise.pdf", sep = "")),
                            file.path(GAPIT_manhattan_plot_save_path, paste("GAPIT.", model,".", colnames(dat)[i], ".Manhattan.Plot.Genomewise.png", sep = "")), sep = " "))
        }

        if(file.exists(file.path(paste("GAPIT.", model,".", colnames(dat)[i], ".QQ-Plot.pdf", sep = "")))){
            system(paste("convert", file.path(paste("GAPIT.", model,".", colnames(dat)[i], ".QQ-Plot.pdf", sep = "")),
                            file.path(GAPIT_qq_plot_save_path, paste("GAPIT.", model,".", colnames(dat)[i], ".QQ-Plot.png", sep = "")), sep = " "))
        }

        gwas_result <- as.data.frame(gapit_result$GWAS, stringsAsFactors = FALSE)

        # Show the GWAS results to user
        cat(rep("\n", 2))
        print(head(gwas_result))

        if(!("FDR_Adjusted_P-values" %in% colnames(gwas_result))){
            gwas_result <- gwas_result %>%
                            add_column(`FDR_Adjusted_P-values` = p.adjust(gwas_result[,4], method = "fdr"), .after = 4) %>%
                            as.data.frame(stringsAsFactors = FALSE)
        }

        # gwas_result$GWAS is formatted in columns below:
        # SNP, Chromosome, Position, P.value, maf, nobs
        # Rsquare.of.Model.without.SNP, Rsquare.of.Model.with.SNP
        # effect, FDR_Adjusted_P-values

        # Remove all the NAs
        # Prevent chromosome column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,2]),]
        # Prevent position column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,3]),]
        # Prevent P.value column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,4]),]
        # Prevent FDR_Adjusted_P-values column has NA
        gwas_result <- gwas_result[!is.na(gwas_result$`FDR_Adjusted_P-values`),]

        # Filter P-values column data base on the p_value_threshold
        if(!is.na(p_value_threshold)){
            gwas_result <- gwas_result[gwas_result[,4] <= p_value_threshold,]
        }

        # Filter FDR_Adjusted_P-values column data base on the p_value_fdr_threshold
        if(!is.na(p_value_fdr_threshold)){
            gwas_result <- gwas_result[gwas_result$`FDR_Adjusted_P-values` <= p_value_fdr_threshold,]
        }

        # Remove all the NAs
        # Prevent chromosome column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,2]),]
        # Prevent position column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,3]),]
        # Prevent P.value column has NA
        gwas_result <- gwas_result[!is.na(gwas_result[,4]),]
        # Prevent FDR_Adjusted_P-values column has NA
        gwas_result <- gwas_result[!is.na(gwas_result$`FDR_Adjusted_P-values`),]

        # If GWAS result has zero row, add a row of NAs
        if(nrow(gwas_result) == 0){
            gwas_result[1,1] <- NA
        }

        # Add trait and method columns
        gwas_result$Trait <- colnames(dat)[i]
        gwas_result$Method <- model[1]

        # Convert chromosome, position, and p.value columns to numeric
        gwas_result[, 2] <- as.numeric(as.character(gwas_result[, 2]))
        gwas_result[, 3] <- as.double(as.character(gwas_result[, 3]))
        gwas_result[, 4] <- as.double(as.character(gwas_result[, 4]))

        # Get LD number and calculate LD start and stop
        gwas_result$LD_number <- as.numeric(ld_number)
        gwas_result$LD_start <- ifelse(gwas_result[,3] - as.double(ld_number) >= 0, gwas_result[,3] - as.double(ld_number), 0)
        gwas_result$LD_end <- gwas_result[,3] + as.double(ld_number)

        # Re-order the table
        gwas_result <- gwas_result[order(gwas_result$LD_end), ]
        gwas_result <- gwas_result[order(gwas_result$LD_start), ]
        # Order the table base on chromosome number which is on column 2 
        gwas_result <- gwas_result[order(gwas_result[,2]), ]

        # Show the GWAS results to user
        cat(rep("\n", 2))
        print(head(gwas_result))

        return(gwas_result)
    }

    #######################################################################
    ## Reset the working directory
    #######################################################################

    # Set back the working directory
    setwd(current_working_directory)

    #######################################################################
    ## Save all GWAS Results
    #######################################################################

    gwas_result_list = split(x = combined_gwas_result, f = combined_gwas_result$Trait)
    for (i in 1:length(gwas_result_list)) {
        gwas_result_filename <- paste("GAPIT.", names(gwas_result_list)[i], ".GWAS.Results.csv", sep = "")
        write.csv(gwas_result_list[[i]], file.path(GAPIT_significant_save_path, gwas_result_filename), row.names = FALSE)
    }

    combined_gwas_result <- as.data.frame(combined_gwas_result, stringsAsFactors = FALSE)
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_end), ]
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_start), ]
    # Order the table base on chromosome number which is on column 2 
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,2]), ]
    gwas_result_filename <- paste("GAPIT.combined.GWAS.Results.csv", sep = "")
    write.csv(combined_gwas_result, file.path(output_path, gwas_result_filename), row.names = FALSE)

    #######################################################################
    ## Everything is done, return combined_gwas_result
    #######################################################################

    if(!is.null(combined_gwas_result) & is.data.frame(combined_gwas_result) & ncol(combined_gwas_result) >= 10){
        return(combined_gwas_result)
    } else{
        return(NULL)
    }

}
