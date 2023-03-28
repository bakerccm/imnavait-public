# perform differential abundance testing using DeSeq2

# load packages

    library("here")
    library("tidyverse")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    # BiocManager::install("DESeq2")
    library("DESeq2")

    library("writexl")

# get data

    # loads imnavait (original data) and imnavait_relative (normalized data)
    load(here("out", "combined", "imnavait_normalized.rdata"))

# set up color scale for plots

    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }

# get rid of dashes in metadata to avoid DeSeq2 warning

    sample_data(imnavait$`16S`)$Depth_range <- sample_data(imnavait$`16S`)$Depth_range %>%
        as.character() %>% gsub(pattern = "-", replacement = "_", .) %>% as.factor()
    sample_data(imnavait$ITS)$Depth_range <- sample_data(imnavait$ITS)$Depth_range %>%
        as.character() %>% gsub(pattern = "-", replacement = "_", .) %>% as.factor()

# run DESeq
    
    # subset original (un-normalized) phyloseq object to depth and months of interest
    # note: can't use depth_selection and month_selection objects here for some reason
    #    imnavait.subset <- subset_samples(imnavait[[dataset_selection]], Depth_range == depth_selection) %>% subset_samples(., Month %in% month_selection)
    #    imnavait.subset <- subset_samples(imnavait[[dataset_selection]], Depth_range == "80-100") %>% subset_samples(., Month %in% c(month1 = "June", month2 = "October"))
    
    # desired significance level (for use with *adjusted* p values)

        alpha = 0.05

    # convert phyloseq objects for each dataset to DESeqDataSet

        imnavait.DESeqDataSet <- lapply(imnavait, phyloseq_to_deseq2, design = ~ Month + Depth_range + Month:Depth_range)
        
    # add single group factor to use in place of interaction

        for (i in names(imnavait.DESeqDataSet)) {
            imnavait.DESeqDataSet[[i]]$group <- factor(paste(imnavait.DESeqDataSet[[i]]$Month, imnavait.DESeqDataSet[[i]]$Depth_range, sep = "."))
            design(imnavait.DESeqDataSet[[i]]) <- ~ group
        }
        
    # run deseq2
    # Note use of sfType = "poscounts" avoids errors caused by every ASV having at least one zero ("ratio")
    # and with non-convergence ("iterate"). It is shown to work well for zero-inflated NB data in this preprint
    # https://www.biorxiv.org/content/10.1101/157982v1 (see also https://support.bioconductor.org/p/100201/)

        imnavait.deseq2 <- lapply(imnavait.DESeqDataSet, DESeq, test="Wald", fitType="parametric", sfType = "poscounts")
        
    # get names to use with results()

        resultsNames(imnavait.deseq2$`16S`)
        resultsNames(imnavait.deseq2$`ITS`)

    # get pairwise DeSeq2 results for each depth
    # note: creates *adjusted* p-values that are by default generated using BH
        
        imnavait.results <- list()
        
        for (i in c("0_20", "20_40", "40_60", "60_80", "80_100")) {
            imnavait.results[[i]] <- list(
                "OctAug" = lapply(imnavait.deseq2, results, cooksCutoff = FALSE, contrast = c("group", paste0("October.", i), paste0("August.", i))),
                "OctJun" = lapply(imnavait.deseq2, results, cooksCutoff = FALSE, contrast = c("group", paste0("October.", i), paste0("June.", i))),
                "AugJun" = lapply(imnavait.deseq2, results, cooksCutoff = FALSE, contrast = c("group", paste0("August.", i), paste0("June.", i)))
            )
        }

    # get names of those ASVs that are significant in any pairwise comparison (by depth)

        imnavait.sigASVs <- list()
        
        for (i in c("0_20", "20_40", "40_60", "60_80", "80_100")) {
            imnavait.sigASVs[[i]] <- list(
                "OctAug" = lapply(imnavait.results[[i]]$OctAug, FUN = function (X) rownames(X)[which(X$padj < alpha)] ),
                "OctJun" = lapply(imnavait.results[[i]]$OctJun, FUN = function (X) rownames(X)[which(X$padj < alpha)] ),
                "AugJun" = lapply(imnavait.results[[i]]$AugJun, FUN = function (X) rownames(X)[which(X$padj < alpha)] )
            )
            imnavait.sigASVs[[i]]$all = mapply(c, imnavait.sigASVs[[i]]$OctAug, imnavait.sigASVs[[i]]$OctJun, imnavait.sigASVs[[i]]$AugJun)
            imnavait.sigASVs[[i]]$all <- lapply(imnavait.sigASVs[[i]]$all, unique)
            imnavait.sigASVs[[i]]$all <- lapply(imnavait.sigASVs[[i]]$all, sort)
            # imnavait.sigASVs[[i]]$OctAug <- NULL
            # imnavait.sigASVs[[i]]$OctJun <- NULL
            # imnavait.sigASVs[[i]]$AugJun <- NULL
        }


    # filter results to those that are in the significant list
        
        imnavait.sigresults <- list()
        
        for (i in c("0_20", "20_40", "40_60", "60_80", "80_100")) {
            imnavait.sigresults[[i]] <- list()
            for (j in c("OctAug", "OctJun", "AugJun")) {
                imnavait.sigresults[[i]][[j]] <- list()
                for (k in c("16S", "ITS")) {
                    imnavait.sigresults[[i]][[j]][[k]] <- imnavait.results[[i]][[j]][[k]][rownames(imnavait.results[[i]][[j]][[k]]) %in% imnavait.sigASVs[[i]][["all"]][[k]], ] %>%
                        as.data.frame()
                    names(imnavait.sigresults[[i]][[j]][[k]])[names(imnavait.sigresults[[i]][[j]][[k]]) != "baseMean"] = paste(names(imnavait.sigresults[[i]][[j]][[k]])[names(imnavait.sigresults[[i]][[j]][[k]]) != "baseMean"], j, sep = ".")
                    imnavait.sigresults[[i]][[j]][[k]] <- imnavait.sigresults[[i]][[j]][[k]] %>% rownames_to_column("ASV")
                }
            }
            for (k in c("16S", "ITS")) {
                imnavait.sigresults[[i]][[k]] <- imnavait.sigresults[[i]][["OctAug"]][[k]] %>%
                    full_join(imnavait.sigresults[[i]][["OctJun"]][[k]], by = c("ASV","baseMean")) %>%
                    full_join(imnavait.sigresults[[i]][["AugJun"]][[k]], by = c("ASV","baseMean"))
                imnavait.sigresults[[i]][[k]] <- imnavait.sigresults[[i]][[k]] %>%
                    cbind(as(tax_table(imnavait[[k]])[imnavait.sigresults[[i]][[k]]$ASV, ], "matrix"))
            }
        }

# write results to Excel file for supplementary material

    write_xlsx(
        x = list(
            "16S 0-20" = imnavait.sigresults$`0_20`$`16S`,
            "16S 20-40" = imnavait.sigresults$`20_40`$`16S`,
            "16S 40-60" = imnavait.sigresults$`40_60`$`16S`,
            "16S 60-80" = imnavait.sigresults$`60_80`$`16S`,
            "16S 80-100" = imnavait.sigresults$`80_100`$`16S`,
            "ITS 0-20" = imnavait.sigresults$`0_20`$`ITS`,
            "ITS 20-40" = imnavait.sigresults$`20_40`$`ITS`,
            "ITS 40-60" = imnavait.sigresults$`40_60`$`ITS`,
            "ITS 60-80" = imnavait.sigresults$`60_80`$`ITS`,
            "ITS 80-100" = imnavait.sigresults$`80_100`$`ITS`
        ), path = here("out", "combined", "deseq2.xlsx")
    )

# save results to rds file for re-use

    saveRDS(deseq2.signif, file = here("out","combined","deseq2.rds"))
