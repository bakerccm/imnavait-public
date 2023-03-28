# analyze beta diversity

# load packages

    library("here")
    library("tidyverse")

    library("cowplot") # to align multiple plots

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    library("vegan")

# Bright colorblind-safe palette from Paul Tol
# https://personal.sron.nl/~pault/
# Figure 1: Bright qualitative colour scheme that is colour-blind safe. The main scheme for lines and their labels.

    colors.tol.bright <- c(
        'blue' = '#4477AA',
        'cyan' = '#EE6677',
        'green' = '#228833',
        'yellow' = '#CCBB44',
        'red' = '#66CCEE',
        'purple' = '#AA3377',
        'grey' = '#BBBBBB'
    )

# get data

    # loads imnavait (original data) and imnavait_relative (normalized data)
    load(here("out", "combined", "imnavait_normalized.rdata"))

# add factor levels to month in metadata

    sample_data(imnavait$`16S`)$Month <- factor(sample_data(imnavait$`16S`)$Month, levels = c("June", "August", "October"))
    sample_data(imnavait$`ITS`)$Month <- factor(sample_data(imnavait$`ITS`)$Month, levels = c("June", "August", "October"))

    sample_data(imnavait_relative$`16S_ASVs`)$Month <- factor(sample_data(imnavait_relative$`16S_ASVs`)$Month, levels = c("June", "August", "October"))
    sample_data(imnavait_relative$`ITS_ASVs`)$Month <- factor(sample_data(imnavait_relative$`ITS_ASVs`)$Month, levels = c("June", "August", "October"))

# add frozen/thawed labels to metadata based on recorded condition (from frost probe)

    # these are almost certainly the same for the two datasets, but do it separately for each in case the ordering of the metadata differs

    for (dataset in c("16S_ASVs", "ITS_ASVs")) {

        sample.condition <- rep("Frozen", nrow(sample_data(imnavait_relative[[dataset]])))

        # June 0-20 = thawed
        sample.condition <- ifelse((sample_data(imnavait_relative[[dataset]])$Month == "June") &
                    (sample_data(imnavait_relative[[dataset]])$Depth_range == "0-20"), "Thawed", sample.condition)

        # Aug 0-20, 20-40, 40-60 = thawed
        sample.condition <- ifelse((sample_data(imnavait_relative[[dataset]])$Month == "August") &
                    (sample_data(imnavait_relative[[dataset]])$Depth_range %in% c("0-20", "20-40", "40-60")), "Thawed", sample.condition)

        # Oct 0-20, 20-40, 40-60 = thawed
        sample.condition <- ifelse((sample_data(imnavait_relative[[dataset]])$Month == "October") &
                    (sample_data(imnavait_relative[[dataset]])$Depth_range %in% c("0-20", "20-40", "40-60")), "Thawed", sample.condition)

        sample_data(imnavait_relative[[dataset]])$sample.condition <- sample.condition

        rm(sample.condition)

    }

# calculate ordinations

    set.seed(123)

    # Bray-Curtis ordinations

        bray_curtis_NMDS <- lapply(imnavait_relative[c('16S_ASVs','ITS_ASVs')], ordinate, method = "NMDS", distance = "bray", trymax = 200, trace = 0)

    # weighted UniFrac ordinations

        # ordinations
        wunifrac_NMDS <- imnavait_relative$`16S_ASVs` %>% ordinate(method = "NMDS", distance = "wunifrac")

        # just the unifrac distances for permanovas etc later
        wunifrac_dists_16S_ASVs <- UniFrac(imnavait_relative$`16S_ASVs`, weighted=TRUE, normalized=TRUE, parallel=FALSE)

        # this gives the same distances as
        #   distance(imnavait_relative$`16S_ASVs`, method="wunifrac", type="samples")
        # and the same distances as used in the ordinations above which are stored in
        #   wunifrac_NMDS$diss

# plot Bray-Curtis NMDS ordinations

    # generate basic plots

        NMDS_plots <- list()

        # Bray-Curtis distances

            for (dataset in c("16S_ASVs", "ITS_ASVs")) {
                NMDS_plots[[dataset]] <- plot_ordination(imnavait_relative[[dataset]], bray_curtis_NMDS[[dataset]],
                        color = "Depth_range", shape = "Month", type = "samples") +
                    theme_bw() +
                    geom_point() + ggtitle(paste0(gsub("_", " ", dataset), " (Bray-Curtis distances)")) +
                    labs(color = "Depth range (cm)", shape = "Sampling month", fill = "Depth range (cm)")
            }


        # weighted UniFrac

            NMDS_plots$wunifrac <- plot_ordination(imnavait_relative$`16S_ASVs`, wunifrac_NMDS, color = "Depth_range", shape = "Month", type="samples") +
                theme_bw() + geom_point() + ggtitle("16S ASVs (weighted UniFrac)") +
                geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                labs(color = "Depth range (cm)", fill = "Depth range (cm)", shape = "Sampling month")

    # generate plots colored by depth and grouped by month/depth

        nmds_monthdepth <- list()

        # 16S Bray-Curtis

            nmds_monthdepth$`16S_ASVs` <- NMDS_plots$`16S_ASVs`  +
                    # use selected colorblind-safe colors
                        scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                        scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                annotate("text", x = 0.5, y = -1.5, label = paste0("stress = ", round(bray_curtis_NMDS$`16S_ASVs`$stress, 2)))

        # ITS Bray-Curtis

            nmds_monthdepth$`ITS_ASVs` <- NMDS_plots$`ITS_ASVs`  +
                    # use selected colorblind-safe colors
                        scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                        scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                annotate("text", x = -0.25, y = -0.7, label = paste0("stress = ", round(bray_curtis_NMDS$`ITS_ASVs`$stress, 2)))

        # 16S weighted UniFrac

            nmds_monthdepth$wunifrac <- NMDS_plots$wunifrac +
                geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                annotate("text", x = -0.03, y = -0.125, label = paste0("stress = ", round(wunifrac_NMDS$stress, 2))) +
                    # use selected colorblind-safe colors
                        scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                        scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                            values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')]))

        # generate composite plots (so that plot areas line up - requires cowplot)

            # with both Bray-Curtis plots
                nmds_monthdepth$composite <- plot_grid(nmds_monthdepth$`16S_ASVs`, nmds_monthdepth$`ITS_ASVs`, align = "v", axis = "lr", ncol=1, rel_heights = c(1,1))

            # with UniFrac for 16S + Bray-Curtis for ITS
            #   nmds_monthdepth$composite <- plot_grid(nmds_monthdepth$`wunifrac`, nmds_monthdepth$`ITS_ASVs`, align = "v", axis = "lr", ncol=1, rel_heights = c(1,1))

        # save plots to file

            # ggsave(plot = nmds_monthdepth$`16S_ASVs`, filename = here("figures","Figure_6a NMDS 16S Bray-Curtis.pdf"), width = 6, height = 4.5, useDingbats = FALSE)
            # ggsave(plot = nmds_monthdepth$`ITS_ASVs`, filename = here("figures","Figure_6b NMDS ITS Bray-Curtis.pdf"), width = 6, height = 4.5, useDingbats = FALSE)

            ggsave(plot = nmds_monthdepth$composite, filename = here("figures","Figure_6ab NMDS plots.pdf"), width = 6, height = 9, useDingbats = FALSE)

            ggsave(plot = nmds_monthdepth$wunifrac, filename = here("figures","Figure_S3 NMDS 16S UniFrac.pdf"), width = 6, height = 4.5, useDingbats = FALSE)

    # generate plots colored by sample condition (from recorded thaw depth) and grouped by month/depth

        nmds_condition <- list()

        # 16S Bray-Curtis

            nmds_condition$`16S_ASVs` <- NMDS_plots$`16S_ASVs` +
                aes(color = sample.condition, fill = sample.condition) +
                geom_polygon(aes(group = Month.Depth_range, fill = sample.condition), alpha = 0.5, color = NA) +
                annotate("text", x = 0.5, y = -1.5, label = paste0("stress = ", round(bray_curtis_NMDS$`16S_ASVs`$stress, 2))) +
                labs(color = "Sample condition", fill = "Sample condition")

        # ITS Bray-Curtis

            nmds_condition$`ITS_ASVs` <- NMDS_plots$`ITS_ASVs` +
                aes(color = sample.condition, fill = sample.condition) +
                geom_polygon(aes(group = Month.Depth_range, fill = sample.condition), alpha = 0.5, color = NA) +
                annotate("text", x = -0.25, y = -0.7, label = paste0("stress = ", round(bray_curtis_NMDS$`ITS_ASVs`$stress, 2))) +
                labs(color = "Sample condition", fill = "Sample condition")

        # 16S weighted UniFrac

            nmds_condition$wunifrac <- NMDS_plots$wunifrac +
                aes(color = sample.condition, fill = sample.condition) +
                geom_polygon(aes(group = Month.Depth_range, fill = sample.condition), alpha = 0.5, color = NA) +
                annotate("text", x = -0.03, y = -0.125, label = paste0("stress = ", round(wunifrac_NMDS$stress, 2))) +
                labs(color = "Sample condition", fill = "Sample condition")

        # generate composite plots (so that plot areas line up - requires cowplot)

            nmds_condition$composite <- plot_grid(nmds_condition$`16S_ASVs`, nmds_condition$`ITS_ASVs`, align = "v", axis = "lr", ncol=1, rel_heights = c(1,1))

            ggsave(plot = nmds_condition$composite, filename = here("figures","Figure_S4ab NMDS plots sample condition.pdf"), width = 6, height = 9, useDingbats = FALSE)

# statistical tests of beta diversity in overall dataset
# -- see https://github.com/rthanis/animal-cdiff/blob/master/analysis/alpha-beta.Rmd

    # 16S Bray-Curtis permanova

        # overall significance

            adonis2(as(otu_table(imnavait_relative$`16S_ASVs`), "matrix") ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"), by = NULL)

            #          Df SumOfSqs      R2      F Pr(>F)
            # Model     6   5.7492 0.48391 5.9384  0.001 ***
            # Residual 38   6.1316 0.51609
            # Total    44  11.8808 1.00000

        # separate lines for month and depth (order doesn't matter, because design is balanced)

            adonis2(as(otu_table(imnavait_relative$`16S_ASVs`), "matrix") ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"))

            #             Df SumOfSqs      R2      F Pr(>F)
            # Month        2   1.0349 0.08711 3.2069  0.002 **
            # Depth_range  4   4.7143 0.39680 7.3042  0.001 ***
            # Residual    38   6.1316 0.51609
            # Total       44  11.8808 1.00000

            adonis2(as(otu_table(imnavait_relative$`16S_ASVs`), "matrix") ~ Depth_range + Month,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"))

            #             Df SumOfSqs      R2      F Pr(>F)    
            # Depth_range  4   4.7143 0.39680 7.3042  0.001 ***
            # Month        2   1.0349 0.08711 3.2069  0.001 ***
            # Residual    38   6.1316 0.51609                  
            # Total       44  11.8808 1.00000                  

        # frozen vs thawed

            adonis2(as(otu_table(imnavait_relative$`16S_ASVs`), "matrix") ~ sample.condition,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"), by = NULL)

            #          Df SumOfSqs      R2      F Pr(>F)
            # Model     1   1.5957 0.13431 6.6712  0.001 ***
            # Residual 43  10.2851 0.86569
            # Total    44  11.8808 1.00000

    # 16S UniFrac permanova

        # overall significance

            adonis2(wunifrac_dists_16S_ASVs ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"), by = NULL)

            #          Df SumOfSqs      R2      F Pr(>F)
            # Model     6  1.02586 0.56926 8.3701  0.001 ***
            # Residual 38  0.77623 0.43074
            # Total    44  1.80209 1.00000

        # separate lines for month and depth (order doesn't matter, because design is balanced)

            adonis2(wunifrac_dists_16S_ASVs ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`16S_ASVs`), "data.frame"))

            #             Df SumOfSqs      R2       F Pr(>F)
            # Month        2  0.15617 0.08666  3.8226  0.002 **
            # Depth_range  4  0.86970 0.48260 10.6439  0.001 ***
            # Residual    38  0.77623 0.43074
            # Total       44  1.80209 1.00000

    # ITS Bray-Curtis permanova

        # overall significance

            adonis2(as(otu_table(imnavait_relative$`ITS_ASVs`), "matrix") ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`ITS_ASVs`), "data.frame"), by = NULL)

            #          Df SumOfSqs     R2      F Pr(>F)
            # Model     6    5.218 0.3247 3.0452  0.001 ***
            # Residual 38   10.852 0.6753
            # Total    44   16.070 1.0000

        # separate lines for month and depth (order doesn't matter, because design is balanced)

            adonis2(as(otu_table(imnavait_relative$`ITS_ASVs`), "matrix") ~ Month + Depth_range,
                data = as(sample_data(imnavait_relative$`ITS_ASVs`), "data.frame"))

            #             Df SumOfSqs     R2      F Pr(>F)
            # Month        2   2.3463 0.1460 4.1078  0.001 ***
            # Depth_range  4   2.8718 0.1787 2.5139  0.002 **
            # Residual    38  10.8523 0.6753
            # Total       44  16.0703 1.0000

        # frozen vs thawed

            adonis2(as(otu_table(imnavait_relative$`ITS_ASVs`), "matrix") ~ sample.condition,
                data = as(sample_data(imnavait_relative$`ITS_ASVs`), "data.frame"), by = NULL)

            #          Df SumOfSqs      R2      F Pr(>F)
            # Model     1    0.939 0.05843 2.6684  0.001 ***
            # Residual 43   15.131 0.94157
            # Total    44   16.070 1.00000

# check beta dispersion assumptions for permanova tests

    # 16S dataset

        bcdists_16S_ASVs <- as(otu_table(imnavait_relative$`16S_ASVs`), "matrix") %>% vegdist(method = "bray", binary = FALSE)

        betadisper_16S_ASVs <- betadisper(bcdists_16S_ASVs,
            group = imnavait_relative$`16S_ASVs` %>% sample_data() %>% pull("Month.Depth_range"))

        permutest(betadisper_16S_ASVs)

        # Permutation test for homogeneity of multivariate dispersions
        # Permutation: free
        # Number of permutations: 999
        # 
        # Response: Distances
        #           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
        # Groups    14 0.24317 0.017369 0.806    999  0.641
        # Residuals 30 0.64654 0.021551                    

    # ITS dataset

        bcdists_ITS_ASVs <- as(otu_table(imnavait_relative$`ITS_ASVs`), "matrix") %>% vegdist(method = "bray", binary = FALSE)

        betadisper_ITS_ASVs <- betadisper(bcdists_ITS_ASVs,
            group = imnavait_relative$`ITS_ASVs` %>% sample_data() %>% pull("Month.Depth_range"))

        permutest(betadisper_ITS_ASVs)

        # Permutation test for homogeneity of multivariate dispersions
        # Permutation: free
        # Number of permutations: 999
        # 
        # Response: Distances
        #           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
        # Groups    14 0.45199 0.032285 1.3188    999  0.252
        # Residuals 30 0.73441 0.024480                     

# permanova results for differences between months at each sampling depth

    # test effect of month at each depth in each dataset

        imnavait_relative_depth <- list()
        imnavait_relative_depth_bray <- list()
        imnavait_relative_depth_permanova <- list()

        for (dataset in c("16S_ASVs", "ITS_ASVs")) {

            imnavait_relative_depth[[dataset]] <- list()
            imnavait_relative_depth_bray[[dataset]] <- list()
            imnavait_relative_depth_permanova[[dataset]] <- list()

            for (depth_range in c("0-20", "20-40", "40-60", "60-80", "80-100")) {

                # subset data
                imnavait_relative_depth[[dataset]][[depth_range]] <- imnavait_relative[[dataset]] %>%
                    subset_samples(Depth_range == depth_range) %>%
                    filter_taxa(function (x) sum(x) > 0, prune = TRUE)

                # calculate pairwise Bray-Curtis dissimilarities among samples
                imnavait_relative_depth_bray[[dataset]][[depth_range]] <- phyloseq::distance(imnavait_relative_depth[[dataset]][[depth_range]], method = "bray")

                # permanova
                set.seed(1234)
                imnavait_relative_depth_permanova[[dataset]][[depth_range]] <- adonis2(imnavait_relative_depth_bray[[dataset]][[depth_range]] ~ Month,
                    data = data.frame(sample_data(imnavait_relative_depth[[dataset]][[depth_range]])),
                    permutations = 9999)

            }

        }

    # summarize test results

        # enframe converts named atomic vectors or lists to one- or two-column data frames
        permanova_summary <- bind_rows(
            "16S" = imnavait_relative_depth_permanova$`16S_ASVs` %>% enframe(name = "Depth_range", value = "permanova"),
            "ITS" = imnavait_relative_depth_permanova$ITS_ASVs %>% enframe(name = "Depth_range", value = "permanova"),
            .id = "dataset"
        )

        # reformat results data frame
        permanova_summary <- permanova_summary %>%
            mutate("SumOfSqs_Month" = map_dbl(permanova, ~.[["SumOfSqs"]][1])) %>%
            mutate("SumOfSqs_Residual" = map_dbl(permanova, ~.[["SumOfSqs"]][2])) %>%
            mutate("SumOfSqs_Total" = map_dbl(permanova, ~.[["SumOfSqs"]][3])) %>%
            mutate("F" = map_dbl(permanova, ~.[["F"]][1])) %>%
            mutate("df1" = map_dbl(permanova, ~.[["Df"]][1])) %>%
            mutate("df2" = map_dbl(permanova, ~.[["Df"]][2])) %>%
            mutate("p.value" = map_dbl(permanova, ~.[["Pr(>F)"]][1])) %>%
            # the false discovery rate corrections don't make a lot of difference with only 5 tests in each dataset
            group_by(dataset) %>%
            mutate("p.value.fdr" = p.adjust(p.value, method = "fdr")) %>% # adjust separately for tests in each dataset
            ungroup() %>%
            mutate("signif" = ifelse(p.value.fdr < 0.05, "*", ""))

    # draw figure summarizing per-depth permanova results

        ss_plot <- permanova_summary %>%
            select(dataset, Depth_range, SumOfSqs_Month, SumOfSqs_Total) %>%
            pivot_longer(starts_with("SumOfSqs"), values_to = "value", names_to = "Sum of Squares") %>%
            mutate(`Sum of Squares` = gsub("SumOfSqs_","", `Sum of Squares`)) %>%
            mutate(`Sum of Squares` = gsub("Month","Sampling month", `Sum of Squares`)) %>%
            rename("Depth range (cm)" = Depth_range) %>%
            ggplot(aes(x = `Depth range (cm)`, y = value, color = `Sum of Squares`)) +
                theme_bw() + geom_point() + geom_path(aes(group = `Sum of Squares`), linetype = "dotted") +
                facet_wrap(~ dataset, nrow = 2) +
                expand_limits(y = 0) + scale_x_discrete(limits=rev) +
                labs(y = "Sums of squares", color = element_blank()) + coord_flip() + theme(legend.position = "bottom")

        r2_plot <- permanova_summary %>%
            mutate(R2 = SumOfSqs_Month / SumOfSqs_Total) %>%
            select(dataset, Depth_range, R2) %>%
            rename("Depth range (cm)" = Depth_range) %>%
            ggplot(aes(x = `Depth range (cm)`, y=R2)) +
                theme_bw() + geom_point() + geom_path(aes(group = dataset), linetype = "dotted") + facet_wrap(~dataset, nrow = 2, strip.position = "right") + expand_limits(y = c(0,0.7)) + scale_x_discrete(limits=rev) +
                labs(y = bquote(Sampling~month~R^2)) + coord_flip()

        # default plot margins are 5.5 pts in all directions
        # theme_grey()$plot.margin
        plot_grid(
            ss_plot + theme(strip.background = element_blank(), strip.text.x = element_blank()),
            r2_plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points")),
            align = "h", axis = "bt", rel_widths = c(1, 0.5))
        ggsave(here("figures","Figure_6c sum of squares.pdf"), width = 6, height = 4.5, useDingbats = FALSE)

# explore alternative distance measures

    # in vegdist(), distance is:
    # Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord".

    # try: canberra, morisita, jaccard, robust.aitchison, chisq

    dist_16S <- otu_table(imnavait_relative$`16S_ASVs`) %>% as("matrix") %>% vegdist(distance = "canberra")

    nmds16S <- metaMDS(dist_16S)

    nmds16S_metadata <- sample_data(imnavait_relative$`16S_ASVs`) %>% as("data.frame")

    identical(rownames(nmds16S$points), rownames(nmds16S_metadata)) # TRUE

    cbind(nmds16S$points, nmds16S_metadata %>% select(Month,Depth_range)) %>%
        ggplot(aes(x=MDS1, y= MDS2, col = Depth_range, shape = Month)) + geom_point()

# mantel tests

    # get bray curtis distances either from ordination object calculated earlier or calculate them again here
    #     bray_curtis_NMDS$`16S_ASVs`$diss
    #     bray_curtis_NMDS$`ITS_ASVs`$diss

        bray_dists <- list(
            `16S_ASVs` = distance(imnavait_relative$`16S_ASVs`, method = "bray", type = "samples"),
            `ITS_ASVs` = distance(imnavait_relative$`ITS_ASVs`, method = "bray", type = "samples")
        )

    # plot distances in each dataset against each other

        plot(bray_dists$`16S_ASVs`, bray_dists$`ITS_ASVs`)

    # mantel test is not significant

        mantel(xdis = bray_dists$`16S_ASVs`, ydis = bray_dists$ITS_ASVs, method = "pearson")

        # Mantel statistic r: 0.001013
        #       Significance: 0.481
