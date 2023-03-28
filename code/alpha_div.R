# analyze alpha diversity

# load packages

    library("here")
    library("tidyverse")
    library("ggbeeswarm") # for geom_quasirandom

    library("cowplot") # for plot_grid

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    # if (!require("BiocManager", quietly = TRUE))
    #     install.packages("BiocManager")
    # BiocManager::install("DESeq2")
    library("DESeq2")

    #install.packages('lme4')
    library('lme4')

# colors from Wong
# https://www.nature.com/articles/nmeth.1618/figures/2

    colors.wong <- c(
        'Black' = rgb(0, 0, 0, maxColorValue = 255),
        'Orange' = rgb(230, 159, 0, maxColorValue = 255),
        'SkyBlue' = rgb(86, 180, 233, maxColorValue = 255),
        'BluishGreen' = rgb(0, 158, 115, maxColorValue = 255),
        'Yellow' = rgb(240, 228, 66, maxColorValue = 255),
        'Blue' = rgb(0, 114, 178, maxColorValue = 255),
        'Vermillion' = rgb(213, 94, 0, maxColorValue = 255),
        'ReddishPurple' = rgb(204, 121, 167, maxColorValue = 255)
    )

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

# calculate diversity measures for the normalized (=relative abundance) ASV datasets
# and add calculated values back to sample_data for phyloseq object

    # Shannon diversity
    sample_data(imnavait_relative$`16S_ASVs`)$shannon.ASVs <- vegan::diversity(otu_table(imnavait_relative$`16S_ASVs`, index = "shannon", MARGIN = 1)
    sample_data(imnavait_relative$`ITS_ASVs`)$shannon.ASVs <- vegan::diversity(otu_table(imnavait_relative$`ITS_ASVs`, index = "shannon", MARGIN = 1)

    # inverse Simpson diversity
    sample_data(imnavait_relative$`16S_ASVs`)$invsimp.ASVs <- vegan::diversity(otu_table(imnavait_relative$`16S_ASVs`, index = "invsimpson", MARGIN = 1)
    sample_data(imnavait_relative$`ITS_ASVs`)$invsimp.ASVs <- vegan::diversity(otu_table(imnavait_relative$`ITS_ASVs`, index = "invsimpson", MARGIN = 1)

# add read depth to raw and normalized datasets

        sample_data(imnavait$`16S`)$total.reads <- sample_sums(imnavait$`16S`)
        sample_data(imnavait$`ITS`)$total.reads <- sample_sums(imnavait$`ITS`)

        sample_data(imnavait_relative$`16S_ASVs`)$total.reads <- sample_sums(imnavait$`16S`)
        sample_data(imnavait_relative$`ITS_ASVs`)$total.reads <- sample_sums(imnavait$`ITS`)

# extract sample data and diversity results, and convert to tidy dataframe

    shannon_sample_data <- lapply(imnavait_relative[c("16S_ASVs","ITS_ASVs")], function (X) {
        sample_data(X) %>% as_tibble() %>%
            select(-starts_with("invsimp")) %>%
            pivot_longer(starts_with("shannon"), names_to = "Index.Grouping", values_to = "Shannon") %>%
            separate("Index.Grouping", c("Index", "Grouping"), sep = "\\.") %>% select(-Index)
    })

    invsimp_sample_data <- lapply(imnavait_relative[c("16S_ASVs","ITS_ASVs")], function (X) {
        sample_data(X) %>% as_tibble() %>%
            select(-starts_with("shannon")) %>%
            pivot_longer(starts_with("invsimp"), names_to = "Index.Grouping", values_to = "Inverse Simpson") %>%
            separate("Index.Grouping", c("Index", "Grouping"), sep = "\\.") %>% select(-Index)
    })

# Is library size correlated with design variables (sampling month and sampling depth in cm)?
# --> plot reads per sample vs sampling month and reads per sample vs sampling depth in cm}

    # plot with sampling month on x axis

        bind_rows(invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            select(Dataset, Month, total.reads, Depth_range) %>%
            ggplot(aes(x= Month, y = total.reads / 1000, color = Depth_range)) +
                theme_bw() +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                geom_quasirandom(width=0.1) +
                facet_wrap(~ Dataset, nrow=2, scales = "free", strip.position="right") +
                expand_limits(y = 0) +
                theme(legend.position="right") +
                labs(color = "Depth range (cm)", y= "Total reads / 1000", x= "Sampling month")

        ggsave(here("figures", "Figure_A total_reads_vs_month_and_depth.pdf"), width=6, height = 4, useDingbats = FALSE)

    # plot with sampling depth on x axis

        bind_rows(invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            select(Dataset, Month, total.reads, Depth_range) %>%
            ggplot(aes(x= Depth_range, y = total.reads / 1000, color = Month)) +
                theme_bw() +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('June', 'August', 'October'),
                        values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                geom_quasirandom(width=0.1) +
                facet_wrap(~ Dataset, nrow=2, scales = "free", strip.position="right") +
                expand_limits(y = 0) +
                theme(legend.position="right") +
                labs(color = "Sampling month", y= "Total reads / 1000", x= "Depth range (cm)")

        ggsave(here("figures", "Figure_B total_reads_vs_depth_and_month.pdf"), width=6, height = 4, useDingbats = FALSE)

# How does read depth per sample compare between the 16S and ITS datasets?

    read_counts_wide <- bind_rows(invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            pivot_wider(c(Month, Depth_range, Replicate), values_from = total.reads, names_from = "Dataset")

    read_counts_wide %>%
        ggplot(aes(x=`16S`/1000, y=`ITS`/1000, color = Depth_range, shape = Month)) +
            theme_bw() +
            geom_point() +
            labs(x = "16S total reads / 1000", y = "ITS total reads / 1000", color = "Sampling depth (cm)", shape = "Sampling month")

    # correlation between datasets is borderline significant (r = 0.2840934, t = 1.943, df = 43, p-value = 0.05858)
    cor.test(read_counts_wide$`16S`, read_counts_wide$ITS)

# Is alpha diversity being driven by differences in library size (despite total sum scaling)?
# --> plot alpha diversity vs total reads

    # Shannon diversity

        bind_rows(shannon_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            ggplot(aes(total.reads / 1000, `Shannon`, shape = Month, color = Depth_range)) +
                theme_bw() +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                geom_point() +
                facet_wrap(~ Dataset, nrow=2, scales = "free", strip.position="right") +
                expand_limits(y = 0) +
                theme(legend.position="right") +
                labs(color = "Depth range (cm)", x= "Total reads / 1000", y = "Shannon diversity", shape = "Sampling month")

        ggsave(here("figures", "Figure_C alpha_diversity_vs_total_reads_shannon.pdf"), width=6, height = 4, useDingbats = FALSE)

        # plots without Oct 0-20, for comparison

            shannon_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
                filter(!(Month == "October" & Depth_range == "0-20")) %>%
                ggplot(aes(x = total.reads, y = `Shannon`)) + geom_point()

            shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
                ggplot(aes(x = total.reads, y = `Shannon`)) + geom_point()

        # stats without Oct 0-20, for comparison

            shannon_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
                filter(!(Month == "October" & Depth_range == "0-20")) %>%
                lm(Shannon ~ total.reads, data = .) %>%
                summary()

            shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
                lm(Shannon ~ total.reads, data = .) %>%
                summary()

    # inverse Simpson diversity

        bind_rows(invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            ggplot(aes(total.reads / 1000, `Inverse Simpson`, shape = Month, color = Depth_range)) +
                theme_bw() +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                geom_point() +
                facet_wrap(~ Dataset, nrow=2, scales = "free", strip.position="right") +
                expand_limits(y = 0) +
                theme(legend.position="right") +
                labs(color = "Depth range (cm)", x= "Total reads / 1000", y = "Inverse Simpson diversity", shape = "Sampling month")

        ggsave(here("figures", "Figure_D alpha_diversity_vs_total_reads_invsimp.pdf"), width=6, height = 4, useDingbats = FALSE)

        # plots without Oct 0-20, for comparison

            invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
                filter(!(Month == "October" & Depth_range == "0-20")) %>%
                ggplot(aes(x = total.reads, y = `Inverse Simpson`)) + geom_point()

            invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
                ggplot(aes(x = total.reads, y = `Inverse Simpson`)) + geom_point()

        # stats without Oct 0-20, for comparison

            invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
                filter(!(Month == "October" & Depth_range == "0-20")) %>%
                rename(Simpson = `Inverse Simpson`) %>%
                lm(Simpson ~ total.reads, data = .) %>%
                summary()

            invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
                rename(Simpson = `Inverse Simpson`) %>%
                lm(Simpson ~ total.reads, data = .) %>%
                summary()

# test out variance stabilizing transformation (VST) from DeSeq2 as an alternative to total sum scaling

    # prepare data and calculate VST

        imnavait_otutable <- lapply(imnavait, otu_table)
        imnavait_otutable_mtx <- lapply(imnavait_otutable, as, "matrix")
        imnavait_otutable_mtx_t <- lapply(imnavait_otutable_mtx, t)
        imnavait_otutable_mtx_t_1 <- lapply(imnavait_otutable_mtx_t, FUN = function (X) X + 1)
        imnavait_vst <- lapply(imnavait_otutable_mtx_t_1, varianceStabilizingTransformation, fitType="local")
        # add minimum to ensure all values are non-negative after transformation
        imnavait_vst_pos <- list(
            `16S` = sweep(imnavait_vst$`16S`, MARGIN = 2, STATS = apply(imnavait_vst$`16S`, MARGIN = 2, FUN = min), FUN = "-"),
            `ITS` = sweep(imnavait_vst$`ITS`, MARGIN = 2, STATS = apply(imnavait_vst$`ITS`, MARGIN = 2, FUN = min), FUN = "-")
        )
        rm(imnavait_otutable, imnavait_otutable_mtx, imnavait_otutable_mtx_t, imnavait_otutable_mtx_t_1, imnavait_vst)

    # calculate diversity from VST data

        shannon.results.vst <- lapply(imnavait_vst_pos, function (X) {
            vegan::diversity(t(X), index = "shannon", MARGIN = 1)}
        )

        invsimp.results.vst <- lapply(imnavait_vst_pos, function (X) {
            vegan::diversity(t(X), index = "invsimpson", MARGIN = 1)}
        )

    # package up diversity results with per-sample read counts for plotting

        vst_results <- bind_rows(
            full_join(
                data.frame(
                     inverse.simpson = invsimp.results.vst$`16S`,
                     shannon = shannon.results.vst$`16S`) %>% rownames_to_column("SampleID"),
                data.frame(total.reads = sample_sums(imnavait$`16S`)) %>% rownames_to_column("SampleID"),
                by = "SampleID"
            ) %>% full_join(sample_data(imnavait$`16S`) %>% as_tibble() %>% select(-total.reads), by = "SampleID"),
            full_join(
                data.frame(
                     inverse.simpson = invsimp.results.vst$`ITS`,
                     shannon = shannon.results.vst$`ITS`) %>% rownames_to_column("SampleID"),
                data.frame(total.reads = sample_sums(imnavait$`ITS`)) %>% rownames_to_column("SampleID"),
                by = "SampleID"
            ) %>% full_join(sample_data(imnavait$`ITS`) %>% as_tibble() %>% select(-total.reads), by = "SampleID")
        )

    # plot diversity against per-sample read counts

        plot_grid(
            vst_results %>% filter(Dataset == "16S") %>% 
                ggplot(aes(x = total.reads/1000, y = shannon)) +
                theme_bw() + geom_point() +
                labs(x = "Total reads / 1000", y= "Shannon diversity"),
            vst_results %>% filter(Dataset == "16S") %>% 
                ggplot(aes(x = total.reads/1000, y = inverse.simpson)) +
                theme_bw() + geom_point() +
                labs(x = "Total reads / 1000", y= "Inverse Simpson diversity"),
            align = "v"
        )

    ggsave(here("figures", "Figure_E diversity_vst_vs_reads.pdf"), width=6, height = 4, useDingbats = FALSE)

# alpha diversity plots vs sampling depth with regression lines

    # Shannon diversity

        bind_rows(shannon_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            ggplot(aes(x = Depth_min, y = `Shannon`, color = Month)) +
                # change background to white but retain gridlines
                    theme_bw() +
                # draw regression lines with jittered data points
                    geom_smooth(method = "lm", alpha = 0.3, aes(color= Month, fill=Month), size = 0.75) +
                    geom_quasirandom(width = 2) +
                # show sampling depth on vertical axis
                    coord_flip() + 
                    scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20), labels = rev(c("0-20","20-40","40-60","60-80","80-100"))) +
                    expand_limits(y = 0) +
                    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
                # add facets and adjust labelling 
                    labs(color = "Sampling month", fill = "Sampling month") +
                    facet_wrap(~ Dataset, scales = "free", ncol=2, strip.position = "top") +
                    scale_fill_manual(breaks = c('June', 'August', 'October'),
                        values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                    theme(legend.position="bottom", legend.justification=c(1,0))

        ggsave(here("figures", "Figure_5a_shannon_lm.pdf"), width=4, height = 4, useDingbats = FALSE)

    # inverse Simpson diversity

        bind_rows(invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs"),
                  invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs")) %>%
            ggplot(aes(x = Depth_min, y = `Inverse Simpson`, color = Month)) +
                # change background to white but retain gridlines
                    theme_bw() +
                # draw regression lines with jittered data points
                    geom_smooth(method = "lm", alpha = 0.3, aes(color= Month, fill=Month), size = 0.75) +
                    geom_quasirandom(width = 2) +
                # show sampling depth on vertical axis
                    coord_flip() +
                    scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20), labels = rev(c("0-20","20-40","40-60","60-80","80-100"))) +
                    expand_limits(y = 0) +
                    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
                # add facets and adjust labelling 
                    labs(x = "Depth range (cm)", color = "Sampling month", fill = "Sampling month") +
                    facet_wrap(~ Dataset, scales = "free", ncol=2, strip.position = "top") +
                    scale_fill_manual(breaks = c('June', 'August', 'October'),
                        values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                    theme(legend.position="bottom", legend.justification=c(1,0))

        ggsave(here("figures", "Figure_5b_invsimp_lm.pdf"), width=4, height = 4, useDingbats = FALSE)

# ********
# Notes about the linear mixed models below:
#  - Month is a factor, but Depth_min is integer.
#  - For the fixed effects, lmer() fits a slope to Depth_min (potentially a different slope for each month, depending on whether the model is fit with an interaction or not)
#  - For the random effects, (1|Month:Depth_min) tries to fit a different random intercept for each combination of Month and Depth_min, i.e. Depth_min is treated as a factor even if it is supplied as integer. Indeed, if you use as.numeric() to ensure that Depth_min is numeric, lmer() complains that it needs to be a factor.
#  - Note that if you wanted to fit random slopes on Depth_min by month, the syntax would be something like (1 + Depth_min|Month).
# ********

# mixed models for 16S Shannon diversity
# (alpha diversity vs depth and month)

    # prepare model data

        model_data_16S_ASVs <- shannon_sample_data$`16S` %>% filter(Grouping == "ASVs")

    # LR test for month and depth interaction is borderline significant

        model_16S_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        model_16S_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 6.0531  2    0.04848 *

    # LR test for month and depth is significant    

        model_16S_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        model_16S_ASVs_mm_reduced <- lmer(Shannon ~ (1|Month:Depth_min), data = model_data_16S_ASVs)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 18.475  5   0.002406 **

    # LR test for depth is significant

        model_16S_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        model_16S_ASVs_mm_reduced <- lmer(Shannon ~ Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq) 
        # 11.964  3   0.007509 **

    # LR test for month is significant

        model_16S_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        model_16S_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + (1|Month:Depth_min), data = model_data_16S_ASVs)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 14.921  4   0.004869 **

# mixed models for 16S inverse Simpson diversity
# (alpha diversity vs depth and month)

    # prepare model data

        model_data_16S_ASVs <- invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
            rename(invSimp = `Inverse Simpson`) # lmer() doesn't like spaces in name

    # LR test for month and depth interaction is (unsurprisingly) significant

        model_16S_ASVs_mm_full <- lmer(invSimp ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        model_16S_ASVs_mm_reduced <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_16S_ASVs)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 7.9368  2     0.0189 *

    # But is this all just driven by the high diversity in October 0-20 cm samples?

    # remove October 0-20 cm samples and refit models

        model_data_16S_ASVs_noOct020 <- invsimp_sample_data$`16S` %>% filter(Grouping == "ASVs") %>%
            rename(invSimp = `Inverse Simpson`) %>%
            filter(!(Month == "October" & Depth_range == "0-20"))

    # model either with or without interaction gives singular fit
    #  -- without Oct 0-20 there is little variation in across months (or depths)

        # with interaction
        model_16S_ASVs_mm_full <- lmer(invSimp ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_16S_ASVs_noOct020)

        # without interaction
        model_16S_ASVs_mm_full <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_16S_ASVs_noOct020)

    # model with just month also gives singular fit, presumably again reflecting minimal variation among months

        model_16S_ASVs_mm_full <- lmer(invSimp ~ Month + (1|Month:Depth_min), data = model_data_16S_ASVs_noOct020)

    # model with just depth fits, and LR test of depth is non-significant

        model_16S_ASVs_mm_full <- lmer(invSimp ~ Depth_min + (1|Month:Depth_min), data = model_data_16S_ASVs_noOct020)
        model_16S_ASVs_mm_reduced <- lmer(invSimp ~ (1|Month:Depth_min), data = model_data_16S_ASVs_noOct020)
        anova(model_16S_ASVs_mm_full, model_16S_ASVs_mm_reduced)

        # Chisq Df Pr(>Chisq)
        # 0.628  1     0.4281

# mixed models for ITS Shannon diversity
# (alpha diversity vs depth and month)

    # prepare model data

        model_data_ITS_ASVs <- shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs")

    # LR test for month and depth interaction is significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 9.1648  2    0.01023 *

    # LR test for month and depth is significant                                                                                   

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 15.43  5   0.008676 **

    # LR test for depth is significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq) 
        # 13.383  3   0.003877 **

    # LR test for month is significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 11.82  4    0.01874 *

    # But is depth slope still significant without the August samples? ... nope!

    # prepare model data without August

        model_data_ITS_ASVs_noAug <- shannon_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
            filter(Month != "August")

    # LR test for month and depth interaction is non-significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 0.0708  1     0.7901

    # LR test for month and depth is non-significant                                                                                   

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 1.7332  2     0.4204

    # LR test for depth is non-significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq) 
        # 0.1698  1     0.6803

    # LR test for month is non-significant

        model_ITS_ASVs_mm_full <- lmer(Shannon ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(Shannon ~ Depth_min + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 1.5882  1     0.2076

# mixed models for ITS inverse Simpson diversity
# (alpha diversity vs depth and month)

    # prepare model data

        model_data_ITS_ASVs <- invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
            rename(invSimp = `Inverse Simpson`) # lmer() doesn't like spaces in name

    # model with month and depth interaction has singular fit

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)

    # but model without interaction fits OK

    # LR test for month and depth is significant

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(invSimp ~ (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 12.786  3   0.005122 **

    # LR test for depth is significant

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(invSimp ~ Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq) 
        # 8.8983  1   0.002854 **

    # LR test for month is (borderline) significant

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        model_ITS_ASVs_mm_reduced <- lmer(invSimp ~ Depth_min + (1|Month:Depth_min), data = model_data_ITS_ASVs)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 6.4329  2     0.0401 *

    # Is the depth effect just driven by the August samples?

    # prepare model data without August

        model_data_ITS_ASVs_noAug <- invsimp_sample_data$`ITS` %>% filter(Grouping == "ASVs") %>%
            rename(invSimp = `Inverse Simpson`) %>% # lmer() doesn't like spaces in name
            filter(Month != "August")

    # model with month and depth has singular fit, either with or without interaction

        # with interaction
        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min * Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)

        # without interaction
        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min + Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)

    # LR test for month is non-significant (but borderline)

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Month + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(invSimp ~ (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 3.6655  1    0.05555 .

    # LR test for depth is non-significant

        model_ITS_ASVs_mm_full <- lmer(invSimp ~ Depth_min + (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        model_ITS_ASVs_mm_reduced <- lmer(invSimp ~ (1|Month:Depth_min), data = model_data_ITS_ASVs_noAug)
        anova(model_ITS_ASVs_mm_full, model_ITS_ASVs_mm_reduced)

        #  Chisq Df Pr(>Chisq)
        # 1.7411  1      0.187
