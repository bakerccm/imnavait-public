# analyze soil nutrients and filtered porewater chemistry data

# load packages

    library("here")

    #install.packages('lme4')
    library('lme4')

    library("vegan")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    library("tidyverse")

    library("readxl")
    library("scales")
    library("purrr")

    library("ggrepel")
    library("cowplot")

# replacements for months

    SP.months <- c("SP1" = "June", "SP2" = "August", "SP3" = "October")

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

    # loads soil_nutrients
    load(here("data", "chemistry", "soil_nutrients.rdata"))

    # loads porewater_chemistry_data and porewater_technical_data
    load(here("data", "chemistry", "porewater_chemistry.rdata"))

# clean up and reformat porewater chemistry data

    # attach detection and reporting limit information

        porewater_chemistry_data <- porewater_chemistry_data %>% left_join(porewater_technical_data, by = "element")

    # convert ND's to NA's but save original ND info in separate column

        porewater_chemistry_data <- porewater_chemistry_data %>%
            mutate(ND = ifelse(concentration=="ND", TRUE, FALSE)) %>%
            mutate(concentration = ifelse(ND, NA, concentration)) %>%
            mutate(concentration = as.numeric(concentration))

    # extract Depth_min from Depth_range

        porewater_chemistry_data <- porewater_chemistry_data %>%
            separate(Depth_range, c("Depth_min", "Depth_max"), sep = "-", remove = FALSE) %>%
            mutate(Depth_min = as.integer(Depth_min), Depth_max = as.integer(Depth_max))

    # reorder columns

        porewater_chemistry_data <- porewater_chemistry_data %>%
            select(Month, Depth_min, Depth_range, element, element_symbol, concentration, ND, everything())

# inspect porewater chemistry data

    # 23 metals were assayed, but some depth/month combinations have missing data
    
        porewater_chemistry_data %>% select(Depth_range, Month) %>% table()

    # no reported concentrations are less than the detection limit (as expected, since presumably these 
    # are shown as ND), but a small number are below the reporting limit (8, vs 184 above the reporting limit)

        porewater_chemistry_data %>%
            # restrict to data not shown as ND
                filter(!(ND)) %>%
            # are any measurements less than either detection or reporting limit?
                mutate(`conc<DL` = ifelse(concentration < detection_limit, TRUE, FALSE)) %>%
                mutate(`conc<RL` = ifelse(concentration < reporting_limit, TRUE, FALSE)) %>%
                select(-units) %>%
            # summarize results
                select(`conc<DL`, `conc<RL`) %>%
                table()

    # which elements have lots of ND measurements?
    # each element may have up to 12 measurements (see table above)

        porewater_chemistry_data %>% select(element, ND) %>% table()

    # how many measurements above reporting limit does each element have?

        porewater_chemistry_data %>%
            # which measurements are not 'ND' and also above reporting limit?
                mutate(`conc>RL` = !(ND) & (concentration > reporting_limit)) %>%
            # summarize results
                select(`element`, `conc>RL`) %>%
                table()
    
    # which elements have all measurements above reporting limit?
    
        elements.above.reporting.limit <- porewater_chemistry_data %>%
            # which measurements are not 'ND' and also above reporting limit?
                mutate(`conc>RL` = !(ND) & (concentration > reporting_limit)) %>%
            # are all measurements above reporting limit?
                group_by(`element`, element_symbol) %>%
                summarize(`all.conc>RL` = (sum(`conc>RL`) == n()), .groups = "drop") %>%
            # remove elements with measurements below reporting limit
                filter(`all.conc>RL`)
        
        elements.above.reporting.limit

# plot filtered porewater data by month and depth
# N.B. only uses elements with all measurements above reporting limit,
# so all data point plotted are reliable (i.e. above reporting limit) in this graph

    # prepare plot data

        porewater_plot_data <- porewater_chemistry_data %>%
            # remove elements with measurements below reporting limit
                filter(element_symbol %in% elements.above.reporting.limit$element_symbol) %>%
            # clean up element labels
                separate(element, into = "element", sep = " ", extra = "drop") %>%
            # clean up sampling month
                mutate(Month = factor(Month, levels = c("June","August","October"))) %>%
                rename(`Sampling month` = Month)

    # draw plot and save to file
        
        porewater_plot_data %>%
            ggplot(aes(x = Depth_min, y = concentration, group = `Sampling month`, col = `Sampling month`)) +
                # change background to white but retain gridlines
                    theme_bw() +
                geom_point() + geom_line(alpha = 0.3) +
                facet_wrap(~element, scales = "free_x", ncol = 3) +
                labs(x = "Depth range (cm)", y = "Concentration (mg/L)", color = "Sampling month") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                scale_color_manual(breaks = c('June', 'August', 'October'),
                    values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                theme(legend.position = "bottom") + 
                # ggtitle("Porewater chemistry (filtered samples)") +
                coord_flip() +
                scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20), labels = rev(c("0-20","20-40","40-60","60-80","80-100")))

        ggsave(here::here("figures", "Figure_4 pore water chemistry.pdf"), width = 6, height = 9, useDingbats = FALSE)

# plot soil nutrient data by depth and month

    # draw panel with C, N and LOI (these are measured in percent mass)
    # gives a warning because October 0-20 cm only has one LOI measurement
    
        soil.nutrients.plot.CNLOI <- soil_nutrients %>%
            filter(measurement %in% c("C", "N", "LOI")) %>%
            ggplot(aes(x = Depth_min, y = value, group = Month, col = Month)) +
                theme_bw() +
                stat_summary(fun = "mean", geom = "line") +
                geom_point() +
                facet_wrap(~measurement, scales = "free_x", nrow = 1) +
                labs(x = "Depth range (cm)", y = "% mass", color = "Sampling month") +
                scale_color_manual(breaks = c('June', 'August', 'October'),
                    values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(legend.position = "bottom") +
                expand_limits(y = 0) +
                coord_flip() +
                scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20),
                    labels = rev(c("0-20","20-40","40-60","60-80","80-100")))

    # draw panel with pH (unitless measure)
    
        soil.nutrients.plot.pH <- soil_nutrients %>%
            filter(measurement %in% c("pH")) %>%
            ggplot(aes(x = Depth_min, y = value, group = Month, col = Month)) +
                theme_bw() +
                stat_summary(fun = "mean", geom = "line") +
                geom_point() + facet_wrap(~measurement, scales = "free_x", nrow = 1) +
                labs(x = "Depth range (cm)", y = "pH", color = "Sampling month") +
                scale_color_manual(breaks = c('June', 'August', 'October'),
                    values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(legend.position = "bottom") +
                coord_flip() +
                scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20),
                    labels = rev(c("0-20","20-40","40-60","60-80","80-100")))

    # draw panel with NH4-N (ppm)
    
        soil.nutrients.plot.NH4 <- soil_nutrients %>%
            filter(measurement %in% c("NH4-N")) %>%
            ggplot(aes(x = Depth_min, y = value, group = Month, col = Month)) +
                theme_bw() +
                stat_summary(fun = "mean", geom = "line") +
                geom_point() + facet_wrap(~measurement, scales = "free_x", nrow = 1) +
                labs(x = "Depth range (cm)", y = "ppm", color = "Sampling month") +
                scale_color_manual(breaks = c('June', 'August', 'October'),
                    values = as.character(colors.wong[c('Orange', 'SkyBlue', 'BluishGreen')])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(legend.position = "bottom") +
                coord_flip() +
                scale_x_reverse(name = "Depth range (cm)", breaks = seq(80, 0, -20),
                    labels = rev(c("0-20","20-40","40-60","60-80","80-100")))

    # arrange panels

        soil.nutrients.top_row <- soil.nutrients.plot.CNLOI + theme(legend.position = "none")

        soil.nutrients.bottom_row <- plot_grid(
            soil.nutrients.plot.pH + theme(legend.position = "none"),
            soil.nutrients.plot.NH4 + theme(legend.position = "none"),
            labels = c('(b)', '(c)'), label_size = 12, align = "h"
        )

        soil.nutrients.legend <- get_legend(
          soil.nutrients.plot.CNLOI + 
            guides(color = guide_legend(nrow = 1)) +
            theme(legend.position = "bottom")
        )

        soil.nutrients.plot.grid <- plot_grid(
            soil.nutrients.top_row,
            soil.nutrients.bottom_row, labels = c('(a)', ''), label_size = 12, ncol = 1
        )

        soil.nutrients.plot <- plot_grid(soil.nutrients.plot.grid, soil.nutrients.legend, ncol = 1, rel_heights = c(1, .1))

    # save combined plot to file
    
        ggsave(plot = soil.nutrients.plot, here("figures", "Figure_3 soil CN with points.pdf"), width=6, height = 8, useDingbats = FALSE)

# test for relationships between soil nutrients and depth/month

    # Notes:
    # - For the random effects, (1|Month:Depth_min) tries to fit a different random intercept for each combination of
    #   Month and Depth_min, i.e. Depth_min is treated as a factor even if it is supplied as integer. Indeed, if you
    #   use as.numeric() to ensure that Depth_min is numeric, lmer() complains that it needs to be a factor.
    # - If you wanted to fit random slopes on Depth_min by month, the syntax would be something like (1 + Depth_min|Month).

    # carbon

        model_C_full <- lmer(value ~ Depth_min + Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "C"))
            
        model_C_reduced <- lmer(value ~ Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "C"))
            
        # full model is better, i.e. negative slope on Depth_min is significant
        #     chisq = 4.8734, df = 1, p = 0.02727
        
            anova(model_C_full, model_C_reduced)
            summary(model_C_full)

    # nitrogen
    
        model_N_full <- lmer(value ~ Depth_min + Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "N"))
            
        model_N_reduced <- lmer(value ~ Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "N"))
            
        # full model is not quite significant at alpha = 0.05 but it is borderline, i.e. negative slope on Depth_min is borderline significant
        #     chisq = 3.5859, df = 1, p = 0.05827 
        
            anova(model_N_full, model_N_reduced)
            summary(model_N_full)

    # LOI
    
        model_LOI_full <- lmer(value ~ Depth_min + Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "LOI"))
            
        model_LOI_reduced <- lmer(value ~ Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "LOI"))
            
        # full model is better, i.e. negative slope on Depth_min is significant
        #     chisq = 4.0548, df = 1, p = 0.04405
        
            anova(model_LOI_full, model_LOI_reduced)
            summary(model_LOI_full)

    # pH
    
        model_pH_full <- lmer(value ~ Depth_min + Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "pH"))
            
        model_pH_reduced <- lmer(value ~ Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "pH"))
            
        # full model is clearly not significantly better than reduced model. i.e. no support for negative slope on Depth_min
        #     chisq = 0.0577, df = 1, p = 0.8102
        
            anova(model_pH_full, model_pH_reduced)
            summary(model_pH_full)

    # NH4-N
    
        model_NH4N_full <- lmer(value ~ Depth_min + Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "NH4-N") %>% filter(value < 2000)) # note remove outlier
            
        model_NH4N_reduced <- lmer(value ~ Month + (1|Month:Depth_min),
            data = soil_nutrients %>% filter(measurement == "NH4-N") %>% filter(value < 2000)) # note remove outlier
            
        # full model is significantly better than reduced model. i.e. support for POSITIVE slope on Depth_min
        #     chisq = 8.6343 , df = 1, p = 0.003299
        
            anova(model_NH4N_full, model_NH4N_reduced)
            summary(model_NH4N_full)

# test for correlations among soil nutrients

    cor.test.data <- list()

    # correlations among C, N and LOI are all significant

        # prepare data
        cor.test.data$CNLOI <- soil_nutrients %>%
            filter(measurement %in% c("C","N", "LOI")) %>%
            pivot_wider(names_from = "measurement", values_from = "value")

        # C vs N
        cor.test.data$CNLOI %>%
            summarize(pearson = cor(C,N)) # r = 0.8813631
        cor.test(x=cor.test.data$CNLOI$C, y=cor.test.data$CNLOI$N) # r = 0.8813631, t = 9.8716, df = 28, p-value = 1.279e-10
            
        # C vs LOI
        cor.test.data$CNLOI %>% 
            select(C,LOI) %>% filter(complete.cases(.)) %>%
            summarize(pearson = cor(C,LOI)) # r = 0.7770625
        cor.test(x=cor.test.data$CNLOI$C, y=cor.test.data$CNLOI$LOI, na.rm = TRUE) # r = 0.7770625, t = 6.415, df = 27, p-value = 7.143e-07
            
        # N vs LOI
        cor.test.data$CNLOI %>%
            select(N,LOI) %>% filter(complete.cases(.)) %>%
            summarize(pearson = cor(N,LOI)) # r = 0.7440147
        cor.test(x=cor.test.data$CNLOI$N, y=cor.test.data$CNLOI$LOI, na.rm = TRUE) # r = 0.7440147, t = 5.786, df = 27, p-value = 3.725e-06

    # correlation between pH and NH4-N is borderline without outlier

        # prepare data
        cor.test.data$pH.NH4N <- soil_nutrients %>%
            filter(measurement %in% c("pH","NH4-N")) %>%
            pivot_wider(names_from = "measurement", values_from = "value")
        
        # pH vs NH4-N (full dataset)
        cor.test.data$pH.NH4N %>%
            summarize(pearson = cor(`pH`,`NH4-N`))
        # r = 0.6977814, t = 5.1546, df = 28, p-value = 1.817e-05
        cor.test(x=cor.test.data$pH.NH4N$pH, y=cor.test.data$pH.NH4N$`NH4-N`)

        # but looks like it's all driven by one technical replicate
        cor.test.data$pH.NH4N %>%
            ggplot(aes(x=pH, y=`NH4-N`)) + geom_point()

        # pH vs NH4-N (without outlier)
        cor.test.data$pH.NH4N.nooutlier <- cor.test.data$pH.NH4N %>%
            filter(!(Month == "August" & Depth_range == "0-20"))
        # r = 0.3617501, t = 1.9786, df = 26, p-value = 0.05855
        cor.test(x=cor.test.data$pH.NH4N.nooutlier$pH, y=cor.test.data$pH.NH4N.nooutlier$`NH4-N`) 

        cor.test.data$pH.NH4N.nooutlier %>%
            ggplot(aes(x=pH, y=`NH4-N`)) + geom_point()

# dbRDA plots


# clean up and reformat porewater chemistry data

    porewater_df <- 


# reformat soil nutrient data as wide

    # prepare data for dbRDA plots
    # -- scale soil nutrient and porewater chemistry data

        # soil nutrient data 

            soil_nutrients_scaled <- soil_nutrients %>%
                # reformat as wide
                    pivot_wider(id_cols = c(Month, Depth_range, Depth_min, Depth_max, replicate), names_from = measurement, values_from = value) %>%
                # calculate mean of technical replicates
                    select(-replicate) %>%
                    group_by(Month, Depth_range, Depth_min, Depth_max) %>%
                    dplyr::summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
                # scale and center data for each column
                    mutate(across(c(-Month, -Depth_range, -Depth_min, -Depth_max), scale)) %>%
                # remove OM since it's just calculated as OM = (2 x C)
                    select(-OM)

        # porewater chemistry data

            porewater_scaled <- porewater_chemistry_data %>%
                # only keep elements with reliable data (see analysis above)
                    filter(element_symbol %in% elements.above.reporting.limit$element_symbol) %>%
                # now that concentration column is only numbers, convert to double
                    mutate(concentration = as.double(concentration)) %>%
                # reformat as wide
                    pivot_wider(id_cols = c(Month, Depth_range), names_from = element_symbol, values_from = concentration) %>%
                # reformat month as factor
                    mutate(Month = factor(Month, levels = SP.months)) %>%
                # scale and center data for each column
                    mutate(across(c(-Month, -Depth_range), scale))

        # merge scaled environmental data 

            env_data_scaled <- full_join(soil_nutrients_scaled, porewater_scaled, by = c("Month", "Depth_range"))

    # prepare data for dbRDA plots
    # -- add scaled environmental data to metadata of new phyloseq object

        # extract 16S and ITS phyloseq objects
        
            imnavait_rel_env <- imnavait_relative[c("16S_ASVs","ITS_ASVs")]

        # add metadata to 16S phyloseq object
        
            sample_data(imnavait_rel_env$`16S_ASVs`) <- sample_data(imnavait_rel_env$`16S_ASVs`) %>% as("data.frame") %>%
                # add environmental data
                    full_join(env_data_scaled, by = c("Month", "Depth_range", "Depth_min", "Depth_max")) %>%
                # add rownames so it can be re-imported to phyloseq
                    column_to_rownames("SampleID") %>% mutate(SampleID = rownames(.)) %>%
                # reformat ready for import to phyloseq
                    select(SampleID, everything()) %>% sample_data()

        # add metadata to ITS phyloseq object
        
            sample_data(imnavait_rel_env$`ITS_ASVs`) <- sample_data(imnavait_rel_env$`ITS_ASVs`) %>% as("data.frame") %>%
                # add environmental data
                    full_join(env_data_scaled, by = c("Month", "Depth_range", "Depth_min", "Depth_max")) %>%
                # add rownames so it can be re-imported to phyloseq
                    column_to_rownames("SampleID") %>% mutate(SampleID = rownames(.)) %>%
                # reformat ready for import to phyloseq
                    select(SampleID, everything()) %>% sample_data()

    # prepare dbRDA ordinations
    # -- note that we can do two analyses:
    #    (i) soil nutrients only, with all depths/months
    #   (ii) soil nutrients + porewater chemistry, but limited to depths/months with porewater data

        dbRDAs <- list(
            # all months
                `16S_ASVs_CN` = imnavait_rel_env$`16S_ASVs` %>% ordinate(method = "CAP", distance = "bray", formula = ~ C + N + LOI + pH + NH4.N),
                `ITS_ASVs_CN` = imnavait_rel_env$`ITS_ASVs` %>% ordinate(method = "CAP", distance = "bray", formula = ~ C + N + LOI + pH + NH4.N),
            # some depths/months are missing porewater data
                `16S_ASVs_all` = imnavait_rel_env$`16S_ASVs` %>%
                    subset_samples(complete.cases(sample_data(imnavait_rel_env$`16S_ASVs`)
                        [,c("C", "N", "LOI", "pH", "NH4.N", "Fe", "Al", "As", "Ba", "Cr", "Cu", "Mn", "Ni", "Sr", "Ti", "V", "Zn")])) %>%
                    ordinate(method = "CAP", distance = "bray", formula = ~ C + N + LOI + pH + NH4.N + Fe + Al + As + Ba + Cr + Cu + Mn + Ni + Sr + Ti + V + Zn),
                `ITS_ASVs_all` = imnavait_rel_env$`ITS_ASVs` %>%
                    subset_samples(complete.cases(sample_data(imnavait_rel_env$`ITS_ASVs`)
                        [,c("C", "N", "LOI", "pH", "NH4.N", "Fe", "Al", "As", "Ba", "Cr", "Cu", "Mn", "Ni", "Sr", "Ti", "V", "Zn")])) %>%
                    ordinate(method = "CAP", distance = "bray", formula = ~ C + N + LOI + pH + NH4.N + Fe + Al + As + Ba + Cr + Cu + Mn + Ni + Sr + Ti + V + Zn)
        )

        dbRDA.ordiplots <- lapply(dbRDAs, ordiplot, type = 'n')

        dbRDA.ordiplot.arrows <- lapply(dbRDA.ordiplots, function (X) {
            X[["biplot"]] %>% as.data.frame() %>% rownames_to_column("element")
        })

    # draw dbRDA plots

        dbRDA.ggplots <- list()

        # 16S : soil nutrients only (all depths/months)

            dbRDA.ggplots$`16S_ASVs_CN` <- plot_ordination(imnavait_rel_env$`16S_ASVs`, dbRDAs$`16S_ASVs_CN`, color = "Depth_range", shape = "Month", type="sites") +
                # change background to white but retain gridlines
                    theme_bw() +
                ggtitle("16S ASVs: dbRDA") + geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                labs(color = "Depth range (cm)", fill = "Depth range (cm)", shape = "Sampling month") +
                geom_segment(data = dbRDA.ordiplot.arrows$`16S_ASVs_CN`, aes(x = 0, xend = CAP1*1.5, y=0, yend = CAP2*1.5, shape = NULL),
                    col = "gray60",  arrow = arrow(length = unit(0.1,"in"))) +
                geom_text(data = dbRDA.ordiplot.arrows$`16S_ASVs_CN`, aes(x = CAP1*2, y = CAP2*2, shape = NULL, label = element),
                    col = "gray30", size = 2)  +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                    scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')]))

            ## Figure S5a ##
            ggsave(here("figures","Figure_S5a dbRDA_16S_ASVs_CN.pdf"), plot = dbRDA.ggplots$`16S_ASVs_CN`, width = 6, height = 4.5, useDingbats = FALSE)

            dbRDAs$`16S_ASVs_CN` # constrained axes = 24.5989% of inertia
            anova(dbRDAs$`16S_ASVs_CN`) # F = 2.5174, df = 5,39, p = 0.001

        # ITS : soil nutrients only (all depths/months)

            dbRDA.ggplots$`ITS_ASVs_CN` <- plot_ordination(imnavait_rel_env$`ITS_ASVs`, dbRDAs$`ITS_ASVs_CN`, color = "Depth_range", shape = "Month", type="sites") +
                # change background to white but retain gridlines
                    theme_bw() +
                ggtitle("ITS ASVs: dbRDA") + geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                labs(color = "Depth range (cm)", fill = "Depth range (cm)", shape = "Sampling month") +
                geom_segment(data = dbRDA.ordiplot.arrows$`ITS_ASVs_CN`, aes(x = 0, xend = CAP1*1.5, y=0, yend = CAP2*1.5, shape = NULL),
                    col = "gray60",  arrow = arrow(length = unit(0.1,"in"))) +
                geom_text(data = dbRDA.ordiplot.arrows$`ITS_ASVs_CN`, aes(x = CAP1*2, y = CAP2*2, shape = NULL, label = element),
                    col = "gray30", size = 2)  +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                    scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')]))

            ## Figure S5b ##
            ggsave(here("figures","Figure_S5b dbRDA_ITS_ASVs_CN.pdf"), plot = dbRDA.ggplots$`ITS_ASVs_CN`, width = 6, height = 4.5, useDingbats = FALSE)

            dbRDAs$`ITS_ASVs_CN` # constrained axes = 22.965% of inertia
            anova(dbRDAs$`ITS_ASVs_CN`) # F = 2.2926, df = 5,39, p = 0.001

        # 16S : soil nutrients + porewater chemistry (only some depths/months)

            dbRDA.ggplots$`16S_ASVs_all` <- plot_ordination(imnavait_rel_env$`16S_ASVs`, dbRDAs$`16S_ASVs_all`, color = "Depth_range", shape = "Month", type="sites") +
                # change background to white but retain gridlines
                    theme_bw() +
                ggtitle("16S ASVs: dbRDA") + geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                labs(color = "Depth range (cm)", fill = "Depth range (cm)", shape = "Sampling month") +
                geom_segment(data = dbRDA.ordiplot.arrows$`16S_ASVs_all`, aes(x = 0, xend = CAP1*1.5, y=0, yend = CAP2*1.5, shape = NULL),
                    col = "gray60",  arrow = arrow(length = unit(0.1,"in"))) +
                geom_text(data = dbRDA.ordiplot.arrows$`16S_ASVs_all`, aes(x = CAP1*2, y = CAP2*2, shape = NULL, label = element),
                    col = "gray30", size = 2)  +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                    scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')]))

            ## Figure 7a ##
            ggsave(here("figures","Figure_7a dbRDA_16S_ASVs_all.pdf"), plot = dbRDA.ggplots$`16S_ASVs_all`, width = 6, height = 4.5, useDingbats = FALSE)

            dbRDAs$`16S_ASVs_all` # constrained axes = 68.4309% of inertia
            anova(dbRDAs$`16S_ASVs_all`) # F = 4.6462, df = 11,24, p = 0.001

        # ITS : soil nutrients + porewater chemistry (only some depths/months)

            dbRDA.ggplots$`ITS_ASVs_all` <- plot_ordination(imnavait_rel_env$`ITS_ASVs`, dbRDAs$`ITS_ASVs_all`, color = "Depth_range", shape = "Month", type="sites") +
                # change background to white but retain gridlines
                    theme_bw() +
                ggtitle("ITS ASVs: dbRDA") + geom_polygon(aes(group = Month.Depth_range, fill = Depth_range), alpha = 0.5, color = NA) +
                #xlim(-4,2) +
                labs(color = "Depth range (cm)", fill = "Depth range (cm)", shape = "Sampling month") +
                geom_segment(data = dbRDA.ordiplot.arrows$`ITS_ASVs_all`, aes(x = 0, xend = CAP1*1.5, y=0, yend = CAP2*1.5, shape = NULL),
                    col = "gray60",  arrow = arrow(length = unit(0.1,"in"))) +
                geom_text(data = dbRDA.ordiplot.arrows$`ITS_ASVs_all`, aes(x = CAP1*2, y = CAP2*2, shape = NULL, label = element),
                    col = "gray30", size = 2)  +
                # use selected colorblind-safe colors
                    scale_color_manual(breaks = c('0-20', '20-40', '40-60', '60-80', '80-100'),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')])) +
                    scale_fill_manual(breaks = as.character(c('0-20', '20-40', '40-60', '60-80', '80-100')),
                        values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red')]))

            ## Figure 7b ##
            ggsave(here("figures","Figure_7b dbRDA_ITS_ASVs_all.pdf"), plot = dbRDA.ggplots$`ITS_ASVs_all`, width = 6, height = 4.5, useDingbats = FALSE)

            dbRDAs$`ITS_ASVs_all` # constrained axes = 53.1964% of inertia
            anova(dbRDAs$`ITS_ASVs_all`) # F = 2.4369, df = 11,24, p = 0.001
