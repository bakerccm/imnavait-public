# code to analyze thermistor data
# - Figure 2ab thermistor data

# load packages

    library("here")
    library("tidyverse")
    library("cowplot") # to combine panels into a single figure

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

# load data

    # loads imnavait.temps
    load(here("data", "temperature", "imnavait_temps.rdata"))

# plot panels showing thermistor data

    thermistor.panel.data <- imnavait.temps %>% 
        pivot_wider(everything(), names_from = "Replicate", values_from = "Temperature (°C)", names_prefix = "Temperature (°C) #") %>%
        # add mean temperature at each depth
            mutate(`Temperature (°C)` = (`Temperature (°C) #1` + `Temperature (°C) #2`) / 2)

    # panel a : all data
    thermistor.panel.a <- thermistor.panel.data %>%
        ggplot(aes(x = Datetime)) +
            # change background to white but retain gridlines
                theme_bw() +
            # axes
                scale_x_datetime(date_labels = "%e %b %Y", limits = c(as.POSIXct('2019/06/06'), as.POSIXct('2020/06/06'))) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            # rectangle to show location of panel b
                geom_rect(aes(xmin = as.POSIXct('2019/06/06'), xmax = as.POSIXct('2019/10/20'), ymin = -2, ymax = 9),
                    fill = NA, color = "darkgray") +
            # plot data
                geom_ribbon(aes(ymin = `Temperature (°C) #1`, ymax = `Temperature (°C) #2`, fill = `Depth (cm)`),
                    alpha = 0.3, color = NA) +
                geom_line(aes(y = `Temperature (°C)`, color = `Depth (cm)`)) +
                labs(x="Date", y = "Soil temperature (°C)") +
            # use selected colorblind-safe colors
                scale_color_manual(breaks = as.character(seq(20,120,20)),
                    values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red', 'purple')])) +
                scale_fill_manual(breaks = as.character(seq(20,120,20)),
                    values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red', 'purple')]))

    # panel b : zoomed in to show sampling period
    sampling.dates.POSIXct <- as.POSIXct(c('2019/06/06', '2019/08/05', '2019/10/08'))
    sampling.dates.df <- data.frame(x=sampling.dates.POSIXct, xend=sampling.dates.POSIXct, y=c(-2,-1.5,-1), yend=c(0.5, 7.0, 0.5))
    thermistor.panel.b <- thermistor.panel.data %>%
        ggplot(aes(x = Datetime)) +
            # change background to white but retain gridlines
                theme_bw() +
            # axes
                scale_x_datetime(date_labels = "%e %b", limits = c(as.POSIXct('2019/06/06'), as.POSIXct('2019/10/20'))) +
                ylim(-2, 9) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            # lines and arrows to show sampling dates
                geom_segment(data = sampling.dates.df, aes(x = x, xend = xend, y = y, yend = yend), linetype = 3, color = "darkgray") +
                geom_segment(data = sampling.dates.df, aes(x = x, xend = xend, y = yend + 1.1, yend = yend + 0.1), arrow = arrow(length = unit(0.2,"cm")), color = "darkgray") +
            # plot data
                geom_ribbon(aes(ymin = `Temperature (°C) #1`, ymax = `Temperature (°C) #2`, fill = `Depth (cm)`), alpha = 0.3, color = NA) +
                geom_line(aes(y = `Temperature (°C)`, color = `Depth (cm)`)) +
                labs(x="Date", y = "Soil temperature (°C)") +
            # use selected colorblind-safe colors
                scale_color_manual(breaks = as.character(seq(20,120,20)),
                    values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red', 'purple')])) +
                scale_fill_manual(breaks = as.character(seq(20,120,20)),
                    values = as.character(colors.tol.bright[c('blue', 'cyan', 'green', 'yellow', 'red', 'purple')]))

    # combine panels into a single plot
    plot_grid(
        thermistor.panel.a +
            theme(axis.title.x = element_blank(), legend.position = "none"),
        thermistor.panel.b +
            theme(axis.title.x = element_blank(), legend.position = "bottom") +
            guides(color = guide_legend(nrow=1), fill = guide_legend(nrow=1)),
        align = "v", axis = 'lr', ncol = 1, rel_heights = c(0.475, 0.525))

    # save plots to file
    ggsave(here("figures", "Figure_2ab soil temperature.pdf"), width=5, height = 8)

# look for temperature inversion

    # by the end of November, both the 0-20 and 20-40 cm depths have plunged in temperature so that the highest temperature layers are now at 40-60 and 60-80 cm
    thermistor.panel.data %>%
        ggplot(aes(x = Datetime)) +
            geom_ribbon(aes(ymin = `Temperature (°C) #1`, ymax = `Temperature (°C) #2`, fill = `Depth (cm)`), alpha = 0.3, color = NA) +
            geom_line(aes(y = `Temperature (°C)`, color = `Depth (cm)`)) +
            labs(x="Date", y = "Soil temperature (°C)") +
            # axes
                xlim(as.POSIXct('2019/11/01'), as.POSIXct('2019/12/01')) + ylim(-1, 1)

# get max average temperatures for main text

    imnavait.temps %>%
        # only including 12 months of data in this paper
            filter(Datetime >= as.POSIXct('2019/06/06') & Datetime < as.POSIXct('2020/06/06')) %>%
        # average thermistors at each depth
            group_by(Datetime, `Depth (cm)`) %>%
            summarize(`mean Temperature (°C)` = mean(`Temperature (°C)`), .groups = 'drop') %>%
        # find max at each depth
            group_by(`Depth (cm)`) %>%
            summarize(`maximum mean Temperature (°C)` = max(`mean Temperature (°C)`), .groups = 'drop')

    # 20     8.717
    # 40     4.204
    # 60    -0.018
    # 80    -0.115
    # 100   -0.283
    # 120   -0.353
