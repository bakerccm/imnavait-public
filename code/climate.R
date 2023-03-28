# code for climate data from nearby NOAA weather stations
# - Figure 1bc climate graphs

# load packages

    library("here")
    library("tidyverse")
    library("scales")
    library("cowplot")

    # devtools::install_github("wilkelab/ungeviz")
    library("ungeviz") # for geom_hpline()

# get data

    # 30 year temperature normals

        toolik.lake.30y <- read.csv(file = here("data", "climate", "USW00096409_ToolikLake_1991_2020_monthly_normals.csv"))

    # 30 year precipitation normals

        imnavait.creek.30y <- read.csv(file = here("data", "climate", "USS0049T01S_ImnavaitCreek_1991_2020_monthly_normals.csv"))

    # 2019 precipitation data

        imnavait.creek.2019 <- read.delim(file = here("data", "climate", "Imnavait_Creek_2019_precipitation.txt"),
            header = TRUE, sep = "\t", comment.char = "#")

# get colors to match ggplot colors in other figures

    gg_cols <- hue_pal()(2)
    names(gg_cols) <- c("red", "blue")

# 30-year temperature normals for Toolik Lake station

    temperature.plot <- toolik.lake.30y %>%
            mutate(month_name = month.name[DATE], month_abb = month.abb[DATE]) %>%
            mutate(month_name = factor(month_name, levels = month.name), month_abb = factor(month_abb, levels = month.abb)) %>%
            # convert °F to °C
                mutate(MLY.TMIN.NORMAL.C = (MLY.TMIN.NORMAL - 32) * 5/9, MLY.TMAX.NORMAL.C = (MLY.TMAX.NORMAL - 32) * 5/9) %>%
        ggplot(aes(x= month_abb)) +
            # change background to white but retain gridlines
                theme_bw() +
            # add temperature data
                geom_point(aes(y = `MLY.TMIN.NORMAL.C`), col = gg_cols["blue"]) +
                geom_line(aes(y = `MLY.TMIN.NORMAL.C`, group = 1), col = gg_cols["blue"]) +
                geom_point(aes(y = `MLY.TMAX.NORMAL.C`), col = gg_cols["red"]) +
                geom_line(aes(y = `MLY.TMAX.NORMAL.C`, group = 1), col = gg_cols["red"]) +
            # add zero degree C line
                geom_hline(yintercept = 0, linetype = "dotted") +
            # add labels
                labs(x = "Month", y = "Air temperature (°C)") +
                annotate(
                    geom = "text", x = 2, y = -2,
                    label = "Monthly\nmax", hjust = 0.5, vjust = 1, size = 4, col = gg_cols["red"]
                ) +
                annotate(
                    geom = "text", x = 5.2, y = -23,
                    label = "Monthly\nmin", hjust = 0.5, vjust = 0.5, size = 4, col = gg_cols["blue"]
                ) +
                annotate(
                    geom = "text", x = 0.8, y = 15.5,
                    label = "Toolik Lake", hjust = 0, vjust = 1, size = 4
                ) +
                annotate(
                    geom = "text", x = 11.2 , y = 0,
                    label = "0°C", hjust = 0.5, vjust = -0.5, size = 4
                )

# 30-year precipitation normals and 2019 precipitation values for Imnavait Creek station

    precip.data.30y <- imnavait.creek.30y %>%
            # generate month names and abbreviations from month numbers
                mutate(month_name = month.name[month], month_abb = month.abb[month]) %>%
                mutate(month_name = factor(month_name, levels = month.name), month_abb = factor(month_abb, levels = month.abb)) %>%
            # convert inches to mm
                mutate(MLY.PRCP.NORMAL.mm = MLY.PRCP.NORMAL * 25.4)

    precip.data.2019 <- imnavait.creek.2019 %>%
            # months are already provided as abbreviated names
                rename(month_abb = Month) %>%
                mutate(month_abb = factor(month_abb, levels = month.abb))

    precipitation.plot <- precip.data.30y %>%
        ggplot(aes(x= month_abb)) +
            # change background to white but retain gridlines
                theme_bw() +
            # add 30 year precipitation normals
                geom_col(aes(y = MLY.PRCP.NORMAL.mm)) +
            # add 2019 values
                geom_hpline(data = precip.data.2019, aes(y = Precipitation_Increment_mm), col = gg_cols["red"], size = 0.75) +
            # add labels
                labs(x = "Month", y = "Precipitation (mm)") +
                annotate(
                    geom = "text", x = 0.8, y = 150,
                    label = "Imnavait Creek", hjust = 0, vjust = 1, size = 4
                )

# combine temperature and precipiation plots

    plot_grid(
        temperature.plot + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()),
        precipitation.plot, align = "v", axis = "lr", ncol=1, rel_heights = c(3, 2)
    )

    ggsave(here("figures","Figure_1bc_climate.pdf"), width = 4, height = 4, useDingbats = FALSE)
