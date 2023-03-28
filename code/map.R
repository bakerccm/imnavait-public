# code to draw map
# - Fig 1a map of field sampling location

# see https://r-spatial.org/r/2018/10/25/ggplot2-sf.html for guidance

# load packages

    library("here")
    library("scales")
    library("tidyverse")
    library("sf")
    library("rnaturalearth")
    library("rnaturalearthdata")
    library("rgeos")
    library("ggspatial")

# get colors to match ggplot colors in other figures

	gg_cols <- hue_pal()(2)
	names(gg_cols) <- c("red", "blue")

# coords of Imnavait field site

	imnavait <- list(label = "Imnavait Creek\nsampling location", Latitude = 68.611949, Longitude = -149.315718)

# get map

	world <- ne_countries(scale = "medium", returnclass = "sf")

# plot map

	ggplot(data = world) +
		# change background to white but retain gridlines
			theme_bw() +
		# plot map
			geom_sf() +
		# north arrow
			annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
		# customize x axis
			scale_x_continuous(breaks = c(-170, -160, -150, -140)) +
		# map limits
			coord_sf(xlim = imnavait$Longitude + c(-30,10), ylim = imnavait$Latitude + c(-20,5), expand = FALSE) +
		# show sampling field site with label
			geom_point(data = as.data.frame(imnavait), aes(x=Longitude, y= Latitude), color = gg_cols[1]) +
			geom_text(data = as.data.frame(imnavait), aes(x=Longitude, y= Latitude, label = label), color = gg_cols[1], nudge_x = -1.85, nudge_y = -2.1)

	ggsave(here("figures", "Figure_1a_alaska_map.pdf"), width=3.5, height = 4, useDingbats = FALSE)
