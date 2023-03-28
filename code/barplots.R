# draw taxonomy barplots

# load packages

    library("here")
    library("tidyverse")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    library("cowplot") # for plot_grid

    #install.packages("RColorBrewer")
    library("RColorBrewer")
    # display.brewer.all() # show all the color palettes

# get data

    # loads imnavait (original data) and imnavait_relative (normalized data)
    load(here("out", "combined", "imnavait_normalized.rdata"))

# depth ranges (in cm) used in sampling

    depth_range_labels <- c("0-20","20-40","40-60","60-80","80-100")

# number of taxonomic groups (phyla, classes etc) to display

    num.taxa <- 10
    
    # plots below based on top num.taxa groups by abundance in normalized ASV tables
    
    # any further taxa are removed (so the difference between the plotted column and 1
    # is effectively 'other')

# top num.taxa Silva phyla

    # plot_bar() produces a basic plot but doesn't give a lot of formatting flexibility, so do it manually here ...

    # filter to top num.taxa phyla
        topn.names <- names(sort(taxa_sums(imnavait_relative$`16S_silva_phylum`), decreasing=TRUE))[1:num.taxa]
        topn.phyloseq <- prune_taxa(topn.names, imnavait_relative$`16S_silva_phylum`)

    # reorganize data for plotting into data.frame
        topn.phyloseq_otu <- otu_table(topn.phyloseq) %>% as.data.frame() %>% rownames_to_column("SampleID")
        topn.phyloseq_tax <- tax_table(topn.phyloseq) %>% as.data.frame() %>% rownames_to_column("ASV") %>% select(ASV, silva_phylum)
        topn.phyloseq_sample <- sample_data(topn.phyloseq) %>% as("data.frame") %>% select(SampleID, Month, Depth_range, Replicate)
        topn.df <- full_join(topn.phyloseq_sample, topn.phyloseq_otu, by = "SampleID") %>%
            pivot_longer(-c(SampleID, Month, Depth_range, Replicate), names_to = "ASV", values_to = "relative.abundance") %>%
            left_join(topn.phyloseq_tax, by = "ASV") %>% select(-ASV) %>%
            mutate(Depth_range = factor(Depth_range, levels = depth_range_labels))

    # create plot
        barplot_silva_phylum <- topn.df %>%
            ggplot(aes(x = as.numeric(Depth_range) + 0.25 * Replicate - 0.5, y= relative.abundance, fill = silva_phylum)) +
                theme_bw() +
                geom_col(width=0.2, position = "stack", col = "black", size = 0.3) + 
                scale_fill_brewer(palette = "Set3") +
                facet_grid(~ Month) +
                labs (x = "Depth range (cm)", y= "Relative abundance", fill = "Silva phylum") +
                scale_x_continuous(breaks = 1:5, labels = depth_range_labels) + 
                theme(axis.text.x = element_text(size=8))

# top num.taxa unite classes

    # plot_bar() produces a basic plot but doesn't give a lot of formatting flexibility, so do it manually here ...

    # filter to top num.taxa classes
        topn.names <- names(sort(taxa_sums(imnavait_relative$ITS_unite_class), decreasing=TRUE))[1:num.taxa]
        topn.phyloseq <- prune_taxa(topn.names, imnavait_relative$ITS_unite_class)

    # reorganize data for plotting into data.frame
        topn.phyloseq_otu <- otu_table(topn.phyloseq) %>% as.data.frame() %>% rownames_to_column("SampleID")
        topn.phyloseq_tax <- tax_table(topn.phyloseq) %>% as.data.frame() %>% rownames_to_column("ASV") %>% select(ASV, unite_class)
        topn.phyloseq_sample <- sample_data(topn.phyloseq) %>% as("data.frame") %>% select(SampleID, Month, Depth_range, Replicate)
        topn.df <- full_join(topn.phyloseq_sample, topn.phyloseq_otu, by = "SampleID") %>%
            pivot_longer(-c(SampleID, Month, Depth_range, Replicate), names_to = "ASV", values_to = "relative.abundance") %>%
            left_join(topn.phyloseq_tax, by = "ASV") %>% select(-ASV) %>%
            mutate(Depth_range = factor(Depth_range, levels = depth_range_labels))

    # create plot
        barplot_unite_class <- topn.df %>%
            mutate(unite_class = gsub("c__", "", unite_class)) %>%
            mutate(unite_class = gsub("_cls_Incertae_sedis", " class Incertae sedis", unite_class)) %>%
            ggplot(aes(x = as.numeric(Depth_range) + 0.25 * Replicate - 0.5, y= relative.abundance, fill = unite_class)) +
                theme_bw() +
                geom_col(width=0.2, position = "stack", col = "black", size = 0.3) + facet_grid(~ Month) +
                scale_fill_brewer(palette = "Set3") + 
                labs (x = "Depth range (cm)", y= "Relative abundance", fill = "Unite class") +
                scale_x_continuous(breaks = 1:5, labels = depth_range_labels) + 
                theme(axis.text.x = element_text(size=8))

# create composite barplot figure (so the plot areas line up)

    plot_grid(
        barplot_silva_phylum + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
        barplot_unite_class + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
        align = "v", axis = "lr", ncol=1, rel_heights = c(1,1))

    ggsave(here("figures","Figure_S2ab taxon barplots.pdf"), width = 10, height = 6, useDingbats = FALSE)
