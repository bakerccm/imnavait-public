# code to normalize reads by total sum scaling, add SEPP tree, and
# package up data from 16S and ITS datasets into a single rds file

# load packages

    library("here")
    library("tidyverse")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

# get file names

    args = commandArgs(trailingOnly=TRUE)
    
	filenames <- list(
		"16S" = args[1], # .rds file containing 16S phyloseq object with chloroplasts etc removed
		"16S_sepp_tree" = args[2], # 16S placement tree from SEPP
		"ITS" = args[3], # .rds file containing ITS phyloseq object
		"output" = args[4] # .rdata file containing original and normalized phyloseq objects for 16S and ITS
	)

# get input data

    imnavait <- lapply(filenames[c("16S", "ITS")], readRDS)

# attach SEPP tree for 16S data

    imnavait_16S_sepp_tree <- read_tree_greengenes(filenames$`16S_sepp_tree`)

    imnavait$`16S` <- merge_phyloseq(imnavait$`16S`, phy_tree(imnavait_16S_sepp_tree))

# Normalize sequence tables to relative read counts, and aggregate to genus and family
# (see, e.g., code at https://github.com/rthanis/animal-cdiff/blob/master/analysis/alpha-beta.Rmd)

# N.B. tax_glom removes ASVs without the relevant taxonomic assignment
# so this code calculates total sum scaling once (at ASV level, so without tax_glom)
# and then uses those relative abundances for the aggregated versions, so the
# latter will likely have relative abundances that do not sum to 1 (per sample).

imnavait_relative_ASVs <- list(
    `16S_ASVs` = imnavait$`16S` %>% transform_sample_counts(function (x) x / sum(x)),
    `ITS_ASVs` = imnavait$`ITS` %>% transform_sample_counts(function (x) x / sum(x))
)

# do each database separately, since the taxon ranks to the left of the column used for tax_glom should be only the hierarchy for the database being glommed

    tax_table_16S <- tax_table(imnavait_relative_ASVs$`16S_ASVs`) %>% as("matrix") %>% as.data.frame(stringsAsFactors = FALSE)
    rownames(tax_table_16S) <- rownames(tax_table(imnavait_relative_ASVs$`16S_ASVs`))

    tax_table_ITS <- tax_table(imnavait_relative_ASVs$`ITS_ASVs`) %>% as("matrix") %>% as.data.frame(stringsAsFactors = FALSE)
    rownames(tax_table_ITS) <- rownames(tax_table(imnavait_relative_ASVs$`ITS_ASVs`))

    # silva
    tax_table(imnavait_relative_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("silva")) %>% as.matrix() %>% tax_table()
    imnavait_relative_silva <- list(
        `16S_silva_genus` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("silva_genus"),
        `16S_silva_family` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("silva_family"),
        `16S_silva_order` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("silva_order"),
        `16S_silva_class` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("silva_class"),
        `16S_silva_phylum` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("silva_phylum")
    )

    # rdp
    tax_table(imnavait_relative_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("rdp")) %>% as.matrix() %>% tax_table()
    imnavait_relative_rdp <- list(
        `16S_rdp_genus` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("rdp_genus"),
        `16S_rdp_family` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("rdp_family"),
        `16S_rdp_order` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("rdp_order"),
        `16S_rdp_class` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("rdp_class"),
        `16S_rdp_phylum` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("rdp_phylum")
    )
    
    # decipher
    tax_table(imnavait_relative_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("decipher")) %>% as.matrix() %>% tax_table()
    imnavait_relative_decipher <- list(
        `16S_decipher_genus` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("decipher_genus"),
        `16S_decipher_family` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("decipher_family"),
        `16S_decipher_order` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("decipher_order"),
        `16S_decipher_class` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("decipher_class"),
        `16S_decipher_phylum` = imnavait_relative_ASVs$`16S_ASVs` %>% tax_glom("decipher_phylum")
    )

    # unite
    tax_table(imnavait_relative_ASVs$`ITS_ASVs`) <- tax_table_ITS %>% select(starts_with("unite")) %>% as.matrix() %>% tax_table()
    imnavait_relative_unite <- list(
        `ITS_unite_genus` = imnavait_relative_ASVs$`ITS_ASVs` %>% tax_glom("unite_genus"),
        `ITS_unite_family` = imnavait_relative_ASVs$`ITS_ASVs` %>% tax_glom("unite_family"),
        `ITS_unite_order` = imnavait_relative_ASVs$`ITS_ASVs` %>% tax_glom("unite_order"),
        `ITS_unite_class` = imnavait_relative_ASVs$`ITS_ASVs` %>% tax_glom("unite_class"),
        `ITS_unite_phylum` = imnavait_relative_ASVs$`ITS_ASVs` %>% tax_glom("unite_phylum")
    )
    
    # restore tax_tables
    tax_table(imnavait_relative_ASVs$`16S_ASVs`) <- tax_table_16S %>% as.matrix() %>% tax_table()
    tax_table(imnavait_relative_ASVs$`ITS_ASVs`) <- tax_table_ITS %>% as.matrix() %>% tax_table()

    # concatenate results
    imnavait_relative <- c(imnavait_relative_ASVs, imnavait_relative_silva, imnavait_relative_rdp, imnavait_relative_decipher, imnavait_relative_unite)

# save original (imnavait) and normalized (imnavait_relative) phyloseq objects

    save(imnavait, imnavait_relative, file = filenames$output)
