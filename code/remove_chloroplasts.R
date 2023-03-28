# code to filter undesired taxa from 16S dataset

# load packages

    library("here")
    library("phyloseq")
    library("Biostrings")
    library("tidyverse")

# get file names

    args = commandArgs(trailingOnly=TRUE)
    
    filenames <- list(
        input_phyloseq = args[1], # .rds file containing 16S phyloseq object
        output_phyloseq = args[2], # .rds file for 16S phyloseq object with chloroplasts, mitochondria etc filtered out
        output_fasta = args[3] # .fasta file for 16S sequences with chloroplasts, mitochondria etc filtered out
    )

# get input data

    phy.16S <- readRDS(filenames$input_phyloseq)

# examine taxonomic assignments

    # extract taxon data as tibbles
    taxa.16S <- tax_table(phy.16S) %>% as("matrix") %>% as_tibble()

    taxa.16S %>% colnames()
    #  [1] "silva_kingdom"          "silva_phylum"           "silva_class"            "silva_order"            "silva_family"           "silva_genus"            "silva_species_single"   "silva_species_multiple" "rdp_kingdom"
    # [10] "rdp_phylum"             "rdp_class"              "rdp_order"              "rdp_family"             "rdp_genus"              "rdp_species_single"     "rdp_species_multiple"   "decipher_domain"        "decipher_phylum"
    # [19] "decipher_class"         "decipher_order"         "decipher_family"        "decipher_genus"         "decipher_species"

    # Eukaryotes
    taxa.16S %>% filter(silva_kingdom == "Eukaryota") # unclear what this is, lots of NAs
    taxa.16S %>% filter(rdp_kingdom == "Eukaryota") # 5 rows of mitochondria, all included in silva_family == "Mitochondria"
    taxa.16S %>% filter(decipher_domain == "Eukaryota") # none

    # Archaea
    taxa.16S %>% filter(silva_kingdom == "Archaea") # 62 rows
    taxa.16S %>% filter(rdp_kingdom == "Archaea") # 63 rows
    taxa.16S %>% filter(decipher_domain == "Archaea") # 35 rows

    # mitochondria
    taxa.16S %>% filter(silva_family == "Mitochondria") # 245 rows
    taxa.16S %>% filter(decipher_family == "Mitochondria") # 7 rows, all included in silva_family == "Mitochondria"

    # NAs
    taxa.16S %>% filter(is.na(silva_kingdom)) # 2 rows
    taxa.16S %>% filter(is.na(rdp_kingdom)) #31 rows incl some mitochondria (according to silva) and some other stuff, often questionably identified
    taxa.16S %>% filter(is.na(decipher_domain)) #2,928 rows, many of which seem adequately identified by silva or rdp

    # chloroplasts
    taxa.16S %>% filter(silva_order == "Chloroplast" | silva_order == "Cyanobacteria/Chloroplast") # 39 rows
    taxa.16S %>% filter(rdp_phylum == "Chloroplast" | rdp_phylum == "Cyanobacteria/Chloroplast") # 43 rows
    taxa.16S %>% filter(rdp_class == "Chloroplast" | rdp_class == "Cyanobacteria/Chloroplast") # 27 rows
    taxa.16S %>% filter(rdp_order == "Chloroplast" | rdp_order == "Cyanobacteria/Chloroplast") # 27 rows
    taxa.16S %>% filter(decipher_order == "Chloroplast" | decipher_order == "Cyanobacteria/Chloroplast") # 13 rows

# filter out undesired taxa
    
    ## note that subset_taxa removes NAs ##

    phy.16S.cleaned <- phy.16S %>%
        # Eukaryotes and NAs
            subset_taxa(silva_kingdom != "Eukaryota") %>% # also removes NAs in silva_kingdom
            subset_taxa(rdp_kingdom != "Eukaryota") %>%  # also removes NAs in rdp_kingdom
            subset_taxa(decipher_domain != "Eukaryota" | is.na(decipher_domain)) %>% # many of the decipher_domain NAs seem to be assigned just fine with silva or RDP
        # Archaea
            subset_taxa(silva_kingdom != "Archaea" | is.na(silva_kingdom)) %>%  ## fix these
            subset_taxa(rdp_kingdom != "Archaea" | is.na(rdp_kingdom)) %>%
            subset_taxa(decipher_domain != "Archaea" | is.na(decipher_domain)) %>%
        # mitochondria
            subset_taxa(silva_family != "Mitochondria" | is.na(silva_family)) %>%
            subset_taxa(decipher_family != "Mitochondria" | is.na(decipher_family)) %>%
        # chloroplasts
            subset_taxa(silva_order != "Chloroplast" | is.na(silva_order)) %>%
            subset_taxa(silva_order != "Cyanobacteria/Chloroplast" | is.na(silva_order)) %>%
            subset_taxa(rdp_phylum != "Chloroplast" | is.na(rdp_phylum)) %>%
            subset_taxa(rdp_phylum != "Cyanobacteria/Chloroplast" | is.na(rdp_phylum)) %>%
            subset_taxa(rdp_class != "Chloroplast" | is.na(rdp_class)) %>%
            subset_taxa(rdp_class != "Cyanobacteria/Chloroplast" | is.na(rdp_class)) %>%
            subset_taxa(rdp_order != "Chloroplast" | is.na(rdp_order)) %>%
            subset_taxa(rdp_order != "Cyanobacteria/Chloroplast" | is.na(rdp_order)) %>%
            subset_taxa(decipher_order != "Chloroplast" | is.na(decipher_order)) %>%
            subset_taxa(decipher_order != "Cyanobacteria/Chloroplast" | is.na(decipher_order))

# look at 16S lengths
    
    # in the unfiltered dataset

        refseq(phy.16S) %>% lengths() %>% table() # lots of long sequences
        
        original_longseqs <- prune_taxa((refseq(phy.16S) %>% lengths()) >= 260, phy.16S)
        tax_table(original_longseqs) # lots of them are mitochondria

    # in the filtered dataset
        
        refseq(phy.16S.cleaned) %>% lengths() %>% table() # filtering doesn't totally clean up the length variation in 16S but it is pretty good
        
        filtered_longseqs <- prune_taxa((refseq(phy.16S.cleaned) %>% lengths()) >= 260, phy.16S.cleaned)
        tax_table(filtered_longseqs)
        # filtered_longseqs %>% refseq() %>% writeXStringSet(filepath = "long_16S_asvs.fasta")

# filter 16S taxa based on fragment lengths
        
    phy.16S.cleaned <- prune_taxa((refseq(phy.16S.cleaned) %>% lengths()) %in% 251:257, phy.16S.cleaned)

# save cleaned and filtered phyloseq object

    saveRDS(phy.16S.cleaned, file = filenames$output_phyloseq)

# save sequences to FASTA file

    writeXStringSet(refseq(phy.16S.cleaned), filepath = filenames$output_fasta)
