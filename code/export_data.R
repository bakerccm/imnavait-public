# exports data for supplementary files

# load packages

    library("readxl")
    library("writexl")
    library("phyloseq")
    library("tidyverse")

# get file names

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        "input_rdata" = args[1], # .rdata file containing original and normalized phyloseq objects for 16S and ITS
        "input_soil_nutrients" = args[2], # rdata file containing soil nutrients data
        "input_porewater_chemistry" = args[3], # rdata file containing porewater chemistry data
        "input_temperatures" = args[4], # rdata file containing thermistor data
        "output_soil_nutrients" = args[5], # CSV file for soil nutrients data
        "output_porewater_chemistry" = args[6], # CSV file for porewater chemistry data
        "output_porewater_technical" = args[7], # CSV file for porewater technical data
        "output_temperatures" = args[8], # CSV file for thermistor data
        "output_asv_16S" = args[9], # filename for 16S ASV table output (text file)
        "output_asv_ITS" = args[10], # filename for ITS ASV table output (text file)
        "output_taxon_16S" = args[11], # filename for 16S taxon table output (text file)
        "output_taxon_ITS" = args[12] # filename for ITS taxon table output (text file)
    )

# get input data

    # loads imnavait (original data) and imnavait_relative (normalized data)
    load(filenames$input_rdata)

    # loads soil_nutrients
    load(filenames$input_soil_nutrients)

    # loads porewater_chemistry_data and porewater_technical_data
    load(filenames$input_porewater_chemistry)

    # loads imnavait.temps
    load(filenames$input_soil_nutrients)

# save soil nutrients data as CSV file

    soil_nutrients %>%
        # reformat as wide
            pivot_wider(id_cols = c(Month, Depth_range, Depth_min, Depth_max, replicate),
                names_from = measurement, values_from = value) %>%
        # clean up
            arrange(Month, Depth_range, replicate) %>%
            select(-Depth_min, -Depth_max) %>%
        # write to file
            write_csv(., file = filenames$output_soil_nutrients)

# save porewater chemistry data as CSV file

    # porewater chemistry data
            
        porewater_chemistry_data %>%
            select(Month, Depth_range, element, element_symbol, concentration) %>%
            write_csv(., file = filenames$output_porewater_chemistry)
            
    # porewater technical data

        write_csv(porewater_technical_data, file = filenames$output_porewater_technical)

# export temperature data to file

    imnavait.temps %>% 
        rename("Depth_cm" = "Depth (cm)",   "Temperature_degC" = "Temperature (°C)") %>%
        mutate(Temperature_degC = format(Temperature_degC, nsmall = 3)) %>%
        arrange(Datetime, Depth_cm, Replicate) %>%
        write_csv(., file = filenames$output_temperatures)

# export ASV tables

    export.asv.table <- list()

    for (dataset in c("16S","ITS")) {
        export.asv.table[[dataset]] <- imnavait[[dataset]] %>% otu_table() %>% as.data.frame()
        export.asv.table[[dataset]] <- export.asv.table[[dataset]][ ,order(colnames(export.asv.table[[dataset]]))] %>%
            rownames_to_column("SampleID") %>% arrange(SampleID)
        write_delim(export.asv.table[[dataset]], file = filenames[[paste0("output_asv_", dataset)]])
    }

    # view samples of exported tables
    # export.asv.table[["16S"]][1:5,1:5]
    # export.asv.table[["ITS"]][1:5,1:5]

# export taxon tables

    export.taxon.table <- list()

    for (dataset in c("16S","ITS")) {
        export.taxon.table[[dataset]] <- imnavait[[dataset]] %>% tax_table() %>% as.data.frame() %>% rownames_to_column("ASV") %>% arrange(ASV)
        write_delim(export.taxon.table[[dataset]], file = filenames[[paste0("output_taxon_", dataset)]])
    }

    # view samples of exported tables
    # export.taxon.table[["16S"]][1:5,1:5]
    # export.taxon.table[["ITS"]][1:5,1:5]
