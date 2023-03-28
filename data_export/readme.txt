SUPPLEMENTARY FILE S12
README FOR SUPPLEMENTARY DATASET FILES S4 THROUGH S11

------------------------------------------------------
File S4 imnavait_thermistor_data.csv

Table of data from thermistors in soil

The table contains the following fields:

Datetime: Date and time in ISO 8601 format
Depth_cm: Measurement depth in cm.
Replicate: Identifier for the thermistor, i.e. 1 or 2 at each depth
Temperature_degC: Measured temperature in degrees Celsius.

------------------------------------------------------
File S5 imnavait_soil_nutrient_data.csv

Table of data on soil nutrients and pH

The table contains the following fields:

Month: Sampling month
Depth_range: Sampling depth in cm.
replicate: Identifier for the replicate, i.e. a or b
C, N: Percentage carbon and nitrogen content respectively measured with dry combustion and direct measurement of total nutrients with Elementar Macro Cube
OM: Organic matter calculated using total organic carbon x 2 as per review by Pribyl (2010) (Geoderma)
LOI: Loss on ignition for organic matter determination
pH: Measured in 1:1 soil:water ratio on Hanna HI5522 benchtop meter
NH4-N: 2M KCl extraction measured on Lachat QuikChem 8500 Series 2 flow injection analyzer

Not shown in table:

NO3-N: 2M KCl extraction measured using Griess reagents on VWR V1200 spectrophotometer. All measurements for NO3-N were below quantifiable limits, i.e. < 0.2 ppm NO3-N detected in 2M KCl extract

------------------------------------------------------
File S6 imnavait_porewater_chemistry_data_filtered.csv
File S7 imnavait_porewater_technical_data.csv

Tables of data on porewater metal concentrations and associated detection/reporting limits.

The tables contain the following fields:

Month: Sampling month
Depth_range: Sampling depth in cm.
element: Metal name
element_symbol: Metal symbol
concentration: Concentration in mg/L
detection_limit: Detection limit (i.e. lowest concentration that can reliably be distinguished from zero, in mg/L)
reporting_limit: Reporting limit (i.e. lowest concentration that can reliably be measured, in mg/L)
units: Units used for concentration, detection_limit and reporting_limit (i.e. mg/L)

------------------------------------------------------
File S8 imnavait_ASV_table_16S.txt
File S9 imnavait_ASV_table_ITS.txt

Tables of raw (not normalized) read counts broken down by sample and ASV at the end of our bioinformatic pipeline.

The tables contain the following fields:

SampleID: Sample identifying label. Each microbial sample has a unique identifier with format <Month>.<MinDepth>.<MaxDepth>.<Rx> where Month is the sampling month, MinDepth and MaxDepth describe the sampling depth in cm (e.g. 20 and 40 respectively for samples from the 20-40 cm depth range), and Rx is an identifier for the replicate (i.e. R1, R2 or R3).

ASV_16S_xxxxx, ASV_ITS_xxxxx: ASV identifying label. Each ASV has a unique identifier with format ASV_<Dataset>_<Number> where Dataset identifies whether the ASV is from the 16S or ITS dataset, and Number is an identifier for the ASV within each dataset (note that two ASVs in different datasets may have the same Number).

------------------------------------------------------
File S10 imnavait_taxon_table_16S.txt
File S11 imnavait_taxon_table_ITS.txt

Tables of taxonomy assignments for all ASVs. Each ASV is shown with taxonomy assigned from Silva, RDP or DECIPHER (16S) , or UNITE (ITS). Empty cells indicate no taxonomy assignment could be made at that level.

The tables contain the following fields:

ASV: ASV identifying label. Each ASV has a unique identifier with format ASV_<Dataset>_<Number> where Dataset identifies whether the ASV is from the 16S or ITS dataset, and Number is an identifier for the ASV within each dataset (note that two ASVs in different datasets may have the same Number).

silva_kingdom: [16S table only] Silva database kingdom
silva_phylum: [16S table only] Silva database phylum
silva_class: [16S table only] Silva database class
silva_order: [16S table only] Silva database order
silva_family: [16S table only] Silva database family
silva_genus: [16S table only] Silva database genus
silva_species_single: [16S table only] Silva database species, if uniquely determined
silva_species_multiple: [16S table only] Silva database species list, potentially including multiple possible species
	
rdp_kingdom: [16S table only] RDP database kingdom
rdp_phylum: [16S table only] RDP database phylum
rdp_class: [16S table only] RDP database class
rdp_order: [16S table only] RDP database order
rdp_family: [16S table only] RDP database family
rdp_genus: [16S table only] RDP database genus
rdp_species_single: [16S table only] RDP database species, if uniquely determined
rdp_species_multiple: [16S table only] RDP database species list, potentially including multiple possible species
	
decipher_domain: [16S table only] DECIPHER kingdom
decipher_phylum: [16S table only] DECIPHER phylum
decipher_class: [16S table only] DECIPHER class
decipher_order: [16S table only] DECIPHER order
decipher_family: [16S table only] DECIPHER family
decipher_genus: [16S table only] DECIPHER genus
decipher_species: [16S table only] DECIPHER species
	
unite_kingdom: [ITS table only] UNITE kingdom
unite_phylum: [ITS table only] UNITE phylum
unite_class: [ITS table only] UNITE class
unite_order: [ITS table only] UNITE order
unite_family: [ITS table only] UNITE family
unite_genus: [ITS table only] UNITE genus
unite_species: [ITS table only] UNITE species

------------------------------------------------------
