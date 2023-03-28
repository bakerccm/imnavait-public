# snakefile for Imnavait manuscript

# Chris Baker
# https://github.com/bakerccm

from snakemake.utils import min_version
min_version("6.4.1")

################################
## get config file

configfile: 'config/config.yaml'

################################
## wildcard constraints

wildcard_constraints:
    sample = '[^_/]+',
    file = '[^/\.]+'

################################
## get sample metadata

import pandas as pd
METADATA = {
    '16S': pd.read_csv(config['sample_metadata']['16S'], sep = '\t', index_col = 'SampleID'),
    'ITS': pd.read_csv(config['sample_metadata']['ITS'], sep = '\t', index_col = 'SampleID')
}
SAMPLES = {
    '16S': list(METADATA['16S'].index),
    'ITS': list(METADATA['ITS'].index)
}

################################
## default rule

rule all:
    input:
        # fastqc on raw data
            expand('out/16S/fastqc/{sample}-{read}_fastqc.{ext}', sample = SAMPLES['16S'], read = ["R1","R2"], ext = ["html","zip"]),
            expand('out/ITS/fastqc/{sample}-{read}_fastqc.{ext}', sample = SAMPLES['ITS'], read = ["R1","R2"], ext = ["html","zip"]),
        # read quality profile plots
            expand('out/16S/read_quality_profiles/{sample}.pdf', sample = SAMPLES['16S']),
            expand('out/ITS/read_quality_profiles/{sample}.pdf', sample = SAMPLES['ITS']),
        # output from normalize_phyloseq
            expand("out/combined/imnavait_normalized.rdata"),
        # figures and other output for manuscript
            # from export_data.R
                "data_export/imnavait_soil_nutrient_data.csv",
                "data_export/imnavait_porewater_chemistry_data_filtered.csv",
                "data_export/imnavait_porewater_technical_data.csv",
                "data_export/imnavait_thermistor_data.csv",
                expand("data_export/imnavait_{type}_table_{dataset}.txt", type = ["asv","taxon"], dataset = ["16S","ITS"]),
            # from map.R
                "figures/Figure_1a_alaska_map.pdf",
            # from climate.R
                "figures/Figure_1bc_climate.pdf",
            # from soil_temperature.R
                "figures/Figure_2ab soil temperature.pdf",
            # from soil_nutrients_and_porewater_chemistry.R
                "figures/Figure_3 soil CN with points.pdf",
                "figures/Figure_4 pore water chemistry.pdf",
                "figures/Figure_7a dbRDA_16S_ASVs_all.pdf",
                "figures/Figure_7b dbRDA_ITS_ASVs_all.pdf",
                "figures/Figure_S5a dbRDA_16S_ASVs_CN.pdf",
                "figures/Figure_S5b dbRDA_ITS_ASVs_CN.pdf",
            # from alpha_div.R
                "figures/Figure_5a_shannon_lm.pdf",
                "figures/Figure_5b_invsimp_lm.pdf",
                "figures/Figure_A total_reads_vs_month_and_depth.pdf",
                "figures/Figure_B total_reads_vs_depth_and_month.pdf",
                "figures/Figure_C alpha_diversity_vs_total_reads_shannon.pdf",
                "figures/Figure_D alpha_diversity_vs_total_reads_invsimp.pdf",
                "figures/Figure_E diversity_vst_vs_reads.pdf",
            # from beta_div.R
                "figures/Figure_6ab NMDS plots.pdf",
                "figures/Figure_6c sum of squares.pdf",
                "figures/Figure_S3 NMDS 16S UniFrac.pdf",
                "figures/Figure_S4ab NMDS plots sample condition.pdf",
            # from barplots.R
                "figures/Figure_S2ab taxon barplots.pdf",
            # from deseq2.R
                expand("out/combined/deseq2.{type}", type = ["xlsx","rds"])

################################
## check quality of raw sequence files

rule fastqc:
    input:
        lambda wildcards: config['raw_data'][wildcards.dataset] + '/' + wildcards.sample + '-' + wildcards.read + '.fastq.gz'
    output:
        'out/{dataset}/fastqc/{sample}-{read}_fastqc.html',
        'out/{dataset}/fastqc/{sample}-{read}_fastqc.zip'
    params:
        output_dir='out/{dataset}/fastqc'
    conda:
        'envs/fastqc-0.11.9.yaml'
    shell:
        'fastqc -o {params.output_dir} {input}'

################################
## remove Ns and check primer orientation in preparation for running cutadapt

rule filter_Ns:
    input:
        read1=lambda wildcards: config['raw_data'][wildcards.dataset] + '/' + wildcards.sample + '-R1' + '.fastq.gz',
        read2=lambda wildcards: config['raw_data'][wildcards.dataset] + '/' + wildcards.sample + '-R2' + '.fastq.gz'
    output:
        read1='out/{dataset}/demultiplexed_no_Ns/{sample}-R1.fastq.gz',
        read2='out/{dataset}/demultiplexed_no_Ns/{sample}-R2.fastq.gz',
        summary='out/{dataset}/demultiplexed_no_Ns/{sample}-summary.txt'
    params:
        primer_fwd=lambda wildcards: config['primers'][wildcards.dataset]['forward'],
        primer_rev=lambda wildcards: config['primers'][wildcards.dataset]['reverse']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        Rscript code/filter_Ns_and_check_primers.R \
            {input.read1} {input.read2} \
            {params.primer_fwd} {params.primer_rev} \
            {output.read1} {output.read2} {output.summary}
        '''

################################
## use cutadapt to remove primers

# note don't length trim ITS sequences
rule cutadapt:
    input:
        read1='out/{dataset}/demultiplexed_no_Ns/{sample}-R1.fastq.gz',
        read2='out/{dataset}/demultiplexed_no_Ns/{sample}-R2.fastq.gz'
    output:
        read1='out/{dataset}/cutadapt/{sample}-R1.fastq.gz',
        read2='out/{dataset}/cutadapt/{sample}-R2.fastq.gz'
    log:
        'out/{dataset}/cutadapt/{sample}.log'
    params:
        primer_fwd = lambda wildcards: config['primers'][wildcards.dataset]['forward'],
        primer_rev = lambda wildcards: config['primers'][wildcards.dataset]['reverse'],
        min_length = config['cutadapt']['min_length']
    conda:
        'envs/cutadapt-3.5.yaml'
    shell:
        '''
        # these are the primer sequences
            primer_fwd={params.primer_fwd}
            primer_rev={params.primer_rev}
        # get reverse complements
            primer_fwd_rc=`echo ${{primer_fwd}} | tr ACGTMRWSYKVHDBNacgtmrwsykvhdbn TGCAKYWSRMBDHVNtgcakywsrmbdhvn | rev`
            primer_rev_rc=`echo ${{primer_rev}} | tr ACGTMRWSYKVHDBNacgtmrwsykvhdbn TGCAKYWSRMBDHVNtgcakywsrmbdhvn | rev`
        cutadapt -g ${{primer_fwd}} -a ${{primer_rev_rc}} -G ${{primer_rev}} -A ${{primer_fwd_rc}} -n 2 -m {params.min_length} \
            -o {output.read1} -p {output.read2} {input.read1} {input.read2} >{log}
        '''

################################
## plot read quality profiles

rule read_quality_profiles:
    input:
        read1='out/{dataset}/cutadapt/{sample}-R1.fastq.gz',
        read2='out/{dataset}/cutadapt/{sample}-R2.fastq.gz'
    output:
        'out/{dataset}/read_quality_profiles/{sample}.pdf'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        'Rscript code/plot_read_profile.R {input.read1} {input.read2} {output}'

################################
## filter and trim in dada2

# note use of different rules for 16S and ITS
# - 16S workflow uses truncLen to force fixed read lengths
# - ITS workflow does not do this (since ITS can vary in length), but does enforce minLen to remove spurious short reads

rule filter_and_trim_16S:
    input:
        read1 = 'out/16S/cutadapt/{sample}-R1.fastq.gz',
        read2 = 'out/16S/cutadapt/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/16S/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/16S/filterAndTrim/{sample}-R2.fastq.gz'
    params:
        maxN = config['filterAndTrim']['16S']['maxN'],
        truncQ = config['filterAndTrim']['16S']['truncQ'],
        maxEE_read1 = config['filterAndTrim']['16S']['maxEE_read1'],
        maxEE_read2 = config['filterAndTrim']['16S']['maxEE_read2'],
        truncLen_read1 = config['filterAndTrim']['16S']['truncLen_read1'],
        truncLen_read2 = config['filterAndTrim']['16S']['truncLen_read2']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        Rscript code/filter_and_trim_16S.R \
            {input.read1} {input.read2} {output.read1} {output.read2} \
            {params.maxN} {params.truncQ} {params.maxEE_read1} {params.maxEE_read2} {params.truncLen_read1} {params.truncLen_read2}
        '''

rule filter_and_trim_ITS:
    input:
        read1 = 'out/ITS/cutadapt/{sample}-R1.fastq.gz',
        read2 = 'out/ITS/cutadapt/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/ITS/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/ITS/filterAndTrim/{sample}-R2.fastq.gz'
    params:
        maxN = config['filterAndTrim']['ITS']['maxN'],
        truncQ = config['filterAndTrim']['ITS']['truncQ'],
        maxEE_read1 = config['filterAndTrim']['ITS']['maxEE_read1'],
        maxEE_read2 = config['filterAndTrim']['ITS']['maxEE_read2'],
        minLen = config['filterAndTrim']['ITS']['minLen']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        Rscript code/filter_and_trim_ITS.R \
            {input.read1} {input.read2} {output.read1} {output.read2} \
            {params.maxN} {params.truncQ} {params.maxEE_read1} {params.maxEE_read2} {params.minLen}
        '''

################################
## learn error rates in dada2

# arrange input data in folders by read direction to facilitate learnerrors rule
rule learnerrors_inputs:
    input:
        read1 = 'out/{dataset}/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/{dataset}/filterAndTrim/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/{dataset}/learnerrors/R1/{sample}-R1.fastq.gz',
        read2 = 'out/{dataset}/learnerrors/R2/{sample}-R2.fastq.gz'
    shell:
        '''
        ln -s ../../../../{input.read1} {output.read1}
        ln -s ../../../../{input.read2} {output.read2}
        '''

rule learnerrors:
    input:
        # all the fastq files for the specified read direction (R1 or R2) for the specified dataset (16S or ITS)
        lambda wildcards: ['out/' + wildcards.dataset + '/learnerrors/' + wildcards.read + '/' + sample + '-' + wildcards.read + '.fastq.gz' for sample in SAMPLES[wildcards.dataset]]
    output:
        rds = 'out/{dataset}/learnerrors/{read}_error_profile.rds',
        pdf = 'out/{dataset}/learnerrors/{read}_error_profile.pdf'
    log:
        'out/{dataset}/learnerrors/{read}_error_profile.log'
    params:
        input_dir = 'out/{dataset}/learnerrors/{read}',
        nbases = config['learnErrors']['nbases']
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: 1
    shell:
        'Rscript code/learn_errors.R {params.input_dir} {output.rds} {output.pdf} {params.nbases} 1>{log}'

################################
## sequence inference for each sample

# consider running dada in pooled mode to pick up rare sequences
# see https://benjjneb.github.io/dada2/tutorial.html

rule dada:
    input:
        fastq_read1 = 'out/{dataset}/learnerrors/R1/{sample}-R1.fastq.gz',
        fastq_read2 = 'out/{dataset}/learnerrors/R2/{sample}-R2.fastq.gz',
        error_profile_read1 = 'out/{dataset}/learnerrors/R1_error_profile.rds',
        error_profile_read2 = 'out/{dataset}/learnerrors/R2_error_profile.rds'
    output:
        'out/{dataset}/dada/{sample}.rds'
    log:
        'out/{dataset}/dada/{sample}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        'Rscript code/dada.R {input.fastq_read1} {input.fastq_read2} {input.error_profile_read1} {input.error_profile_read2} {output} >{log}'

rule dada_merge_samples:
    input:
        lambda wildcards: ['out/' + wildcards.dataset + '/dada/' + sample + '.rds' for sample in SAMPLES[wildcards.dataset]]
    output:
        sample_list = 'out/{dataset}/dada/sample_list.txt',
        sequence_table = 'out/{dataset}/dada/sequence_table.rds'
    log:
        'out/{dataset}/dada/sequence_table.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # write list of input files
            rm -rf {output.sample_list}; for i in {input};do echo $i>>{output.sample_list}; done
        # make sequence table
            Rscript code/merge_samples.R {output.sample_list} {output.sequence_table} >{log}
        '''

################################
## remove chimeras

rule remove_chimeras:
    input:
        'out/{dataset}/dada/sequence_table.rds'
    output:
        'out/{dataset}/remove_chimeras/sequence_table.rds'
    log:
        'out/{dataset}/remove_chimeras/sequence_table.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        'Rscript code/remove_chimeras.R {input} {output} >{log}'

################################
## assign taxonomy

# uses dada2::assignTaxonomy() and dada2::assignSpecies()
# returns a file of unambiguous matches and a file of multiple matches
rule assign_taxonomy_16S:
    input:
        sequence_table = 'out/16S/remove_chimeras/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database],
        assignSpecies_refFasta = lambda wildcards: config['databases']['assignSpecies'][wildcards.database]
    output:
        single_match = 'out/16S/assign_taxonomy/{database}_single.rds',
        multiple_matches = 'out/16S/assign_taxonomy/{database}_multiple.rds'
    log:
        'out/16S/assign_taxonomy/{database}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: 2
    shell:
        '''
        Rscript code/assign_taxonomy_16S.R {input.sequence_table} \
            {input.assignTaxonomy_refFasta} {input.assignSpecies_refFasta} \
            {output.single_match} {output.multiple_matches} >{log}
        '''

# uses dada2::assignTaxonomy() only
rule assign_taxonomy_ITS:
    input:
        sequence_table = 'out/ITS/remove_chimeras/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database]
    output:
        'out/ITS/assign_taxonomy/{database}.rds'
    log:
        'out/ITS/assign_taxonomy/{database}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: 2
    shell:
        'Rscript code/assign_taxonomy_ITS.R {input.sequence_table} {input.assignTaxonomy_refFasta} {output} >{log}'

# uses DECIPHER::IdTaxa()
rule decipher:
    input:
        sequence_table = 'out/16S/remove_chimeras/sequence_table.rds',
        decipher_ref = lambda wildcards: config['databases']['DECIPHER'][wildcards.database]
    output:
        'out/16S/decipher/{database}.rds'
    log:
        'out/16S/decipher/{database}.log'
    conda:
        'envs/dada2-1.18.0-decipher.yaml'
    threads: 2
    shell:
        'Rscript code/decipher.R {input.sequence_table} {input.decipher_ref} {output} >{log}'

################################
## export results as phyloseq object

rule export_to_phyloseq_16S:
    input:
        sequence_table = "out/16S/remove_chimeras/sequence_table.rds",
        silva_single = "out/16S/assign_taxonomy/silva_single.rds",
        silva_multiple = "out/16S/assign_taxonomy/silva_multiple.rds",
        rdp_single = "out/16S/assign_taxonomy/rdp_single.rds",
        rdp_multiple = "out/16S/assign_taxonomy/rdp_multiple.rds",
        decipher = "out/16S/decipher/silva.rds",
        metadata = config['sample_metadata']['16S']
    output:
        phyloseq="out/16S/phyloseq/phyloseq.rds",
        fasta="out/16S/phyloseq/phyloseq.fasta"
    log:
        "out/16S/phyloseq/phyloseq.log"
    conda:
        "envs/phyloseq.yaml"
    shell:
        '''
        Rscript code/export_to_phyloseq_16S.R {input.sequence_table} \
        {input.silva_single} {input.silva_multiple} {input.rdp_single} {input.rdp_multiple} {input.decipher} \
        {input.metadata} {output.phyloseq} {output.fasta} >{log}
        '''

rule export_to_phyloseq_ITS:
    input:
        sequence_table = "out/ITS/remove_chimeras/sequence_table.rds",
        unite = "out/ITS/assign_taxonomy/unite.rds",
        metadata = config['sample_metadata']['ITS']
    output:
        "out/ITS/phyloseq/phyloseq.rds"
    log:
        "out/ITS/phyloseq/phyloseq.log"
    conda:
        "envs/phyloseq.yaml"
    shell:
        "Rscript code/export_to_phyloseq_ITS.R {input.sequence_table} {input.unite} {input.metadata} {output} >{log}"

################################
## remove chloroplasts, mitochondria etc

rule remove_chloroplasts:
    input:
        "out/16S/phyloseq/phyloseq.rds"
    output:
        rds = "out/16S/phyloseq/phyloseq_cleaned.rds",
        fasta = "out/16S/phyloseq/phyloseq_cleaned.fasta"
    conda:
        "envs/phyloseq-tidyverse.yaml"
    shell:
        "Rscript code/remove_chloroplasts.R {input} {output.rds} {output.fasta}"

################################
## phylogenetic placement

# note: SEPP needs to be installed separately (not through conda) in the 'software' folder,
# following instructions at https://github.com/smirarab/sepp/tree/master/sepp-package
# -- the sepp_env conda environment provides a suitable runtime environment but does not actually install SEPP
rule run_sepp:
    input:
        "out/16S/phyloseq/phyloseq_cleaned.fasta"
    output:
        "out/16S/phyloseq/phyloseq_cleaned_placement.tog.tre"
    conda:
        "envs/sepp_env.yaml"
    params:
        output_prefix="phyloseq_cleaned"
    threads:
        12
    shell:
        # running sepp by cd'ing into this directory seems to be the way to do it but it requires some convoluted relative paths to get it to work here
        '''
        cd software
        ./sepp-package/run-sepp.sh ../{input} {params.output_prefix} -x {threads}
        mv *{params.output_prefix}* ../out/16S/phyloseq
        '''

################################
## normalize results by total-sum scaling

# also adds SEPP placement to phyloseq object
rule normalize_phyloseq:
    input:
        phyloseq_16S = "out/16S/phyloseq/phyloseq_cleaned.rds",
        sepp_tree_16S = "out/16S/phyloseq/phyloseq_cleaned_placement.tog.tre",
        phyloseq_ITS = "out/ITS/phyloseq/phyloseq.rds"
    output:
        "out/combined/imnavait_normalized.rdata"
    log:
        "out/combined/imnavait_normalized.log"
    conda:
        "envs/phyloseq-tidyverse.yaml"
    shell:
        "Rscript code/normalize_phyloseq.R {input.phyloseq_16S} {input.sepp_tree_16S} {input.phyloseq_ITS} {output} >{log}"

################################
## export data and results to tables

rule export_data:
    input:
        amplicons="out/combined/imnavait_normalized.rdata",
        soil_nutrients="data/chemistry/soil_nutrients.rdata",
        porewater_chemistry="data/chemistry/porewater_chemistry.rdata",
        temperatures="data/temperature/imnavait_temps.rdata"
    output:
        soil_nutrients="data_export/imnavait_soil_nutrient_data.csv",
        porewater_chemistry="data_export/imnavait_porewater_chemistry_data_filtered.csv",
        porewater_technical="data_export/imnavait_porewater_technical_data.csv",
        temperatures="data_export/imnavait_thermistor_data.csv",
        asv_16S = "data_export/imnavait_asv_table_16S.txt",
        asv_ITS = "data_export/imnavait_asv_table_ITS.txt",
        taxon_16S = "data_export/imnavait_taxon_table_16S.txt",
        taxon_ITS = "data_export/imnavait_taxon_table_ITS.txt"
    conda:
        "envs/phyloseq-tidyverse.yaml"
    shell:
        '''
        Rscript code/export_data.R {input.amplicons} {input.soil_nutrients} {input.porewater_chemistry} {input.temperatures} \
        {output.soil_nutrients} {output.porewater_chemistry} {output.porewater_technical} {output.temperatures} \
        {output.asv_16S} {output.asv_ITS} {output.taxon_16S} {output.taxon_ITS}
        '''

################################
## code for manuscript ##

# We have added the following rules, for our statistical analysis, to this snakefile to document
# dependencies and output, and so that the steps display in the rulegraph -- however, we didn't 
# actually run these in snakemake, and most of this code is probably more usefully run interactively

rule map:
    output:
        "figures/Figure_1a_alaska_map.pdf"
    shell:
        "Rscript code/map.R"

rule climate:
    input:
        "data/climate/USW00096409_ToolikLake_1991_2020_monthly_normals.csv",
        "data/climate/USS0049T01S_ImnavaitCreek_1991_2020_monthly_normals.csv",
        "data/climate/Imnavait_Creek_2019_precipitation.txt"
    output:
        "figures/Figure_1bc_climate.pdf"
    shell:
        "Rscript code/climate.R"

rule soil_temperature:
    input:
        "data/temperature/imnavait_temps.rdata"
    output:
        "figures/Figure_2ab soil temperature.pdf"
    shell:
        "Rscript code/soil_temperature.R"

rule soil_nutrients_and_porewater_chemistry:
    input:
        "out/combined/imnavait_normalized.rdata",
        "data/chemistry/soil_nutrients.rdata",
        "data/chemistry/porewater_chemistry.rdata"
    output:
        "figures/Figure_3 soil CN with points.pdf",
        "figures/Figure_4 pore water chemistry.pdf",
        "figures/Figure_7a dbRDA_16S_ASVs_all.pdf",
        "figures/Figure_7b dbRDA_ITS_ASVs_all.pdf",
        "figures/Figure_S5a dbRDA_16S_ASVs_CN.pdf",
        "figures/Figure_S5b dbRDA_ITS_ASVs_CN.pdf"
    shell:
        "Rscript code/soil_nutrients_and_porewater_chemistry.R"

rule alpha_div:
    input:
        "out/combined/imnavait_normalized.rdata"
    output:
        "figures/Figure_5a_shannon_lm.pdf",
        "figures/Figure_5b_invsimp_lm.pdf",
        "figures/Figure_A total_reads_vs_month_and_depth.pdf",
        "figures/Figure_B total_reads_vs_depth_and_month.pdf",
        "figures/Figure_C alpha_diversity_vs_total_reads_shannon.pdf",
        "figures/Figure_D alpha_diversity_vs_total_reads_invsimp.pdf",
        "figures/Figure_E diversity_vst_vs_reads.pdf"
    shell:
        "Rscript code/alpha_div.R"

rule beta_div:
    input:
        "out/combined/imnavait_normalized.rdata"
    output:
        "figures/Figure_6ab NMDS plots.pdf",
        "figures/Figure_6c sum of squares.pdf",
        "figures/Figure_S3 NMDS 16S UniFrac.pdf",
        "figures/Figure_S4ab NMDS plots sample condition.pdf"
    shell:
        "Rscript code/beta_div.R"

rule barplots:
    input:
        "out/combined/imnavait_normalized.rdata"
    output:
        "figures/Figure_S2ab taxon barplots.pdf"
    shell:
        "Rscript code/barplots.R"

rule deseq2:
    input:
        "out/combined/imnavait_normalized.rdata"
    output:
        "out/combined/deseq2.xlsx",
        "out/combined/deseq2.rds"
    shell:
        "Rscript code/deseq2.R"

################################
