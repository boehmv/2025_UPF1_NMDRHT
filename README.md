# Computational analyses for: Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome

------------------------------------------------------------------------

This repository contains the analyses scripts for the project: <br/> [**Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome**](https://doi.org/10.1101/2024.03.04.583328) <br/> (available as bioRxiv preprint)

## Scope

This repository primarily aims to provide transparent insight into the high-throughput analysis steps used in the study of UPF1 depletion via conditional degron tags (CDT) in human cells. :warning: **NOTE:** The complete pipeline is currently not optimized to run on different computing infrastructures in a standardized/portable manner. This means that all required packages have to be installed manually and configured accordingly to reproduce the results.

## Features / Requirements

-   Complete analysis of multiple RNA-Seq datasets (provided in FASTQ format; see [here](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/Metadata/UPF1_NMDRHT_datasets_forGitHub.csv) for dataset overview and [here](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/Metadata/UPF1_NMDRHT_datasets_experiments_forGitHub.csv) for individual sample identification), mapped to [Gencode v42](https://www.gencodegenes.org/human/release_42.html) / GRCh38.primary_assembly supplemented with SIRVomeERCCome (from Lexogen; [download](https://www.lexogen.com/wp-content/uploads/2018/08/SIRV_Set3_Sequences_170612a-ZIP.zip)) using STAR (short-read) or minimap2 (long-read).

-   Transcript quantification was performed using Salmon in mapping-based mode with a decoy-aware transcriptome index (either GENCODE.v42 or consolidated NMDRHT.v1.2) and the options --numGibbsSamples 30 --useVBOpt --gcBias --seqBias

-   Differential gene expression (DGE) was analyzed via DESeq2 or swish, differential transcript expression (DTE) via edgeR, differential transcript usage (DTU) via IsoformSwitchAnalyzeR and alternative splicing (AS) via LeafCutter.

-   The main pre-revision [CRSA_V009.sh](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/CRSA_V009.sh) or post-revision Bash script [CRSA_V010.sh](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/CRSA_V010.sh) runs the complete pipeline for standard short-read RNA-Seq data or individual modules using the options (see CRSA_V010.sh -h) and requires a design file specifying the following:

    -   reference type (gencode.v42.SIRVomeERCCome was used in this study)
    -   sequencing design (single- or paired-end reads)
    -   study name
    -   folder locations (srvdir for raw file locations, mydir for analyses output)
    -   location of the experiment file which specifies sample IDs and condition

-   Please see the provided [design.txt](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/Example/design.txt) file example for more information concerning this design file. An example for the tab-delimited [experiment.txt](https://github.com/boehmv/2025_UPF1_NMDRHT/blob/main/Analyses/Example/experiment.txt) file is provided as well. Please see the comments in CRSA_V009.sh or CRSA_V010.sh for further instructions

-   To run/reproduce the complete analysis script, many modules require specific tools. Please make sure you have the following tools installed and configured if required:

    -   [STAR](https://github.com/alexdobin/STAR) - version 2.7.10b was used for the analyses - with genome indices generated using GRCh38.primary.SIRVomeERCCome.fa and gencode.v42.SIRVomeERCCome.annotation.gtf (both reference files can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). The following code was used for genome index generation:

    ```         
    STAR   --runMode genomeGenerate   --runThreadN 15   --genomeDir /home/volker/reference/gencode.v42.SIRVomeERCCome  --genomeFastaFiles /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa --sjdbGTFfile /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf   --sjdbOverhang 99
    ```

    -   [Alfred](https://github.com/tobiasrausch/alfred) - version v0.2.6 was used for the analyses
    -   [samtools](http://www.htslib.org/) - version 1.16.1 (using htslib 1.16) was used for the analyses
    -   [IGV tools](http://software.broadinstitute.org/software/igv/download) - version 2.14.1 or 2.17.2 was used for the analyses - make sure you have the gencode.v42.SIRVomeERCCome.chrom.sizes file (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)) located in /PATH/TO/IGV/lib/genomes
    -   [Salmon](https://github.com/COMBINE-lab/salmon) - version v1.9.0 was used for the analyses - with an index generated using gentrome.v42.SIRV.ERCC.fa.gz and decoys.txt (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). A separate conda environment was created for Salmon. The following code was used for index generation:

    ```         
    salmon index -t /home/volker/reference/Gencode/gentrome.v42.SIRV.ERCC.fa.gz -d /home/volker/reference/Gencode/decoys.txt -p 12 -i /home/volker/reference/Transcriptome/gencode.v42.SIRVomeERCCome --gencode
    ```

    -   [DESeq2](https://github.com/mikelove/DESeq2) - version 1.40.1 was used for the analyses. The tx2gene file used for the analyses can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)
    -   [IsoformSwitchAnalyzeR](https://github.com/kvittingseerup/IsoformSwitchAnalyzeR) - version 1.18.0 was used for the analyses.
    -   [LeafCutter](https://github.com/davidaknowles/leafcutter) - version v0.2.9 was used for the analyse. :memo: **NOTE:** small changes in the /scripts of LeafCutter maintained gene IDs from Gencode (changed in [gtf_to_exons_vb.R](https://github.com/boehmv/2024-UPF1-degron/blob/main/scripts/LeafCutter/scripts/gtf_to_exons.R) and [leafcutter_ds.R](https://github.com/boehmv/2024-UPF1-degron/blob/main/scripts/LeafCutter/scripts/leafcutter_ds.R))
    -   [FastQC](https://github.com/s-andrews/FastQC) - version 0.11.9 was used for the analyses
    -   [MultiQC](https://github.com/ewels/MultiQC) - version v1.14 was used for the analyses

-   Additionally, many analyses were run using a plethora of R packages (including swish, edgeR, ...), please see the session info for the individual R scripts for more information.

-   All analyses were performed on a 16-core (2x Intel(R) Xeon(R) CPU E5-2687W v2 \@ 3.40GHz) workstation with 128 GB RAM running Ubuntu 22.04.2 LTS

-   Please make sure to change installation and file paths in the respective scripts to match your local environment

## Individual scripts

The specialized scripts called by the main CRSA_V009.sh (pre-revision) or CRSA_V010.sh (post-revision) script can be found [here](https://github.com/boehmv/2025_UPF1_NMDRHT/tree/main/Analyses). 

Config files and scripts for analyzing Ribo-Seq and long-read RNA-Seq data are provided in the folder as well.

## Feedback / Questions

Feedback is welcome! For any question, please email: [boehmv\@uni.koeln.de](mailto:boehmv@uni.koeln.de){.email} or [create an issue](https://github.com/boehmv/2025_UPF1_NMDRHT/issues)

## Citation

### Journal article

TBD

### bioRxiv preprint

Volker Boehm, Damaris Wallmeroth, Paul O. Wulf, Luiz Gustavo Teixeira Alves, Oliver Popp, Maximilian Riedel, Emanuel Wyler, Marek Franitza, Jennifer V. Gerbracht, Kerstin Becker, Karina Polkovnychenko, Simone Del Giudice, Nouhad Benlasfer, Philipp Mertins, Markus Landthaler, Niels H. Gehring (2024) **Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated transcriptome in human cells**. bioRxiv 2024.03.04.583328; doi: [https://doi.org/10.1101/2020.07.07.191437](https://doi.org/10.1101/2024.03.04.583328)
