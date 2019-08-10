## Evolutionary selective constraints acting on the stop codon across land plants

#### **Table of Contents**
+ [Motivation](#motivation)  
+ [Installation](#install)
+ [File Descriptions](#filedesc)
+ [Usage & Flow](#usage)
+ [Output](#output)

### Motivation
<a name="motivation"/>

All genomes are under an evolutionary pressure and struggle to keep the functional portion
of the DNA. The rate at which favorable genes are retained and deleterious ones are lost is exerted by a
parameter which is the ratio of synonymous($dS$) to non synonymous($dN$) mutation rates. Substitutions do
not alter the coded amino acid while mutations do. This makes substitutions helpful and mutations harmful.
When $\frac{dN}{dS}$ < 1, substitution rate is greater than mutation rate and the gene is said to be under purifying selection. Purifying selection favors synonymous substitutions than non-synonymous mutations thereby preventing change of an amino acid residue at a give position. In conventional models of substitution only the sense(non-stop) codons are accounted for while the non-sense codons are omitted because they do not contribute to amino acid changes. Since stop codons function with varying efficiencies, they can be read-through and have the ability to alter the final protein products. When combined with other
mechanisms like ribosome stalling and mRNA regulation, stop codons can indirectly modulate protein
synthesis. This gives meaning to stop codon preservation and substitution, thereby creating the need
to include them in standard models of substitution. Stop codons have a low probability of undergoing
mutations but the pressure acting on their rate of substitution is only vaguely addressed. [*Seioghe et al.*](https://github.com/cseoighe/StopEvol) have introduced a new model which incorporates stop codons into the general Muse & Gaut substitution model. The model is constructed based on the assumption that stop codons are also under selection pressure and co-evolve with the genes. The extended Muse & Gaut model, casually called extMG model, has been applied on mammalian orthologous sequences and 50% of the genes were found to be under purifying selection. In this study the extMG model of substitution, for all 64 codons, is used to estimate $\phi$ (rate of substitution between stop codons) for plant ortholog families under the **Viridiplantae** clade. 

### Installation
<a name="install"/> 

+ ##### Requires
    + [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
    + gawk
    + [seqkit](https://bioinf.shenwei.me/seqkit/download/)
    + [codonphyml](https://sourceforge.net/projects/codonphyml/)


1. [**Download the archive as a zip**](https://github.com/vizkidd/stop_codon_plants/archive/master.zip) and extract it
2. These flat files from OrthoDB are required to be in same directory as the extracted files
    + [levels.tab](https://v100.orthodb.org/download/odb10v0_levels.tab.gz)
    + [species.tab](https://v100.orthodb.org/download/odb10v0_species.tab.gz)
    + [levels2species.tab](https://v100.orthodb.org/download/odb10v0_level2species.tab.gz)
    + [OG2genes.tab](https://v100.orthodb.org/download/odb10v0_OG2genes.tab.gz)
    + [OGs.tab](https://v100.orthodb.org/download/odb10v0_OGs.tab.gz)

> **NOTE:**

 + The files are from OrthoDB version 10. But any other versions which are compatible and have the same file structure can be used.
 + Make sure to extract and rename the flat files to the names provided above.

### File Descriptions
<a name="filedesc"/> 

##### ***Flat Files***

  + [levels.tab](https://v100.orthodb.org/download/odb10v0_levels.tab.gz) - File with NCBI taxonomic nodes [*first column*] and node names [*second column*].
      ![][levels] 
  + [species.tab](https://v100.orthodb.org/download/odb10v0_species.tab.gz) - Contains organism IDs [*second column*] and organism names [*third column*].
      ![][species]
  + [OGs.tab](https://v100.orthodb.org/download/odb10v0_OGs.tab.gz) - Contains information about orthologous groups (***OG***), OG ID [*first column*] & OG name [*last column*]. OG ID has a format of [*cluster_ID*]at[*taxa_node*]
      ![][OGs] 
  + [levels2species.tab](https://v100.orthodb.org/download/odb10v0_level2species.tab.gz) - Connects NCBI taxa ID [*first column*] to organism IDs (same as species IDs) [*second column*] and also provides information about number of hops & NCBI taxonomic levels [*last column*]
      ![][levels2species] 
  + [OG2genes.tab](https://v100.orthodb.org/download/odb10v0_OG2genes.tab.gz) - Connects OGs to genes. Contains OG IDs [*first column*] & gene IDs [*last column*]. Each gene ID is of the format [*organism_ID*]:[*gene_ID*]. (This can be used to accumulate clusters based on organisms or genes.)
      ![][OG2genes]

##### ***Scripts***

Main scripts are

+ [START.sh](START.sh) - Starts the pipeline. Selects clusters based on the organism_ID and cluster_ID provided. Passes the clusters one at a time to [download_data.sh](download_data.sh).
+ [download_data.sh](download_data.sh) - Download the CDS of gene_IDs in each cluster and selects one sequence for each organism.
+ [process_data.sh](process_data.sh) - Process the data. Applies the extended Muse & Gaut model.

Run these after the clusters are downloaded (or if cluster download is stopped)

+ [MIXMOD_BOOTSTRAP.sh](MIXMOD_BOOTSTRAP.sh) - Applies the mixture model and performs bootstrapping

The overall flow including the misc scripts goes like this


+ [START.sh](START.sh) &#8595;
    + [extract_org_ids.py](extract_org_ids.py)
    + [get_cluster_counts.py](get_cluster_counts.py)
    + [names2dict.py](names2dict.py)
    + [download_data.sh](download_data.sh) &#8595;
        + [orthodb_api.py](orthodb_api.py)
        + [cds_from_ncbi.py](cds_from_ncbi.py) **
        + [remove_empty_sequences.py](remove_empty_sequences.py)
        + [select_sequences.py](select_sequences.py)
        + [process_data.sh](process_data.sh) &#8595;
            + [trim_fasta_header.py](trim_fasta_header.py)
            + [align_seqs_init.R](align_seqs_init.R)
            + [stop_codon_locations.py](stop_codon_locations.py)
            + [extract_cluster_sequences.py](extract_cluster_sequences.py)
            + [align_seqs.R](align_seqs.R)
            + [stopcodon.R](stopcodon.R)
            + [cut_stops.py](cut_stops.py)
            + [fasta2relaxedPhylip.pl](fasta2relaxedPhylip.pl)
            + [parameter_formatter.py](parameter_formatter.py) .
+ [MIXMOD_BOOTSTRAP.sh](MIXMOD_BOOTSTRAP.sh) &#8595;
    + [mixture_model.rscript](mixture_model.rscript)
    + [bootstrap.r](bootstrap.r)
    + [plot.r](plot.r)

> **NOTE:** Scripts were written for *slurm*, if you don't have slurm modify the sbatch lines in [START.sh](START.sh), [download_data.sh](download_data.sh), [process_data.sh](process_data.sh) and just call the command without sbatch.

> **NOTE:** **Scripts for Phytomine and ensemble have also been provided, check Misc

### Usage & Flow
<a name="usage"/> 

![Fig. *Taxonomic level coverage of example orthologous groups. <span style="color: #ff0000">Organisms</span> and the <span style="color: #00ff00">clusters</span> in which they part.*][taxa]

Pipeline can be started using the START.sh script,

```{bash}
sh START.sh <organism_ID> <upper_limit_scale> <cluster_ID> <min_orgs>
```

eg,

>sh START.sh <span style="color: #ff0000">3193</span> 10 <span style="color: #00ff00">33090</span> 20

Both the <span style="color: #ff0000">organism</span> and the <span style="color: #00ff00">cluster</span> IDs can be selected from the first column of the levels.tab file. For this study, organisms in <span style="color: #ff0000">Embryophyta</span> and the clusters in <span style="color: #00ff00">Viridiplantae</span> are selected.

The pipeline can

+ For each cluster
    + Downloads CDS (coding sequences) from NCBI for the cluster.
    + Selects one sequence for each organism in the cluster.
    + Applies the extended Muse & Gaut model which estimates the $\kappa$, $\omega$ ,$\phi$ & *treescale* for the cluster.
+ For all clusters
    + Applies the mixture model
    + Performs bootstrapping

*On an average, each sequence takes 4~5 seconds to download in order to respect the NCBI query laws.*

[taxa]: Figures/taxa.png "Coverage of different NCBI Taxanomic Levels"
[levels]: Figures/levels.png "OrthoDB levels.tab file"
[levels2species]: Figures/levels2species.png "OrthoDB levels2species.tab file"
[OGs]: Figures/OGs.png "OrthoDB OGs.tab file"
[species]: Figures/species.png "OrthoDB species.tab file"
[OG2genes]: Figures/OG2genes.png "OrthoDB OG2genes.tab file"

### Output
<a name="output"/> 

*Under Construction*