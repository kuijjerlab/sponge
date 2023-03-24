# SPONGE

## Simple PriOr Network GEnerator

This repository contains the SPONGE package, which allows the generation 
of human prior gene regulatory networks based mainly on the data from 
the JASPAR database. It also uses HomoloGene to find the human analogs 
of vertebrate transcription factors, Ensemble to collect all the 
promoter regions in the human genome, UniProt for symbol matching, and
STRING to retrieve protein-protein interactions between transcription
factors.

Currently, it is restricted only to human networks and relatively narrow
range of settings, but further improvements are planned. A simple 
workflow would be as follows:

``` python
# Import the class definition
from sponge.sponge import Sponge
# Create the SPONGE object
sponge_obj = Sponge()
# Select the vertebrate transcription factors from JASPAR
sponge_obj.select_tfs()
# Find human homologs for the TFs if possible
sponge_obj.find_human_homologs()
# Filter the matches of the JASPAR bigbed file to the ones in the
# promoters of human transcripts
sponge_obj.filter_matches()
# Retrieve the protein-protein interactions between the transcription
# factors from the STRING database
sponge_obj.retrieve_ppi()
# Write the PPI prior to a file
sponge_obj.write_ppi_prior()
# Aggregate the filtered matches on promoters to genes
sponge_obj.aggregate_matches()
# Write the final motif prior to a file
sponge_obj.write_motif_prior()
```

SPONGE will attempt to download the files it needs into a temporary 
directory (`.sponge_temp` by default). Paths can be provided if these
files were downloaded in advance. The JASPAR bigbed file required for
filtering is huge (> 100 GB), so the download might take some time. Make
sure you're running SPONGE somewhere that has enough space!
