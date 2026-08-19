# Acinetobacter baumannii Antibacterial Drug Target Identification

Computational identification and prioritization of potential antibacterial protein targets from the Acinetobacter baumannii genome.

## Project Overview

This project focuses on systematically narrowing down the A. baumannii genome to a focused set of potential antibacterial drug targets.

The analysis combines essentiality, conservation, human homology exclusion, functional annotation, and KEGG-based functional prioritization to identify proteins that are biologically relevant and suitable for further antibacterial research.

## Target Prioritization

The final prioritization was focused specifically on core antibacterial functions:

* DNA / replication / repair
* Cell wall / envelope
* Protein synthesis / ribosome

Redundant proteins associated with the same KEGG Orthology (KO) were reduced to representative candidates before the final target prioritization.

Proteins associated with multiple core antibacterial functions were given higher priority.

## Final Result

The complete workflow resulted in 15 shortlisted antibacterial target candidates.

These candidates represent the final output of the target-identification and prioritization stage and can be taken forward for subsequent computer-aided drug discovery (CADD) studies.

[View the final 15 targets]
(final_15_targets.csv)

[View the analysis code]
(target_prioritization.py)

## Project Status

Completed

The target identification and prioritization workflow has been completed, with 15 candidate antibacterial targets selected for further evaluation.

## Author

Abdul Mateen
MS Biotechnology
COMSATS University Abbottabad
