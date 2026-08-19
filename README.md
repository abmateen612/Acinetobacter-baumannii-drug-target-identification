# Acinetobacter baumannii Antibacterial Drug Target Identification

## Overview

Antimicrobial resistance is making infections caused by *Acinetobacter baumannii* increasingly difficult to treat. One approach to addressing this problem is to identify bacterial proteins that are essential for survival and sufficiently different from human proteins to serve as potential antibacterial drug targets.

This project presents a computational workflow for identifying and prioritizing potential antibacterial protein targets from the *A. baumannii* genome.

The workflow combines biological relevance, conservation, human homology exclusion, functional annotation, and pathway-level information to progressively narrow the genome down to a focused set of candidate targets.

## Research Objective

The main objective was to answer:

> **Which proteins in the *A. baumannii* genome are promising candidates for further antibacterial drug discovery?**

The analysis was designed to prioritize proteins that are biologically important to the bacterium while reducing candidates that are less suitable for antibacterial targeting.

## Target Prioritization Strategy

The workflow considered several factors during target identification and prioritization:

- Protein essentiality
- Conservation
- Human homology exclusion
- Functional annotation
- KEGG-based functional information
- Biological relevance to antibacterial drug discovery

After the initial filtering and annotation steps, the remaining proteins were evaluated based on their involvement in important bacterial functions.

## Priority Functional Categories

The final prioritization focused on three major antibacterial target categories:

### DNA / Replication / Repair

Proteins involved in maintaining and replicating bacterial genetic material.

### Cell Wall / Cell Envelope

Proteins associated with bacterial cell wall or envelope functions, which can provide attractive opportunities for selective antibacterial targeting.

### Protein Synthesis / Ribosome

Proteins involved in translation and bacterial protein synthesis.

## Redundancy Reduction

Multiple proteins can be associated with the same KEGG Orthology (KO) and may represent functionally redundant candidates.

To reduce redundancy, proteins associated with the same KO were narrowed down to representative candidates before the final prioritization.

Proteins associated with multiple core antibacterial functional categories were given higher priority because of their potential biological importance.

## Final Result

The complete workflow produced:

**15 shortlisted antibacterial protein targets**

These 15 candidates represent the final output of the target-identification and prioritization stage.

The shortlisted targets can be taken forward for subsequent computational drug discovery studies, including:

- Protein structure investigation
- Binding-site analysis
- Ligand screening
- Molecular docking
- Molecular dynamics
- Experimental validation

## Project Workflow

```text
A. baumannii Genome
        ↓
Essentiality Analysis
        ↓
Conservation Analysis
        ↓
Human Homology Exclusion
        ↓
Functional Annotation
        ↓
KEGG Functional Analysis
        ↓
Core Antibacterial Functions
        ↓
KO-based Redundancy Reduction
        ↓
Target Prioritization
        ↓
15 Shortlisted Targets
        ↓
Future CADD Studies
