# =============================================================================
# # Phase 2: Genome annotation
# =============================================================================

import pyrodigal
from Bio import SeqIO
import pandas as pd

# ============================================================
# 1. INPUT GENOME
# ============================================================

genome_path = r"E:\Genomics\A.baumanii_genome\GCF_009035845.1\GCF_009035845.1.fna"

# Output files
protein_file = r"E:\Genomics\A.baumanii_genome\A_baumannii_proteins.faa"
annotation_file = r"E:\Genomics\A.baumanii_genome\A_baumannii_gene_annotation.csv"


# ============================================================
# 2. READ THE COMPLETE GENOME
# ============================================================

records = list(SeqIO.parse(genome_path, "fasta"))

print("Number of contigs:", len(records))

total_genome_length = sum(len(record.seq) for record in records)

print("Total genome length:", total_genome_length, "bp")


# ============================================================
# 3. GENE PREDICTION
# ============================================================

finder = pyrodigal.GeneFinder(meta=True)

all_genes = []

for record in records:

    sequence = str(record.seq)

    genes = finder.find_genes(sequence)

    for gene_number, gene in enumerate(genes, 1):

        all_genes.append({
            "contig": record.id,
            "gene_number": gene_number,
            "gene": gene
        })


print("Total predicted CDS:", len(all_genes))


# ============================================================
# 4. SAVE PROTEIN SEQUENCES
# ============================================================

with open(protein_file, "w") as f:

    for i, item in enumerate(all_genes, 1):

        gene = item["gene"]

        protein = gene.translate()

        protein_id = f"AB_gene_{i}"

        f.write(f">{protein_id}\n")
        f.write(f"{protein}\n")


# ============================================================
# 5. CREATE ANNOTATION TABLE
# ============================================================

annotation = []

for i, item in enumerate(all_genes, 1):

    gene = item["gene"]

    protein = gene.translate()

    protein_id = f"AB_gene_{i}"

    # Pyrodigal coordinates are 0-based internally
    start = gene.begin + 1
    end = gene.end

    strand = "+" if gene.strand == 1 else "-"

    annotation.append({
        "protein_id": protein_id,
        "contig": item["contig"],
        "gene_number": item["gene_number"],
        "start": start,
        "end": end,
        "strand": strand,
        "CDS_length_bp": end - start + 1,
        "protein_length_aa": len(protein)
    })


df = pd.DataFrame(annotation)

df.to_csv(annotation_file, index=False)


# ============================================================
# 6. BASIC GENOME/ANNOTATION SUMMARY
# ============================================================

print("\n========== ANNOTATION SUMMARY ==========")

print("Genome length:", total_genome_length, "bp")
print("Number of contigs:", len(records))
print("Predicted CDS:", len(df))
print("Average protein length:",
      round(df["protein_length_aa"].mean(), 2), "aa")

print("\nFiles created:")
print(protein_file)
print(annotation_file)

print("\nFirst 5 genes:")
print(df.head())





# =============================================================================
# ## Phase: 3 :Gene Essential analysis 
# =============================================================================
import os

# ============================================================
# STEP 1 — CHECK INPUT FILES
# ============================================================

# CHANGE THIS:
# Complete path to A. baumannii protein FASTA file
protein_file = r"E:\Genomics\A.baumanii_genome\A_baumannii_proteins.faa"

# CHANGE THIS:
# Complete path to DEG protein FASTA file
deg_file = r"E:\Genomics\A.baumanii_genome\DEG\deg_proteins.fasta"


# ============================================================
# CHECK FILES
# ============================================================

print("A. baumannii protein file:")
print(protein_file)

print("\nDEG protein file:")
print(deg_file)

print("\nFILE CHECK")
print("--------------------------------")

print(
    "A. baumannii proteins:",
    os.path.exists(protein_file)
)

print(
    "DEG proteins:",
    os.path.exists(deg_file)
)


# ============================================================
# STOP IF FILES ARE MISSING
# ============================================================

if not os.path.exists(protein_file):

    raise FileNotFoundError(
        "A. baumannii protein file not found. "
        "Check protein_file."
    )

if not os.path.exists(deg_file):

    raise FileNotFoundError(
        "DEG protein file not found. "
        "Check deg_file."
    )


print("\nBoth protein files are ready.")




# ============================================================
# PHASE 3 — ESSENTIALITY ANALYSIS
# STEP 2 — BLASTp
# ============================================================

import os
import subprocess


# ============================================================
# CHANGE THESE LOCATIONS ONLY IF NEEDED
# ============================================================

# A. baumannii protein FASTA
protein_file = r"E:\Genomics\A.baumanii_genome\A_baumannii_proteins.faa"

# DEG protein FASTA
deg_file = r"E:\Genomics\A.baumanii_genome\DEG\deg_proteins.fasta"

# BLASTp executable
blastp = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastp.exe"

# BLAST result file
blast_output = r"E:\Genomics\A.baumanii_genome\A_baumannii_vs_DEG.tsv"


# ============================================================
# CHECK FILES
# ============================================================

print("FILE CHECK")
print("--------------------------------")

print("A. baumannii proteins:", os.path.exists(protein_file))
print("DEG proteins:", os.path.exists(deg_file))
print("BLASTp:", os.path.exists(blastp))


# ============================================================
# STOP IF SOMETHING IS MISSING
# ============================================================

if not os.path.exists(protein_file):
    raise FileNotFoundError(
        "A. baumannii protein file not found."
    )

if not os.path.exists(deg_file):
    raise FileNotFoundError(
        "DEG protein file not found."
    )

if not os.path.exists(blastp):
    raise FileNotFoundError(
        "BLASTp executable not found."
    )


# ============================================================
# RUN BLASTp
# ============================================================

print("\nRunning BLASTp...")
print("Wait Babu Bhaiya Wait.")


command = [
    blastp,

    "-query",
    protein_file,

    "-subject",
    deg_file,

    "-out",
    blast_output,

    "-outfmt",
    "6 qseqid sseqid pident length "
    "mismatch gapopen qstart qend "
    "sstart send evalue bitscore",

    "-evalue",
    "1e-5",

    "-max_target_seqs",
    "10"
]


result = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True
)


# ============================================================
# RESULT
# ============================================================

if result.stdout:
    print(result.stdout)

if result.stderr:
    print(result.stderr)


if result.returncode == 0:

    print("\nBLASTp completed successfully.")
    print("Results saved to:")
    print(blast_output)

else:

    print("\nBLASTp failed.")
    print("Return code:", result.returncode)



# ============================================================
# PHASE 3 — ESSENTIALITY ANALYSIS
# STEP 3 — ANALYZE BLASTp RESULTS
# ============================================================

import pandas as pd


# ============================================================
# CHANGE THIS:
# BLASTp result file from Step 2
# ============================================================

blast_output = r"E:\Genomics\A.baumanii_genome\A_baumannii_vs_DEG.tsv"


# ============================================================
# BLAST RESULT COLUMN NAMES
# ============================================================

columns = [
    "query_id",
    "subject_id",
    "identity",
    "alignment_length",
    "mismatch",
    "gap_open",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore"
]





# ============================================================
# READ BLAST RESULTS
# ============================================================

blast_results = pd.read_csv(
    blast_output,
    sep="\t",
    names=columns
)


# ============================================================
# SUMMARY
# ============================================================

print("BLAST results loaded successfully.")
print("--------------------------------")

print(
    "Total BLAST matches:",
    len(blast_results)
)

print(
    "A. baumannii proteins with matches:",
    blast_results["query_id"].nunique()
)


# ============================================================
# SHOW FIRST 5 RESULTS
# ============================================================

print("\nFirst 5 BLAST matches:")

print(
    blast_results.head()
)



# ============================================================
# PHASE 3 — ESSENTIALITY ANALYSIS
# STEP 4 — FILTER STRONG DEG MATCHES
# ============================================================

from Bio import SeqIO


# ============================================================
# GET LENGTH OF EACH A. BAUMANNII PROTEIN
# ============================================================

protein_lengths = {}

for record in SeqIO.parse(protein_file, "fasta"):
    protein_lengths[record.id] = len(record.seq)


# ============================================================
# ADD PROTEIN LENGTH
# ============================================================

blast_results["query_length"] = (
    blast_results["query_id"].map(protein_lengths)
)


# ============================================================
# CALCULATE ALIGNMENT LENGTH
# ============================================================

blast_results["query_aligned_length"] = (
    blast_results["query_end"]
    - blast_results["query_start"]
    + 1
)


# ============================================================
# CALCULATE QUERY COVERAGE
# ============================================================

blast_results["query_coverage"] = (
    blast_results["query_aligned_length"]
    / blast_results["query_length"]
) * 100


# ============================================================
# FILTER STRONG MATCHES
# ============================================================

essential_candidates = blast_results[
    (blast_results["identity"] >= 30) &
    (blast_results["query_coverage"] >= 70) &
    (blast_results["evalue"] <= 1e-5)
].copy()


# ============================================================
# KEEP BEST DEG MATCH FOR EACH PROTEIN
# ============================================================

essential_candidates = (
    essential_candidates
    .sort_values(
        ["query_id", "bitscore"],
        ascending=[True, False]
    )
    .drop_duplicates(
        subset="query_id"
    )
)


# ============================================================
# SUMMARY
# ============================================================

print("ESSENTIALITY FILTERING")
print("--------------------------------")

print(
    "Initial BLAST matches:",
    len(blast_results)
)

print(
    "Proteins with BLAST matches:",
    blast_results["query_id"].nunique()
)

print(
    "Essentiality candidates:",
    len(essential_candidates)
)


# ============================================================
# SHOW TOP 20
# ============================================================

print("\nTop candidates:")

print(
    essential_candidates[
        [
            "query_id",
            "subject_id",
            "identity",
            "query_coverage",
            "evalue",
            "bitscore"
        ]
    ].head(20)
)



# ============================================================
# PHASE 3 — ESSENTIALITY ANALYSIS
# STEP 5 — SAVE ESSENTIALITY RESULTS
# ============================================================

# CHANGE THIS:
essential_output = r"E:\Genomics\A.baumanii_genome\essential_protein_candidates.csv"


# ============================================================
# SAVE
# ============================================================

essential_candidates.to_csv(
    essential_output,
    index=False
)


# ============================================================
# CONFIRM
# ============================================================

print("Essentiality analysis completed.")
print("--------------------------------")
print(
    "Essential protein candidates:",
    len(essential_candidates)
)
print("Results saved to:")
print(essential_output)




# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 1 — FIND PROTEIN FASTA FILES
# ============================================================

# CHANGE THIS ONLY IF YOUR FOLDER LOCATION CHANGES
strain_folder = r"E:\Genomics\A.baumanii_genome\Other_strains\ncbi_dataset"


# ============================================================
# FIND ALL protein.faa FILES
# ============================================================

protein_files = []

for root, dirs, files in os.walk(strain_folder):

    for file in files:

        if file.lower() == "protein.faa":

            protein_files.append(
                os.path.join(root, file)
            )


# ============================================================
# DISPLAY RESULTS
# ============================================================

print("PROTEIN FASTA FILES FOUND")
print("--------------------------------")

print(
    "Number of protein.faa files:",
    len(protein_files)
)

for i, file in enumerate(protein_files, 1):

    print(
        f"{i}. {file}"
    )


# ============================================================
# CHECK
# ============================================================

if len(protein_files) == 0:

    raise FileNotFoundError(
        "No protein.faa files were found."
    )

print("\nProtein FASTA files successfully located.")



# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 2 — COMBINE 10 STRAIN PROTEIN FASTA FILES
# ============================================================

# CHANGE THIS ONLY IF NEEDED
combined_proteins = r"E:\Genomics\A.baumanii_genome\Other_strains\A_baumannii_10_strains_combined.faa"


# ============================================================
# COMBINE PROTEINS
# ============================================================

total_proteins = 0

with open(combined_proteins, "w") as output:

    for protein_file_path in protein_files:

        # Get GCF/GCA accession from the folder name
        accession = os.path.basename(
            os.path.dirname(protein_file_path)
        )

        # Read proteins from this strain
        for record in SeqIO.parse(protein_file_path, "fasta"):

            # Add strain accession to protein ID
            record.id = accession + "|" + record.id
            record.description = ""

            SeqIO.write(
                record,
                output,
                "fasta"
            )

            total_proteins += 1


# ============================================================
# CHECK RESULT
# ============================================================

print("10-strain protein database created.")
print("--------------------------------")
print(
    "Strains combined:",
    len(protein_files)
)

print(
    "Total protein sequences:",
    total_proteins
)

print("Saved to:")
print(combined_proteins)




# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 2.5 — CREATE ESSENTIAL CANDIDATE PROTEIN FASTA
# ============================================================

candidate_proteins = r"E:\Genomics\A.baumanii_genome\essential_candidate_proteins.faa"


# ============================================================
# GET CANDIDATE PROTEIN IDs
# ============================================================

candidate_ids = set(
    essential_candidates["query_id"]
    .astype(str)
)


print("Candidate proteins to extract:", len(candidate_ids))


# ============================================================
# LOAD ORIGINAL PROTEIN FASTA
# ============================================================

protein_records = SeqIO.to_dict(
    SeqIO.parse(protein_file, "fasta")
)


print("Proteins in original FASTA:", len(protein_records))


# ============================================================
# EXTRACT CANDIDATE PROTEINS
# ============================================================

found_ids = []
missing_ids = []

with open(candidate_proteins, "w") as output:

    for protein_id in candidate_ids:

        if protein_id in protein_records:

            SeqIO.write(
                protein_records[protein_id],
                output,
                "fasta"
            )

            found_ids.append(protein_id)

        else:

            missing_ids.append(protein_id)


# ============================================================
# CHECK RESULT
# ============================================================

print("\nCANDIDATE FASTA CREATION")
print("--------------------------------")

print("Candidates requested:", len(candidate_ids))
print("Sequences found:", len(found_ids))
print("Sequences missing:", len(missing_ids))

print("\nSaved to:")
print(candidate_proteins)


if len(missing_ids) == 0:

    print("\n✓ All candidate proteins successfully extracted.")

else:

    print("\n⚠ Missing protein IDs:")
    print(missing_ids[:20])



# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 3 — BLASTp AGAINST 10 OTHER STRAINS
# ============================================================

candidate_proteins = r"E:\Genomics\A.baumanii_genome\essential_candidate_proteins.faa"

combined_proteins = r"E:\Genomics\A.baumanii_genome\Other_strains\A_baumannii_10_strains_combined.faa"

blastp = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastp.exe"

conservation_blast_output = r"E:\Genomics\A.baumanii_genome\A_baumannii_conservation_BLAST.tsv"


# ============================================================
# CHECK FILES
# ============================================================

print("FILE CHECK")
print("--------------------------------")

print("Candidate proteins:", os.path.exists(candidate_proteins))
print("10-strain proteins:", os.path.exists(combined_proteins))
print("BLASTp:", os.path.exists(blastp))


# ============================================================
# RUN BLASTp
# ============================================================

print("\nRunning conservation BLASTp...")
print("1,607 candidates vs 10-strain protein database")
print("Please wait...")


command = [
    blastp,
    "-query", candidate_proteins,
    "-subject", combined_proteins,
    "-out", conservation_blast_output,
    "-outfmt",
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore",
    "-evalue", "1e-5",
    "-max_target_seqs", "10"
]


result = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True
)


# ============================================================
# CHECK BLAST RESULT
# ============================================================

if result.returncode == 0:

    print("\n✓ Conservation BLASTp completed successfully.")
    print("Results saved to:")
    print(conservation_blast_output)

    print(
        "\nOutput file exists:",
        os.path.exists(conservation_blast_output)
    )

else:

    print("\n✗ Conservation BLASTp failed.")
    print("Return code:", result.returncode)
    print("\nError:")
    print(result.stderr)




# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 4 — ANALYZE CONSERVATION BLASTp RESULTS
# ============================================================

# ============================================================
# BLAST RESULT COLUMN NAMES
# ============================================================

conservation_columns = [
    "query_id",
    "subject_id",
    "identity",
    "alignment_length",
    "mismatch",
    "gap_open",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore"
]


# ============================================================
# READ BLAST RESULTS
# ============================================================

conservation_results = pd.read_csv(
    conservation_blast_output,
    sep="\t",
    names=conservation_columns
)


print("CONSERVATION BLAST RESULTS")
print("--------------------------------")

print("Total BLAST matches:",
      len(conservation_results))

print("Candidate proteins with matches:",
      conservation_results["query_id"].nunique())


# ============================================================
# GET CANDIDATE PROTEIN LENGTHS
# ============================================================

candidate_lengths = {}

for record in SeqIO.parse(candidate_proteins, "fasta"):
    candidate_lengths[record.id] = len(record.seq)


conservation_results["query_length"] = (
    conservation_results["query_id"].map(candidate_lengths)
)


# ============================================================
# CALCULATE QUERY COVERAGE
# ============================================================

conservation_results["query_aligned_length"] = (
    conservation_results["query_end"]
    - conservation_results["query_start"]
    + 1
)

conservation_results["query_coverage"] = (
    conservation_results["query_aligned_length"]
    / conservation_results["query_length"]
) * 100


# ============================================================
# IDENTIFY STRAIN
# ============================================================

# Combined FASTA IDs have:
# accession|original_protein_id

conservation_results["strain"] = (
    conservation_results["subject_id"]
    .astype(str)
    .str.split("|")
    .str[0]
)


# ============================================================
# FILTER STRONG CONSERVATION HITS
# ============================================================

strong_conservation = conservation_results[
    (conservation_results["identity"] >= 30) &
    (conservation_results["query_coverage"] >= 70) &
    (conservation_results["evalue"] <= 1e-5)
].copy()


print("\nSTRONG CONSERVATION HITS")
print("--------------------------------")

print("Strong BLAST matches:",
      len(strong_conservation))

print("Candidate proteins with strong matches:",
      strong_conservation["query_id"].nunique())

print("Strains detected:",
      strong_conservation["strain"].nunique())


# ============================================================
# KEEP BEST HIT FOR EACH CANDIDATE × STRAIN
# ============================================================

best_hits = (
    strong_conservation
    .sort_values(
        ["query_id", "strain", "bitscore"],
        ascending=[True, True, False]
    )
    .drop_duplicates(
        subset=["query_id", "strain"]
    )
)


# ============================================================
# COUNT CONSERVED STRAINS
# ============================================================

strain_counts = (
    best_hits
    .groupby("query_id")["strain"]
    .nunique()
    .reset_index(name="conserved_strains")
)


# ============================================================
# CREATE COMPLETE CONSERVATION SUMMARY
# ============================================================

conservation_summary = (
    pd.DataFrame({
        "query_id": list(candidate_lengths.keys())
    })
    .merge(
        strain_counts,
        on="query_id",
        how="left"
    )
)


conservation_summary["conserved_strains"] = (
    conservation_summary["conserved_strains"]
    .fillna(0)
    .astype(int)
)


# ============================================================
# CORE CONSERVED PROTEINS
# ============================================================

total_strains = len(protein_files)

core_conserved = conservation_summary[
    conservation_summary["conserved_strains"] == total_strains
].copy()


# ============================================================
# SUMMARY
# ============================================================

print("\n========================================")
print("CONSERVATION SUMMARY")
print("========================================")

print("Total candidate proteins:",
      len(conservation_summary))

print("Other strains:",
      total_strains)

print("Candidates conserved in all strains:",
      len(core_conserved))

print("Candidates conserved in at least 8 strains:",
      len(
          conservation_summary[
              conservation_summary["conserved_strains"] >= 8
          ]
      ))

print("Candidates conserved in at least 5 strains:",
      len(
          conservation_summary[
              conservation_summary["conserved_strains"] >= 5
          ]
      ))


# ============================================================
# CONSERVATION DISTRIBUTION
# ============================================================

print("\nConservation distribution:")

print(
    conservation_summary[
        "conserved_strains"
    ].value_counts().sort_index()
)


# ============================================================
# SHOW FIRST 20 CORE CONSERVED PROTEINS
# ============================================================

print("\nFirst 20 core conserved proteins:")

print(core_conserved.head(20))



# ============================================================
# PHASE 4 — CONSERVATION ANALYSIS
# STEP 5 — CREATE CORE CONSERVED PROTEIN FASTA
# ============================================================

core_conserved_proteins = (
    r"E:\Genomics\A.baumanii_genome\conserved_protein\core_conserved_proteins.faa"
)


# ============================================================
# GET CORE CONSERVED PROTEIN IDs
# ============================================================

core_ids = set(
    core_conserved["query_id"]
    .astype(str)
)


print("Core conserved proteins:", len(core_ids))


# ============================================================
# EXTRACT SEQUENCES
# ============================================================

core_found = []
core_missing = []

with open(core_conserved_proteins, "w") as output:

    for protein_id in core_ids:

        if protein_id in candidate_lengths:

            # Retrieve sequence from original candidate FASTA
            record = protein_records[protein_id]

            SeqIO.write(
                record,
                output,
                "fasta"
            )

            core_found.append(protein_id)

        else:

            core_missing.append(protein_id)


# ============================================================
# CHECK RESULT
# ============================================================

print("\nCORE CONSERVED FASTA")
print("--------------------------------")

print("Requested:", len(core_ids))
print("Sequences written:", len(core_found))
print("Missing:", len(core_missing))

print("\nSaved to:")
print(core_conserved_proteins)


if len(core_missing) == 0:

    print("\n✓ All 1,471 core conserved proteins successfully saved.")

else:

    print("\n⚠ Missing proteins:")
    print(core_missing[:20])
    
    
    
    
    # ============================================================
# PHASE 5 — HUMAN HOMOLOGY ANALYSIS
# STEP 1 — FIND HUMAN PROTEIN FASTA
# ============================================================

human_folder = r"E:\Genomics\A.baumanii_genome\human_proteome\ncbi_dataset\GCF_000001405.40"


# ============================================================
# FIND FASTA FILES
# ============================================================

human_fasta_files = []

for root, dirs, files in os.walk(human_folder):

    for file in files:

        if file.lower().endswith((".faa", ".fasta", ".fa")):

            human_fasta_files.append(
                os.path.join(root, file)
            )


# ============================================================
# DISPLAY RESULTS
# ============================================================

print("HUMAN PROTEIN FASTA FILES")
print("--------------------------------")

print(
    "Number of FASTA files found:",
    len(human_fasta_files)
)

for i, file in enumerate(human_fasta_files, 1):

    print(f"{i}. {file}")


# ============================================================
# CHECK
# ============================================================

if len(human_fasta_files) == 0:

    raise FileNotFoundError(
        "No human protein FASTA file found."
    )

print("\n✓ Human protein FASTA successfully located.")





# ============================================================
# PHASE 5 — HUMAN HOMOLOGY ANALYSIS
# STEP 2 — BLASTp AGAINST HUMAN PROTEOME
# ============================================================

human_proteins = r"E:\Genomics\A.baumanii_genome\human_proteome\ncbi_dataset\GCF_000001405.40\protein.faa"

human_blast_output = r"E:\Genomics\A.baumanii_genome\human_homology_BLAST.tsv"


# ============================================================
# CHECK FILES
# ============================================================

print("FILE CHECK")
print("--------------------------------")

print(
    "Core conserved proteins:",
    os.path.exists(core_conserved_proteins)
)

print(
    "Human proteins:",
    os.path.exists(human_proteins)
)

print(
    "BLASTp:",
    os.path.exists(blastp)
)


# ============================================================
# STOP IF FILE MISSING
# ============================================================

if not os.path.exists(core_conserved_proteins):
    raise FileNotFoundError(
        "Core conserved protein FASTA not found."
    )

if not os.path.exists(human_proteins):
    raise FileNotFoundError(
        "Human protein FASTA not found."
    )

if not os.path.exists(blastp):
    raise FileNotFoundError(
        "BLASTp executable not found."
    )


# ============================================================
# RUN BLASTp
# ============================================================

print("\nRunning BLASTp against human proteome...")
print("1,471 A. baumannii proteins vs human proteins")
print("Wait please...")


command = [
    blastp,

    "-query",
    core_conserved_proteins,

    "-subject",
    human_proteins,

    "-out",
    human_blast_output,

    "-outfmt",
    "6 qseqid sseqid pident length "
    "mismatch gapopen qstart qend "
    "sstart send evalue bitscore",

    "-evalue",
    "1e-5",

    "-max_target_seqs",
    "10"
]


result = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True
)


# ============================================================
# RESULT
# ============================================================

if result.returncode == 0:

    print("\n✓ Human homology BLASTp completed successfully.")

    print("Results saved to:")
    print(human_blast_output)

    print(
        "\nOutput file exists:",
        os.path.exists(human_blast_output)
    )

else:

    print("\n✗ Human homology BLASTp failed.")

    print(
        "Return code:",
        result.returncode
    )

    print("\nBLAST error:")
    print(result.stderr)



# ============================================================
# PHASE 5 — HUMAN HOMOLOGY ANALYSIS
# STEP 3 — ANALYZE HUMAN BLASTp RESULTS
# ============================================================

human_columns = [
    "query_id",
    "subject_id",
    "identity",
    "alignment_length",
    "mismatch",
    "gap_open",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore"
]


# ============================================================
# READ HUMAN BLAST RESULTS
# ============================================================

human_results = pd.read_csv(
    human_blast_output,
    sep="\t",
    names=human_columns
)


print("HUMAN BLAST RESULTS")
print("--------------------------------")

print(
    "Total BLAST matches:",
    len(human_results)
)

print(
    "Core conserved proteins with matches:",
    human_results["query_id"].nunique()
)


# ============================================================
# GET QUERY PROTEIN LENGTHS
# ============================================================

human_query_lengths = {}

for record in SeqIO.parse(
    core_conserved_proteins,
    "fasta"
):

    human_query_lengths[record.id] = len(record.seq)


human_results["query_length"] = (
    human_results["query_id"]
    .map(human_query_lengths)
)


# ============================================================
# CALCULATE QUERY COVERAGE
# ============================================================

human_results["query_aligned_length"] = (
    human_results["query_end"]
    - human_results["query_start"]
    + 1
)


human_results["query_coverage"] = (
    human_results["query_aligned_length"]
    / human_results["query_length"]
) * 100


# ============================================================
# FILTER POTENTIAL HUMAN HOMOLOGS
# ============================================================

human_homologs = human_results[
    (human_results["identity"] >= 30) &
    (human_results["query_coverage"] >= 70) &
    (human_results["evalue"] <= 1e-5)
].copy()


# ============================================================
# KEEP BEST HUMAN HIT PER BACTERIAL PROTEIN
# ============================================================

human_homologs = (
    human_homologs
    .sort_values(
        ["query_id", "bitscore"],
        ascending=[True, False]
    )
    .drop_duplicates(
        subset="query_id"
    )
)


# ============================================================
# IDENTIFY BACTERIA-SPECIFIC CANDIDATES
# ============================================================

all_core_ids = set(
    core_conserved["query_id"]
    .astype(str)
)

human_homolog_ids = set(
    human_homologs["query_id"]
    .astype(str)
)


bacteria_specific_ids = (
    all_core_ids - human_homolog_ids
)


# ============================================================
# SUMMARY
# ============================================================

print("\n========================================")
print("HUMAN HOMOLOGY SUMMARY")
print("========================================")

print(
    "Core conserved proteins:",
    len(all_core_ids)
)

print(
    "Potential human homologs:",
    len(human_homolog_ids)
)

print(
    "Bacteria-specific candidates:",
    len(bacteria_specific_ids)
)


# ============================================================
# SHOW HUMAN HOMOLOGS
# ============================================================

print("\nFirst 20 potential human homologs:")

print(
    human_homologs[
        [
            "query_id",
            "subject_id",
            "identity",
            "query_coverage",
            "evalue",
            "bitscore"
        ]
    ].head(20)
)


# ============================================================
# PHASE 5 — HUMAN HOMOLOGY ANALYSIS
# STEP 4 — CREATE BACTERIA-SPECIFIC PROTEIN FASTA
# ============================================================

bacteria_specific_proteins = (
    r"E:\Genomics\A.baumanii_genome\bacteria_specific_proteins.faa"
)


# ============================================================
# EXTRACT BACTERIA-SPECIFIC SEQUENCES
# ============================================================

found_ids = []
missing_ids = []

with open(bacteria_specific_proteins, "w") as output:

    for protein_id in bacteria_specific_ids:

        if protein_id in protein_records:

            SeqIO.write(
                protein_records[protein_id],
                output,
                "fasta"
            )

            found_ids.append(protein_id)

        else:

            missing_ids.append(protein_id)


# ============================================================
# CHECK RESULT
# ============================================================

print("BACTERIA-SPECIFIC FASTA")
print("--------------------------------")

print(
    "Bacteria-specific candidates:",
    len(bacteria_specific_ids)
)

print(
    "Sequences written:",
    len(found_ids)
)

print(
    "Missing:",
    len(missing_ids)
)

print("\nSaved to:")
print(bacteria_specific_proteins)


if len(missing_ids) == 0:

    print(
        "\n✓ All bacteria-specific candidate proteins "
        "successfully saved."
    )

else:

    print("\n⚠ Missing protein IDs:")
    print(missing_ids[:20])
    
    
    
    
    # ============================================================
# PHASE 6 — KEGG ANALYSIS
# STEP 2 — INSPECT BLASTKOALA OUTPUT
# ============================================================

kegg_file = r"E:\Genomics\A.baumanii_genome\KEGG\blast_koala.txt"


print("KEGG FILE CHECK")
print("--------------------------------")

print(
    "BlastKOALA file:",
    os.path.exists(kegg_file)
)


if not os.path.exists(kegg_file):

    raise FileNotFoundError(
        "blast_koala.txt not found. Check the file path."
    )


# ============================================================
# SHOW FIRST 30 LINES
# ============================================================

print("\nBLASTKOALA FILE PREVIEW")
print("--------------------------------")

with open(kegg_file, "r", encoding="utf-8") as f:

    for i, line in enumerate(f):

        print(line.rstrip())

        if i >= 29:
            break
        
        
        
# ============================================================
# PHASE 6 — KEGG ANALYSIS
# STEP 3 — LOAD BLASTKOALA KO ASSIGNMENTS
# ============================================================

# Read BlastKOALA result
kegg_results = pd.read_csv(
    kegg_file,
    sep="\t",
    header=None,
    names=["protein_id", "KO"],
    dtype=str
)


# ============================================================
# CLEAN DATA
# ============================================================

kegg_results["protein_id"] = (
    kegg_results["protein_id"]
    .str.strip()
)

kegg_results["KO"] = (
    kegg_results["KO"]
    .fillna("")
    .str.strip()
)


# ============================================================
# SUMMARY
# ============================================================

total_entries = len(kegg_results)

annotated_entries = (
    kegg_results["KO"] != ""
).sum()

unannotated_entries = (
    kegg_results["KO"] == ""
).sum()


print("KEGG KO ANNOTATION SUMMARY")
print("--------------------------------")

print(
    "Total proteins:",
    total_entries
)

print(
    "Proteins with KO:",
    annotated_entries
)

print(
    "Proteins without KO:",
    unannotated_entries
)

print(
    "Annotation rate:",
    round(
        annotated_entries / total_entries * 100,
        2
    ),
    "%"
)


# ============================================================
# CHECK DUPLICATES
# ============================================================

print(
    "\nUnique protein IDs:",
    kegg_results["protein_id"].nunique()
)

print(
    "Unique KO IDs:",
    kegg_results.loc[
        kegg_results["KO"] != "",
        "KO"
    ].nunique()
)


# ============================================================
# SHOW FIRST 20 ANNOTATED PROTEINS
# ============================================================

print("\nFirst 20 KO assignments:")

print(
    kegg_results[
        kegg_results["KO"] != ""
    ].head(20)
)





# ============================================================
# PHASE 6 — KEGG ANALYSIS
# STEP 4 — KEGG PATHWAY MAPPING
# CHECKPOINT + CORRECTED VERSION
# ============================================================

import os
import requests
import time


# ============================================================
# OUTPUT LOCATION
# ============================================================

kegg_folder = r"E:\Genomics\A.baumanii_genome\KEGG"

os.makedirs(kegg_folder, exist_ok=True)

mapping_file = os.path.join(
    kegg_folder,
    "KO_pathway_mapping.csv"
)


# ============================================================
# GET UNIQUE KO IDs
# ============================================================

ko_ids = (
    kegg_results.loc[
        kegg_results["KO"] != "",
        "KO"
    ]
    .drop_duplicates()
    .tolist()
)


print("KEGG PATHWAY MAPPING")
print("--------------------------------")
print("Unique KOs:", len(ko_ids))


# ============================================================
# LOAD EXISTING CHECKPOINT
# ============================================================

if os.path.exists(mapping_file):

    old_mapping = pd.read_csv(
        mapping_file,
        dtype=str
    )

    # --------------------------------------------------------
    # KEEP ONLY VALID mapXXXXX PATHWAYS
    # --------------------------------------------------------

    old_mapping = old_mapping[
        old_mapping["pathway_id"]
        .astype(str)
        .str.startswith("map")
    ].copy()

    old_mapping = old_mapping.drop_duplicates()

    completed_kos = set(
        old_mapping["KO"].dropna()
    )

    print(
        "Previously mapped KOs:",
        len(completed_kos)
    )

else:

    old_mapping = pd.DataFrame(
        columns=["KO", "pathway_id"]
    )

    completed_kos = set()

    print(
        "No previous checkpoint found."
    )


# ============================================================
# FIND REMAINING KOs
# ============================================================

remaining_kos = [
    ko for ko in ko_ids
    if ko not in completed_kos
]


print(
    "KOs remaining:",
    len(remaining_kos)
)

print("--------------------------------")


# ============================================================
# MAP REMAINING KOs
# ============================================================

new_records = []

for i, ko in enumerate(remaining_kos, 1):

    try:

        response = requests.get(
            f"https://rest.kegg.jp/link/pathway/ko:{ko}",
            timeout=20
        )


        # ====================================================
        # SUCCESSFUL KEGG RESPONSE
        # ====================================================

        if response.status_code == 200:

            lines = response.text.strip().splitlines()

            found_pathway = False

            for line in lines:

                if not line:
                    continue

                parts = line.split("\t")

                if len(parts) != 2:
                    continue

                ko_id = parts[0].replace(
                    "ko:",
                    ""
                )

                pathway_id = parts[1].replace(
                    "path:",
                    ""
                )

                # --------------------------------------------
                # KEEP ONLY mapXXXXX
                # --------------------------------------------

                if pathway_id.startswith("map"):

                    new_records.append({
                        "KO": ko_id,
                        "pathway_id": pathway_id
                    })

                    found_pathway = True


            # ------------------------------------------------
            # Save KO even if it has no pathway
            # ------------------------------------------------

            if not found_pathway:

                new_records.append({
                    "KO": ko,
                    "pathway_id": ""
                })


        else:

            print(
                f"KEGG error for {ko}: "
                f"HTTP {response.status_code}"
            )

            continue


    except Exception as e:

        print(
            f"\nConnection error for {ko}: {e}"
        )

        print(
            "Saving current progress..."
        )

        break


    # ========================================================
    # CHECKPOINT EVERY 25 KOs
    # ========================================================

    if i % 25 == 0 or i == len(remaining_kos):

        if new_records:

            new_df = pd.DataFrame(
                new_records
            )

            old_mapping = pd.concat(
                [
                    old_mapping,
                    new_df
                ],
                ignore_index=True
            )

            old_mapping = (
                old_mapping
                .drop_duplicates()
            )

            old_mapping.to_csv(
                mapping_file,
                index=False
            )

            new_records = []


        print(
            f"Processed {i}/{len(remaining_kos)} "
            f"| Total saved: "
            f"{old_mapping['KO'].nunique()}"
        )


    time.sleep(0.1)


# ============================================================
# FINAL SAVE
# ============================================================

if new_records:

    new_df = pd.DataFrame(
        new_records
    )

    old_mapping = pd.concat(
        [
            old_mapping,
            new_df
        ],
        ignore_index=True
    )

    old_mapping = (
        old_mapping
        .drop_duplicates()
    )

    old_mapping.to_csv(
        mapping_file,
        index=False
    )


# ============================================================
# FINAL DATAFRAME
# ============================================================

ko_pathway_df = old_mapping.copy()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n========================================")
print("KEGG PATHWAY MAPPING SUMMARY")
print("========================================")

print(
    "Unique KOs submitted:",
    len(ko_ids)
)

print(
    "KOs processed:",
    ko_pathway_df["KO"].nunique()
)

print(
    "KO-pathway associations:",
    len(
        ko_pathway_df[
            ko_pathway_df["pathway_id"] != ""
        ]
    )
)

print(
    "KOs mapped to pathways:",
    ko_pathway_df.loc[
        ko_pathway_df["pathway_id"] != "",
        "KO"
    ].nunique()
)

print(
    "KEGG pathways detected:",
    ko_pathway_df.loc[
        ko_pathway_df["pathway_id"] != "",
        "pathway_id"
    ].nunique()
)


print("\nSaved to:")
print(mapping_file)


print("\nFirst 20 valid pathway mappings:")

print(
    ko_pathway_df[
        ko_pathway_df["pathway_id"] != ""
    ].head(20)
)

# =============================================================================
# 
# ## Quality control check
# 
# =============================================================================

print("KEGG MAPPING FILE CHECK")
print("--------------------------------")

print("Total rows:", len(ko_pathway_df))

print(
    "Unique KO IDs:",
    ko_pathway_df["KO"].nunique()
)

print(
    "Unique KO-pathway pairs:",
    ko_pathway_df[
        ko_pathway_df["pathway_id"] != ""
    ][
        ["KO", "pathway_id"]
    ].drop_duplicates().shape[0]
)

print(
    "Unique KEGG pathways:",
    ko_pathway_df.loc[
        ko_pathway_df["pathway_id"] != "",
        "pathway_id"
    ].nunique()
)

print("\nDuplicate rows:")
print(
    ko_pathway_df.duplicated(
        ["KO", "pathway_id"]
    ).sum()
)




# ============================================================
# PHASE 6 — KEGG ANALYSIS
# STEP 5 — CONNECT PROTEINS → KO → PATHWAYS
# ============================================================

# ============================================================
# KEEP ONLY VALID KO ANNOTATIONS
# ============================================================

protein_ko = kegg_results[
    kegg_results["KO"] != ""
].copy()


# ============================================================
# KEEP ONLY VALID KEGG PATHWAYS
# ============================================================

ko_pathways = ko_pathway_df[
    ko_pathway_df["pathway_id"]
    .astype(str)
    .str.startswith("map")
].copy()


# ============================================================
# MERGE PROTEINS WITH PATHWAYS
# ============================================================

protein_pathways = protein_ko.merge(
    ko_pathways,
    on="KO",
    how="left"
)


# ============================================================
# SUMMARY
# ============================================================

print("PROTEIN → KO → PATHWAY MAPPING")
print("--------------------------------")

print(
    "KO-annotated proteins:",
    protein_ko["protein_id"].nunique()
)

print(
    "Proteins with KEGG pathways:",
    protein_pathways.loc[
        protein_pathways["pathway_id"].notna(),
        "protein_id"
    ].nunique()
)

print(
    "Unique pathways:",
    protein_pathways.loc[
        protein_pathways["pathway_id"].notna(),
        "pathway_id"
    ].nunique()
)

print(
    "Protein-pathway associations:",
    protein_pathways[
        protein_pathways["pathway_id"].notna()
    ].shape[0]
)


# ============================================================
# SAVE RESULT
# ============================================================

protein_pathway_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\protein_KO_pathway_mapping.csv"
)

protein_pathways.to_csv(
    protein_pathway_file,
    index=False
)


# ============================================================
# SHOW EXAMPLES
# ============================================================

print("\nFirst 20 protein → KO → pathway mappings:")

print(
    protein_pathways[
        protein_pathways["pathway_id"].notna()
    ].head(20)
)


print("\nSaved to:")
print(protein_pathway_file)



# ============================================================
# PHASE 6 — KEGG ANALYSIS
# STEP 6 — FUNCTIONAL / PATHWAY PRIORITIZATION
# ============================================================

# ============================================================
# GET KEGG PATHWAY NAMES
# One KEGG request for the complete pathway list
# ============================================================

response = requests.get(
    "https://rest.kegg.jp/list/pathway",
    timeout=30
)


if response.status_code != 200:

    raise RuntimeError(
        "KEGG pathway list could not be retrieved."
    )


# ============================================================
# CREATE PATHWAY NAME TABLE
# ============================================================

pathway_names = []

for line in response.text.strip().splitlines():

    parts = line.split("\t")

    if len(parts) >= 2:

        pathway_id = parts[0].replace(
            "path:",
            ""
        )

        pathway_name = parts[1].strip()

        # Keep only map pathways
        if pathway_id.startswith("map"):

            pathway_names.append({
                "pathway_id": pathway_id,
                "pathway_name": pathway_name
            })


pathway_names_df = pd.DataFrame(
    pathway_names
)


# ============================================================
# MERGE PATHWAY NAMES WITH PROTEIN-PATHWAY DATA
# ============================================================

protein_pathways_named = protein_pathways.merge(
    pathway_names_df,
    on="pathway_id",
    how="left"
)


# ============================================================
# CREATE PATHWAY SUMMARY
# ============================================================

pathway_summary = (
    protein_pathways_named[
        protein_pathways_named["pathway_id"].notna()
    ]
    .groupby(
        ["pathway_id", "pathway_name"]
    )
    .agg(
        proteins=(
            "protein_id",
            "nunique"
        ),
        KOs=(
            "KO",
            "nunique"
        )
    )
    .reset_index()
    .sort_values(
        "proteins",
        ascending=False
    )
)


# ============================================================
# DISPLAY PATHWAY SUMMARY
# ============================================================

print("KEGG FUNCTIONAL / PATHWAY SUMMARY")
print("--------------------------------")

print(
    "Unique pathways:",
    pathway_summary["pathway_id"].nunique()
)

print(
    "Proteins represented in pathways:",
    protein_pathways_named[
        protein_pathways_named["pathway_id"].notna()
    ]["protein_id"].nunique()
)


print("\nTop KEGG pathways by number of candidate proteins:")

print(
    pathway_summary.head(30)
)


# ============================================================
# DEFINE BIOLOGICALLY RELEVANT TARGET CATEGORIES
# ============================================================

priority_keywords = {

    "Cell wall / envelope":
        [
            "cell wall",
            "peptidoglycan",
            "lipopolysaccharide",
            "lipid A",
            "cell envelope"
        ],

    "Amino acid biosynthesis":
        [
            "biosynthesis of amino acids",
            "amino acid metabolism"
        ],

    "Nucleotide metabolism":
        [
            "nucleotide metabolism",
            "purine metabolism",
            "pyrimidine metabolism"
        ],

    "Energy metabolism":
        [
            "energy metabolism",
            "oxidative phosphorylation",
            "carbon metabolism"
        ],

    "Cofactor / vitamin metabolism":
        [
            "biosynthesis of cofactors",
            "metabolism of cofactors",
            "vitamin"
        ],

    "DNA / replication / repair":
        [
            "replication",
            "dna repair",
            "recombination"
        ],

    "Protein synthesis":
        [
            "ribosome",
            "translation"
        ]
}


# ============================================================
# ASSIGN FUNCTIONAL CATEGORY
# ============================================================

def assign_category(pathway_name):

    if pd.isna(pathway_name):

        return "Unclassified"


    name = pathway_name.lower()


    for category, keywords in priority_keywords.items():

        for keyword in keywords:

            if keyword in name:

                return category


    return "Other"


protein_pathways_named["functional_category"] = (
    protein_pathways_named["pathway_name"]
    .apply(assign_category)
)


# ============================================================
# COUNT PROTEINS BY FUNCTIONAL CATEGORY
# ============================================================

category_summary = (
    protein_pathways_named[
        protein_pathways_named["pathway_id"].notna()
    ]
    .groupby("functional_category")
    .agg(
        proteins=(
            "protein_id",
            "nunique"
        ),
        pathways=(
            "pathway_id",
            "nunique"
        )
    )
    .reset_index()
    .sort_values(
        "proteins",
        ascending=False
    )
)


# ============================================================
# DISPLAY CATEGORY SUMMARY
# ============================================================

print("\nFUNCTIONAL CATEGORY SUMMARY")
print("--------------------------------")

print(
    category_summary
)


# ============================================================
# CREATE CANDIDATE PRIORITIZATION TABLE
# ============================================================

candidate_priority = (
    protein_pathways_named[
        protein_pathways_named["pathway_id"].notna()
    ]
    .groupby("protein_id")
    .agg(
        KO=(
            "KO",
            "first"
        ),

        pathway_count=(
            "pathway_id",
            "nunique"
        ),

        pathway_names=(
            "pathway_name",
            lambda x: "; ".join(
                sorted(
                    set(x.dropna())
                )
            )
        ),

        functional_categories=(
            "functional_category",
            lambda x: "; ".join(
                sorted(
                    set(x)
                )
            )
        )
    )
    .reset_index()
)


# ============================================================
# FUNCTIONAL PRIORITY FLAG
# ============================================================

high_value_categories = [
    "Cell wall / envelope",
    "Amino acid biosynthesis",
    "Nucleotide metabolism",
    "Energy metabolism",
    "Cofactor / vitamin metabolism",
    "DNA / replication / repair",
    "Protein synthesis"
]


candidate_priority["functional_priority"] = (
    candidate_priority[
        "functional_categories"
    ]
    .apply(
        lambda x:
        "High"
        if any(
            category in x
            for category in high_value_categories
        )
        else "Standard"
    )
)


# ============================================================
# SORT CANDIDATES
# ============================================================

candidate_priority = (
    candidate_priority
    .sort_values(
        [
            "functional_priority",
            "pathway_count"
        ],
        ascending=[
            True,
            False
        ]
    )
)


# ============================================================
# SAVE RESULT
# ============================================================

priority_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\functional_candidate_prioritization.csv"
)

candidate_priority.to_csv(
    priority_file,
    index=False
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n========================================")
print("FUNCTIONAL PRIORITIZATION SUMMARY")
print("========================================")

print(
    "Proteins with pathway information:",
    len(candidate_priority)
)

print(
    "High functional-priority proteins:",
    (
        candidate_priority[
            "functional_priority"
        ] == "High"
    ).sum()
)

print("\nTop candidate proteins:")

print(
    candidate_priority.head(20)
)


print("\nSaved to:")
print(priority_file)




# ============================================================
# PHASE 7 — FINAL TARGET PRIORITIZATION
# STEP 1 — BIOLOGICAL TARGET FILTER
# ============================================================

import pandas as pd


# ============================================================
# LOAD PHASE 6 RESULT
# ============================================================

priority_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\functional_candidate_prioritization.csv"
)

candidate_priority = pd.read_csv(priority_file)


# ============================================================
# STEP 1 — KEEP HIGH-PRIORITY PROTEINS
# ============================================================

candidates = candidate_priority[
    candidate_priority["functional_priority"] == "High"
].copy()


print("PHASE 7 — STEP 1: BIOLOGICAL TARGET FILTER")
print("===========================================")

print("High-priority candidates:", len(candidates))
print("Unique KOs:", candidates["KO"].nunique())


# ============================================================
# STEP 2 — KEEP IMPORTANT ANTIBACTERIAL FUNCTIONS
# ============================================================

strong_target_categories = [
    "Cell wall / envelope",
    "DNA / replication / repair",
    "Protein synthesis",
    "Nucleotide metabolism",
    "Amino acid biosynthesis",
    "Cofactor / vitamin metabolism"
]


target_candidates = candidates[
    candidates["functional_categories"]
    .fillna("")
    .apply(
        lambda x: any(
            category in x
            for category in strong_target_categories
        )
    )
].copy()


# ============================================================
# REMOVE DUPLICATE PROTEINS
# ============================================================

target_candidates = (
    target_candidates
    .drop_duplicates(subset=["protein_id"])
    .reset_index(drop=True)
)


# ============================================================
# SUMMARY
# ============================================================

print("\nBIOLOGICAL TARGET FILTER RESULT")
print("--------------------------------")

print(
    "Potential target candidates:",
    len(target_candidates)
)

print(
    "Unique KOs:",
    target_candidates["KO"].nunique()
)


# ============================================================
# FUNCTIONAL DISTRIBUTION
# ============================================================

print("\nFUNCTIONAL DISTRIBUTION")
print("-----------------------")

for category in strong_target_categories:

    count = target_candidates[
        target_candidates["functional_categories"]
        .fillna("")
        .str.contains(category, regex=False)
    ]["protein_id"].nunique()

    print(f"{category}: {count}")


# ============================================================
# DISPLAY CANDIDATES
# ============================================================

print("\nTARGET CANDIDATES")
print("-----------------")

print(
    target_candidates[
        [
            "protein_id",
            "KO",
            "functional_categories",
            "pathway_names"
        ]
    ].to_string(index=False)
)


# ============================================================
# SAVE
# ============================================================

target_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\preliminary_target_candidates.csv"
)

target_candidates.to_csv(
    target_file,
    index=False
)

print("\nSaved to:")
print(target_file)


# ============================================================
# PHASE 7 — FINAL TARGET PRIORITIZATION
# STEP 2 — KO-LEVEL TARGET REDUNDANCY CHECK
# ============================================================
# LOAD STEP 1 RESULT
# ============================================================

target_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\preliminary_target_candidates.csv"
)

targets = pd.read_csv(target_file)


print("PHASE 7 — STEP 2: TARGET REDUNDANCY CHECK")
print("==========================================")


# ============================================================
# BASIC SUMMARY
# ============================================================

print("Total protein candidates:", len(targets))
print("Unique KOs:", targets["KO"].nunique())


# ============================================================
# CHECK PROTEINS PER KO
# ============================================================

ko_summary = (
    targets
    .groupby("KO")
    .agg(
        protein_count=("protein_id", "nunique"),
        proteins=(
            "protein_id",
            lambda x: "; ".join(sorted(x.astype(str)))
        ),
        functional_categories=(
            "functional_categories",
            lambda x: "; ".join(sorted(set(x.dropna())))
        ),
        pathway_names=(
            "pathway_names",
            lambda x: "; ".join(sorted(set(x.dropna())))
        )
    )
    .reset_index()
)


# ============================================================
# DISPLAY REDUNDANCY
# ============================================================

print("\nKO-LEVEL SUMMARY")
print("----------------")

print(
    ko_summary[
        [
            "KO",
            "protein_count",
            "proteins",
            "functional_categories"
        ]
    ].to_string(index=False)
)


# ============================================================
# IDENTIFY REDUNDANT KOs
# ============================================================

redundant_kos = ko_summary[
    ko_summary["protein_count"] > 1
]


print("\nREDUNDANT KOs")
print("-------------")

if len(redundant_kos) == 0:

    print("No redundant KOs found.")
    print("Each candidate protein represents a unique KO.")

else:

    print(
        redundant_kos[
            [
                "KO",
                "protein_count",
                "proteins"
            ]
        ].to_string(index=False)
    )


# ============================================================
# SELECT ONE REPRESENTATIVE PROTEIN PER KO
# ============================================================

representative_targets = (
    targets
    .sort_values(
        ["KO", "protein_id"]
    )
    .drop_duplicates(
        subset=["KO"],
        keep="first"
    )
    .reset_index(drop=True)
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\nREPRESENTATIVE TARGET SET")
print("-------------------------")

print(
    "Protein candidates before:",
    len(targets)
)

print(
    "Unique target KOs:",
    len(representative_targets)
)

print(
    "Proteins removed as KO-level duplicates:",
    len(targets) - len(representative_targets)
)


# ============================================================
# SAVE
# ============================================================

representative_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\representative_target_candidates.csv"
)

representative_targets.to_csv(
    representative_file,
    index=False
)


print("\nSaved to:")
print(representative_file)



# ============================================================
# PHASE 7 — FINAL TARGET PRIORITIZATION
# STEP 3 — CORE ANTIBACTERIAL FUNCTION SCORE
# ============================================================

import pandas as pd


# ============================================================
# LOAD STEP 2 RESULT
# ============================================================

target_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\representative_target_candidates.csv"
)

targets = pd.read_csv(target_file)


print("PHASE 7 — STEP 3: CORE ANTIBACTERIAL FUNCTION SCORE")
print("===================================================")

print("Total unique targets:", len(targets))


# ============================================================
# CORE ANTIBACTERIAL FUNCTIONS
# EACH FUNCTION = 3 POINTS
# ============================================================

core_function_scores = {

    "Cell wall / envelope": 3,

    "DNA / replication / repair": 3,

    "Protein synthesis": 3
}


# ============================================================
# CALCULATE SCORE
# ============================================================

def calculate_antibacterial_score(categories):

    if pd.isna(categories):
        return 0

    score = 0

    for function, points in core_function_scores.items():

        if function in categories:
            score += points

    return score


targets["antibacterial_score"] = (
    targets["functional_categories"]
    .apply(calculate_antibacterial_score)
)


# ============================================================
# EXCLUDE PROTEINS WITH NONE OF THE THREE FUNCTIONS
# ============================================================

antibacterial_targets = targets[
    targets["antibacterial_score"] > 0
].copy()


# ============================================================
# RANK TARGETS
# ============================================================

antibacterial_targets = (
    antibacterial_targets
    .sort_values(
        "antibacterial_score",
        ascending=False
    )
    .reset_index(drop=True)
)


# ============================================================
# ASSIGN RANK
# ============================================================

antibacterial_targets["target_rank"] = (
    antibacterial_targets.index + 1
)


# ============================================================
# SUMMARY
# ============================================================

print("\nANTIBACTERIAL TARGET SUMMARY")
print("----------------------------")

print(
    "Targets before filtering:",
    len(targets)
)

print(
    "Antibacterial targets retained:",
    len(antibacterial_targets)
)


# ============================================================
# SCORE DISTRIBUTION
# ============================================================

print("\nSCORE DISTRIBUTION")
print("------------------")

print(
    antibacterial_targets[
        "antibacterial_score"
    ]
    .value_counts()
    .sort_index(ascending=False)
)


# ============================================================
# TOP TARGETS
# ============================================================

print("\nTOP ANTIBACTERIAL TARGETS")
print("-------------------------")

print(
    antibacterial_targets[
        [
            "target_rank",
            "protein_id",
            "KO",
            "functional_categories",
            "antibacterial_score"
        ]
    ]
    .head(30)
    .to_string(index=False)
)


# ============================================================
# SAVE
# ============================================================

output_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\antibacterial_target_ranking.csv"
)

antibacterial_targets.to_csv(
    output_file,
    index=False
)


print("\nSaved to:")
print(output_file)



# ============================================================
# PHASE 7 — FINAL TARGET PRIORITIZATION
# STEP 4 — MECHANISM-BASED TARGET SHORTLIST
# ============================================================
# LOAD STEP 3 RESULT
# ============================================================

input_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\antibacterial_target_ranking.csv"
)

targets = pd.read_csv(input_file)


print("PHASE 7 — STEP 4: MECHANISM-BASED TARGET SHORTLIST")
print("===================================================")

print("Antibacterial candidates:", len(targets))


# ============================================================
# DEFINE THREE MAIN ANTIBACTERIAL MECHANISMS
# ============================================================

mechanism_keywords = {

    "DNA replication / repair": [
        "DNA / replication / repair",
        "DNA replication",
        "DNA repair",
        "Homologous recombination",
        "Mismatch repair",
        "Base excision repair",
        "Nucleotide excision repair"
    ],

    "Cell wall / envelope": [
        "Cell wall / envelope",
        "Peptidoglycan biosynthesis",
        "LPS biosynthesis",
        "Lipopolysaccharide",
        "Beta-lactam resistance",
        "Vancomycin resistance"
    ],

    "Protein synthesis / ribosome": [
        "Protein synthesis",
        "Ribosome"
    ]
}


# ============================================================
# ASSIGN MECHANISM
# ============================================================

def assign_mechanism(row):

    text = (
        str(row["functional_categories"])
        + " "
        + str(row["pathway_names"])
    )

    matched = []

    for mechanism, keywords in mechanism_keywords.items():

        for keyword in keywords:

            if keyword.lower() in text.lower():

                matched.append(mechanism)
                break

    if len(matched) == 0:
        return "Other"

    return "; ".join(matched)


targets["mechanism"] = targets.apply(
    assign_mechanism,
    axis=1
)


# ============================================================
# REMOVE TARGETS THAT CANNOT BE ASSIGNED
# ============================================================

targets = targets[
    targets["mechanism"] != "Other"
].copy()


# ============================================================
# COUNT TARGETS BY MECHANISM
# ============================================================

print("\nMECHANISM DISTRIBUTION")
print("----------------------")

print(
    targets["mechanism"]
    .value_counts()
)


# ============================================================
# CREATE A PRIORITY SCORE WITHIN EACH MECHANISM
# ============================================================
#
# Existing antibacterial score:
# 3 = one core antibacterial function
#
# We use pathway information only as a tie-breaker.
# No new biological score is created.
# ============================================================

targets["pathway_count"] = (
    targets["pathway_names"]
    .fillna("")
    .apply(
        lambda x: len(
            set(
                p.strip()
                for p in x.split(";")
                if p.strip()
            )
        )
    )
)


# ============================================================
# RANK WITHIN EACH MECHANISM
# ============================================================

targets = (
    targets
    .sort_values(
        [
            "mechanism",
            "pathway_count"
        ],
        ascending=[
            True,
            False
        ]
    )
    .reset_index(drop=True)
)


# ============================================================
# SELECT REPRESENTATIVE TARGETS
# ============================================================
#
# Maximum number from each mechanism:
# DNA / repair       = 5
# Cell wall/envelope = 5
# Protein synthesis  = 5
#
# Total maximum = 15
# ============================================================

shortlist_per_mechanism = 5

shortlist = (
    targets
    .groupby(
        "mechanism",
        group_keys=False
    )
    .head(shortlist_per_mechanism)
    .copy()
    .reset_index(drop=True)
)


# ============================================================
# FINAL RANK
# ============================================================

shortlist["shortlist_rank"] = (
    shortlist.index + 1
)


# ============================================================
# SUMMARY
# ============================================================

print("\nFINAL SHORTLIST")
print("---------------")

print(
    "Total shortlisted targets:",
    len(shortlist)
)

print(
    "\nTargets per mechanism:"
)

print(
    shortlist["mechanism"]
    .value_counts()
)


# ============================================================
# DISPLAY SHORTLIST
# ============================================================

print("\nSHORTLISTED TARGETS")
print("-------------------")

print(
    shortlist[
        [
            "shortlist_rank",
            "protein_id",
            "KO",
            "mechanism",
            "functional_categories",
            "pathway_names",
            "antibacterial_score",
            "pathway_count"
        ]
    ]
    .to_string(index=False)
)


# ============================================================
# SAVE FINAL SHORTLIST
# ============================================================

output_file = (
    r"E:\Genomics\A.baumanii_genome\KEGG"
    r"\final_antibacterial_target_shortlist.csv"
)

shortlist.to_csv(
    output_file,
    index=False
)


print("\nSaved to:")
print(output_file)


 
 

