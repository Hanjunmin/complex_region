#!/usr/bin/env bash

set -euo pipefail

# phylogenetic analysis of block3
mafft --auto --quiet --thread 32 block3.sequences.fasta > block3.sequences.MSA.fasta

trimal -in block3.sequences.MSA.fasta -out block3.sequences.MSA.trimmed.fasta -gt 0.9

beast -threads 16 block3.sequences.MSA.trimmed.xml


# phylogenetic analysis of block5
mafft --auto --quiet --thread 32 block5.sequences.fasta > block5.sequences.MSA.fasta

trimal -in block5.sequences.MSA.fasta -out block5.sequences.MSA.trimmed.fasta -gt 0.9

iqtree2 -s block5.sequences.MSA.trimmed.fasta -B 1000 -redo -m MFP -T 32 --prefix block5.sequences.MSA.trimmed.phylogeny
