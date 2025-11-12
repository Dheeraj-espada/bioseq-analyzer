#!/bin/bash
# Analyze all downloaded NCBI sequences

echo "Analyzing Human Cancer Genes..."
python sequence_analyzer.py real_data/human_cancer_genes.fasta \
    --codon-usage \
    --find-motif ATG \
    --output-dir results/human_genes

echo ""
echo "Analyzing E. coli Genes..."
python sequence_analyzer.py real_data/ecoli_genes.fasta \
    --codon-usage \
    --find-motif ATG \
    --output-dir results/ecoli_genes

echo ""
echo "Analyzing Viral Genomes..."
python sequence_analyzer.py real_data/viral_genomes.fasta \
    --codon-usage \
    --min-orf 200 \
    --output-dir results/viral_genomes

echo ""
echo "✅ All analyses complete!"
echo "Check the results/ directory for outputs"
