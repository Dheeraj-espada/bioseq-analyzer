#!/usr/bin/env python3
"""
Download real sequences from NCBI GenBank for testing
Uses Biopython's Entrez module to fetch sequences
"""

from Bio import Entrez, SeqIO
import time
import os

# IMPORTANT: Always tell NCBI who you are
Entrez.email = "dheerajbabu133@gmail.com"  # REPLACE WITH YOUR EMAIL

def download_gene_sequences(gene_ids, output_file="real_data/ncbi_sequences.fasta"):
    """
    Download sequences from NCBI GenBank by accession numbers
    
    Args:
        gene_ids (list): List of GenBank accession numbers
        output_file (str): Output FASTA file path
    """
    print(f"📥 Downloading {len(gene_ids)} sequences from NCBI GenBank...")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    sequences = []
    
    for i, gene_id in enumerate(gene_ids, 1):
        try:
            print(f"   [{i}/{len(gene_ids)}] Fetching {gene_id}...", end=" ")
            
            # Fetch sequence from GenBank
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            sequences.append(record)
            print(f"✅ ({len(record)} bp)")
            
            # Be nice to NCBI servers - wait between requests
            time.sleep(0.5)
            
        except Exception as e:
            print(f"❌ Error: {e}")
            continue
    
    # Write all sequences to file
    if sequences:
        SeqIO.write(sequences, output_file, "fasta")
        print(f"\n✅ Saved {len(sequences)} sequences to {output_file}")
        print(f"   Total base pairs: {sum(len(s.seq) for s in sequences):,}")
        return output_file
    else:
        print("\n❌ No sequences downloaded")
        return None


def download_human_genes():
    """Download a collection of important human genes"""
    print("="*60)
    print("Downloading Human Disease-Related Genes")
    print("="*60)
    
    # Famous human genes with their GenBank accessions
    human_genes = {
        'NM_007294.3': 'BRCA1 (Breast cancer susceptibility)',
        'NM_000546.6': 'TP53 (Tumor protein p53)',
        'NM_005228.5': 'EGFR (Epidermal growth factor receptor)',
        'NM_002467.6': 'MYC (MYC proto-oncogene)',
        'NM_033360.4': 'KRAS (Kirsten rat sarcoma)',
        'NM_004333.6': 'BRAF (B-Raf proto-oncogene)',
        'NM_000059.4': 'BRCA2 (Breast cancer susceptibility 2)',
        'NM_004985.5': 'KRAS (Alternative transcript)',
    }
    
    print("\nGenes to download:")
    for acc, name in human_genes.items():
        print(f"   • {acc}: {name}")
    
    print()
    gene_ids = list(human_genes.keys())
    output_file = download_gene_sequences(gene_ids, "real_data/human_cancer_genes.fasta")
    
    return output_file


def download_bacterial_genes():
    """Download E. coli genes for comparison"""
    print("\n" + "="*60)
    print("Downloading E. coli Genes (for codon usage comparison)")
    print("="*60)
    
    ecoli_genes = {
        'NP_414542.1': 'lacZ (Beta-galactosidase)',
        'NP_414543.1': 'lacY (Lactose permease)',
        'NP_414544.1': 'lacA (Thiogalactoside acetyltransferase)',
        'NP_418820.1': 'rpoB (RNA polymerase beta subunit)',
        'NP_417804.1': 'groEL (Chaperonin)',
    }
    
    print("\nGenes to download:")
    for acc, name in ecoli_genes.items():
        print(f"   • {acc}: {name}")
    
    print()
    gene_ids = list(ecoli_genes.keys())
    output_file = download_gene_sequences(gene_ids, "real_data/ecoli_genes.fasta")
    
    return output_file


def download_viral_genomes():
    """Download small viral genomes"""
    print("\n" + "="*60)
    print("Downloading Viral Sequences")
    print("="*60)
    
    viral_seqs = {
        'NC_045512.2': 'SARS-CoV-2 complete genome',
        'NC_001802.1': 'HIV-1 complete genome',
        'NC_001489.1': 'Hepatitis B virus',
    }
    
    print("\nSequences to download:")
    for acc, name in viral_seqs.items():
        print(f"   • {acc}: {name}")
    
    print()
    gene_ids = list(viral_seqs.keys())
    output_file = download_gene_sequences(gene_ids, "real_data/viral_genomes.fasta")
    
    return output_file


def create_sample_script():
    """Create a script to analyze all downloaded data"""
    script_content = """#!/bin/bash
# Analyze all downloaded NCBI sequences

echo "Analyzing Human Cancer Genes..."
python sequence_analyzer.py real_data/human_cancer_genes.fasta \\
    --codon-usage \\
    --find-motif ATG \\
    --output-dir results/human_genes

echo ""
echo "Analyzing E. coli Genes..."
python sequence_analyzer.py real_data/ecoli_genes.fasta \\
    --codon-usage \\
    --find-motif ATG \\
    --output-dir results/ecoli_genes

echo ""
echo "Analyzing Viral Genomes..."
python sequence_analyzer.py real_data/viral_genomes.fasta \\
    --codon-usage \\
    --min-orf 200 \\
    --output-dir results/viral_genomes

echo ""
echo "✅ All analyses complete!"
echo "Check the results/ directory for outputs"
"""
    
    with open('analyze_all_ncbi.sh', 'w') as f:
        f.write(script_content)
    
    os.chmod('analyze_all_ncbi.sh', 0o755)
    print("\n✅ Created analyze_all_ncbi.sh script")


def main():
    """Main function"""
    print("\n🧬 NCBI Sequence Downloader for BioSeq Analyzer")
    print("="*60)
    print("\n⚠️  IMPORTANT: Update Entrez.email in this script first!")
    print("   Set it to your email address (line 13)")
    print()
    
    # Download different datasets
    human_file = download_human_genes()
    ecoli_file = download_bacterial_genes()
    viral_file = download_viral_genomes()
    
    # Create analysis script
    create_sample_script()
    
    # Summary
    print("\n" + "="*60)
    print("📊 Download Summary")
    print("="*60)
    
    if human_file:
        print(f"✅ Human genes: {human_file}")
    if ecoli_file:
        print(f"✅ E. coli genes: {ecoli_file}")
    if viral_file:
        print(f"✅ Viral genomes: {viral_file}")
    
    print("\n📝 Next Steps:")
    print("   1. Run: ./analyze_all_ncbi.sh")
    print("   2. Or analyze individually:")
    print("      python sequence_analyzer.py real_data/human_cancer_genes.fasta --codon-usage")
    print()


if __name__ == "__main__":
    main()
