#!/usr/bin/env python3
"""
Comprehensive Sequence Analysis Tool
Author: Your Name
Date: 2024
Description: Analyzes FASTA sequences with statistics, ORF finding, translation, and visualization
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set style for visualizations
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)


class SequenceAnalyzer:
    """Main class for sequence analysis"""
    
    def __init__(self, fasta_file):
        """
        Initialize the analyzer with a FASTA file
        
        Args:
            fasta_file (str): Path to FASTA file
        """
        self.fasta_file = fasta_file
        self.sequences = {}
        self.results = []
        self.all_orfs = []
        self.codon_usage = None
    
    def read_fasta(self):
        """Read and parse FASTA file"""
        print(f"\U0001F4D6 Reading FASTA file: {self.fasta_file}")
        try:
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                self.sequences[record.id] = record
            print(f"\u2705 Successfully read {len(self.sequences)} sequences")
            return True
        except FileNotFoundError:
            print(f"\u274C Error: File {self.fasta_file} not found")
            return False
        except Exception as e:
            print(f"\u274C Error reading file: {e}")
            return False

    def calculate_statistics(self):
        """Calculate comprehensive statistics for each sequence"""
        print("\n\U0001F4CA Calculating sequence statistics...")
        
        for seq_id, record in self.sequences.items():
            seq = str(record.seq)
            
            # Basic statistics
            length = len(seq)
            gc_content = gc_fraction(record.seq) * 100
            
            # Nucleotide composition
            composition = Counter(seq.upper())
            a_count = composition.get('A', 0)
            t_count = composition.get('T', 0)
            g_count = composition.get('G', 0)
            c_count = composition.get('C', 0)
            n_count = composition.get('N', 0)
            
            # AT/GC ratio
            at_content = ((a_count + t_count) / length * 100) if length > 0 else 0
            
            # Molecular weight
            mol_weight = molecular_weight(record.seq, seq_type='DNA')
            
            # Store results
            result = {
                'Sequence_ID': seq_id,
                'Length_bp': length,
                'GC_Content_%': round(gc_content, 2),
                'AT_Content_%': round(at_content, 2),
                'A_Count': a_count,
                'T_Count': t_count,
                'G_Count': g_count,
                'C_Count': c_count,
                'N_Count': n_count,
                'Molecular_Weight_Da': round(mol_weight, 2)
            }
            
            self.results.append(result)
        
        print(f"\u2705 Statistics calculated for {len(self.results)} sequences")
        return self.results

    def find_orfs(self, min_length=100):
        """
        Find all Open Reading Frames (ORFs) in sequences
        
        Args:
            min_length (int): Minimum ORF length in nucleotides
        """
        print(f"\n\U0001F50D Finding ORFs (minimum length: {min_length} bp)...")
        
        start_codon = 'ATG'
        stop_codons = {'TAA', 'TAG', 'TGA'}
        
        for seq_id, record in self.sequences.items():
            seq = str(record.seq).upper()
            
            # Check all three reading frames on both strands
            for strand, nuc in [(+1, seq), (-1, str(record.seq.reverse_complement()))]:
                for frame in range(3):
                    trans = str(Seq(nuc[frame:]).translate(to_stop=False))
                    trans_len = len(trans)
                    aa_start = 0
                    
                    while aa_start < trans_len:
                        # Find start codon (M in protein)
                        aa_start = trans.find('M', aa_start)
                        if aa_start == -1:
                            break
                        
                        # Find next stop codon
                        aa_end = trans_len
                        for stop in ['*']:
                            stop_pos = trans.find(stop, aa_start)
                            if stop_pos != -1 and stop_pos < aa_end:
                                aa_end = stop_pos
                        
                        # Calculate nucleotide positions
                        nt_start = frame + aa_start * 3
                        nt_end = frame + aa_end * 3 + 3
                        orf_length = nt_end - nt_start
                        
                        # Check if ORF meets minimum length
                        if orf_length >= min_length:
                            orf_seq = nuc[nt_start:nt_end]
                            protein_seq = str(Seq(orf_seq).translate(to_stop=True))
                            
                            orf_data = {
                                'Sequence_ID': seq_id,
                                'Strand': '+' if strand == 1 else '-',
                                'Frame': frame + 1,
                                'Start_Position': nt_start + 1,
                                'End_Position': nt_end,
                                'Length_bp': orf_length,
                                'Length_aa': len(protein_seq),
                                'ORF_Sequence': orf_seq,
                                'Protein_Sequence': protein_seq
                            }
                            self.all_orfs.append(orf_data)
                        
                        aa_start += 1
        
        print(f"\u2705 Found {len(self.all_orfs)} ORFs across all sequences")
        return self.all_orfs

    def codon_usage_analysis(self, output_dir=None, visualize=False):
        """Compute codon usage across all input sequences and optionally generate a heatmap/barplot."""
        print("\n\U0001F9EC Calculating codon usage...")
        codon_counts = defaultdict(int)
        total_codons = 0
        all_codons = [a+b+c for a in 'ATGC' for b in 'ATGC' for c in 'ATGC']
        
        for record in self.sequences.values():
            seq = str(record.seq).upper()
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if len(codon) == 3 and set(codon) <= set('ATGC'):
                    codon_counts[codon] += 1
                    total_codons += 1
        
        frequencies = {codon: codon_counts.get(codon, 0) / total_codons if total_codons else 0 for codon in all_codons}
        self.codon_usage = frequencies

        df_codon = pd.DataFrame({
            'Codon': list(frequencies.keys()),
            'Frequency': list(frequencies.values())
        })
        df_codon = df_codon.sort_values('Codon')

        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            df_codon.to_csv(f"{output_dir}/codon_usage.csv", index=False)
            print(f"\u2705 Codon usage exported to {output_dir}/codon_usage.csv")

        if visualize:
            plt.figure(figsize=(16, 6))
            sns.barplot(data=df_codon, x='Codon', y='Frequency', color='royalblue')
            plt.title('Codon Usage Frequency', fontsize=14, fontweight='bold')
            plt.ylabel('Relative Frequency')
            plt.xlabel('Codon')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/codon_usage_barplot.png", dpi=300, bbox_inches='tight')
            plt.close()
            print(f"\u2705 Codon usage barplot saved to {output_dir}/codon_usage_barplot.png")

    def translate_sequences(self):
        """Translate sequences and analyze proteins"""
        print("\n\U0001F9EC Translating sequences...")
        
        translation_results = []
        
        for seq_id, record in self.sequences.items():
            # Translate in all 6 frames
            for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                    length = 3 * ((len(nuc) - frame) // 3)
                    if length < 3:
                        continue
                    
                    trans = str(nuc[frame:frame+length].translate(to_stop=False))
                    
                    # Analyze protein properties if translation is valid
                    if len(trans) > 0 and '*' not in trans:
                        try:
                            protein_analysis = ProteinAnalysis(trans)
                            
                            translation_results.append({
                                'Sequence_ID': seq_id,
                                'Strand': '+' if strand == 1 else '-',
                                'Frame': frame + 1,
                                'Protein_Length_aa': len(trans),
                                'Molecular_Weight_Da': round(protein_analysis.molecular_weight(), 2),
                                'Isoelectric_Point': round(protein_analysis.isoelectric_point(), 2),
                                'Aromaticity': round(protein_analysis.aromaticity(), 4),
                                'Instability_Index': round(protein_analysis.instability_index(), 2),
                                'Protein_Sequence': trans[:50] + '...' if len(trans) > 50 else trans
                            })
                        except Exception as e:
                            # Skip invalid protein sequences
                            continue
        
        print(f"\u2705 Generated {len(translation_results)} translations")
        return translation_results
<<<<<<< HEAD

    def generate_visualizations(self, output_dir='output/visualizations'):
=======
    
    def visualize_codon_usage(self, codon_data=None, output_dir='output/visualizations'):
        """
        Generate visualizations for codon usage analysis
        
        Args:
            codon_data (list): Codon usage data from analyze_codon_usage()
            output_dir (str): Output directory for plots
        """
        if codon_data is None:
            codon_data = self.analyze_codon_usage()
        
        if not codon_data:
            print("⚠️  No codon data to visualize")
            return
        
        df_codon = pd.DataFrame(codon_data)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        print("📊 Generating codon usage visualizations...")
        
        # 1. Top 20 Most Frequent Codons Bar Chart
        plt.figure(figsize=(14, 8))
        top20 = df_codon.head(20)
        colors = plt.cm.viridis(range(len(top20)))
        
        bars = plt.bar(range(len(top20)), top20['Frequency_%'], color=colors, alpha=0.8, edgecolor='black')
        plt.xticks(range(len(top20)), 
                   [f"{row['Codon']}\n({row['Amino_Acid']})" for _, row in top20.iterrows()],
                   rotation=0, fontsize=10)
        plt.ylabel('Frequency (%)', fontsize=12, fontweight='bold')
        plt.xlabel('Codon (Amino Acid)', fontsize=12, fontweight='bold')
        plt.title('Top 20 Most Frequent Codons', fontsize=14, fontweight='bold', pad=20)
        plt.grid(axis='y', alpha=0.3, linestyle='--')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}%',
                    ha='center', va='bottom', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/codon_frequency_top20.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Codon Usage by Amino Acid (Grouped)
        plt.figure(figsize=(16, 8))
        aa_grouped = df_codon.groupby('Amino_Acid').agg({
            'Count': 'sum',
            'Frequency_%': 'sum'
        }).sort_values('Count', ascending=False)
        
        # Amino acid full names
        aa_names = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
            '*': 'Stop'
        }
        
        colors_aa = plt.cm.Set3(range(len(aa_grouped)))
        bars = plt.bar(range(len(aa_grouped)), aa_grouped['Frequency_%'], 
                      color=colors_aa, alpha=0.8, edgecolor='black')
        
        labels = [f"{aa}\n{aa_names.get(aa, aa)}" for aa in aa_grouped.index]
        plt.xticks(range(len(aa_grouped)), labels, rotation=0, fontsize=10)
        plt.ylabel('Total Frequency (%)', fontsize=12, fontweight='bold')
        plt.xlabel('Amino Acid', fontsize=12, fontweight='bold')
        plt.title('Codon Usage Grouped by Amino Acid', fontsize=14, fontweight='bold', pad=20)
        plt.grid(axis='y', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/codon_usage_by_amino_acid.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Codon Usage Heatmap (by amino acid)
        # Organize codons by amino acid
        aa_order = sorted(df_codon['Amino_Acid'].unique())
        
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Create matrix for heatmap
        codon_matrix = []
        codon_labels = []
        
        for aa in aa_order:
            aa_codons = df_codon[df_codon['Amino_Acid'] == aa].sort_values('Frequency_%', ascending=False)
            for _, row in aa_codons.iterrows():
                codon_matrix.append(row['Frequency_%'])
                codon_labels.append(f"{row['Codon']} ({aa})")
        
        # Create heatmap-style bar chart
        y_pos = range(len(codon_matrix))
        colors_heat = plt.cm.YlOrRd([x/max(codon_matrix) for x in codon_matrix])
        
        plt.barh(y_pos, codon_matrix, color=colors_heat, alpha=0.8, edgecolor='black', linewidth=0.5)
        plt.yticks(y_pos, codon_labels, fontsize=7)
        plt.xlabel('Frequency (%)', fontsize=12, fontweight='bold')
        plt.title('Complete Codon Usage Profile', fontsize=14, fontweight='bold', pad=20)
        plt.grid(axis='x', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/codon_usage_complete.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 4. Pie Chart - Top 10 Codons
        plt.figure(figsize=(10, 10))
        top10 = df_codon.head(10)
        others_freq = df_codon.iloc[10:]['Frequency_%'].sum()
        
        sizes = list(top10['Frequency_%']) + [others_freq]
        labels = [f"{row['Codon']} ({row['Amino_Acid']})" for _, row in top10.iterrows()] + ['Others']
        colors_pie = plt.cm.Set3(range(len(sizes)))
        
        plt.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
               startangle=90, textprops={'fontsize': 10})
        plt.title('Top 10 Codon Distribution', fontsize=14, fontweight='bold', pad=20)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/codon_usage_pie.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Codon usage visualizations saved to {output_dir}/")
    
    def generate_visualizations(self, output_dir='output/visualizations', include_codon=False):
>>>>>>> ad48967 (Save local changes before pulling new updates)
        """Generate comprehensive visualizations"""
        print("\n\U0001F4CA Generating visualizations...")
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        if not self.results:
            print("\u26A0\ufe0f  No results to visualize. Run calculate_statistics() first.")
            return
        
        df = pd.DataFrame(self.results)
        
        # 1. GC Content Distribution
        plt.figure(figsize=(10, 6))
        plt.hist(df['GC_Content_%'], bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        plt.xlabel('GC Content (%)', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title('Distribution of GC Content', fontsize=14, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/gc_content_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Sequence Length Distribution
        plt.figure(figsize=(10, 6))
        plt.hist(df['Length_bp'], bins=20, color='coral', edgecolor='black', alpha=0.7)
        plt.xlabel('Sequence Length (bp)', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title('Distribution of Sequence Lengths', fontsize=14, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/length_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Nucleotide Composition
        nucleotides = ['A_Count', 'T_Count', 'G_Count', 'C_Count']
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
        
        fig, ax = plt.subplots(figsize=(10, 6))
        x = range(len(df))
        bottom = [0] * len(df)
        
        for nuc, color in zip(nucleotides, colors):
            ax.bar(x, df[nuc], bottom=bottom, label=nuc.replace('_Count', ''), 
                   color=color, alpha=0.8)
            bottom = [sum(x) for x in zip(bottom, df[nuc])]
        
        ax.set_xlabel('Sequence Index', fontsize=12)
        ax.set_ylabel('Nucleotide Count', fontsize=12)
        ax.set_title('Nucleotide Composition per Sequence', fontsize=14, fontweight='bold')
        ax.legend(title='Nucleotide', fontsize=10)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/nucleotide_composition.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 4. GC vs AT Content Scatter
        plt.figure(figsize=(10, 6))
        plt.scatter(df['GC_Content_%'], df['AT_Content_%'], 
                   s=100, alpha=0.6, c=df['Length_bp'], cmap='viridis')
        plt.colorbar(label='Sequence Length (bp)')
        plt.xlabel('GC Content (%)', fontsize=12)
        plt.ylabel('AT Content (%)', fontsize=12)
        plt.title('GC vs AT Content', fontsize=14, fontweight='bold')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/gc_vs_at_content.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 5. ORF Analysis (if ORFs were found)
        if self.all_orfs:
            orf_df = pd.DataFrame(self.all_orfs)
            
            plt.figure(figsize=(12, 6))
            
            # ORF length distribution
            plt.subplot(1, 2, 1)
            plt.hist(orf_df['Length_bp'], bins=20, color='purple', edgecolor='black', alpha=0.7)
            plt.xlabel('ORF Length (bp)', fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            plt.title('ORF Length Distribution', fontsize=12, fontweight='bold')
            plt.grid(axis='y', alpha=0.3)
            
            # ORFs per strand
            plt.subplot(1, 2, 2)
            strand_counts = orf_df['Strand'].value_counts()
            plt.bar(strand_counts.index, strand_counts.values, 
                   color=['#FF6B6B', '#4ECDC4'], alpha=0.7)
            plt.xlabel('Strand', fontsize=12)
            plt.ylabel('Number of ORFs', fontsize=12)
            plt.title('ORFs per Strand', fontsize=12, fontweight='bold')
            plt.grid(axis='y', alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f'{output_dir}/orf_analysis.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        print(f"\u2705 Visualizations saved to {output_dir}/")
    
    def export_to_csv(self, output_dir='output/results', include_codon=True, include_motifs=True):
        """Export all results to CSV files"""
        print("\n\U0001F4BE Exporting results to CSV...")
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Export sequence statistics
        if self.results:
            df_stats = pd.DataFrame(self.results)
            stats_file = f'{output_dir}/sequence_statistics.csv'
            df_stats.to_csv(stats_file, index=False)
            print(f"\u2705 Statistics exported to {stats_file}")
        
        # Export ORF data
        if self.all_orfs:
            df_orfs = pd.DataFrame(self.all_orfs)
            orfs_file = f'{output_dir}/orfs_found.csv'
            df_orfs.to_csv(orfs_file, index=False)
            print(f"\u2705 ORFs exported to {orfs_file}")
        
        # Export translation results
        translations = self.translate_sequences()
        if translations:
            df_trans = pd.DataFrame(translations)
            trans_file = f'{output_dir}/translations.csv'
            df_trans.to_csv(trans_file, index=False)
<<<<<<< HEAD
            print(f"\u2705 Translations exported to {trans_file}")
=======
            print(f"✅ Translations exported to {trans_file}")
        
        # Export codon usage
        if include_codon:
            codon_data = self.analyze_codon_usage()
            if codon_data:
                df_codon = pd.DataFrame(codon_data)
                codon_file = f'{output_dir}/codon_usage.csv'
                df_codon.to_csv(codon_file, index=False)
                print(f"✅ Codon usage exported to {codon_file}")
        
        # Export motif findings
        if include_motifs:
            motif_data = self.find_common_motifs()
            if motif_data:
                df_motifs = pd.DataFrame(motif_data)
                motifs_file = f'{output_dir}/common_motifs.csv'
                df_motifs.to_csv(motifs_file, index=False)
                print(f"✅ Common motifs exported to {motifs_file}")
    
    def analyze_codon_usage(self):
        """Analyze codon usage patterns across all sequences"""
        print("\n🧪 Analyzing codon usage...")
        
        from collections import defaultdict
        
        codon_counts = defaultdict(int)
        total_codons = 0
        
        # Genetic code for reference
        genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        for seq_id, record in self.sequences.items():
            seq = str(record.seq).upper()
            
            # Analyze codons in frame 0 (most common reading frame)
            for i in range(0, len(seq)-2, 3):
                codon = seq[i:i+3]
                if len(codon) == 3 and all(n in 'ATGC' for n in codon):
                    codon_counts[codon] += 1
                    total_codons += 1
        
        # Calculate frequencies and organize by amino acid
        codon_data = []
        for codon, count in sorted(codon_counts.items(), key=lambda x: x[1], reverse=True):
            aa = genetic_code.get(codon, 'X')
            frequency = (count / total_codons * 100) if total_codons > 0 else 0
            
            codon_data.append({
                'Codon': codon,
                'Amino_Acid': aa,
                'Count': count,
                'Frequency_%': round(frequency, 2),
                'Per_Thousand': round((count / total_codons * 1000), 2) if total_codons > 0 else 0
            })
        
        print(f"✅ Analyzed {total_codons:,} codons")
        print(f"   Unique codons found: {len(codon_counts)}")
        
        return codon_data
    
    def find_motif(self, motif):
        """
        Find specific motif in all sequences
        
        Args:
            motif (str): DNA motif to search for (e.g., 'TATA', 'CAAT')
        """
        print(f"\n🔎 Searching for motif: {motif}")
        
        results = []
        total_occurrences = 0
        
        for seq_id, record in self.sequences.items():
            seq = str(record.seq).upper()
            motif_upper = motif.upper()
            positions = []
            
            # Find all occurrences
            start = 0
            while True:
                pos = seq.find(motif_upper, start)
                if pos == -1:
                    break
                positions.append(pos + 1)  # 1-based position
                start = pos + 1
            
            if positions:
                total_occurrences += len(positions)
                results.append({
                    'Sequence_ID': seq_id,
                    'Motif': motif,
                    'Count': len(positions),
                    'Positions': ', '.join(map(str, positions[:10])) + ('...' if len(positions) > 10 else ''),
                    'First_Position': positions[0],
                    'Last_Position': positions[-1]
                })
        
        print(f"✅ Found {total_occurrences} occurrences in {len(results)} sequences")
        
        return results
    
    def find_common_motifs(self):
        """Find common regulatory motifs"""
        print("\n🔍 Searching for common regulatory motifs...")
        
        # Common promoter and regulatory motifs
        common_motifs = {
            'TATA_box': 'TATAAA',
            'CAAT_box': 'CCAAT',
            'GC_box': 'GGGCGG',
            'Kozak': 'GCCACCATGG',
            'PolyA_signal': 'AATAAA',
            'Start_codon': 'ATG',
            'Stop_TAA': 'TAA',
            'Stop_TAG': 'TAG',
            'Stop_TGA': 'TGA'
        }
        
        all_motif_results = []
        
        for motif_name, motif_seq in common_motifs.items():
            results = self.find_motif(motif_seq)
            for result in results:
                result['Motif_Name'] = motif_name
                all_motif_results.append(result)
        
        return all_motif_results
>>>>>>> ad48967 (Save local changes before pulling new updates)
    
    def generate_summary_report(self):
        """Generate a summary report"""
        print("\n" + "="*60)
        print("\U0001F4CB SEQUENCE ANALYSIS SUMMARY REPORT")
        print("="*60)
        
        if not self.results:
            print("No analysis results available.")
            return
        
        df = pd.DataFrame(self.results)
        
        print(f"\n\U0001F4CA GENERAL STATISTICS")
        print(f"   Total Sequences: {len(df)}")
        print(f"   Total Base Pairs: {df['Length_bp'].sum():,}")
        print(f"   Average Length: {df['Length_bp'].mean():.2f} bp")
        print(f"   Median Length: {df['Length_bp'].median():.2f} bp")
        print(f"   Min Length: {df['Length_bp'].min()} bp")
        print(f"   Max Length: {df['Length_bp'].max()} bp")
        
        print(f"\n\U0001F9EC GC CONTENT ANALYSIS")
        print(f"   Average GC Content: {df['GC_Content_%'].mean():.2f}%")
        print(f"   Median GC Content: {df['GC_Content_%'].median():.2f}%")
        print(f"   Min GC Content: {df['GC_Content_%'].min():.2f}%")
        print(f"   Max GC Content: {df['GC_Content_%'].max():.2f}%")
        
        print(f"\n\U0001F4E6 NUCLEOTIDE COMPOSITION (Total)")
        print(f"   A: {df['A_Count'].sum():,}")
        print(f"   T: {df['T_Count'].sum():,}")
        print(f"   G: {df['G_Count'].sum():,}")
        print(f"   C: {df['C_Count'].sum():,}")
        print(f"   N: {df['N_Count'].sum():,}")
        
        if self.all_orfs:
            orf_df = pd.DataFrame(self.all_orfs)
            print(f"\n\U0001F50D ORF ANALYSIS")
            print(f"   Total ORFs Found: {len(orf_df)}")
            print(f"   Average ORF Length: {orf_df['Length_bp'].mean():.2f} bp")
            print(f"   Longest ORF: {orf_df['Length_bp'].max()} bp")
            print(f"   ORFs on + strand: {len(orf_df[orf_df['Strand'] == '+'])}")
            print(f"   ORFs on - strand: {len(orf_df[orf_df['Strand'] == '-'])}")
        
        print("\n" + "="*60)


def main():
    """Main function with CLI interface"""
    parser = argparse.ArgumentParser(
        description='Comprehensive Sequence Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python sequence_analyzer.py input.fasta
  python sequence_analyzer.py input.fasta --min-orf 150
  python sequence_analyzer.py input.fasta --output-dir my_results
<<<<<<< HEAD
  python sequence_analyzer.py input.fasta --codon-usage --viz-codon
=======
  python sequence_analyzer.py input.fasta --find-motif TATAAA
  python sequence_analyzer.py input.fasta --codon-usage
>>>>>>> ad48967 (Save local changes before pulling new updates)
        '''
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--min-orf', type=int, default=100,
                       help='Minimum ORF length in bp (default: 100)')
    parser.add_argument('--output-dir', default='output',
                       help='Output directory (default: output)')
    parser.add_argument('--no-viz', action='store_true',
                       help='Skip visualization generation')
<<<<<<< HEAD
    parser.add_argument('--codon-usage', action='store_true',
                       help='Perform codon usage analysis and output CSV')
    parser.add_argument('--viz-codon', action='store_true',
                       help='Generate codon usage frequency plot (use --codon-usage)')
=======
    parser.add_argument('--find-motif', type=str, default=None,
                       help='Search for specific motif (e.g., TATAAA)')
    parser.add_argument('--codon-usage', action='store_true',
                       help='Perform codon usage analysis')
>>>>>>> ad48967 (Save local changes before pulling new updates)
    
    args = parser.parse_args()
    
    analyzer = SequenceAnalyzer(args.fasta_file)
    if not analyzer.read_fasta():
        sys.exit(1)
    
<<<<<<< HEAD
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=args.min_orf)
    
    if args.codon_usage:
        viz_dir = f"{args.output_dir}/visualizations"
        analyzer.codon_usage_analysis(output_dir=viz_dir, visualize=args.viz_codon)
    
=======
    # Run basic analysis
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=args.min_orf)
    
    # Run optional analyses
    if args.find_motif:
        motif_results = analyzer.find_motif(args.find_motif)
        if motif_results:
            print("\n🎯 Motif Search Results:")
            for result in motif_results[:5]:  # Show first 5
                print(f"   {result['Sequence_ID']}: {result['Count']} occurrences")
                print(f"      Positions: {result['Positions']}")
    
    if args.codon_usage:
        codon_data = analyzer.analyze_codon_usage()
        if codon_data:
            print("\n📊 Top 10 Most Frequent Codons:")
            for i, codon in enumerate(codon_data[:10], 1):
                print(f"   {i}. {codon['Codon']} ({codon['Amino_Acid']}): "
                      f"{codon['Count']} ({codon['Frequency_%']:.2f}%)")
    
    # Generate visualizations
>>>>>>> ad48967 (Save local changes before pulling new updates)
    if not args.no_viz:
        viz_dir = f"{args.output_dir}/visualizations"
        analyzer.generate_visualizations(output_dir=viz_dir)
    
<<<<<<< HEAD
    results_dir = f"{args.output_dir}/results"
    analyzer.export_to_csv(output_dir=results_dir)
    analyzer.generate_summary_report()
    
    print("\n\u2705 Analysis complete! Check the output directory for results.")
=======
    # Export results (includes codon and motif analysis if performed)
    results_dir = f"{args.output_dir}/results"
    analyzer.export_to_csv(
        output_dir=results_dir,
        include_codon=args.codon_usage,
        include_motifs=True
    )
    
    # Generate summary
    analyzer.generate_summary_report()
    
    print("\n✅ Analysis complete! Check the output directory for results.")
    print(f"📁 Results: {results_dir}")
    if not args.no_viz:
        print(f"📊 Visualizations: {viz_dir}")

>>>>>>> ad48967 (Save local changes before pulling new updates)

if __name__ == "__main__":
    main()
