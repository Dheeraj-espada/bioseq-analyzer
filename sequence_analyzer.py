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

    def generate_visualizations(self, output_dir='output/visualizations'):
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
    
    def export_to_csv(self, output_dir='output/results'):
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
            print(f"\u2705 Translations exported to {trans_file}")
    
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
  python sequence_analyzer.py input.fasta --codon-usage --viz-codon
        '''
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--min-orf', type=int, default=100,
                       help='Minimum ORF length in bp (default: 100)')
    parser.add_argument('--output-dir', default='output',
                       help='Output directory (default: output)')
    parser.add_argument('--no-viz', action='store_true',
                       help='Skip visualization generation')
    parser.add_argument('--codon-usage', action='store_true',
                       help='Perform codon usage analysis and output CSV')
    parser.add_argument('--viz-codon', action='store_true',
                       help='Generate codon usage frequency plot (use --codon-usage)')
    
    args = parser.parse_args()
    
    analyzer = SequenceAnalyzer(args.fasta_file)
    if not analyzer.read_fasta():
        sys.exit(1)
    
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=args.min_orf)
    
    if args.codon_usage:
        viz_dir = f"{args.output_dir}/visualizations"
        analyzer.codon_usage_analysis(output_dir=viz_dir, visualize=args.viz_codon)
    
    if not args.no_viz:
        viz_dir = f"{args.output_dir}/visualizations"
        analyzer.generate_visualizations(output_dir=viz_dir)
    
    results_dir = f"{args.output_dir}/results"
    analyzer.export_to_csv(output_dir=results_dir)
    analyzer.generate_summary_report()
    
    print("\n\u2705 Analysis complete! Check the output directory for results.")

if __name__ == "__main__":
    main()
