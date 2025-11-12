#!/usr/bin/env python3
"""
Complete Week 4 Pipeline - Integrated Analysis
Run from bioseq-analyzer root directory
"""

import sys
import os
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from enhanced_sequence_parser import EnhancedSequenceParser
from html_report_generator import HTMLReportGenerator
from sequence_analyzer import SequenceAnalyzer
import pandas as pd


def run_complete_pipeline(input_file, output_dir="week4_analysis"):
    """Run complete analysis pipeline"""
    
    print("="*70)
    print("🧬 COMPLETE SEQUENCE ANALYSIS PIPELINE - WEEK 4")
    print("="*70)
    
    # Step 1: Parse input (supports multiple formats)
    print("\n📖 Step 1: Parsing Input File...")
    parser = EnhancedSequenceParser(input_file)
    sequences = parser.parse()
    
    if not sequences:
        return False
    
    # Convert to FASTA for analysis
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fasta_file = f"{output_dir}/converted.fasta"
    parser.convert_format(fasta_file, 'fasta')
    
    # Step 2: Run comprehensive analysis
    print("\n📊 Step 2: Running Analysis...")
    analyzer = SequenceAnalyzer(fasta_file)
    analyzer.read_fasta()
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=100)
    codon_data = analyzer.analyze_codon_usage()
    
    # Step 3: Generate visualizations
    print("\n📈 Step 3: Generating Visualizations...")
    analyzer.generate_visualizations(
        output_dir=f"{output_dir}/visualizations",
        include_codon=True
    )
    
    # Step 4: Export data
    print("\n💾 Step 4: Exporting Data...")
    analyzer.export_to_csv(
        output_dir=f"{output_dir}/results",
        include_codon=True,
        include_motifs=True
    )
    
    # Step 5: Generate HTML report
    print("\n📄 Step 5: Generating HTML Report...")
    report = HTMLReportGenerator("Complete Sequence Analysis - Week 4")
    
    # Add summary
    df_stats = pd.DataFrame(analyzer.results)
    summary = {
        "Total Sequences": str(len(df_stats)),
        "Total Base Pairs": f"{df_stats['Length_bp'].sum():,}",
        "Average GC%": f"{df_stats['GC_Content_%'].mean():.2f}%",
        "ORFs Found": str(len(analyzer.all_orfs))
    }
    report.add_summary(summary)
    
    # Add key findings
    findings = [
        f"Analyzed {len(df_stats)} sequences totaling {df_stats['Length_bp'].sum():,} base pairs",
        f"GC content: {df_stats['GC_Content_%'].mean():.2f}% (range: {df_stats['GC_Content_%'].min():.2f}%-{df_stats['GC_Content_%'].max():.2f}%)",
        f"Identified {len(analyzer.all_orfs)} Open Reading Frames",
        f"Codon usage analysis completed with {len(codon_data)} unique codons" if codon_data else "Analysis complete"
    ]
    report.add_findings(findings)
    
    # Add visualizations
    viz_files = sorted(Path(f"{output_dir}/visualizations").glob("*.png"))
    if viz_files:
        report.add_visualizations("📊 Visualizations", [str(f) for f in viz_files])
    
    # Add statistics table
    report.add_table("📋 Sequence Statistics", df_stats.head(10))
    
    # Generate report
    report.generate(f"{output_dir}/analysis_report.html")
    
    # Summary
    print("\n" + "="*70)
    print("✅ PIPELINE COMPLETE!")
    print("="*70)
    print(f"📁 Results: {output_dir}/results/")
    print(f"📊 Visualizations: {output_dir}/visualizations/")
    print(f"📄 HTML Report: {output_dir}/analysis_report.html")
    print("="*70)
    
    return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python complete_pipeline.py <input_file> [output_dir]")
        print("\nExample:")
        print("  python week4_project/complete_pipeline.py test_sequences.fasta my_analysis")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "week4_analysis"
    
    run_complete_pipeline(input_file, output_dir)
