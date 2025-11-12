#!/usr/bin/env python3
"""Complete Week 4 Pipeline"""

import sys
import os
from pathlib import Path

# Add paths
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
    
    # Step 1: Parse
    print("\n📖 Step 1: Parsing...")
    parser = EnhancedSequenceParser(input_file)
    sequences = parser.parse()
    
    if not sequences:
        return False
    
    # Convert to FASTA
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fasta_file = f"{output_dir}/converted.fasta"
    parser.convert_format(fasta_file, 'fasta')
    
    # Step 2: Analyze
    print("\n📊 Step 2: Analyzing...")
    analyzer = SequenceAnalyzer(fasta_file)
    analyzer.read_fasta()
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=100)
    codon_data = analyzer.analyze_codon_usage()
    
    # Step 3: Visualize
    print("\n📈 Step 3: Visualizing...")
    analyzer.generate_visualizations(
        output_dir=f"{output_dir}/visualizations",
        include_codon=True
    )
    
    # Step 4: Export
    print("\n💾 Step 4: Exporting...")
    analyzer.export_to_csv(
        output_dir=f"{output_dir}/results",
        include_codon=True,
        include_motifs=True
    )
    
    # Step 5: Generate Report
    print("\n📄 Step 5: Generating Report...")
    report = HTMLReportGenerator("Complete Sequence Analysis - Week 4")
    
    df_stats = pd.DataFrame(analyzer.results)
    summary = {
        "Total Sequences": str(len(df_stats)),
        "Total Base Pairs": f"{df_stats['Length_bp'].sum():,}",
        "Average GC%": f"{df_stats['GC_Content_%'].mean():.2f}%",
        "ORFs Found": str(len(analyzer.all_orfs))
    }
    report.add_summary(summary)
    
    findings = [
        f"Analyzed {len(df_stats)} sequences",
        f"Average GC: {df_stats['GC_Content_%'].mean():.2f}%",
        f"Found {len(analyzer.all_orfs)} ORFs",
        "Analysis complete"
    ]
    report.add_findings(findings)
    
    viz_files = sorted(Path(f"{output_dir}/visualizations").glob("*.png"))
    if viz_files:
        report.add_visualizations("📊 Visualizations", [str(f) for f in viz_files])
    
    report.add_table("📋 Statistics", df_stats.head(10))
    report.generate(f"{output_dir}/analysis_report.html")
    
    print("\n" + "="*70)
    print("✅ COMPLETE!")
    print("="*70)
    print(f"📁 Results: {output_dir}/results/")
    print(f"📊 Visualizations: {output_dir}/visualizations/")
    print(f"📄 Report: {output_dir}/analysis_report.html")
    print("="*70)
    
    return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python week4_project/complete_pipeline.py <input_file> [output_dir]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "week4_analysis"
    
    run_complete_pipeline(input_file, output_dir)
