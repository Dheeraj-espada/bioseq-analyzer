#!/usr/bin/env python3
"""
Test script for codon usage and motif finding features
"""

from sequence_analyzer import SequenceAnalyzer
import pandas as pd

def test_codon_usage():
    """Test codon usage analysis"""
    print("="*60)
    print("Testing Codon Usage Analysis")
    print("="*60)
    
    analyzer = SequenceAnalyzer('test_sequences.fasta')
    
    if not analyzer.read_fasta():
        print("Error reading file!")
        return
    
    # Run codon analysis
    codon_data = analyzer.analyze_codon_usage()
    
    if codon_data:
        df = pd.DataFrame(codon_data)
        
        print("\n📊 Top 15 Most Frequent Codons:")
        print(df.head(15).to_string(index=False))
        
        print("\n📊 Codon Usage by Amino Acid:")
        aa_summary = df.groupby('Amino_Acid').agg({
            'Count': 'sum',
            'Codon': 'count'
        }).rename(columns={'Codon': 'Num_Codons'})
        print(aa_summary.sort_values('Count', ascending=False).head(10))
        
        print(f"\n✅ Total unique codons: {len(df)}")
        print(f"✅ Total codons analyzed: {df['Count'].sum():,}")
    else:
        print("❌ No codon data generated")


def test_motif_finding():
    """Test motif finding"""
    print("\n" + "="*60)
    print("Testing Motif Finding")
    print("="*60)
    
    analyzer = SequenceAnalyzer('test_sequences.fasta')
    
    if not analyzer.read_fasta():
        print("Error reading file!")
        return
    
    # Test specific motif
    print("\n🔍 Test 1: Finding ATG (start codon)")
    atg_results = analyzer.find_motif('ATG')
    if atg_results:
        for result in atg_results:
            print(f"\n   Sequence: {result['Sequence_ID']}")
            print(f"   Occurrences: {result['Count']}")
            print(f"   First 10 positions: {result['Positions']}")
    
    # Test TATA box
    print("\n🔍 Test 2: Finding TATAAA (TATA box)")
    tata_results = analyzer.find_motif('TATAAA')
    if tata_results:
        for result in tata_results:
            print(f"\n   Sequence: {result['Sequence_ID']}")
            print(f"   Occurrences: {result['Count']}")
            print(f"   Positions: {result['Positions']}")
    else:
        print("   ⚠️ No TATA box found (normal for coding sequences)")
    
    # Test common motifs
    print("\n🔍 Test 3: Finding all common regulatory motifs")
    common_motifs = analyzer.find_common_motifs()
    if common_motifs:
        df_motifs = pd.DataFrame(common_motifs)
        print(f"\n   Total motif occurrences: {len(df_motifs)}")
        
        # Summary by motif type
        motif_summary = df_motifs.groupby('Motif_Name')['Count'].sum().sort_values(ascending=False)
        print("\n   📊 Motif Summary:")
        for motif_name, count in motif_summary.items():
            print(f"      {motif_name}: {count} occurrences")


def test_custom_motif():
    """Test with user-defined motif"""
    print("\n" + "="*60)
    print("Testing Custom Motif Search")
    print("="*60)
    
    analyzer = SequenceAnalyzer('test_sequences.fasta')
    
    if not analyzer.read_fasta():
        print("Error reading file!")
        return
    
    # Test various custom motifs
    custom_motifs = ['GC', 'CG', 'GGGG', 'CCCC', 'GCCGC']
    
    for motif in custom_motifs:
        print(f"\n🔍 Searching for: {motif}")
        results = analyzer.find_motif(motif)
        if results:
            total_count = sum(r['Count'] for r in results)
            print(f"   ✅ Found in {len(results)} sequences")
            print(f"   Total occurrences: {total_count}")
        else:
            print(f"   ❌ Not found")


def test_full_pipeline():
    """Test complete pipeline with new features"""
    print("\n" + "="*60)
    print("Testing Full Analysis Pipeline")
    print("="*60)
    
    analyzer = SequenceAnalyzer('test_sequences.fasta')
    
    # Read file
    if not analyzer.read_fasta():
        return
    
    # Basic analysis
    analyzer.calculate_statistics()
    analyzer.find_orfs(min_length=100)
    
    # New features
    codon_data = analyzer.analyze_codon_usage()
    motif_data = analyzer.find_common_motifs()
    
    # Export everything
    analyzer.export_to_csv(
        output_dir='test_output/results',
        include_codon=True,
        include_motifs=True
    )
    
    print("\n✅ All files exported to test_output/results/")
    print("   - sequence_statistics.csv")
    print("   - orfs_found.csv")
    print("   - translations.csv")
    print("   - codon_usage.csv")
    print("   - common_motifs.csv")


if __name__ == "__main__":
    print("🧬 Testing New Features in BioSeq Analyzer\n")
    
    # Run all tests
    test_codon_usage()
    test_motif_finding()
    test_custom_motif()
    test_full_pipeline()
    
    print("\n" + "="*60)
    print("✅ All tests completed!")
    print("="*60)
