#!/usr/bin/env python3
"""
Enhanced Sequence Parser - Multi-format Support
Supports: FASTA, GenBank, EMBL
"""

from Bio import SeqIO
from pathlib import Path
import sys
import pandas as pd


class EnhancedSequenceParser:
    """Parse multiple sequence file formats"""
    
    SUPPORTED_FORMATS = {
        '.fasta': 'fasta', '.fa': 'fasta', '.fna': 'fasta',
        '.ffn': 'fasta', '.faa': 'fasta', '.frn': 'fasta',
        '.gb': 'genbank', '.gbk': 'genbank', '.genbank': 'genbank',
        '.embl': 'embl', '.emb': 'embl'
    }
    
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.format = self._detect_format()
        self.sequences = {}
        
    def _detect_format(self):
        """Auto-detect file format from extension"""
        suffix = self.filepath.suffix.lower()
        if suffix in self.SUPPORTED_FORMATS:
            return self.SUPPORTED_FORMATS[suffix]
        return self._detect_from_content()
    
    def _detect_from_content(self):
        """Detect format by reading first line"""
        try:
            with open(self.filepath, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('>'):
                    return 'fasta'
                elif first_line.startswith('LOCUS'):
                    return 'genbank'
                elif first_line.startswith('ID'):
                    return 'embl'
                else:
                    print(f"⚠️  Unknown format, trying FASTA...")
                    return 'fasta'
        except Exception as e:
            print(f"❌ Error reading file: {e}")
            sys.exit(1)
    
    def parse(self):
        """Parse the sequence file"""
        print(f"📖 Reading {self.format.upper()} file: {self.filepath.name}")
        
        try:
            for record in SeqIO.parse(self.filepath, self.format):
                self.sequences[record.id] = record
            
            if len(self.sequences) == 0:
                print(f"⚠️  No sequences found")
                return None
            
            print(f"✅ Parsed {len(self.sequences)} sequences")
            return self.sequences
            
        except Exception as e:
            print(f"❌ Error: {e}")
            # Try alternative formats
            for fmt in ['fasta', 'genbank', 'embl']:
                if fmt != self.format:
                    try:
                        self.sequences = {}
                        for record in SeqIO.parse(self.filepath, fmt):
                            self.sequences[record.id] = record
                        if len(self.sequences) > 0:
                            print(f"✅ Parsed as {fmt.upper()}")
                            self.format = fmt
                            return self.sequences
                    except:
                        continue
            print(f"❌ Could not parse file")
            return None
    
    def get_metadata(self):
        """Extract metadata from sequences"""
        metadata = []
        for seq_id, record in self.sequences.items():
            meta = {
                'ID': seq_id,
                'Name': record.name,
                'Description': record.description,
                'Length': len(record.seq),
                'Format': self.format
            }
            
            if self.format == 'genbank':
                annotations = record.annotations
                meta['Organism'] = annotations.get('organism', 'N/A')
                meta['Taxonomy'] = ', '.join(annotations.get('taxonomy', []))[:100]
                meta['Features'] = len(record.features)
            
            metadata.append(meta)
        return metadata
    
    def convert_format(self, output_file, output_format='fasta'):
        """Convert sequences to different format"""
        if not self.sequences:
            print("⚠️  No sequences to convert")
            return False
        
        try:
            SeqIO.write(self.sequences.values(), output_file, output_format)
            print(f"✅ Converted to {output_format.upper()}: {output_file}")
            return True
        except Exception as e:
            print(f"❌ Conversion failed: {e}")
            return False


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhanced Sequence Parser')
    parser.add_argument('input_file', help='Input sequence file')
    parser.add_argument('--show-metadata', action='store_true',
                       help='Show metadata')
    parser.add_argument('--convert', help='Convert to format')
    parser.add_argument('--output', help='Output file')
    
    args = parser.parse_args()
    
    # Parse file
    seq_parser = EnhancedSequenceParser(args.input_file)
    sequences = seq_parser.parse()
    
    if not sequences:
        sys.exit(1)
    
    # Show metadata
    if args.show_metadata:
        metadata = seq_parser.get_metadata()
        df = pd.DataFrame(metadata)
        print("\n📋 Metadata:")
        print(df.to_string(index=False))
    
    # Convert format
    if args.convert and args.output:
        seq_parser.convert_format(args.output, args.convert)
