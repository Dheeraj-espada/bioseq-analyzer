# BioSeq Analyzer 🧬

A comprehensive Python-based bioinformatics tool for analyzing biological sequences with advanced features including statistics calculation, ORF detection, translation, codon usage analysis, motif finding, and publication-quality visualizations.

![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![BioPython](https://img.shields.io/badge/BioPython-1.81%2B-orange)

## 🌟 Features

### Core Analysis
- **📖 FASTA File Parsing**: Efficient parsing of single and multi-sequence FASTA files
- **📊 Sequence Statistics**: GC%, AT%, nucleotide composition, molecular weight
- **🔍 ORF Detection**: Find Open Reading Frames in all 6 reading frames
- **🧬 Translation**: Translate sequences and analyze protein properties

### Advanced Features
- **🧪 Codon Usage Analysis**: Calculate codon frequency and usage patterns
- **🔎 Motif Finding**: Search for DNA motifs and regulatory elements
- **📋 Regulatory Elements**: Auto-detect TATA box, CAAT box, Kozak, PolyA signals
- **📈 Visualizations**: 10+ publication-quality plots (300 DPI)

### Output
- **💾 CSV Export**: All results in structured format
- **📊 High-Resolution Plots**: PNG images ready for publication

## 🚀 Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Setup
```bash
# Clone repository
git clone https://github.com/Dheeraj-espada/bioseq-analyzer.git
cd bioseq-analyzer

# Install dependencies
pip install -r requirements.txt
```

### Dependencies
```
biopython>=1.81
pandas>=2.0.0
matplotlib>=3.7.0
seaborn>=0.12.0
numpy>=1.24.0
```

## 📖 Quick Start

### Basic Usage
```bash
# Simple analysis
python sequence_analyzer.py sequences.fasta

# With codon usage analysis
python sequence_analyzer.py sequences.fasta --codon-usage

# Search for specific motif
python sequence_analyzer.py sequences.fasta --find-motif TATAAA

# Complete analysis with custom output
python sequence_analyzer.py sequences.fasta \
    --codon-usage \
    --find-motif ATG \
    --min-orf 150 \
    --output-dir my_results
```

### Command Line Arguments
```
positional arguments:
  fasta_file            Input FASTA file

optional arguments:
  -h, --help            Show help message
  --min-orf MIN_ORF     Minimum ORF length in bp (default: 100)
  --output-dir DIR      Output directory (default: output)
  --no-viz              Skip visualization generation
  --codon-usage         Perform codon usage analysis
  --find-motif MOTIF    Search for specific motif
```

## 📁 Output Structure
```
output/
├── results/
│   ├── sequence_statistics.csv      # Basic sequence metrics
│   ├── orfs_found.csv               # All detected ORFs
│   ├── translations.csv             # Translation results
│   ├── codon_usage.csv              # Codon frequency (optional)
│   └── common_motifs.csv            # Regulatory motifs
└── visualizations/
    ├── gc_content_distribution.png
    ├── length_distribution.png
    ├── nucleotide_composition.png
    ├── gc_vs_at_content.png
    └── orf_analysis.png
```

## 📊 Example Analysis

### Terminal Output
```
📖 Reading FASTA file: sequences.fasta
✅ Successfully read 5 sequences

📊 Calculating sequence statistics...
✅ Statistics calculated for 5 sequences

🔍 Finding ORFs (minimum length: 100 bp)...
✅ Found 12 ORFs across all sequences

🧪 Analyzing codon usage...
✅ Analyzed 1,234 codons

📊 Top 10 Most Frequent Codons:
   1. CTG (L): 145 (3.66%)
   2. GAG (E): 132 (3.40%)
   ...
```

## 🔬 Analysis Capabilities

### Sequence Statistics
- Length in base pairs
- GC and AT content percentages
- Individual nucleotide counts
- Molecular weight calculation

### ORF Detection
- Scans all 6 reading frames
- Detects start codons (ATG)
- Identifies stop codons (TAA, TAG, TGA)
- Reports positions and sequences

### Codon Usage
- Frequency of all 64 codons
- Per-thousand normalization
- Amino acid grouping
- Preferred codon identification

### Motif Finding
- Exact sequence matching
- Position reporting
- Custom motif search
- Common regulatory elements

## 📚 Example Datasets

### Test Data
Sample test file included: `test_sequences.fasta`

### Real Data
Download from NCBI using `download_ncbi_data.py`:
```bash
python download_ncbi_data.py
```

## 🛠️ Development

### Project Structure
```
bioseq-analyzer/
├── sequence_analyzer.py         # Main analysis tool
├── download_ncbi_data.py         # NCBI data downloader
├── test_new_features.py          # Test suite
├── requirements.txt              # Dependencies
├── README.md                     # This file
└── test_sequences.fasta          # Sample data
```

### Testing
```bash
# Run test suite
python test_new_features.py

# Verify installation
./verify_week1.sh
```

## 📊 Performance

- **Speed**: Analyzes ~10 Mbp in under 30 seconds
- **Memory**: Efficient handling of large FASTA files
- **Output**: High-resolution (300 DPI) publication-ready figures

## 🤝 Contributing

Contributions welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License.

## 👤 Author

**Dheeraj**
- GitHub: [@Dheeraj-espada](https://github.com/Dheeraj-espada)

## 🙏 Acknowledgments

- BioPython community for excellent documentation
- Bioinformatics best practices

## 📞 Contact

For questions or collaboration:
- Open an issue on GitHub
- Email: dheerajbabu133@gmail.com

---

⭐ If you found this tool helpful, please give it a star on GitHub!
