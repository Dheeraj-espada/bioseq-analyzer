# Usage Examples

## Basic Analysis
```bash
# Analyze a FASTA file
python sequence_analyzer.py mysequences.fasta
```

## Advanced Usage

### Codon Usage Analysis
```bash
python sequence_analyzer.py sequences.fasta --codon-usage
```

### Motif Finding
```bash
# Search for TATA box
python sequence_analyzer.py promoters.fasta --find-motif TATAAA

# Search for start codon
python sequence_analyzer.py genes.fasta --find-motif ATG
```

### Custom Output Directory
```bash
python sequence_analyzer.py data.fasta --output-dir my_analysis
```

### Combined Analysis
```bash
python sequence_analyzer.py sequences.fasta \
    --codon-usage \
    --find-motif ATG \
    --min-orf 150 \
    --output-dir comprehensive_results
```

## Working with Real Data

### Download from NCBI
```bash
# Edit email in download_ncbi_data.py first!
python download_ncbi_data.py

# Analyze downloaded data
python sequence_analyzer.py real_data/human_cancer_genes.fasta --codon-usage
```

## Output Files

All results are saved in the output directory:

- `results/sequence_statistics.csv` - Basic sequence metrics
- `results/orfs_found.csv` - All detected ORFs
- `results/codon_usage.csv` - Codon frequency data
- `results/common_motifs.csv` - Found motifs
- `visualizations/*.png` - All plots (300 DPI)

## Tips

1. Use `--min-orf` to adjust ORF detection sensitivity
2. Use `--no-viz` to skip visualizations for faster processing
3. Check `output/visualizations/` for all generated plots
4. Open CSV files in Excel/LibreOffice for easy viewing
