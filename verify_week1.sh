#!/bin/bash
echo "🔍 Verifying Week 1 Project..."
echo ""

# Check main script
if [ -f "sequence_analyzer.py" ]; then
    echo "✅ sequence_analyzer.py exists"
else
    echo "❌ sequence_analyzer.py missing"
fi

# Check it runs
python sequence_analyzer.py test_sequences.fasta --output-dir verify_test > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Script runs without errors"
else
    echo "❌ Script has errors"
fi

# Check output files
if [ -f "verify_test/results/sequence_statistics.csv" ]; then
    echo "✅ CSV export works"
else
    echo "❌ CSV export failed"
fi

if [ -f "verify_test/visualizations/gc_content_distribution.png" ]; then
    echo "✅ Visualizations work"
else
    echo "❌ Visualizations failed"
fi

# Check codon usage
python sequence_analyzer.py test_sequences.fasta --codon-usage --output-dir verify_test > /dev/null 2>&1
if [ -f "verify_test/results/codon_usage.csv" ]; then
    echo "✅ Codon usage analysis works"
else
    echo "❌ Codon usage failed"
fi

# Check Git
git status > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Git repository initialized"
else
    echo "❌ Git not initialized"
fi

# Check remote
git remote -v | grep -q "github.com"
if [ $? -eq 0 ]; then
    echo "✅ GitHub remote configured"
else
    echo "❌ GitHub remote not configured"
fi

# Clean up
rm -rf verify_test

echo ""
echo "📋 Verification complete!"
