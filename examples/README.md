# spaGAPA Examples

This directory contains example scripts and Jupyter notebooks demonstrating how to use spaGAPA.

## Prerequisites

Before running these examples, make sure you have:

1. **Installed spaGAPA**:
   ```bash
   cd spaGAPA
   conda activate spagapa
   pip install -e ".[dev]"
   ```

2. **Installed R and scAPAtrap**:
   ```R
   # In R console
   install.packages('devtools')
   devtools::install_github("BMILAB/scAPAtrap")
   ```

3. **Prepared your data**:
   - BAM file (sorted and indexed)
   - Spatial coordinates (CSV format)
   - Reference genome (FASTA)
   - Gene annotation (GTF)

## Examples

### 01_scapatrap_usage.py

Basic usage of scAPAtrap wrapper to identify APA sites from spatial transcriptomics data.

```bash
python 01_scapatrap_usage.py
```

**What it does**:
- Runs scAPAtrap on a BAM file
- Parses the output (peaks and counts)
- Applies basic spatial filtering

**Input**:
- `data/sample_spatial.bam` - Aligned reads
- `data/spatial_coordinates.csv` - Spot coordinates
- `data/genome.fa` - Reference genome
- `data/genes.gtf` - Gene annotation

**Output**:
- `results/scapatrap/peaks_meta.csv` - Peak metadata
- `results/scapatrap/peaks_counts.csv` - Count matrix

### Coming Soon

- `02_spatial_validation.ipynb` - Spatial validation of APA sites
- `03_gp_imputation.ipynb` - Gaussian process imputation
- `04_differential_apa.ipynb` - Differential APA analysis
- `05_visualization.ipynb` - Visualization examples
- `06_full_pipeline.ipynb` - Complete analysis pipeline

## Data Format

### Spatial Coordinates CSV

```csv
barcode,x,y
SPOT_1,100.5,200.3
SPOT_2,105.2,198.7
SPOT_3,110.8,195.1
```

### BAM File Requirements

- Must be sorted: `samtools sort input.bam -o sorted.bam`
- Must be indexed: `samtools index sorted.bam`
- Should contain cell/spot barcodes in tags (e.g., CB tag for 10x data)

## Troubleshooting

### R not found
```bash
# Check R installation
which Rscript
R --version

# If not installed, install R:
conda install r-base -c conda-forge
```

### scAPAtrap not found
```R
# In R console
if (!require("devtools")) install.packages("devtools")
devtools::install_github("BMILAB/scAPAtrap")

# Check installation
library(scAPAtrap)
```

### Memory issues
If you encounter memory issues with large datasets:
- Reduce `n_cores` parameter
- Process chromosomes separately
- Use a machine with more RAM

## Getting Help

- Check the [main documentation](../docs/)
- Open an issue on [GitHub](https://github.com/yourusername/spaGAPA/issues)
- Contact: your.email@example.com
