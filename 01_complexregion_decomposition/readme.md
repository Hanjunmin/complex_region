## Complex Region Decomposition

***

A Snakemake-based automated workflow for genome-wide identification and characterization of complex genomic regions. This framework classifies genomic regions into simple regions and seven complex region types (CC1-CC7) through sliding window analysis.

### 🚀 Quick Start

### Installation & Usage
1. Clone Repository
```
git clone https://github.com/your-username/complex-region-decomposition.git
cd complex-region-decomposition
```
2. Setup Environment
```
conda env create -f envs/decomposition.yaml
conda activate complex-region-decomposition
```
3. Configuration
```
# Copy and edit configuration files
cp config/config.yaml.example config/config.yaml
cp config/samples.tsv.example config/samples.tsv

# Edit samples.tsv with your sample information
```
4. Run Pipeline
```
# Dry-run (check workflow)
snakemake -n --use-conda

# Execute workflow (with 8 cores)
snakemake --cores 8 --use-conda
```

### ⚙️ Configuration
The config/config.yaml file allows you to define all paths and parameters. Here's a comprehensive example:

```
# config/config.yaml

# ========== PATH CONFIGURATION ==========
paths:
  # Reference genomes and annotations
  reference:
    genome: "resources/genomes/{genome}.fa"
    genome_index: "resources/genomes/{genome}.fa.fai"
    annotation: "resources/annotations/{genome}.gtf"
    chromosome_sizes: "resources/genomes/{genome}.chrom.sizes"
  
  # Input data directories
  input:
    samples: "config/samples.tsv"
    bam_directory: "data/alignments/"
    bigwig_directory: "data/bigwigs/"
    bed_directory: "data/bed_files/"
    
  # Output directories
  output:
    base_dir: "results"
    intermediate: "results/intermediate"
    final: "results/final"
    clusters: "results/clusters"
    reports: "results/reports"
    visualization: "results/visualization"
    logs: "logs"
    
  # External resources and databases
  resources:
    repeat_masker: "resources/repeats/{genome}_rmsk.bed"
    conservation: "resources/conservation/{genome}_phastCons.bw"
    blacklist: "resources/blacklists/{genome}_blacklist.bed"
    gc_content: "resources/gc_content/{genome}_gc.bw"
    
  # Script paths
  scripts:
    feature_calculation: "scripts/calculate_features.py"
    clustering_analysis: "scripts/decompose_regions.py"
    visualization: "scripts/plot_results.R"
    utilities: "scripts/utils/"

# ========== ANALYSIS PARAMETERS ==========
sliding_window:
  window_size: 1000    # 1 kbp windows
  step_size: 500       # 500bp step size
  min_region_size: 50  # minimum region size to keep

clustering:
  method: "hdbscan"
  min_cluster_size: 50
  min_samples: 5
  n_clusters: 8        # Simple + CC1-CC7
  random_state: 42

features:
  # Feature calculation settings
  gc_content: true
  sequence_complexity: true
  repeat_content: true
  conservation_scores: true
  cpg_content: true
  
  # Normalization
  normalize: true
  scale_features: true

output:
  formats: ["bed", "tsv", "bigwig", "json"]
  create_reports: true
  include_visualizations: true
  compression: false
```

