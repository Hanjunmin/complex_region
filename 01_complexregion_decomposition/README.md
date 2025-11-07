## Complex Region Decomposition

***

A Snakemake-based automated workflow for genome-wide identification and characterization of complex genomic regions. This framework classifies genomic regions into simple regions and seven complex region types (CC1-CC7) through sliding window analysis.

### 🚀 Quick Start

### Installation & Usage
1. Clone Repository
```
git clone https://github.com/Hanjunmin/complex_region.git
cd 01_complexregion_decomposition
```
2. Setup Environment
```
conda env create -f envs/decomposition.yaml
conda activate complex-region-decomposition
```
3. Prepare a New Run Directory and Configuration
```
# Copy and edit configuration files
mkdir -p run
cd run
cp ../complex_decom/config.json ./
```
4. Run Pipeline
```
# Execute workflow (with 20 jobs now)
snakemake -n ../complex_decom/complex_decom.smk -j 20
```

### ⚙️ Configuration
The `config.json` file allows you to define all paths and parameters. Here's a comprehensive example:

```
# config.json

{
  "threads": 30,
  "script_dir": "01_complexregion_decomposition/complex_decom/",
  "base_dir": "01_complexregion_decomposition/run/",
  "data_dir": "01_complexregion_decomposition/example/",
  "reference": "/home/chm13v2.0.fa",
  "software": {
    "gfabase": "/home/pangenome/software/gfabase",
    "RM": "/home/Public/software/RepeatMasker-4.1.4/",
    "TRF": "/home/Public/software/trf-4.09.1/trf"
  }
}

```
## Configuration Explanation

- **script_dir**: Directory containing workflow-related scripts, including Snakemake rules, Python scripts, and other auxiliary files. This directory is usually part of the cloned GitHub repository and is used as the main location for workflow execution.


- **base_dir**: Working directory (output files will be generated here)

- **data_dir**: Input/example data directory

- **reference**: Reference genome sequence path

- **software**: Dictionary of software paths
  - **gfabase**: Path to gfabase installation
  - **RM**: Path to RepeatMasker
  - **TRF**: Path to Tandem Repeat Finder

> ⚠️ All paths should be absolute to avoid errors in Snakemake.

## Example Data

The example dataset includes the following files (available on Zenodo: [Download here]):
```
example_dataset/
├── 1p36.bed                 # Genomic regions
├── 1p36.gfa                 # Graph
├── 1p36.vcf                 # Variant calls
├── FA/                      # FASTA sequences directory
│   ├── C001-CHA-E01-Mat.fa  # Sample FASTA files
│   ├── C002-CHA-E02-Pat.fa
│   ├── ... (other samples)
│   └── <sample_id>.fa       # Each file named as <sample_id>.fa
└── samlis.txt               # Sample list
```


Sample List： (`samlis.txt`)

| Sample ID           | Type       |
|--------------------|------------|
| CHM13v2             | reference  |
| C001-CHA-E01-Mat    | .          |
| C001-CHA-E01-Pat    | .          |

> ⚠️ Make sure that the sample IDs in `samlis.txt` match the FASTA filenames in the `FA/` directory (e.g., `C001-CHA-E01-Mat.fa`).

## Output Data
The analysis pipeline generates the following output files in the `run/03_model/` directory:

```
run/03_model/
├── final_anno.complex.png    # Visualization plot (PNG format)
├── final_anno.complex.pdf    # Visualization plot (PDF format)
└── final.anno                # Annotation results file
```




### Annotation Results

**`final_anno.complex.png`**:  
![Complex Region Decomposition Analysis](https://raw.githubusercontent.com/Hanjunmin/complex_region/main/01_complexregion_decomposition/complex_decom/output/final_anno.complex.png)

This plot shows the complexity annotation of genomic regions.  
- Regions are colored by type: **simple** or complex clusters (**CC1–CC7**).  
- The **y-axis** represents the cluster probability.  
- The **x-axis** shows the genomic position of each region.

**Example data from `final.anno`:**
```
chr     start   end     cluster simple  CC1     CC2     CC3     CC4     CC5     CC6     CC7
chr1    12054001        12055000        simple  0.8592633       0.011185256     0.06440703      0.005208004     0.002919506     0.0015025395    0.0013939425    0.054120403
chr1    12054501        12055500        simple  0.85257804      0.011569369     0.06705679      0.0053913686    0.0030226677    0.0015524194    0.0014420691    0.057387255
chr1    12055001        12056000        simple  0.8765505       0.011200878     0.058783863     0.005138364     0.0028743104    0.001511905     0.0013829212    0.04255722
```
**Columns in `final.anno`:**
| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `start` | Start genomic position |
| `end` | End genomic position |
| `cluster` | Cluster classification |
| `simple` | Simple region probability score |
| `CC1` - `CC7` | Complex component probability scores |

## Overview

`complex_decom.smk` workflow for analyzing complex genomic regions, including preprocessing, index construction, and clustering analysis.

The workflow includes three main modules:

1. **Preprocessing (`preproc`)**  
   Generates basic files for regions, including genome windows, sample paths, and segment databases...

2. **Index Construction (`index_construction`)**  
   Builds indices and annotation files for each genome windows.

3. **Clustering Model (`clustering_model`)**  
   Performs clustering analysis on genomic regions using a DeepClustering model.
