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
### Visualization Results

![Complex Region Decomposition Analysis](https://raw.githubusercontent.com/Hanjunmin/complex_region/main/01_complexregion_decomposition/complex_decom/output/final_anno.complex.png)

This plot displays:
- Genomic regions colored by cluster classification
- Probability scores for simple vs complex components
- Spatial distribution of complex features along the chromosome
- Component-wise contribution patterns (CC1-CC7)
