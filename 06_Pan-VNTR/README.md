## Pan-VNTR analysis

***

A Snakemake-based automated workflow for Pan-VNTR analysis based on pangenome graph.

### 🚀 Quick Start

### Installation & Usage
1. Clone Repository
```
git clone https://github.com/Hanjunmin/complex_region.git
cd 06_Pan-VNTR

```
2. Setup Environment
The workflow can run either in a Conda environment or inside a Singularity container.
#### Using Conda
```
conda env create -f environment.yaml
conda activate pan_VNTR
```
#### Using Singularity
```
singularity exec \
  pan_VNTR_v1.0.1.sif \
  snakemake -s ../Pan_VNTR/PanVNTR_anno.smk \
            --configfile config.json \
            -j 50 -p
```
3. Prepare a New Run Directory and Configuration
```
# edit configuration files
mkdir -p run
cd run
touch config.json
```
4. Run Pipeline
```
# Execute workflow (with 20 jobs now)
snakemake -n ../Pan_VNTR/PanVNTR_anno.smk -j 20
```

### ⚙️ Configuration
The `config.json` file allows you to define all paths and parameters. Here's a comprehensive example:

```
# config.json

{
  "threads": 30,
  "script_dir": "06_Pan-VNTR/Pan_VNTR/",
  "base_dir": "06_Pan-VNTR/Pan_VNTR/run/",
  "data_dir": "06_Pan-VNTR/Pan_VNTR/example/"
}

```
## Configuration Explanation

- **`script_dir`**: Directory containing workflow-related scripts, including Snakemake rules, Python scripts, and other auxiliary files. This directory is usually part of the cloned GitHub repository and is used as the main location for workflow execution.


- **`base_dir`**: Working directory (output files will be generated here)

- **`data_dir`**: Input/example data directory

> ⚠️ All paths should be absolute to avoid errors in Snakemake.

## Example Data

The example dataset includes the following files

```
example_dataset/
├── example_mc.gfa                 # Graph
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
The analysis pipeline generates `run/output/VNTR.anno`.

### Annotation Results

**Example data from `VNTR.anno`:**
```
chr	start	end	motif	CHM13v2	GRCh38	C001-CHA-E01-Pat	C002-CHA-E02-Mat
chr1	100000	101500	ACGTGTA,TTGACGA	12	12	10	11
chr5	200500	201800	TGCAAGT	8	8	9	7
chr3	305200	306700	GGATCAG,ATGTCGA	15	15	14	16
chr2	450100	451300	CAGTGTC	5	5	6	5
chr7	510200	511500	AGTTGCA	9	9	8	10

```
**Columns in `VNTR.anno`:**
| Column                                        | Description                                                                                                                                                        |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `chr`                                         | Chromosome.                                    |
| `start`                                       | Start position.                                                                                                           |
| `end`                                         | End position.                                                                                                             |
| `motif`                                       | Repeat motif(s) of the VNTR. Multiple motifs in a VNTR are separated by commas.                                         |
| `CHM13v2`, `GRCh38`, `C001-CHA-E01-Pat`, etc. | Copy number of the VNTR in the corresponding sample.|
