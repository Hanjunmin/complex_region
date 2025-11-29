## Pan-SD analysis

***

A Snakemake-based automated workflow for Pan-SD analysis based on pangenome graph.

**Workflow**:  
<img src="https://github.com/Hanjunmin/complex_region/blob/main/05_Pan-SD/pan_SD_workflow.png" width="600">

### 🚀 Quick Start

### Installation & Usage
Refer to the [complex_region_decomposition](https://github.com/Hanjunmin/complex_region/tree/main/01_complexregion_decomposition) repository for:

- **Environment setup** using Docker or Conda
- **RepeatMasker configuration** with similar database requirements

```
cd 05_Pan-SD/
wget https://synplotter.sjtu.edu.cn/disk2/COMPLEX/com_analysis/inter_kmer_DB/inter_seg_kmer.pkl ./Pan_SD/inter_sim/DB/
```

### Prepare a New Run Directory and Configuration
```
# edit configuration files
mkdir -p run
cd run
touch config.json
```
### Run Pipeline

```
singularity exec \
  --bind /home/miniconda3/envs/com_decom/share/RepeatMasker/Libraries/famdb:/opt/conda/envs/com_decom/share/RepeatMasker/Libraries/famdb \
  pan_SD_v1.0.1.sif \
  snakemake -s ../Pan_SD/PanSD_anno.smk \
            --configfile config.json \
            -j 50 -p
```

### ⚙️ Configuration
The `config.json` file allows you to define all paths and parameters. Here's a comprehensive example:

```
# config.json

{
  "threads": 30,
  "script_dir": "05_Pan-SD/Pan_SD/",
  "base_dir": "05_Pan-SD/Pan_SD/run/",
  "data_dir": "05_Pan-SD/Pan_SD/example/",
}

```
## Configuration Explanation

- **`script_dir`**: Directory containing workflow-related scripts, including Snakemake rules, Python scripts, and other auxiliary files. This directory is usually part of the cloned GitHub repository and is used as the main location for workflow execution.


- **`base_dir`**: Working directory (output files will be generated here)

- **`data_dir`**: Input/example data directory

> ⚠️ All paths should be absolute to avoid errors in Snakemake.

## Example Data

The example dataset includes the following files ([Download here](https://synplotter.sjtu.edu.cn/disk2/COMPLEX/com_analysis/example_2/)):

```
example_dataset/
├── example_mc.gfa                 # Graph
├── example_mc.vcf                 # Variant calls
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
The analysis pipeline generates the following output files in the `run/output/` directory:

```
├── reference_SD.anno
└── nonreference_SD.anno
```

### Annotation Results

**Example data from `*SD.anno`:**
```
chr_A	start_A	end_A	chr_B	start_B	end_B	sample	identity	init_pair	graph_alignA	graph_alignB	graph_align_ratioA	graph_align_ratioB
CHM13v2_chr1	100000	101500	CHM13v2_chr2	500000	501400	CHM13v2	0.98	chr1:100000-101500@chr2:500000-501400 (reference)	1	1	1.0	1.0
C001-CHA-E01-Pat_chr5	200500	201800	C001-CHA-E01-Pat_chr7	150300	151450	C001-CHA-E01-Pat	0.95	C001-CHA-E01-Pat_chr5:200500-201800@chr7:150300-151450 (non-reference)	1	0	1.0	0.0
C002-CHA-E02-Mat_chr3	305200	306700	C002-CHA-E02-Mat_chr10	450100	451250	C002-CHA-E02-Mat	0.93	C002-CHA-E02-Mat_chr3:305200-306700@chr10:450100-451250 (non-reference)	0	1	0.0	0.9

```
**Columns in `*SD.anno`:**
| Column                                     | Description                                                                                                                                                                                |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `chr_A`, `start_A`, `end_A`                | Coordinates of the SD in **region A**, based on the sample’s genome.                                                        |
| `chr_B`, `start_B`, `end_B`                | Coordinates of the SD in **region B**, based on the sample’s genome.                                                                                                                   |
| `sample`                                   | Sample name (e.g., `CHM13v2`, `GRCh38`, `C001-CHA-E01-Pat`).                                                                                                                         |
| `init_pair_A_B`                                | SD pair coordinates originates from the **reference** or **non-reference**. Each init_pair defines the “canonical” SD to compare across samples. |
| `identity`                                 | Sequence identity between regions A and B for this sample.                                                                                                                                 |
| `graph_alignA`, `graph_alignB`             | Binary values (0 or 1) indicating whether the sample aligns to region A or B in the pangenome graph.                                                                                       |
| `graph_align_ratioA`, `graph_align_ratioB` | Fraction of region A or B that aligns in the pangenome graph (continuous value between 0 and 1).                                                                                           |
