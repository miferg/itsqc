# ecmitser

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-managed-green.svg)](https://docs.conda.io/en/latest/miniconda.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`itsqc` is a Snakemake pipeline for the processing of eukaryotic Internal Transcribed Spacer (ITS) sequences.

The pipeline takes one or more FASTA files as input and performs a full workflow, including:
* File merging and sequence quality filtering
* Sequence orientation and dereplication
* ITS region extraction using **ITSx**

---

## Workflow

The pipeline follows these main processing steps:

`Input FASTA (*.fna)` &rarr; `1. Merge Files` &rarr; `2. Clean & Trim` &rarr; `3. Orient (vsearch)` &rarr; `4. Dereplicate (vsearch)` &rarr; `5. Extract ITS (ITSx)` &rarr; `6. Filter (non-ITS)

---

## 🛠️ Installation

`ecmitser` is built on Snakemake and uses Conda to manage all software dependencies automatically.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/miferg/itsqc.git
    cd itsqc
    ```

2.  **Create the Snakemake environment:**
    
    ```bash
    conda create -c conda-forge -c bioconda -n itsqc-env snakemake
    ```

3.  **Activate the environment:**
    ```bash
    conda activate itsqc-env
    ```

All other dependencies (like `vsearch`, `ITSx`, and Python libraries) will be automatically downloaded and installed by Snakemake into isolated environments when you first run the pipeline.

---

## 🚀 Usage

Run the pipeline from within the cloned `itsqc` directory.

### Example Command

```bash
snakemake --cores 4 --use-conda --config querydir="path/to/my_files" outdir="path/to/itsqc_out"
```

### Configuration Parameters
`--cores <int>`: The number of CPU cores to allocate to the pipeline. The minimum requirement is 4 cores.

`--use-conda`: This flag is mandatory. It tells Snakemake to use the .yaml files in the snakes/ directory to create and manage the necessary software environments.

`--config`: This flag is used to pass required paths to the pipeline:

- `querydir="<path>"`: (Required) The path to the directory containing your input FASTA files.

- `outdir="<path>"`: (Required) The path to the directory where all results will be saved.

---

## 📥 Input


Your input files must be:

- In a single directory (specified by `querydir`).

- In valid FASTA format.

- Named with a `.fna` file extension.

---

## 📤 Output


The pipeline will create the specified `outdir` and organize all results into subdirectories.

### Directory Structure

```
ecmitser_out/
├── clean/            # Cleaned, oriented, and dereplicated sequences
├── itsx_out/         # Raw output files from ITSx
├── merged/           # The single merged input FASTA file
|
├── ecmitser.wits.fna # <-- 2. FASTA file of dereplicated sequences with an ITS region
└── ecmitser.wits.tsv # <-- 3. Table of sequences and their ITSx-detected positions
```

### Key Output Files

1. `ecmitser.wits.fna`: The final, filtered set of unique sequences that were confirmed to contain an ITS region. These are the sequences used as queries in the BLAST search.

2. `ecmitser.wits.tsv`: A simple 3-column table (`SeqID`, `ITSx_Start`, `ITSx_Stop`) for sequences that passed ITSx filtering.

---

## 🔬 Pipeline Details

The Snakefile executes the following steps:

1. Download Database: Downloads the Eukaryome ITS database (`release 2.0`).

2. Merge & Clean: All input *.fna files are concatenated. A custom script (`snakes/trim_and_filter.py`) is run to perform initial sequence cleaning.

3. Orient & Dereplicate: Sequences are oriented against the Eukaryome reference using vsearch `--orient` and then dereplicated using vsearch `--derep_fulllength`.

4. Extract ITS: `ITSx` is run on the unique sequences to find and extract fungal ITS regions.

5. Filter: The ITSx output is parsed. Only sequences where an ITS region was successfully detected ("Not found" is excluded) are kept for the next step. A custom script (`snakes/filter_by_its_content.py`) filters the FASTA file based on this list.

---

## 📄 License

This project is licensed under the MIT License.
