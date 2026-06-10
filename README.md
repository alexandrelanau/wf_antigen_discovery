# wf_antigen_discovery

*Personalized neoantigen discovery from immunopeptidomics data*

## What does this workflow do?

This pipeline takes a set of RNA sequences and immunopeptidomics mass spectrometry data, and identifies peptides from  the query RNA .fasta translation are actually presented on that patient's MHC molecules. In plain terms: it helps discover tumor-specific antigens that could be targeted by the immune system.

The main steps are:

- **Translation:** Query RNA .fasta is translated into all possible protein sequences using ORFfinder.
- **Database construction:** The translation result is concatenated with a reference human proteome (UniProt) to build a custom search database.
- **Mass spec search:** Immunopeptidomics data (HLA pulldown + mass spectrometry) is searched against the database using FragPipe/MSFragger.
- **MHC binding prediction:** Identified peptides are scored for binding to the sample specific HLA alleles using NetMHCpan.
- **Postprocessing:** Results are filtered, annotated against the HLA Ligand Atlas, and summarized into a final report.

## Requirements

- Snakemake >= 7
- Apptainer (Singularity) — for containerized execution
- All tool containers are pre-built and located in `containers/image/`

## Input files you need to provide

- **Set of query RNA sequences:** RNA .fasta
- **HLA alleles file:** A plain text file with the patient's HLA alleles, formatted for NetMHCpan (e.g. `HLA-A01:01,HLA-B07:02`).
- **Immunopeptidomics data:** A directory containing `.raw` mass spectrometry files from an HLA pulldown experiment.
- **Reference proteome:** A reviewed UniProt FASTA with PE tags (path set in config).

## Quick start

1. Clone the repository and navigate to it.
2. Build singularity images in containers
   ```
   singulariy build image/fragpipe.22.0.ubuntu.sif  def/fragpipe.22.0.ubuntu.def
   singularity build image/netmhcpan.4.1.ubuntu.sif  def/netmhcpan.4.1.ubuntu.def
   singularity build image/python.3.12.2.debian.sif  def/python.3.12.2.debian.def
   ```
3. Download Fragpipe dependencies at:
   http://msfragger-upgrader.nesvilab.org/upgrader/
   https://msfragger.arsci.com/ionquant/
   https://msfragger-upgrader.nesvilab.org/diatracer/
   
   or run containers/def/fragpipe22.0.dependencies.bash
5. Copy an existing profile folder (e.g. `workflow/profiles/mariia`) and rename it.
6. Edit `config.yaml` in your new profile to set your resource limits and file paths.
7. Edit `todo.py` in your new profile to define your sample(s).
8. Run the workflow:
   ```
   snakemake --profile workflow/profiles/<your_profile>
   ```
### Example

Example data is provided for you to test the pipeline at example. You need to download spectral  files from the link in example/readme.txt and run the workflow
   ```
   snakemake --profile workflow/profiles/mariia
   ```

## Output

All results are written to `_data/antigene_search/<project_name>/`. Key output files:

| File | Description |
|------|-------------|
| `fragpipe_postprocessing/Summary.*.csv` | Final table of non-canonical peptides with MHC binding predictions and HLA Ligand Atlas annotation |
| `fragpipe_postprocessing/*.NetMHCpan_out.xls` | Raw NetMHCpan binding predictions |
| `fragpipe_workdir.*/*/psm.tsv` | FragPipe peptide-spectrum matches |
| `final_msfragger_database.*.fasta` | The custom search database used for this run |

## Repository structure

```
wf_antigen_discovery/
├── workflow/
│   ├── Snakefile          # Main workflow
│   ├── utils.py           # Helper functions (e.g. get_size_mb)
│   ├── shell -> ../shell  # Symlink to shell scripts and tools
│   └── profiles/          # Per-user config and sample definitions
│       └── <profile>/
│           ├── config.yaml
│           └── todo.py
├── shell/                 # Python scripts, ORFfinder, FragPipe tools
├── containers/
│   ├── def/               # Apptainer definition files
│   └── image/             # Pre-built .sif container images
└── _data/                 # All outputs (created at runtime)
```
