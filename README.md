# merged_annotator

A Nextflow pipeline for genome and metagenome annotation. Combines structural annotation
(Bakta / Pyrodigal / Prodigal / Prokka), functional annotation (EggNOG-mapper, KofamScan,
dbCAN, VFDB), and optional ARG screening (AMRFinderPlus, fARGene, RGI, DeepARG, ABRicate).

---

## Quick start

```bash
nextflow run main.nf -profile standard --mode single --manifest samples.csv
```

The manifest is a CSV with at minimum two columns:

```csv
ID,assembly
sample1,/path/to/sample1.fasta
sample2,/path/to/sample2.fasta
```

---

## Modes

The `--mode` parameter controls which annotation tool is used and sets sensible defaults
for the rest of the pipeline. There are two named presets and four direct tool selectors.

### Named presets

#### `--mode single` (default)

For finished or high-quality draft isolate genomes.

- Annotator: **Pyrodigal** (`-p single`, fast ORF prediction)
- Protein clustering: off

```bash
nextflow run main.nf -profile standard --mode single --manifest samples.csv
```

#### `--mode meta`

For metagenome-assembled genomes (MAGs) or any set of assemblies to be treated as a
comparative/metagenomic dataset.

- Annotator: **Pyrodigal** (`-p meta`, fast ORF prediction)
- Protein clustering: on (MMseqs2 `easy-linclust` before functional annotation)

```bash
nextflow run main.nf -profile standard --mode meta --manifest samples.csv
```

### Direct tool selectors

Pass the tool name as `--mode` to bypass the presets entirely. Pipeline-level switches
(`cluster_proteome`) remain at their defaults (off) unless you also set them.

| `--mode`     | Annotator | `-p` flag   |
|--------------|-----------|-------------|
| `bakta`      | Bakta     | n/a         |
| `pyrodigal`  | Pyrodigal | `-p single` |
| `prodigal`   | Prodigal  | `-p single` |
| `prokka`     | Prokka    | n/a         |

```bash
# Prokka with ARG screening
nextflow run main.nf -profile standard --mode prokka --manifest samples.csv --arg_annotate

# Pyrodigal with protein clustering
nextflow run main.nf -profile standard --mode pyrodigal --manifest samples.csv --cluster_proteome
```

---

## Key parameters

| Parameter              | Default   | Description                                                  |
|------------------------|-----------|--------------------------------------------------------------|
| `--manifest`           | —         | Path to input CSV (required)                                 |
| `--outdir`             | `results` | Output directory                                             |
| `--mode`               | `single`  | Run mode (see above)                                         |
| `--arg_annotate`       | `false`   | Run ARG screening subworkflow                                |
| `--func_annotate`      | `false`   | Run Functional screening subworkflow                         |
| `--cluster_proteome`   | `false`   | Cluster proteins with MMseqs2 before functional annotation   |

### Database paths (pre-configured at Sanger)

| Parameter                | Default path                                                                  |
|--------------------------|-------------------------------------------------------------------------------|
| `--annotation_bakta_db`  | `/data/pam/software/bakta/v6.0/`                                               |
| `--arg_abricate_db`      | `/data/pam/software/abricate/db/`                                             |
| `--arg_amrfinderplus_db` | `/data/pam/software/amrfinder/latest/`                                        |
| `--eggnog_data_dir`      | `/data/pam/software/eggnog/v5.0/`                                             |
| `--dbcan_db`             | `/data/pam/software/run_dbcan/5.2.5/db/`                                      |
| `--vfdb_db`              | `/data/pam/team230/sd28/scratch/secondment_162/dbs/vfdb/VFDB_setB_pro.dmnd`   |
| `--kofam_db`             | `/data/pam/software/kofam_scan/`                                              |

### Key thresholds

| Parameter               | Default   | Description                          |
|-------------------------|-----------|--------------------------------------|
| `--module_completeness` | `0.5`     | KEGG module completeness threshold   |
| `--kofamscan_eval`      | `0.00001` | KofamScan e-value cutoff             |
| `--vfdb_identity`       | `50`      | VFDB minimum identity (%)            |
| `--vfdb_coverage`       | `80`      | VFDB minimum subject coverage (%)    |
| `--vfdb_e_value`        | `1e-10`   | VFDB e-value cutoff                  |
| `--cazyme_hmm_eval`     | `default` | dbCAN HMM e-value cutoff             |
| `--cazyme_hmm_cov`      | `default` | dbCAN HMM minimum coverage           |
| `--mmseqs_min_id`       | `0.8`     | MMseqs2 minimum sequence identity    |
| `--mmseqs_min_cov`      | `0.9`     | MMseqs2 minimum coverage             |

---

## Advanced usage

Tool parameters (translation table, HMM thresholds, per-tool ARG flags etc.)
are defined with sensible defaults in `conf/advanced.config`. Override them on the CLI as
needed — no need to edit the file directly.

To see a full list of these params use the --helpFull command

```bash
nextflow run main.nf --helpFull
```

### Running annotation tools in meta mode via direct tool selector

When using a direct tool selector (`--mode pyrodigal`, `--mode prodigal`, `--mode bakta`)
the tools default to single-genome mode. Use the corresponding flag to force meta mode:

| `--mode`    | Override flag                     | Effect    |
|-------------|-----------------------------------|-----------|
| `pyrodigal` | `--annotation_pyrodigal_meta`     | `-p meta` |
| `prodigal`  | `--annotation_prodigal_meta`      | `-p meta` |
| `bakta`     | `--annotation_bakta_meta`         | `--meta`  |

```bash
nextflow run main.nf -profile standard --mode pyrodigal --manifest samples.csv \
    --annotation_pyrodigal_meta

nextflow run main.nf -profile standard --mode bakta --manifest samples.csv \
    --annotation_bakta_meta
```

### Enabling/disabling individual ARG tools

When `--arg_annotate`, all five tools run by default. Skip individual ones:

```bash
nextflow run main.nf -profile standard --mode single --manifest samples.csv \
    --arg_annotate \
    --arg_skip_deeparg \
    --arg_skip_fargene
```

---

## Output structure

```text
results/
  annotation/
    pyrodigal/<sample_id>/       # or bakta/, prodigal/, prokka/
  func/
    kofamscan/
    dbcan/
    vfdb/
    eggnogmapper/
  arg/                       # only when --arg_annotate
    amrfinderplus/<sample_id>/
    abricate/<sample_id>/
    rgi/<sample_id>/
    deeparg/<sample_id>/
    fargene/<sample_id>/
    hamronization/
    argnorm/
  mmseqs/                    # only when --cluster_proteome
  reports/
    hamronization_summarize/
```
