# contaminate_bam.sh

Used to simulate a contamination of a bam_A by bam_B at a selected rate

Dependencies needed:
- picard.jar (used to downsample bam files)
- samtools (to merge the bam files)

Usage:
- script.sh -p <picard.jar file> -c <bam_contaminator> -d <bam_contaminated> -o <output_folder> -r <contamination_rate>


# conta_check.py

Detects cross-contamination between samples in a multi-sample VCF file.

## Principle

For each ordered sample pair (source → target), the script identifies variants that are:
- **Heterozygous (0/1)** in the source sample (potential contaminator)
- **Homozygous reference (0/0)** in the target sample (potentially contaminated)

The Variant Allele Frequency (VAF) of those variants is measured in the target BAM using `bcftools mpileup`. The estimated contamination is `mean_VAF × 2` (factor of 2 accounts for heterozygosity of the source variants). If this value exceeds the threshold, contamination is flagged.

Samples are read directly from the VCF header — no sample sheet required.

## Analysis pipeline (single genotyping pass)

Earlier versions re-measured the VAF of the same positions once per sample pair, which meant genotyping each target's BAM up to N−1 times. `conta_check_v11.py` avoids that redundant work with a 3-pass design:

1. **Genotype matrix** — the VCF is read once (`bcftools query`) to build a full genotype matrix for every sample.
2. **Informative variants** — that matrix is scanned in-memory to determine, for every sample pair, which variants are informative, and for every sample (as target), the full set of positions that need VAF measurement.
3. **Single genotyping pass** — each sample's BAM is genotyped exactly once (one batched `bcftools mpileup` call per sample, parallelized with `--threads`), and the resulting VAFs are reused across every sample-pair comparison that needs them.

## Requirements

- Python 3.7+
- `bcftools` (in PATH)
- `plotly` (optional, for the HTML contamination report)
- `numpy` (optional, for the VAF distribution plot)
- BAM files must already be indexed (`.bai`)

## Usage

```bash
python conta_check_v11.py \
    --vcf   samples.vcf.gz \
    --bam-folder  /path/to/bams/ \
    --reference   /path/to/ref.fa \
    --output      batch_XXX_conta

# Regenerate the report from a previous results file
python conta_check_v11.py \
    --reload-from batch_XXX_conta.results.tsv \
    --output      batch_XXX_conta_replot
```

## Arguments

| Argument | Description | Default |
|---|---|---|
| `--vcf` / `-v` | Multi-sample VCF (.vcf or .vcf.gz) | required |
| `--bam-folder` / `-b` | Folder containing indexed BAM files | required |
| `--reference` / `-r` | Reference genome FASTA (indexed) | required |
| `--output` / `-o` | Output file prefix | required |
| `-t` / `--threshold` | Contamination threshold | 0.01 (1%) |
| `-m` / `--min-variants` | Minimum informative variants required | 500 |
| `-j` / `--threads` | Parallel threads (used for genotyping) | 8 |
| `--reload-from` | Reload existing results TSV, regenerate the report only | — |
| `--vaf-plot-only` | Generate only the VAF distribution plot (diagnostic) | — |

## BAM file naming

The script looks for BAM files in `--bam-folder` matching the sample name:
1. Exact match: `{sample}.sorted.grouped.dedup.BQSR.bam`
2. Any `.bam` file whose name contains the sample ID (must be unique and indexed)

## Output files

| File | Description |
|---|---|
| `{prefix}.results.tsv` | Per-pair results (VAF, estimated contamination, status) |
| `{prefix}.matrix.estimated_conta.tsv` | Square matrix of estimated contamination values |
| `{prefix}.heatmap.html` | Interactive **Contamination Report** (see below) |

### Contamination Report layout

The `{prefix}.heatmap.html` file is a single-page report with three vertical sections:

1. **Header** — run summary: number of samples analyzed, threshold, min variants, number of contaminations detected, number of insufficient-data pairs, generation timestamp, and the analysis parameters used (VCF, BAM folder, reference, threads — or the reloaded results file when using `--reload-from`).
2. **Cross-contamination matrix** — the heatmap (source × target) with its legend on the right. The matrix takes about 2/3 of the page width, the legend the remaining third. Cell size and font sizes automatically shrink for batches with many samples, so the matrix never overflows the page width.
3. **VAF distribution** — density curves of heterozygous SNV VAFs per sample (requires `plotly` + `numpy`, skipped in `--reload-from` mode since it needs the VCF).

## Notes

- Up to 5 000 informative variants are used per pair (random downsampling if exceeded) — the random sample is drawn from the whole genome-wide list, not just the first positions.
- Pairs with fewer than `--min-variants` informative variants are flagged as `INSUFFICIENT_DATA` (shown in orange in the report).

