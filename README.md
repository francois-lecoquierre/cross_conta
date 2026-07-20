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

Samples are read directly from the multi VCF header.

## Requirements

- Python 3.7+
- `bcftools` (in PATH)
- `plotly` (optional, for HTML heatmap)
- `numpy` (optional, for VAF distribution plot)

## Usage

```bash
python conta_check_v9.py \
    --vcf   multi.vcf.gz \
    --bam-folder  /path/to/bams/ \
    --reference   /path/to/ref.fa \
    --output      batch_XXX_conta

# Regenerate plots from a previous results file
python conta_check_v9.py \
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
| `-j` / `--threads` | Parallel threads | 8 |
| `--reload-from` | Reload existing results TSV, regenerate plots only | — |
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
| `{prefix}.heatmap.html` | Interactive heatmap with hover details and VAF plot |

## Notes

- Up to 5 000 informative variants are used per pair (random downsampling if exceeded).
- Pairs with fewer than `--min-variants` informative variants are flagged as `INSUFFICIENT_DATA` (shown in orange in the heatmap).
- The heatmap embeds a VAF distribution plot for heterozygous SNVs per sample (requires plotly + numpy).
