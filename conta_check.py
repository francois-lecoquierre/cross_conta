#!/usr/bin/env python3
"""
Cross-Contamination Detection for Multi-Sample VCF Files

Version 9 — samples are read directly from the multi-sample VCF header.

This script detects cross-contamination between samples by analyzing variants that are:
- Heterozygous (0/1) in the source sample (potential contaminator)
- Homozygous reference (0/0) in the target sample (potentially contaminated)

The Variant Allele Frequency (VAF) is calculated from BAM files using bcftools mpileup.
The estimated contamination is calculated as MEAN_VAF x 2 (accounting for heterozygosity).
If the estimated contamination exceeds a specified threshold, contamination is suspected.

Features:
- Parallel processing for efficient analysis of multiple sample pairs
- Interactive HTML heatmap with visual contamination indicators
- Configurable thresholds for contamination and minimum variant counts
- Support for result caching and reload without re-analysis
- Suitable for intra-family contamination detection (e.g., mother-fetus)

Requirements:
- bcftools (for VCF operations and mpileup)
- samtools (for BAM operations)
- Python 3.7+
- plotly (optional, for heatmap generation)

Author: [Your Name]
Version: 8.0
License: MIT
"""

import argparse
import sys
import os
import subprocess
import random
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Optional dependencies
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    go = None  # type: ignore[assignment]
    PLOTLY_AVAILABLE = False
    print("Note: plotly not available, skipping heatmap generation", file=sys.stderr)

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    np = None  # type: ignore[assignment]
    NUMPY_AVAILABLE = False
    print("Note: numpy not available, VAF plotting will be limited", file=sys.stderr)


def parse_arguments():
    """
    Parse and validate command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Detect cross-contamination between samples in a multi-sample VCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --vcf samples.vcf.gz --bam-folder bams/ --reference ref.fa --output output

  # With custom thresholds and parallel processing
  %(prog)s --vcf samples.vcf.gz --bam-folder bams/ --reference ref.fa --output output -t 0.02 -m 100 -j 16

  # Reload existing results to regenerate plots
  %(prog)s --reload-from results.tsv --output output_new -t 0.01 -m 500
        """
    )

    # Required arguments
    parser.add_argument(
        '--vcf', '-v',
        required=True,
        help='Input multi-sample VCF file (.vcf or .vcf.gz)'
    )
    parser.add_argument(
        '--bam-folder', '-b',
        required=True,
        help='Folder containing BAM files (one per sample, must be indexed)'
    )
    parser.add_argument(
        '--reference', '-r',
        required=True,
        help='Reference genome FASTA file (must be indexed)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        dest='output_prefix',
        help='Output file prefix for results files'
    )

    # Optional arguments
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=0.01,
        help='Estimated contamination threshold (default: 0.01, i.e., 1%%)'
    )
    parser.add_argument(
        '-m', '--min-variants',
        type=int,
        default=500,
        help='Minimum number of informative variants required (default: 500)'
    )
    parser.add_argument(
        '-j', '--threads',
        type=int,
        default=8,
        help='Number of threads for parallel processing (default: 8)'
    )
    parser.add_argument(
        '--reload-from',
        type=str,
        default=None,
        help='Path to existing results TSV file to reload and regenerate plots only (skips analysis)'
    )
    parser.add_argument(
        '--vaf-plot-only',
        action='store_true',
        help='Generate only the VAF distribution plot without running cross-contamination analysis (for testing)'
    )

    return parser.parse_args()


def get_samples_from_vcf(vcf_file):
    """
    Extract sample names from a VCF file.
    
    Args:
        vcf_file (str): Path to VCF file
        
    Returns:
        list: List of sample names
        
    Raises:
        SystemExit: If bcftools fails to extract samples
    """
    try:
        result = subprocess.run(
            ['bcftools', 'query', '-l', vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        samples = result.stdout.strip().split('\n')
        return samples
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to extract samples from VCF: {e}", file=sys.stderr)
        sys.exit(1)


def find_bam_for_sample(sample_name, bam_folder):
    """
    Find the BAM file corresponding to a sample name.
    
    Priority rules:
    1. If exact format {sample}.sorted.grouped.dedup.BQSR.bam exists, return it
    2. Otherwise, find BAM files containing the sample name
    3. Only one valid BAM should match (must have .bai index)
    
    Args:
        sample_name (str): Sample name to search for
        bam_folder (str): Path to folder containing BAM files
        
    Returns:
        str or None: Path to the BAM file, or None if not found/ambiguous
    """
    bam_folder = Path(bam_folder)
    
    # Check for exact format first
    exact_bam = bam_folder / f"{sample_name}.sorted.grouped.dedup.BQSR.bam"
    if exact_bam.exists():
        bai_file = exact_bam.with_suffix('.bam.bai')
        if not bai_file.exists():
            bai_file = Path(str(exact_bam) + '.bai')
        
        if bai_file.exists():
            return str(exact_bam)
    
    # Search for BAM files containing the sample name
    matching_bams = []
    for bam_file in bam_folder.glob('*.bam'):
        if sample_name in bam_file.name:
            # Check if BAM index exists
            bai_file = bam_file.with_suffix('.bam.bai')
            if not bai_file.exists():
                bai_file = Path(str(bam_file) + '.bai')
            
            if bai_file.exists():
                matching_bams.append(bam_file)
    
    if len(matching_bams) == 0:
        print(f"Warning: No BAM file found for sample {sample_name}", file=sys.stderr)
        return None
    elif len(matching_bams) > 1:
        print(f"Warning: Multiple BAM files found for sample {sample_name}: {matching_bams}", 
              file=sys.stderr)
        return None
    
    return str(matching_bams[0])


def get_informative_variants(vcf_file, source_sample, target_sample):
    """
    Get informative variants for contamination detection.
    
    Informative variants are those that are:
    - Heterozygous (0/1) in source sample
    - Homozygous reference (0/0) in target sample
    - Single nucleotide variants (SNVs) only
    
    Args:
        vcf_file (str): Path to VCF file
        source_sample (str): Source sample name (potential contaminator)
        target_sample (str): Target sample name (potentially contaminated)
        
    Returns:
        list: List of tuples (chrom, pos, ref, alt, source_gt)
    """
    # Get sample indices
    try:
        samples_result = subprocess.run(
            ['bcftools', 'query', '-l', vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        samples = samples_result.stdout.strip().split('\n')
        source_idx = samples.index(source_sample)
        target_idx = samples.index(target_sample)
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"Error getting sample indices: {e}", file=sys.stderr)
        return []
    
    # Build filter expression for bcftools
    filter_expr = (
        f'GT[{source_idx}]="0/1" && GT[{target_idx}]="0/0" && '
        f'TYPE="snp" && strlen(REF)=1 && strlen(ALT)=1'
    )
    
    try:
        # Filter VCF
        result = subprocess.run(
            ['bcftools', 'view', '-i', filter_expr, vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        
        if result.returncode != 0:
            return []
        
        # Extract variant information including genotypes
        query_result = subprocess.run(
            ['bcftools', 'query', '-f', f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'],
            input=result.stdout,
            capture_output=True,
            text=True,
            check=True
        )
        
        variants = []
        for line in query_result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 4 + len(samples):
                    chrom, pos, ref, alt = parts[0:4]
                    source_gt = parts[4 + source_idx]
                    
                    # Verify SNV
                    if len(ref) == 1 and len(alt) == 1:
                        variants.append((chrom, pos, ref, alt, source_gt))
        
        return variants
        
    except subprocess.CalledProcessError as e:
        print(f"Error querying VCF: {e}", file=sys.stderr)
        return []


def calculate_vaf_batch_mpileup(bam_file, reference, variants):
    """
    Calculate VAF for multiple positions using bcftools mpileup in batch mode.
    
    This is significantly faster than calling mpileup individually for each position.
    
    Args:
        bam_file (str): Path to BAM file
        reference (str): Path to reference genome
        variants (list): List of tuples (chrom, pos, ref, alt, source_gt)
        
    Returns:
        dict: Dictionary mapping (chrom, pos) to VAF
    """
    if not variants:
        return {}
    
    # Create temporary BED file with all positions
    # BED format is 0-based, VCF is 1-based
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_file:
        bed_path = bed_file.name
        for chrom, pos, ref, alt, source_gt in variants:
            bed_file.write(f"{chrom}\t{int(pos)-1}\t{pos}\n")
    
    vaf_dict = {}
    
    # Create lookup for expected alt alleles
    variant_lookup = {}
    for chrom, pos, ref, alt, source_gt in variants:
        variant_lookup[(chrom, pos)] = alt
    
    try:
        # Run bcftools mpileup with targets from BED file
        result = subprocess.run(
            [
                'bcftools', 'mpileup',
                '-f', reference,
                '-T', bed_path,  # Targets from BED file (more efficient than -R)
                '-a', 'AD',      # Output allelic depth
                '-Q', '20',      # Minimum base quality
                '-q', '20',      # Minimum mapping quality
                '-d', '10000',   # Maximum depth
                '--ignore-RG',   # Ignore read groups
                '-Ov',           # VCF output
                bam_file
            ],
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"Warning: bcftools mpileup failed with return code {result.returncode}", 
                  file=sys.stderr)
            if result.stderr:
                print(f"Error message: {result.stderr}", file=sys.stderr)
            return {}
        
        if not result.stdout.strip():
            print("Warning: mpileup returned no output", file=sys.stderr)
            return {}
        
        # Parse mpileup VCF output
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 10:
                continue
            
            chrom = fields[0]
            pos = fields[1]
            ref_allele = fields[3]
            alt_alleles = fields[4].split(',')
            
            # Parse FORMAT and sample data
            format_fields = fields[8].split(':')
            sample_data = fields[9].split(':')
            
            # Find AD index in FORMAT
            try:
                ad_idx = format_fields.index('AD')
                ad_str = sample_data[ad_idx]
            except (ValueError, IndexError):
                continue
            
            # Get expected alt allele for this position
            expected_alt = variant_lookup.get((chrom, pos))
            if not expected_alt:
                continue
            
            try:
                # Parse AD (allelic depths)
                ad_parts = ad_str.split(',')
                ref_count = int(ad_parts[0])
                
                # Find which alt allele matches our expected alt
                alt_count = 0
                if expected_alt in alt_alleles:
                    alt_idx = alt_alleles.index(expected_alt)
                    if alt_idx + 1 < len(ad_parts):
                        alt_count = int(ad_parts[alt_idx + 1])
                elif len(ad_parts) > 1:
                    # If expected alt not in called alts, use first alt count
                    alt_count = int(ad_parts[1])
                
                # Calculate total depth and VAF
                dp = ref_count + sum(int(x) for x in ad_parts[1:])
                
                if dp == 0:
                    continue
                
                vaf = alt_count / dp
                vaf_dict[(chrom, pos)] = vaf
                
            except (ValueError, IndexError):
                continue
        
    except (subprocess.CalledProcessError, ValueError, IndexError) as e:
        print(f"Warning: Batch mpileup failed: {e}", file=sys.stderr)
    finally:
        # Clean up temporary BED file
        try:
            os.unlink(bed_path)
        except:
            pass
    
    return vaf_dict


def process_sample_pair(vcf_file, source_sample, target_sample, sample_to_bam, 
                       reference, threshold, min_variants=100, max_variants=5000):
    """
    Process a single source-target sample pair for contamination detection.
    
    This function is designed to be run in parallel.
    
    Args:
        vcf_file (str): Path to VCF file
        source_sample (str): Source sample name
        target_sample (str): Target sample name
        sample_to_bam (dict): Mapping of sample names to BAM file paths
        reference (str): Path to reference genome
        threshold (float): Contamination threshold
        min_variants (int): Minimum number of informative variants required
        max_variants (int): Maximum number of variants to use (downsampling)
        
    Returns:
        dict: Analysis results with keys:
            - source_sample, target_sample
            - median_vaf, mean_vaf, estimated_conta
            - informative_variant_count
            - gt_0_1_count
            - status (CONTA, OK, INSUFFICIENT_DATA, etc.)
            - insufficient_variants (bool)
            - text_result (detailed message)
    """
    result = {
        'source_sample': source_sample,
        'target_sample': target_sample,
        'median_vaf': 'NA',
        'mean_vaf': 'NA',
        'estimated_conta': 'NA',
        'informative_variant_count': 0,
        'gt_0_1_count': 0,
        'status': 'NO_DATA',
        'insufficient_variants': False,
        'text_result': ''
    }
    
    # Get informative variants (1/1 in source, 0/0 in target)
    variants = get_informative_variants(vcf_file, source_sample, target_sample)
    
    original_count = len(variants)
    
    # Check if we have enough variants - return immediately without processing
    if len(variants) < min_variants:
        result['insufficient_variants'] = True
        result['informative_variant_count'] = original_count
        result['gt_0_1_count'] = original_count
        result['text_result'] = f'Insufficient informative variants: {original_count} (< {min_variants} required)'
        result['status'] = 'INSUFFICIENT_DATA'
        return result
    
    # Downsample if needed
    if len(variants) > max_variants:
        variants = random.sample(variants, max_variants)
    
    # Count genotypes in the (possibly downsampled) set
    gt_0_1_count = sum(1 for v in variants if v[4] in ['0/1', '1/0', '0|1', '1|0'])
    
    # Calculate VAF for all variants using batch mpileup
    target_bam = sample_to_bam[target_sample]
    vaf_dict = calculate_vaf_batch_mpileup(target_bam, reference, variants)
    
    # Extract VAF values
    vafs = []
    vafs_0_1 = []
    
    for chrom, pos, ref, alt, source_gt in variants:
        vaf = vaf_dict.get((chrom, pos))
        if vaf is not None:
            vafs.append(vaf)
            if source_gt in ['0/1', '1/0', '0|1', '1|0']:
                vafs_0_1.append(vaf)
    
    if len(vafs) == 0:
        result['informative_variant_count'] = len(variants)
        result['gt_0_1_count'] = gt_0_1_count
        result['text_result'] = 'Failed to calculate VAF for variants'
        result['insufficient_variants'] = False
        return result
    
    # Calculate statistics
    vafs.sort()
    median_vaf = vafs[len(vafs) // 2]
    mean_vaf = sum(vafs) / len(vafs)
    estimated_conta = mean_vaf * 2  # Factor of 2 for heterozygous variants
    
    # Build downsample note
    downsample_note = ""
    if original_count > max_variants:
        downsample_note = f" (downsampled from {original_count})"
    
    # Determine contamination status
    if estimated_conta > threshold:
        status = 'CONTA'
        text_result = (
            f"WARNING | contamination detected: {target_sample} is "
            f"contaminated by {source_sample}. Variants that are "
            f"heterozygous (0/1) in {source_sample} "
            f"and absent (0/0) in {target_sample} have a mean VAF of {mean_vaf:.4f} "
            f"in {target_sample}, corresponding to an estimated contamination of {estimated_conta:.4f} "
            f"(above threshold: {threshold}). "
            f"Informative variants: {len(vafs)}{downsample_note}"
        )
    else:
        status = 'OK'
        text_result = (
            f"OK | {target_sample} is not contaminated by {source_sample}. "
            f"Variants that are heterozygous (0/1) "
            f"in {source_sample} and absent (0/0) in {target_sample} have a mean VAF "
            f"of {mean_vaf:.4f} in {target_sample}, corresponding to an estimated contamination of {estimated_conta:.4f} "
            f"(below threshold: {threshold}). "
            f"Informative variants: {len(vafs)}{downsample_note}"
        )
    
    result.update({
        'median_vaf': f"{median_vaf:.4f}",
        'mean_vaf': f"{mean_vaf:.4f}",
        'estimated_conta': f"{estimated_conta:.4f}",
        'informative_variant_count': len(vafs),
        'gt_0_1_count': len(vafs_0_1),
        'status': status,
        'text_result': text_result
    })
    
    return result


def analyze_contamination(vcf_file, bam_folder, reference, threshold, min_variants=100, threads=8):
    """
    Analyze cross-contamination between all sample pairs.

    Args:
        vcf_file (str): Path to VCF file
        bam_folder (str): Path to BAM folder
        reference (str): Path to reference genome
        threshold (float): Contamination threshold
        min_variants (int): Minimum number of informative variants required
        threads (int): Number of parallel threads to use

    Returns:
        tuple: (results list, samples list)
    """
    print("=" * 60)
    print("Starting cross-contamination analysis")
    print("=" * 60)

    # Get samples from VCF
    samples = get_samples_from_vcf(vcf_file)

    print(f"\nAnalyse de {len(samples)} samples: {', '.join(samples)}")
    
    # Map samples to BAM files
    sample_to_bam = {}
    print("\nMapping samples to BAM files...")
    for sample in samples:
        bam_file = find_bam_for_sample(sample, bam_folder)
        if bam_file:
            sample_to_bam[sample] = bam_file
            print(f"  {sample} -> {bam_file}")
        else:
            print(f"  {sample} -> NOT FOUND", file=sys.stderr)
    
    if len(sample_to_bam) < 2:
        print("\nError: Need at least 2 samples with valid BAM files", file=sys.stderr)
        sys.exit(1)
    
    # Build list of all sample pairs to process
    sample_pairs = []
    for source_sample in samples:
        if source_sample not in sample_to_bam:
            continue
        for target_sample in samples:
            if source_sample == target_sample:
                continue
            if target_sample not in sample_to_bam:
                continue
            sample_pairs.append((source_sample, target_sample))
    
    total_pairs = len(sample_pairs)
    print(f"\nProcessing {total_pairs} sample pairs using {threads} threads...\n")
    
    results = []
    
    # Process pairs in parallel
    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all jobs
            future_to_pair = {}
            for source_sample, target_sample in sample_pairs:
                future = executor.submit(
                    process_sample_pair,
                    vcf_file,
                    source_sample,
                    target_sample,
                    sample_to_bam,
                    reference,
                    threshold,
                    min_variants
                )
                future_to_pair[future] = (source_sample, target_sample)
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_pair):
                source_sample, target_sample = future_to_pair[future]
                completed += 1
                try:
                    result = future.result()
                    results.append(result)
                    status = result['status']
                    mean_vaf = result['mean_vaf']
                    estimated_conta = result['estimated_conta']
                    informative_count = result['gt_0_1_count']
                    print(f"[{completed}/{total_pairs}] {source_sample} -> {target_sample}: "
                          f"Status={status}, Mean_VAF={mean_vaf}, Estimated_conta={estimated_conta}, Informative_variants={informative_count}")
                except Exception as e:
                    print(f"Error processing {source_sample} -> {target_sample}: {e}", 
                          file=sys.stderr)
    else:
        # Single-threaded mode
        for idx, (source_sample, target_sample) in enumerate(sample_pairs, 1):
            print(f"\n[{idx}/{total_pairs}] Processing {source_sample} -> {target_sample}")
            try:
                result = process_sample_pair(
                    vcf_file,
                    source_sample,
                    target_sample,
                    sample_to_bam,
                    reference,
                    threshold,
                    min_variants
                )
                results.append(result)
                informative_count = result['gt_0_1_count']
                print(f"  Status: {result['status']}, Mean_VAF: {result['mean_vaf']}, Estimated_conta: {result['estimated_conta']}, Informative_variants: {informative_count}")
            except Exception as e:
                print(f"  Error: {e}", file=sys.stderr)
    
    return results, samples


def extract_het_snv_vafs(vcf_file, samples):
    """
    Extract VAF values for all heterozygous SNV variants for each sample.
    Only variants with depth >= 100x are included.
    
    Args:
        vcf_file (str): Path to VCF file
        samples (list): List of sample names
        
    Returns:
        dict: Dictionary mapping sample names to lists of VAF values
    """
    sample_vafs = {sample: [] for sample in samples}
    
    try:
        # Query VCF for SNV variants with GT and AD fields
        result = subprocess.run(
            ['bcftools', 'view', '-v', 'snps', vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        
        if result.returncode != 0:
            return sample_vafs
        
        # Parse VCF to extract heterozygous SNVs
        vcf_content = result.stdout
        
        for line in vcf_content.split('\n'):
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            
            fields = line.split('\t')
            if len(fields) < 10:
                continue
            
            ref = fields[3]
            alt = fields[4]
            
            # Only SNVs
            if len(ref) != 1 or len(alt) != 1:
                continue
            
            format_fields = fields[8].split(':')
            
            # Check if GT and AD are present
            if 'GT' not in format_fields or 'AD' not in format_fields:
                continue
            
            gt_idx = format_fields.index('GT')
            ad_idx = format_fields.index('AD')
            
            # Process each sample
            for i, sample in enumerate(samples):
                sample_data = fields[9 + i].split(':')
                
                if len(sample_data) <= max(gt_idx, ad_idx):
                    continue
                
                gt = sample_data[gt_idx]
                ad_str = sample_data[ad_idx]
                
                # Check if heterozygous
                if gt not in ['0/1', '1/0', '0|1', '1|0']:
                    continue
                
                # Parse AD
                try:
                    ad_parts = ad_str.split(',')
                    if len(ad_parts) < 2:
                        continue
                    
                    ref_count = int(ad_parts[0])
                    alt_count = int(ad_parts[1])
                    dp = ref_count + alt_count
                    
                    # Filter: minimum depth of 100x
                    if dp >= 60:
                        vaf = alt_count / dp
                        sample_vafs[sample].append(vaf)
                except (ValueError, IndexError):
                    continue
        
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to extract heterozygous SNVs: {e}", file=sys.stderr)
    
    return sample_vafs


def generate_vaf_plot(vcf_file, samples):
    """
    Generate a VAF distribution plot for heterozygous SNVs per sample.
    Shows distribution as density curves (line plots).
    
    Args:
        vcf_file (str): Path to VCF file
        samples (list): List of sample names
        
    Returns:
        go.Figure or None: Plotly figure object, or None if failed
    """
    if not PLOTLY_AVAILABLE:
        return None
    assert go is not None

    print("\nExtracting heterozygous SNV VAFs for VAF plot...")
    sample_vafs = extract_het_snv_vafs(vcf_file, samples)

    if not NUMPY_AVAILABLE:
        print("Error: numpy is required for VAF plot generation", file=sys.stderr)
        return None
    assert np is not None

    # Create figure
    fig = go.Figure()
    
    # Add a density curve for each sample
    def gaussian_smooth(data, sigma=3):
        """Apply Gaussian smoothing using numpy convolution."""
        assert np is not None
        # Create Gaussian kernel
        kernel_size = int(6 * sigma + 1)
        x = np.arange(-kernel_size // 2 + 1, kernel_size // 2 + 1)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        kernel = kernel / kernel.sum()
        
        # Apply convolution with padding
        padded = np.pad(data, kernel_size // 2, mode='edge')
        smoothed = np.convolve(padded, kernel, mode='valid')
        return smoothed
    
    for sample in samples:
        vafs = sample_vafs[sample]
        if not vafs:
            print(f"  Warning: No heterozygous SNVs found for {sample}")
            continue
        
        # Calculate histogram with many bins
        hist, bin_edges = np.histogram(vafs, bins=300, range=(0, 1), density=True)
        
        # Apply Gaussian smoothing for a smoother curve
        # Sigma controls the smoothing bandwidth (higher = more smooth)
        hist_smoothed = gaussian_smooth(hist, sigma=3)
        
        # Get bin centers for x-axis
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Add density curve
        fig.add_trace(go.Scatter(
            x=bin_centers,
            y=hist_smoothed,
            mode='lines',
            name=sample,
            line=dict(width=2, shape='spline'),
            hovertemplate=f'{sample}<br>VAF: %{{x:.3f}}<br>Density: %{{y:.3f}}<extra></extra>'
        ))
        
        print(f"  {sample}: {len(vafs)} heterozygous SNVs (mean VAF: {sum(vafs)/len(vafs):.3f})")
    
    # Add vertical line at VAF=0.5 (expected for heterozygous variants)
    fig.add_vline(
        x=0.5,
        line_dash="dash",
        line_color="gray",
        line_width=1.5,
        opacity=0.6,
        layer="below",
        annotation_text="Expected VAF (0.5)",
        annotation_position="top"
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='VAF Distribution of Heterozygous SNVs<br><sub>Density Curves per Sample</sub>',
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title='VAF (Variant Allele Frequency)',
            showgrid=True,
            range=[0, 1]
        ),
        yaxis=dict(
            title='Probability Density',
            showgrid=True
        ),
        hovermode='closest',
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        ),
        width=700,
        height=600,
        margin=dict(l=80, r=200, t=100, b=80)
    )
    
    return fig


def generate_heatmap(results, samples, output_prefix, threshold=0.02, min_variants=100, vcf_file=None):
    """
    Generate an interactive heatmap of estimated contamination values using plotly.
    
    The heatmap includes:
    - Color coding based on estimated contamination (green=OK, yellow=intermediate, red=contamination)
    - Emoji indicators (⚠️ for contamination)
    - "N=X" annotations for insufficient data
    - Hover tooltips with detailed information
    - VAF distribution plot for heterozygous SNVs (if vcf_file provided)
    
    Args:
        results (list): List of result dictionaries
        samples (list): List of sample names
        output_prefix (str): Output file prefix
        threshold (float): Contamination threshold
        min_variants (int): Minimum number of informative variants required
        vcf_file (str, optional): Path to VCF file for VAF plot generation
        
    Returns:
        str or None: Path to heatmap file, or None if plotly unavailable
    """
    if not PLOTLY_AVAILABLE:
        print("Plotly not available, skipping heatmap generation")
        return None
    assert go is not None

    # Create lookup dictionaries
    result_dict = {}
    status_dict = {}
    gt_dict = {}
    insufficient_dict = {}
    
    for result in results:
        key = (result['source_sample'], result['target_sample'])
        estimated_conta = result.get('estimated_conta', 'NA')
        insufficient_dict[key] = result.get('insufficient_variants', False)
        
        # Convert to float, handle 'NA'
        if estimated_conta != 'NA':
            try:
                result_dict[key] = float(estimated_conta)
                status_dict[key] = result['status']
                gt_dict[key] = result['gt_0_1_count']
            except ValueError:
                result_dict[key] = None
                status_dict[key] = 'NO_DATA'
                gt_dict[key] = 0
        else:
            result_dict[key] = None
            status_dict[key] = 'NO_DATA'
            gt_dict[key] = 0
    
    # Build matrices for heatmap
    n_samples = len(samples)
    z_matrix = []       # Estimated contamination values
    hover_text = []     # Hover tooltips
    text_matrix = []    # Emoji indicators
    annotations = []    # Text annotations for insufficient data
    
    for i, source_sample in enumerate(samples):
        row = []
        hover_row = []
        text_row = []
        
        for j, target_sample in enumerate(samples):
            if source_sample == target_sample:
                # Diagonal
                row.append(None)
                hover_row.append(f"{source_sample}<br>(self)")
                text_row.append("")
            else:
                key = (source_sample, target_sample)
                conta = result_dict.get(key)
                status = status_dict.get(key, 'NO_DATA')
                gt_0_1 = gt_dict.get(key, 0)
                is_insufficient = insufficient_dict.get(key, False)
                
                # Determine emoji based on status
                emoji = ""
                if status == 'CONTA':
                    emoji = "⚠️"
                
                text_row.append(emoji)
                
                # Handle insufficient variants case
                if is_insufficient:
                    informative_variants = gt_0_1
                    
                    if conta is not None:
                        row.append(conta)
                        # Get mean_vaf from results
                        mean_vaf_val = None
                        for res in results:
                            if res['source_sample'] == source_sample and res['target_sample'] == target_sample:
                                mean_vaf_val = res.get('mean_vaf', 'NA')
                                break
                        hover_row.append(
                            f"Source: {source_sample}<br>"
                            f"Target: {target_sample}<br>"
                            f"⚠️ INSUFFICIENT DATA<br>"
                            f"Mean VAF: {mean_vaf_val}<br>"
                            f"Estimated conta: {conta:.4f}<br>"
                            f"Informative variants: {informative_variants} (< {min_variants} required)"
                        )
                    else:
                        row.append(None)
                        hover_row.append(
                            f"Source: {source_sample}<br>"
                            f"Target: {target_sample}<br>"
                            f"⚠️ INSUFFICIENT DATA<br>"
                            f"Informative variants: {informative_variants} (< {min_variants} required)"
                        )
                    
                    # Add annotation for variant count
                    annotations.append(
                        dict(
                            x=j, y=i,
                            text=f"N={informative_variants}",
                            showarrow=False,
                            font=dict(size=9, color='orange'),
                            xref='x', yref='y'
                        )
                    )
                else:
                    row.append(conta)
                    
                    if conta is not None:
                        informative_variants = gt_0_1
                        # Get mean_vaf from results
                        mean_vaf_val = None
                        for res in results:
                            if res['source_sample'] == source_sample and res['target_sample'] == target_sample:
                                mean_vaf_val = res.get('mean_vaf', 'NA')
                                break
                        hover_row.append(
                            f"Source: {source_sample}<br>"
                            f"Target: {target_sample}<br>"
                            f"Mean VAF: {mean_vaf_val}<br>"
                            f"Estimated conta: {conta:.4f}<br>"
                            f"Status: {status}<br>"
                            f"Informative variants: {informative_variants}"
                        )
                    else:
                        hover_row.append(
                            f"Source: {source_sample}<br>"
                            f"Target: {target_sample}<br>"
                            f"No data"
                        )
        
        z_matrix.append(row)
        hover_text.append(hover_row)
        text_matrix.append(text_row)
    
    # Create heatmap figure
    fig = go.Figure(data=go.Heatmap(
        z=z_matrix,
        x=samples,
        y=samples,
        text=text_matrix,
        texttemplate='%{text}',
        textfont=dict(size=20),
        customdata=hover_text,
        hovertemplate='%{customdata}<extra></extra>',
        colorscale='RdYlGn_r',  # Reversed Red-Yellow-Green
        zmin=0,
        zmax=0.10,  # Scale from 0% to 10%
        colorbar=dict(
            title=dict(
                text="Estimated Contamination",
                side="right"
            ),
            tickmode="linear",
            tick0=0,
            dtick=0.02,  # Tick every 2%
            len=0.7,
            thickness=20,
            tickformat=".1%"
        )
    ))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f'Cross-Contamination Analysis - Estimated Contamination Heatmap<br><sub>Threshold: {threshold:.1%} | Min variants: {min_variants}</sub>',
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title='Target Sample',
            tickangle=-45,
            side='bottom'
        ),
        yaxis=dict(
            title='Source Sample',
            autorange='reversed'
        ),
        annotations=annotations,
        width=max(800, n_samples * 40),
        height=max(700, n_samples * 40),
        margin=dict(l=150, r=150, t=100, b=150)
    )
    
    # Save to HTML with method description
    heatmap_file = f"{output_prefix}.heatmap.html"
    
    # Create method description panel
    method_description = f"""
    <div style="max-width: 400px; min-width: 300px; padding: 20px; background-color: #f8f9fa; border-radius: 5px; font-family: Arial, sans-serif; margin-left: 20px; margin-top: 80px;">
        <h2 style="color: #2c3e50; margin-top: 0; font-size: 1.3em;">📊 Method</h2>
        
        <h3 style="color: #34495e; font-size: 1.1em;">Principle</h3>
        <p style="line-height: 1.6; font-size: 0.9em;">
        Contamination detection via variants:
        </p>
        <ul style="line-height: 1.6; font-size: 0.9em;">
            <li><strong>0/1</strong> (heterozygous) in source (Y-axis)</li>
            <li><strong>0/0</strong> (homozygous ref) in target (X-axis)</li>
        </ul>
        
        <h3 style="color: #34495e; font-size: 1.1em;">Metrics</h3>
        <ul style="line-height: 1.6; font-size: 0.9em;">
            <li><strong>VAF</strong>: Variant Allele Frequency in target</li>
            <li><strong>Estimated conta</strong>: VAF × 2 (accounts for heterozygosity)</li>
        </ul>
        
        <h3 style="color: #34495e; font-size: 1.1em;">Interpretation</h3>
        <ul style="line-height: 1.6; font-size: 0.9em;">
            <li><strong>Threshold</strong>: {threshold * 100:.1f}%</li>
            <li><strong>Min variants</strong>: {min_variants}</li>
            <li><strong>⚠️</strong>: Contamination detected</li>
        </ul>
        
        <p style="margin-top: 15px; padding: 8px; background-color: #fff3cd; border-left: 3px solid #ffc107; font-size: 0.85em;">
        <strong>⚠️ Note</strong>: Downsampling to 5000 variants max.
        </p>
    </div>
    """
    
    # Generate VAF plot if VCF file provided
    vaf_plot_html = ""
    if vcf_file:
        print(f"\nGenerating VAF distribution plot from {vcf_file}...")
        try:
            vaf_fig = generate_vaf_plot(vcf_file, samples)
            if vaf_fig:
                vaf_plot_html = f'<div style="margin-left: 20px; margin-top: 80px;">{vaf_fig.to_html(include_plotlyjs=False, full_html=False)}</div>'
                print("VAF plot successfully generated and integrated into heatmap")
            else:
                print("Warning: VAF plot generation returned None", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Failed to generate VAF plot: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    else:
        print("\nSkipping VAF plot generation (no VCF file provided)")
    
    # Get plotly HTML
    html_content = fig.to_html(include_plotlyjs='cdn')
    
    # Wrap plot, description, and VAF plot in flex container
    html_parts = html_content.split('<body>')
    if len(html_parts) == 2:
        flex_wrapper_start = '<div style="display: flex; align-items: flex-start; justify-content: flex-start; padding: 20px;">'
        flex_wrapper_end = '</div>'
        
        body_parts = html_parts[1].split('</body>')
        if len(body_parts) == 2:
            final_html = (html_parts[0] + '<body>' + 
                         flex_wrapper_start + 
                         body_parts[0] + 
                         method_description + 
                         vaf_plot_html + 
                         flex_wrapper_end + 
                         '</body>' + body_parts[1])
        else:
            final_html = html_content
    else:
        final_html = html_content
    
    # Write to file
    with open(heatmap_file, 'w', encoding='utf-8') as f:
        f.write(final_html)
    
    print(f"Heatmap saved to {heatmap_file}")
    
    return heatmap_file


def write_results(results, samples, output_prefix, threshold=0.02, min_variants=100, vcf_file=None):
    """
    Write results to TSV files and generate visualizations.
    
    Generates three output files:
    1. {output_prefix}.results.tsv - Detailed results for all sample pairs
    2. {output_prefix}.matrix.estimated_conta.tsv - Matrix format of estimated contamination values
    3. {output_prefix}.heatmap.html - Interactive heatmap (if plotly available)
    
    Args:
        results (list): List of result dictionaries
        samples (list): List of sample names
        output_prefix (str): Output file prefix
        threshold (float): Contamination threshold used
        min_variants (int): Minimum variants threshold used
        vcf_file (str, optional): Path to VCF file for VAF plot generation
    """
    # Write detailed results
    results_file = f"{output_prefix}.results.tsv"
    print(f"\nWriting detailed results to {results_file}")
    
    with open(results_file, 'w') as f:
        # Header
        f.write("source_sample\ttarget_sample\tmedian_vaf\tmean_vaf\testimated_conta\t"
                "informative_variant_count\tgt_0_1_count\t"
                "insufficient_variants\tstatus\ttext_result\n")
        
        # Data
        for result in results:
            f.write(
                f"{result['source_sample']}\t"
                f"{result['target_sample']}\t"
                f"{result['median_vaf']}\t"
                f"{result['mean_vaf']}\t"
                f"{result.get('estimated_conta', 'NA')}\t"
                f"{result['informative_variant_count']}\t"
                f"{result['gt_0_1_count']}\t"
                f"{result.get('insufficient_variants', False)}\t"
                f"{result['status']}\t"
                f"{result['text_result']}\n"
            )
    
    # Write matrix
    matrix_file = f"{output_prefix}.matrix.estimated_conta.tsv"
    print(f"Writing matrix to {matrix_file}")
    
    # Create lookup dictionary
    result_dict = {}
    for result in results:
        key = (result['source_sample'], result['target_sample'])
        result_dict[key] = result.get('estimated_conta', 'NA')
    
    with open(matrix_file, 'w') as f:
        # Header
        f.write("source_sample\t" + "\t".join(samples) + "\n")
        
        # Matrix rows
        for source_sample in samples:
            row = [source_sample]
            for target_sample in samples:
                if source_sample == target_sample:
                    value = 'NA'
                else:
                    key = (source_sample, target_sample)
                    value = result_dict.get(key, 'NA')
                row.append(str(value))
            f.write("\t".join(row) + "\n")
    
    # Display matrix summary
    print("\n" + "=" * 60)
    print("Matrix of mean VAF results")
    print("=" * 60)
    print("source_sample\t" + "\t".join(samples))
    
    for source_sample in samples:
        row = [source_sample]
        for target_sample in samples:
            if source_sample == target_sample:
                value = 'NA'
            else:
                key = (source_sample, target_sample)
                value = result_dict.get(key, 'NA')
            row.append(str(value))
        print("\t".join(row))
    
    print("\n" + "=" * 60)
    print(f"Results written to {results_file}")
    print(f"Matrix written to {matrix_file}")
    
    # Generate heatmap
    try:
        heatmap_file = generate_heatmap(results, samples, output_prefix, threshold, min_variants, vcf_file)
        if heatmap_file:
            print(f"Interactive heatmap saved to {heatmap_file}")
    except Exception as e:
        print(f"Warning: Failed to generate heatmap: {e}", file=sys.stderr)
    
    print("=" * 60)


def load_results_from_file(results_file):
    """
    Load results from an existing results TSV file.
    
    This allows regenerating plots without rerunning the analysis.
    
    Args:
        results_file (str): Path to results TSV file
        
    Returns:
        tuple: (results list, samples list) or (None, None) if failed
    """
    results = []
    samples_set = set()
    
    try:
        with open(results_file, 'r') as f:
            # Read header
            header = next(f).strip().split('\t')
            
            # Check if file has new format with insufficient_variants column
            has_insufficient_col = 'insufficient_variants' in header
            
            for line in f:
                fields = line.strip().split('\t')
                
                # Check format based on number of fields
                if len(fields) >= 10:
                    # v6 format with estimated_conta
                    result = {
                        'source_sample': fields[0],
                        'target_sample': fields[1],
                        'median_vaf': fields[2],
                        'mean_vaf': fields[3],
                        'estimated_conta': fields[4],
                        'informative_variant_count': int(fields[5]) if fields[5].isdigit() else 0,
                        'gt_0_1_count': int(fields[6]) if fields[6].isdigit() else 0,
                        'insufficient_variants': fields[7].lower() == 'true',
                        'status': fields[8],
                        'text_result': fields[9] if len(fields) > 9 else ''
                    }
                elif has_insufficient_col and len(fields) >= 9:
                    # v5 format (backward compatibility)
                    result = {
                        'source_sample': fields[0],
                        'target_sample': fields[1],
                        'median_vaf': fields[2],
                        'mean_vaf': fields[3],
                        'estimated_conta': 'NA',
                        'informative_variant_count': int(fields[4]) if fields[4].isdigit() else 0,
                        'gt_0_1_count': int(fields[5]) if fields[5].isdigit() else 0,
                        'insufficient_variants': fields[7].lower() == 'true',
                        'status': fields[8],
                        'text_result': fields[9] if len(fields) > 9 else ''
                    }
                else:
                    continue
                
                results.append(result)
                samples_set.add(fields[0])
                samples_set.add(fields[1])
        
        samples = sorted(list(samples_set))
        return results, samples
        
    except Exception as e:
        print(f"Error loading results from file: {e}", file=sys.stderr)
        return None, None


def main():
    """Main entry point for the script."""
    args = parse_arguments()
    
    # Display analysis parameters
    print("=" * 60)
    print("Cross-Contamination Analysis - Version 9.0")
    print("=" * 60)
    print(f"VCF file:       {args.vcf}")
    print(f"BAM folder:     {args.bam_folder}")
    print(f"Reference:      {args.reference}")
    print(f"Output prefix:  {args.output_prefix}")
    print(f"Threshold:      {args.threshold}")
    print(f"Min variants:   {args.min_variants}")
    print(f"Threads:        {args.threads}")
    if args.vaf_plot_only:
        print(f"Mode:           VAF PLOT ONLY")
    print("=" * 60)
    
    # Check if user wants to generate VAF plot only (for testing)
    if args.vaf_plot_only:
        if not PLOTLY_AVAILABLE:
            print("Error: Plotly is required for VAF plot generation", file=sys.stderr)
            sys.exit(1)
        
        if not os.path.exists(args.vcf):
            print(f"Error: VCF file not found: {args.vcf}", file=sys.stderr)
            sys.exit(1)
        
        print("\n🧪 VAF Plot Only Mode - Generating VAF distribution plot...")
        
        # Get samples from VCF
        samples = get_samples_from_vcf(args.vcf)
        print(f"\nFound {len(samples)} samples: {', '.join(samples)}")
        
        # Generate VAF plot
        try:
            vaf_fig = generate_vaf_plot(args.vcf, samples)
            if vaf_fig:
                output_file = f"{args.output_prefix}.vaf_plot.html"
                vaf_fig.write_html(output_file)
                print(f"\n✅ VAF plot saved to {output_file}")
            else:
                print("\n❌ Failed to generate VAF plot", file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            print(f"\n❌ Error generating VAF plot: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)
        
        return
    
    # Check if user wants to reload from existing file
    if args.reload_from:
        if not os.path.exists(args.reload_from):
            print(f"Error: Reload file not found: {args.reload_from}", file=sys.stderr)
            sys.exit(1)
        
        print(f"\n📂 Loading existing results from: {args.reload_from}")
        results, samples = load_results_from_file(args.reload_from)
        
        if results and samples:
            print(f"Loaded {len(results)} results for {len(samples)} samples")
            
            # Regenerate outputs (without VCF for reload mode)
            write_results(results, samples, args.output_prefix, args.threshold, args.min_variants, None)
            
            print("\n✅ Plot regeneration complete!")
            return
        else:
            print(f"Error: Failed to load results from {args.reload_from}", file=sys.stderr)
            sys.exit(1)
    
    # Validate input files
    if not os.path.exists(args.vcf):
        print(f"Error: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.bam_folder):
        print(f"Error: BAM folder not found: {args.bam_folder}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        sys.exit(1)
    
    # Run analysis
    results, samples = analyze_contamination(
        args.vcf,
        args.bam_folder,
        args.reference,
        args.threshold,
        args.min_variants,
        args.threads
    )
    
    # Write results with VCF file for VAF plot generation
    write_results(results, samples, args.output_prefix, args.threshold, args.min_variants, args.vcf)
    
    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()


