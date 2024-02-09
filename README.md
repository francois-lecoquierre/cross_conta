# cross_conta
Scripts for NGS contamination

## contaminate_bam.sh

Used to simulate a contamination of a bam_A by bam_B at a selected rate

Dependencies needed:
- picard.jar (used to downsample bam files)
- samtools (to merge the bam files)

Usage:
- script.sh -p <picard.jar file> -c <bam_contaminator> -d <bam_contaminated> -o <output_folder> -r <contamination_rate>


## conta_check.sh

Takes a multi-vcf as an input. For each sample pair, detects the contamination of sample_A (source) on sample_B (target). 
The script first finds informative genotypes : 
- Heterozygous on source
- Wild type on target
And computes the mean allelic fraction on sample B.
If the mean variant allelic fraction (VAF) is over the threshold, the target sample is considered contaminated by the source sample.
In summary, the output value is the mean VAF in sample_B of variant that are het in sample_A and wt in sample_B. Therefore the contamination rate can be seen as twice as this value.

Dependencies needed:
- bcftools

Usage:
- conta_check.sh -v vcf_file -o output_prefix -t threshold
