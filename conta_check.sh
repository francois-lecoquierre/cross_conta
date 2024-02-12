


# this bash script is used to the cross contamination between samples in a vcf
# for each source_sample in the vcf, it will check if the source_sample is present in all other samples (target_samples)
# to do so, we use the variants that are 0/1 in the source_sample and 0/0 in the target_samples
# for each target sample we compute the median VAF of the variants that are 0/0 in the target sample and 0/1 in the source sample
# if the median VAF is above a certain threshold, we consider that the source_sample is contaminated by the target_sample
# the threshold is set to 0.1 by default but can be changed with the -t option
# the output is a tabulated file with the interaction between all samples

# usage: bash conta_check.sh -v vcf_file -o output_file -t threshold


# parse options
while getopts "v:o:t:" option
do
    case $option in
        v) vcf_file=$OPTARG;;
        o) output_prefix=$OPTARG;;
        t) threshold=$OPTARG;;
    esac
done

# check that all the options are set
if [ -z $vcf_file ] || [ -z $output_prefix ] || [ -z $threshold ]; then
    echo "Synopsys: conta_check.sh -v vcf_file -o output_prefix -t threshold"
    exit 1
fi

# if the separator is a comma, we replace it by a point
threshold=$(echo "$threshold" | tr ',' '.')
# display the options
echo '******************'
echo 'Options:'
echo "vcf_file: $vcf_file"
echo "output_prefix: $output_prefix"
echo "threshold: $threshold"
echo '******************'


# check if the vcf file exists
if [ ! -f $vcf_file ]; then
    echo "vcf file not found"
    exit 1
fi



# get the list of samples in the vcf
samples=$(bcftools query -l $vcf_file)
sample_count=$(echo $samples | wc -w)

# create the header of the output file
output_result=$output_prefix.results.tsv
echo -e "source_sample\ttarget_sample\tmedian_vaf\tmean_vaf\tinformative_variant_count\tstatus\ttext_result" > $output_result

# create the temporary folder if it does not exist
if [ ! -d tmp ]; then
    mkdir tmp
fi

# for each source sample
for source_sample in $samples
do
    # we count the number of source samples processed for user feedback
    source_sample_count=$(($source_sample_count + 1))
    # we get the index of the source sample in the vcf
    source_sample_index=$(bcftools query -l $vcf_file | grep -nx $source_sample | cut -d: -f1)
    source_sample_index=$(($source_sample_index - 1))
    # display status
    echo "***"
    echo "***"
    echo "***"
    echo "source sample #$source_sample_count/$sample_count: $source_sample ($source_sample_index in vcf)"
    # for each target sample
    for target_sample in $samples
        do
        target_sample_index=$(bcftools query -l $vcf_file | grep -nx $target_sample | cut -d: -f1)
        target_sample_index=$(($target_sample_index - 1)) # can fail if two sample match the grep...
        # we skip the source sample
        if [ $source_sample == $target_sample ]; then
            continue
        fi
        # we get the variants that are 0/1 in the source sample and 0/0 in the target sample
        # we write the result in a temporary file as chrom  pos  ref  alt  and the genotype info of the target sample that are needed to compute the VAF : dp[0]  ad
        bcftools view -i "GT[$source_sample_index] == \"0/1\" && GT[$target_sample_index] == \"0/0\"" $vcf_file | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%AD{1}\t%DP]\n" -s $target_sample > tmp/$source_sample.$target_sample.tmp
        # we remove the lines where the depth ($6) is 0 and we add the VAF to each line as field $7 (which is computed as $5/$6)
        awk -F"\t" '$6 != 0 {print $0"\t"$5/$6}' tmp/$source_sample.$target_sample.tmp > tmp/$source_sample.$target_sample.vaf
        # we count the informative variants
        informative_variant_count=$(wc -l tmp/$source_sample.$target_sample.vaf | awk '{print $1}')
        # we compute the median and mean VAF
        median_vaf=$(awk '{print $7}' tmp/$source_sample.$target_sample.vaf | sort -n | awk '{data[NR]=$1} END{print data[int(NR/2)+1]}')
        median_vaf=$(echo "$median_vaf" | tr ',' '.')
        mean_vaf=$(awk '{print $7}' tmp/$source_sample.$target_sample.vaf | awk '{sum+=$1} END{print sum/NR}')
        mean_vaf=$(echo "$mean_vaf" | tr ',' '.')
        # we check if the mean VAF is above the threshold using awk
        if [ $(awk -v mean_vaf="$mean_vaf" -v threshold="$threshold" 'BEGIN {print (mean_vaf > threshold)}') -eq 1 ]; then
            result_summary="CONTA"
            text_result="WARNING | contamination detected: $target_sample is contaminated by $source_sample: heterozygous variants (GT=0/1) of $source_sample that are absent (GT=0/0) in $target_sample have a mean VAF of $mean_vaf in $target_sample (which is above the threshold: $threshold). Number of informative variants used: $informative_variant_count"
        else
            result_summary="OK"
            text_result="OK | $target_sample is not contaminated by $source_sample: heterozygous variants (GT=0/1) of $source_sample that are absent (GT=0/0) in $target_sample have a mean VAF of $mean_vaf in $target_sample (which is below the threshold: $threshold). Number of informative variants used: $informative_variant_count"
        fi
        # we append the result to the output file
        line="$source_sample\t$target_sample\t$median_vaf\t$mean_vaf\t$informative_variant_count\t$result_summary\t$text_result"
        echo -e $line
        echo -e $line >> $output_result
    
    done
    #
done



# Display the results as a matrix
echo "***************************"
echo "Matrix of the results"
echo "***************************"

output_matrix_file=$output_prefix.matrix.mean.tsv
header="source_sample"
for sample in $samples; do
    header="$header\t$sample"
done
echo -e "$header" > $output_matrix_file
echo -e "$header"

# Créer la matrice en parcourant chaque échantillon source
for source_sample in $samples; do
    line="$source_sample"
    for target_sample in $samples; do
        value=$(awk -v source="$source_sample" -v target="$target_sample" '$1 == source && $2 == target {print $4}' "$output_result")
        # Si la valeur n'est pas trouvée, utiliser 0
        if [ -z "$value" ]; then
            value='NA'
        fi
        line="$line\t$value"
    done

    echo -e "$line"
    echo -e "$line" >> $output_matrix_file
done

echo "***************************"
echo "Results written in $output_result"
echo "Matrix written in $output_matrix_file"
echo "***************************"

# remove the temporary folder
# rm -r tmp


