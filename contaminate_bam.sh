
# This script is used to contaminate a bam file with another bam file
# synopsys : script.sh -p <picard.jar file> -c <bam_contaminator> -d <bam_contaminated> -o <output_folder> -r <contamination_rate>
# parse args
while getopts p:c:d:o:r: flag
do
    case "${flag}" in
        p) picard=${OPTARG};;
        c) bam_contaminator=${OPTARG};;
        d) bam_contaminated=${OPTARG};;
        o) output_folder=${OPTARG};;
        r) contamination_rate=${OPTARG};;
    esac
done

# if one arg is missing, display usage and exit
if [ -z "$picard" ] || [ -z "$bam_contaminator" ] || [ -z "$bam_contaminated" ] || [ -z "$output_folder" ] || [ -z "$contamination_rate" ]
then
    echo "Usage: $0 -p <picard.jar file> -c <bam_contaminator> -d <bam_contaminated> -o <output_folder> -r <contamination_rate>"
    exit 1
fi

# if the output folder does not exist, create it
if [ ! -d "$output_folder" ]
then
    mkdir -p $output_folder
fi

# tmp folder is ./tmp
tmp_folder="./tmp"
if [ ! -d "$tmp_folder" ]
then
    mkdir -p $tmp_folder
fi

# display the args
echo "************************************************************************"
echo "picard: $picard"
echo "bam_contaminator: $bam_contaminator"
echo "bam_contaminated: $bam_contaminated"
bam_contaminated_prefix=$(basename $bam_contaminated)
echo "output_folder: $output_folder"
echo "contamination_rate: $contamination_rate"
# compute the contaminated rate which is 1 - contamination_rate
if $(awk -v n="$contamination_rate" 'BEGIN {exit !(n >= 0 && n <= 1)}');
then
    echo "contamination_rate must be between 0 and 1"
    exit 1
fi
contaminated_rate=$(awk "BEGIN { print 1 - $contamination_rate }")
echo "contaminated_rate: $contaminated_rate"
echo "************************************************************************"

# use picard to downsample the bam files
echo "Downsampling bam files"
java -jar $picard DownsampleSam I=$bam_contaminator O=$tmp_folder/contaminator.$contamination_rate.bam P=$contamination_rate
java -jar $picard DownsampleSam I=$bam_contaminated O=$tmp_folder/contaminated.$contamination_rate.bam P=$contaminated_rate
echo "Downsampling done"
echo "contaminator.$contamination_rate.bam created"
echo "contaminated.$contamination_rate.bam created"

# merge the sub_bam files
echo "Merging the sub_bam files"
samtools merge $output_folder/${bam_contaminated_prefix}.contaminated.$contamination_rate.bam $tmp_folder/contaminated.$contamination_rate.bam $tmp_folder/contaminator.$contamination_rate.bam
echo "Merge done"

# index the merged bam file
echo "Indexing the merged bam file"
samtools index $output_folder/${bam_contaminated_prefix}.contaminated.$contamination_rate.bam
echo "Index done"

# sort the merged bam file
# echo "Sorting the merged bam file"
# samtools sort $output_folder/${bam_contaminated_prefix}.contaminated.$contamination_rate.bam -o $output_folder/${bam_contaminated_prefix}.contaminated.$contamination_rate.sorted.bam
# echo "Sort done"

# index the sorted merged bam file
# echo "Indexing the sorted merged bam file"
# samtools index $output_folder/${bam_contaminated_prefix}.contaminated.$contamination_rate.sorted.bam
# echo "Index done"

# delete tmp files
rm $tmp_folder/contaminator.$contamination_rate.bam
rm $tmp_folder/contaminated.$contamination_rate.bam



