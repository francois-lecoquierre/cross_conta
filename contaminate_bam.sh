
# This script is used to contaminate a bam file with another bam file
# synopsys : script.sh -p <picard.jar file> -c <bam_contaminator> -d <bam_contaminated> -o <output_folder> -r <contamination_rate>
# parse args

# parse options
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
echo "*** Options:"
echo "picard: $picard"
echo "bam_contaminator: $bam_contaminator"
echo "bam_contaminated: $bam_contaminated"
bam_contaminated_prefix=$(basename $bam_contaminated)
echo "output_folder: $output_folder"
# compute the contaminated rate which is 1 - contamination_rate
if awk -v rate="$contamination_rate" 'BEGIN { if (rate >= 0 && rate <= 1) exit 0; else exit 1 }'; then
    echo "contamination_rate: $contamination_rate"
else
    echo "contamination_rate must be between 0 and 1"
    echo "contamination_rate: $contamination_rate"
    exit 1
fi
contaminated_rate=$(awk "BEGIN { print 1 - $contamination_rate }")
echo "contaminated_rate: $contaminated_rate"
echo "************************************************************************"

# compute the contaminator_extraction_rate and the contaminated_extraction_rate
# the bam that we want to generate will be of the size of the smallest, so that it can be generated for any contamination rate between 0 and 1
# example : if the contaminator has 1000 reads and the contaminated has 2000 reads: the output bam will be of set to 1000 reads
# if the conta rate is 0.1, we will extract 100 reads from the contaminator and 900 reads from the contaminated
# therefore in this example
# contamination rate = 0.1
# contaminator_extraction_rate = 0.1
# contaminated rate = 0.9
# contaminated_extraction_rate = 0.9 * (1000 / 2000) = 0.45

# get the number of reads in the contaminator and the contaminated
echo "*** Adjusting on differential number of reads between bam files"
echo "Counting the number of reads in the bam files"
contaminator_read_count=$(samtools view -c $bam_contaminator)
contaminated_read_count=$(samtools view -c $bam_contaminated)


# compute the extraction rates
if [ $contaminator_read_count -lt $contaminated_read_count ]
then
    contaminator_extraction_rate=$contamination_rate
    contaminated_extraction_rate=$(awk "BEGIN { print $contaminated_rate * ($contaminator_read_count / $contaminated_read_count) }")
    predicted_output_read_count=$contaminator_read_count
else
    contaminated_extraction_rate=$contaminated_rate
    contaminator_extraction_rate=$(awk "BEGIN { print $contamination_rate * ($contaminated_read_count / $contaminator_read_count) }")
    predicted_output_read_count=$contaminated_read_count
fi

# display the sample names
echo "Extracting the sample names from the bam files"
contaminator_sample_name=$(samtools view -H $bam_contaminator | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
contaminated_sample_name=$(samtools view -H $bam_contaminated | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
echo "contaminator_sample_name: $contaminator_sample_name"
echo "contaminated_sample_name: $contaminated_sample_name"
echo "************************************************************************"


# compute the extraction rates
if [ $contaminator_read_count -lt $contaminated_read_count ]
then
    contaminator_extraction_rate=$contamination_rate
    contaminated_extraction_rate=$(awk "BEGIN { print $contaminated_rate * ($contaminator_read_count / $contaminated_read_count) }")
    predicted_output_read_count=$contaminator_read_count
else
    contaminated_extraction_rate=$contaminated_rate
    contaminator_extraction_rate=$(awk "BEGIN { print $contamination_rate * ($contaminated_read_count / $contaminator_read_count) }")
    predicted_output_read_count=$contaminated_read_count
fi


# display the extraction rates
echo "contaminator_read_count: $contaminator_read_count"
echo "contaminated_read_count: $contaminated_read_count"
echo "contaminator_extraction_rate: $contaminator_extraction_rate"
echo "contaminated_extraction_rate: $contaminated_extraction_rate"
echo "output bam will have approx. $predicted_output_read_count reads"
echo "************************************************************************"


# use picard to downsample the bam files
echo "Extracting random reads from both bam files"
contamiator_bam=$tmp_folder/contaminator.${contaminator_sample_name}.$contamination_rate.bam
contaminated_bam=$tmp_folder/contaminated.${contaminated_sample_name}.$contamination_rate.bam
java -jar $picard DownsampleSam I=$bam_contaminator O=$contamiator_bam P=$contaminator_extraction_rate
java -jar $picard DownsampleSam I=$bam_contaminated O=$contaminated_bam P=$contaminated_extraction_rate
echo "Downsampling done"
echo "$contamiator_bam created"
echo "$contaminated_bam created"
echo "************************************************************************"

# merge the sub_bam files
echo "Merging the sub_bam files"
merged_sample_name="$contaminated_sample_name.$contamination_rate"
merged_bam=$tmp_folder/${contaminated_sample_name}.contaminated.$contamination_rate.bam
output_bam=$output_folder/${contaminated_sample_name}.contaminated.$contamination_rate.newheader.bam
samtools merge $merged_bam $contamiator_bam $contaminated_bam
echo "Merging done"
echo "$merged_bam created"
echo "************************************************************************"


# update apply the contaminated sample name to the read groups of the merged bam file
echo "Updating the header of the merged bam file"
echo "Sample name that will be used in the merged bam file: $merged_sample_name"
samtools view -H $merged_bam | sed "s/SM:.*\t/SM:$merged_sample_name\t/" > $merged_bam.newheader.sam
samtools reheader $merged_bam.newheader.sam $merged_bam > $output_bam
# mv $merged_bam.reheader.bam $merged_bam
echo "************************************************************************"

# index the merged bam file
echo "Indexing the merged bam file"
samtools index $output_bam
echo "Index done"
echo "************************************************************************"

# delete tmp files
# echo "Cleaning up"
# rm $tmp_folder/header.$contamination_rate.sam
# rm $contamiator_bam
# rm $contaminated_bam
# rm $merged_bam
# echo "Clean up done"
# echo "************************************************************************"



