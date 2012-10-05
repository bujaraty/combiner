#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/"

#define default values
STEP_SIZE_DEFAULT="0.0002"
ITERATION_DEFAULT="10000"
MIN_FIRST_HIDDEN_NODES_DEFAULT="3"
MAX_FIRST_HIDDEN_NODES_DEFAULT="5"
MIN_SECOND_HIDDEN_NODES_DEFAULT="3"
MAX_SECOND_HIDDEN_NODES_DEFAULT="5"
OUTPUT_DIR_DEFAULT=$scriptdir/result

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-N FILE    specify path of featured neutral variants file (required)
-P FILE    specify path of featured pathogenic variants file (required)
-s VALUE   set step size (default is $STEP_SIZE_DEFAULT)
-i VALUE   set number of iteration (default is $ITERATION_DEFAULT)
-j VALUE   set number of minimum hidden nodes in the first hidden layer (default is $MIN_FIRST_HIDDEN_NODES_DEFAULT)
-k VALUE   set number of maximum hidden nodes in the first hidden layer (default is $MAX_FIRST_HIDDEN_NODES_DEFAULT)
-m VALUE   set number of minimum hidden nodes in the second hidden layer (default is $MIN_SECOND_HIDDEN_NODES_DEFAULT)
-n VALUE   set number of maximum hidden nodes in the second hidden layer (default is $MAX_SECOND_HIDDEN_NODES_DEFAULT)
-O DIR     set output directory (default is $OUTPUT_DIR_DEFAULT). This directory will contain figures and parameter files
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#get file
while getopts "N:P:i:j:k:m:n:s:O:" OPTION; do
  case "$OPTION" in
    N)
      featured_neutral_data="$OPTARG"
      ;;
    P)
      featured_pathogenic_data="$OPTARG"
      ;;
    i)
      iteration="$OPTARG"
      ;;
    j)
      min_first_hidden_nodes="$OPTARG"
      ;;
    k)
      max_first_hidden_nodes="$OPTARG"
      ;;
    m)
      min_second_hidden_nodes="$OPTARG"
      ;;
    n)
      max_second_hidden_nodes="$OPTARG"
      ;;
    s)
      step_size="$OPTARG"
      ;;
    O)
      output_dir="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

[ ! -z $featured_neutral_data ] || die "No neutral variants file provided"
[ ! -z $featured_pathogenic_data ] || die "No pathogenic variants file provided"
[ -f $featured_neutral_data ] || die "$featured_neutral_data is not a valid file name"
[ -f $featured_pathogenic_data ] || die "$featured_pathogenic_data is not a valid file name"

#setting default values:
: ${step_size=$STEP_SIZE_DEFAULT}
: ${iteration=$ITERATION_DEFAULT}
: ${min_first_hidden_nodes=$MIN_FIRST_HIDDEN_NODES_DEFAULT}
: ${max_first_hidden_nodes=$MAX_FIRST_HIDDEN_NODES_DEFAULT}
: ${min_second_hidden_nodes=$MIN_SECOND_HIDDEN_NODES_DEFAULT}
: ${max_second_hidden_nodes=$MAX_SECOND_HIDDEN_NODES_DEFAULT}
: ${output_dir=$OUTPUT_DIR_DEFAULT}

#show the values as read in by the flags
cat <<EOF
training configuration:
step size               : $step_size
iteration               : $iteration
min first hidden nodes  : $min_first_hidden_nodes
max first hidden nodes  : $max_first_hidden_nodes
min second hidden nodes : $min_second_hidden_nodes
max second hidden nodes : $max_second_hidden_nodes
neutral data            : $featured_neutral_data
pathogenic data         : $featured_pathogenic_data
output directory        : $output_dir
EOF

if [ ! -e $output_dir ]
then
	mkdir $output_dir
fi

if [ ! -e $scriptdir/tmp ]
then
	mkdir $scriptdir/tmp
fi

shuffled_featured_neutral_data=$scriptdir/tmp/shuffled_featured_neutral_data
shuffled_featured_pathogenic_data=$scriptdir/tmp/shuffled_featured_pathogenic_data

tmp_training_dataset=$scriptdir/tmp/tmp_training_dataset
tmp_validation_dataset=$scriptdir/tmp/tmp_validation_dataset
tmp_test_dataset=$scriptdir/tmp/tmp_test_dataset

training_dataset=$scriptdir/tmp/training_dataset
validation_dataset=$scriptdir/tmp/validation_dataset
test_dataset=$scriptdir/tmp/test_dataset

function separate_training_dataset {
	cat $1 > tmp.txt
	record_count=$( cat tmp.txt | wc -l  )

	training_count=$( printf "%.0f\n" $( echo "scale=2;"$record_count"*70/100" | bc ) )
	validating_count=$( printf "%.0f\n" $( echo "scale=2;"$record_count"*15/100" | bc ) )

	sed -n "1,$training_count"p tmp.txt >> $tmp_training_dataset
	sed -n "$[$training_count+1],$[$training_count+$validating_count]"p tmp.txt >> $tmp_validation_dataset
	sed -n "$[$training_count+$validating_count+1],$"p tmp.txt >> $tmp_test_dataset

	rm tmp.txt
}  

if [ -e $tmp_training_dataset ]
then
	rm $tmp_training_dataset
fi

if [ -e $tmp_validation_dataset ]
then
	rm $tmp_validation_dataset
fi

if [ -e $tmp_test_dataset ]
then
	rm $tmp_test_dataset
fi

shuf $featured_neutral_data > $shuffled_featured_neutral_data
shuf $featured_pathogenic_data > $shuffled_featured_pathogenic_data

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
	grep "^"$i"|" $shuffled_featured_neutral_data | separate_training_dataset
	grep "^"$i"|" $shuffled_featured_pathogenic_data | separate_training_dataset
done

sort -k1 $tmp_training_dataset > $training_dataset
sort -k1 $tmp_validation_dataset > $validation_dataset
sort -k1 $tmp_test_dataset > $test_dataset

#module load matlab/7.13
matlab -nosplash -nodesktop -r "cd "$scriptdir"; try, calibrate_perceptron($step_size, $iteration, $min_first_hidden_nodes:$max_first_hidden_nodes, $min_second_hidden_nodes:$max_second_hidden_nodes, '$training_dataset', '$validation_dataset', '$test_dataset', '$output_dir'); end; quit"


