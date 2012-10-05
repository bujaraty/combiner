#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_prediction

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/"

#define default values
OUTPUT_FILE_DEFAULT=$scriptdir/result/predicted_probability.txt

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-P FILE    specify the parameters file (required)
-f FILE    specify path of featured variants file (required)
-o FILE    set output file (default is $OUTPUT_FILE_DEFAULT).
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#get file
while getopts "f:o:P:" OPTION; do
  case "$OPTION" in
    f)
      featured_input="$OPTARG"
      ;;
    o)
      output_file="$OPTARG"
      ;;
    P)
      params_file="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

[ ! -z $featured_input ] || die "No input file provided"
[ ! -z $params_file ] || die "No parameters file provided"
[ -f $featured_input ] || die "$featured_input is not a valid file name"
[ -f $params_file ] || die "$params_file is not a valid file name"

#setting default values:
: ${output_file=$OUTPUT_FILE_DEFAULT}

#show the values as read in by the flags
cat <<EOF
training configuration:
input file      : $featured_input
output file     : $output_file
parameters file : $params_file
EOF

#module load matlab/7.13
matlab -nosplash -nodesktop -r "cd "$scriptdir"; try, predict_pathogenic_probability('$params_file', '$featured_input', '$output_file'); end; quit"


