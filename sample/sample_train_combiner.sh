#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_sample_train_combiner

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/sample/"

module load matlab/7.13
$scriptdir/../train_combiner.sh -N $scriptdir/new_featured_sorted_probabilistic_Neutral_SNP -P $scriptdir/new_featured_sorted_probabilistic_Pathogenic_SNP -i 500 -s 0.007 -m 4 -n 5 -j 4 -k 5


