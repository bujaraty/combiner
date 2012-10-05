#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_sample_prediction

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/sample/"

#module load matlab/7.13
$scriptdir/../predict_pathogenic_probability.sh -P $scriptdir/combiner_params_0.3849.mat -f $scriptdir/featured_input -o $scriptdir/predict_pathogenic_probability.out

stty sane
