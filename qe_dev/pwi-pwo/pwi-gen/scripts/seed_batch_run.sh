#!/usr/bin/env bash

#SBATCH --job-name LaCoO3
#SBATCH -o slurm.out
#SBATCH -N 2
#SBATCH -n 96
#SBATCH --cpus-per-task=1
#SBATCH -p short
#SBATCH --signal=USR1@300
#SBATCH --time 24:00:00

#USAGE
#sbatch <path to this script (seed_batch_run.sh)> <Directory to run in>
#It is safest to cd into the directory you will be running in
#This will run calculations in every directory in the provided directory

module load slurm
module load borah-base borah-libraries
module load intel/mkl oneapi/2022.1.0
module load gcc/12.1.0
module load openmpi/4.1.3/oneapi
module load python/3.9.7

PREFIX="$HOME/q-e"
BIN_DIR="$PREFIX/bin"
PSEUDO_DIR="$PREFIX/pseudo"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export ESPRESSO_TMPDIR="$HOME/scratch/"
export ESPRESSO_PSEUDO="$HOME/q-e/pseudo"

PARA_PREFIX="srun --mpi=pmix -N $SLURM_NNODES -n $SLURM_NTASKS"
PARA_POSTFIX="-nk 2"
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"

cd $1 || exit 1

touch checkpoint.log
for dir in $(ls -d ./*/); do
	cd $dir
	outfile="espresso"
	echo "Currently running in $dir"
	if test ! -f *.pwi; then
		echo "No .pwi input file found in $dir" >> /dev/stderr && exit 1
	fi
	if ! grep -q "^$dir$" ../checkpoint.log; then
		if ! grep -q '^!' "$outfile".pwo; then
			$PW_COMMAND -i *.pwi 1>"$outfile".pwo 2>"$outfile".error
			test ! -s "$outfile".error && rm "$outfile".error
			wait
			echo "$dir" >> ../checkpoint.log
		else
			echo "$dir" >> ../checkpoint.log
		fi
	fi
	cd ../
done