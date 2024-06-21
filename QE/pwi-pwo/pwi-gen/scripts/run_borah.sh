#!/usr/bin/bash

#SBATCH --job-name NAMEHERE
#SBATCH -o slurm.out
#SBATCH -N 2
#SBATCH -n 96
#SBATCH --cpus-per-task=1
#SBATCH --signal=USR1@300
#SBATCH --time 24:00:00

#USAGE
#sbatch <path to this script (seed_batch_run.sh)> <Directory to run in>
#It is safest to cd into the directory you will be running in
#This will run calculations in every directory in the provided directory

module load slurm
module load borah-base
module load borah-libraries
module load intel/mkl oneapi/2022.1.0
module load gcc/12.1.0
module load openmpi/4.1.3/oneapi
module load python/3.9.7

PREFIX="$HOME/q-e"
BIN_DIR="$PREFIX/bin"
PSEUDO_DIR="$PREFIX/pseudo"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export ESPRESSO_TMPDIR="$HOME/scratch/"
export ESPRESSO_PSEUDO="$HOME/q-e/pseudo"

PARA_PREFIX="srun --mpi=pmix -N $SLURM_NNODES -n $SLURM_NTASKS"
PARA_POSTFIX="-nk 2"
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX"

if test -z "$1"; then 
    DIR=./
else
    DIR="$1/"
fi
echo "Running in $DIR"

#cd $DIR || echo "$DIR not found" && exit 1
cd $DIR
echo "Running in $DIR"
#cd  "0" || echo "$DIR not found" && exit 1
cd 0
echo "In directory "0""

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

for i in $(ls -d ./*/); do
    cd "$i"
    if test -e "espresso.pwo" && grep -Pq "CPU.*WALL" "espresso.pwo"; then
        if grep -Pq "BFGS Geometry Optimization" "espresso.pwo"; then
            if grep -Pq "bfgs converged in.*" "espresso.pwo"; then
                run=False
            else
                run=True
            fi
        else
            grep -Pq "total energy.*=.*" "espresso.pwo"; then
                run=False
            else
                run=True
            fi
        fi
    else
        run=True
    fi
    if $run; then
        if test -e "espresso.in"; then
            echo -E "$1"
            $PW_COMMAND --in "espresso.pwi" 1> "espresso.pwo" 2> "espresso.error"
            status=$?
            if $? ; then
                echo "SUCCEEDED\n"
            else
                echo "FAILED\n"
            fi
        else
            echo "$1 NO INPUT FILE SKIPPING"
        fi
    else
        echo "$1 DONE SKIPPING"
    fi
    cd ../
done

cd  "../1" || echo "$DIR not found" && exit 1

for i in $(ls -d ./*/); do
    cd "$i"
    if test -e "espresso.pwo" && grep -Pq "CPU.*WALL" "espresso.pwo"; then
        if grep -Pq "BFGS Geometry Optimization" "espresso.pwo"; then
            if grep -Pq "bfgs converged in.*" "espresso.pwo"; then
                run=False
            else
                run=True
            fi
        else
            grep -Pq "total energy.*=.*" "espresso.pwo"; then
                run=False
            else
                run=True
            fi
        fi
    else
        run=True
    fi
    if $run; then
        if test -e "espresso.in"; then
            echo -E "$1"
            $PW_COMMAND --in "espresso.pwi" 1> "espresso.pwo" 2> "espresso.error"
            status=$?
            if $? ; then
                echo "SUCCEEDED\n"
            else
                echo "FAILED\n"
            fi
        else
            echo "$1 NO INPUT FILE SKIPPING"
        fi
    else
        echo "$1 DONE SKIPPING"
    fi
    cd ../
done