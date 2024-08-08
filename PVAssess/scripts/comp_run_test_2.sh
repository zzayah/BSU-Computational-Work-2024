#!/usr/bin/bash

#SBATCH --job-name=NAMEHERE
#SBATCH -o slurm.out
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --cpus-per-task=1
#SBATCH --signal=USR1@300
#SBATCH --time=72:00:00

# USAGE
# sbatch <path to this script (seed_batch_run.sh)> <Directory to run in>
# This script will recursively go through every directory and run calculations where applicable.

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

# Start from the top-level directory provided as input
start_dir=$1

# Process directories recursively
for d1 in "$start_dir"/*/; do
    for d2 in "$d1"*/; do
        for d3 in "$d2"*/; do
            # Correctly handle path construction
            full_path="${d3%/}"
            echo "Entering directory $full_path"
            cd "$full_path" || { echo "Could not enter directory $full_path, skipping"; continue; }

            if test -e "espresso.pwo" && grep -Pq "CPU.*WALL" "espresso.pwo"; then
                if grep -Pq "BFGS Geometry Optimization" "espresso.pwo"; then
                    if grep -Pq "bfgs converged in.*" "espresso.pwo"; then
                        run=False
                    else
                        run=True
                    fi
                else
                    if grep -Pq "total energy.*=.*" "espresso.pwo"; then
                        run=False
                    else
                        run=True
                    fi
                fi
            else
                run=True
            fi

            if [ "$run" = True ]; then
                if test -e "espresso.pwi"; then
                    echo "Running calculation in $full_path"
                    $PW_COMMAND --in "espresso.pwi" 1> "espresso.pwo" 2> "espresso.error"
                    status=$?
                    if [ $status -eq 0 ]; then
                        echo "SUCCEEDED in $full_path\n"
                    else
                        echo "FAILED in $full_path\n"
                    fi
                else
                    echo "NO INPUT FILE in $full_path, SKIPPING"
                fi
            else
                echo "Calculation already completed in $full_path, SKIPPING"
            fi

            cd ../../../ || exit 1  # Go back three levels to the original directory
        done
    done
done
