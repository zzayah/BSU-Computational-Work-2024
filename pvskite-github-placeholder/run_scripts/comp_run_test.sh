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
# It is safest to cd into the directory you will be running in
# This will run calculations in every directory in the provided directory

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

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Function to process directories
process_directory() {
    local dir=$1
    cd $dir
    for sub_dir in $(ls -d */); do
        cd "$sub_dir"
        if [[ -e "espresso.pwi" ]]; then
            if [[ -e "espresso.pwo" ]] && grep -Pq "CPU.*WALL" "espresso.pwo"; then
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
            
            if [[ $run == True ]]; then
                sed -i "s/pseudo_dir.*/pseudo_dir       =  '$ESPRESSO_PSEUDO'/g" espresso.pwi
                echo -E "$dir/$sub_dir"
                $PW_COMMAND -i "espresso.pwi" 1> "espresso.pwo" 2> "espresso.error"
                status=$?
                if [[ $status -eq 0 ]]; then
                    echo "SUCCEEDED\n"
                else
                    echo "FAILED\n"
                fi
            else
                echo "$dir/$sub_dir DONE SKIPPING"
            fi
        else
            echo "$dir/$sub_dir NO INPUT FILE SKIPPING"
        fi
        cd ..
    done
    cd ..
}

# Main script
for parent_dir in $(ls -d */); do
    cd $parent_dir
    for sub_parent_dir in $(ls -d */); do
        cd $sub_parent_dir
        if [[ -d "0" ]]; then
            echo "Processing $parent_dir/$sub_parent_dir/0"
            process_directory 0
        fi
        if [[ -d "1" ]]; then
            echo "Processing $parent_dir/$sub_parent_dir/1"
            process_directory 1
        fi
        cd ..
    done
    cd ..
done
