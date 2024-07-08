#!/usr/bin/bash
#SBATCH --job-name=NAMEHERE
#SBATCH -o slurm.out
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --cpus-per-task=1
#SBATCH --signal=USR1@300
#SBATCH --time=24:00:00

#USAGE
#sbatch <path to this script (run_espresso.sh)>

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

# Change to the directory where the script is located

# Check if the input file exists
if test -e "espresso.pwi"; then
    # Update the pseudo_dir in the input file
    # sed -i "s/pseudo_dir.*/pseudo_dir       =  "\'$ESPRESSO_PSEUDO\'"/g" espresso.pwi
    echo "Running espresso.pwi"
    $PW_COMMAND -i "espresso.pwi" 1> "espresso.pwo" 2> "espresso.error"
    status=$?
    if [ $status -eq 0 ]; then
        echo "SUCCEEDED"
    else
        echo "FAILED"
    fi
else
    echo "NO INPUT FILE: espresso.pwi NOT FOUND"
fi
