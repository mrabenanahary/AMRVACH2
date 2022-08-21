#!/bin/bash
#SBATCH --job-name=Lee2001_us_300yrs
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=7200
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --output=slurm-Lee2001_us_300yrs-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mem=64gb
#SBATCH --tmp=200gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir ${SLURM_JOB_NAME}

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/${SLURM_JOB_NAME}.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more ${SLURM_JOB_NAME}.par
mpiexec amrvac -i ${SLURM_JOB_NAME}.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/${SLURM_JOB_NAME}/* /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
