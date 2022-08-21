#!/bin/bash
#SBATCH --job-name=Lee2001_s_300yrs_1e4yrs
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --nodelist=tycho70
#SBATCH --time=7200
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --output=slurm-Lee2001_s_300yrs_1e4yrs-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mem=64gb
#SBATCH --tmp=200gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
if [ -d "/data/$USER/Output/${SLURM_JOB_NAME}" ] 
then
    srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf /data/$USER/Output/${SLURM_JOB_NAME}
fi
mkdir ${SLURM_JOB_NAME}

cd $SCRATCH
#mkdir amrvac
#srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES cp -r $AMRVAC_DIR .
SCRATCH_AMRVAC_PROJECT=$SCRATCH/amrvac/mysimus/HD/HH_YSO_Jet
cd $SCRATCH_AMRVAC_PROJECT


srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p uniform/${SLURM_JOB_NAME}
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .


export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


echo $SLURM_CPUS_ON_NODE


srun --ntasks=$SLURM_JOB_NUM_NODES more Parfiles/Article/${SLURM_JOB_NAME}.par
mpiexec -np $SLURM_CPUS_ON_NODE amrvac -i Parfiles/Article/${SLURM_JOB_NAME}.par
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/.

#srun --ntasks=$SLURM_JOB_NUM_NODES cp Output/${SLURM_JOB_NAME}/Jet_CC_0100.dat .
srun --ntasks=$SLURM_JOB_NUM_NODES more Parfiles/Article/uniform/${SLURM_JOB_NAME}.par
mpiexec -np $SLURM_CPUS_ON_NODE amrvac -i Parfiles/Article/uniform/${SLURM_JOB_NAME}.par

srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/${SLURM_JOB_NAME}/* /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES mv uniform /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/uniform/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
