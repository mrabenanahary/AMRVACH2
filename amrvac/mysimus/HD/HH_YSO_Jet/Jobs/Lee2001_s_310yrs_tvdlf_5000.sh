#!/bin/bash
#SBATCH --job-name=Lee2001_s_310yrs_tvdlf_5000
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=7200
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --output=slurm-Lee2001_s_310yrs_tvdlf_5000-%j.out
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
#srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
#srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p uniform
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p uniform/${SLURM_JOB_NAME}
#srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/${SLURM_JOB_NAME}.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/uniform/${SLURM_JOB_NAME}.par uniform/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jet_CC_0050.dat .

echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
#srun --ntasks=$SLURM_JOB_NUM_NODES more ${SLURM_JOB_NAME}.par
#mpiexec amrvac -i ${SLURM_JOB_NAME}.par
#srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/.

#srun --ntasks=$SLURM_JOB_NUM_NODES cp Output/${SLURM_JOB_NAME}/Jet_CC_0010.dat .
srun --ntasks=$SLURM_JOB_NUM_NODES more uniform/${SLURM_JOB_NAME}/${SLURM_JOB_NAME}.par
mpiexec amrvac -i uniform/${SLURM_JOB_NAME}/${SLURM_JOB_NAME}.par

#srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/${SLURM_JOB_NAME}/* /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES mv uniform /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/uniform/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
