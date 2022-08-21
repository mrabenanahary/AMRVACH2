#!/bin/bash
#SBATCH --job-name=Lee2001_us_115yrs_fgrad
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=7200
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --mem=64gb
#SBATCH --tmp=200gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${HOSTNAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir Lee2001_us_115yrs_fgrad

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/Lee2001_us_115yrs_fgrad
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/Lee2001_us_115yrs_fgrad.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more Lee2001_us_115yrs_fgrad.par
mpiexec amrvac -i Lee2001_us_115yrs_fgrad.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/Lee2001_us_115yrs_fgrad/* /data/$USER/Output/Lee2001_us_115yrs_fgrad/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOBID}.out /data/$USER/Output/Lee2001_us_115yrs_fgrad/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
