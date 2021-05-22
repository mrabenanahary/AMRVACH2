#!/bin/bash
#SBATCH --job-name=Shang06e
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=4320
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --mem=32gb
#SBATCH --tmp=100gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${HOSTNAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir Shang06e

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/Shang06e
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/Shang06/Shang06e.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more Shang06e.par
mpiexec amrvac -i Shang06e.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/Shang06e/Jet_Shang06_0100.dat /data/$USER/Output/Shang06e/.
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/Shang06e/* /data/$USER/Output/Shang06e/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOBID}.out /data/$USER/Output/Shang06e/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
