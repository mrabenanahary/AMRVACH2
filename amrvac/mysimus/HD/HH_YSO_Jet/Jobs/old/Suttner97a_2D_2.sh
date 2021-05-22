#!/bin/bash
#SBATCH --job-name=Suttner97_2D_2
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
mkdir Suttner97_2D_2

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/Suttner97_2D
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/mi_these/Suttner97_2D_2.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more Suttner97_2D_2.par
mpiexec amrvac -i Suttner97_2D_2.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/Suttner97_2D/* /data/$USER/Output/Suttner97_2D_2/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOBID}.out /data/$USER/Output/Suttner97_2D_2/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
