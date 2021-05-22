#!/bin/bash
#SBATCH --job-name=Jet_Zhang
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --time=15
#SBATCH --partition=short
#SBATCH --mail-user=mialy.rabenanahary@obspm.fr
#SBATCH --mem=4gb
#SBATCH --tmp=120gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run
cd $SCRATCH

srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/Zhang2019/Parfiles/zheng19.par .
##srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/Zhang2019/JetI_1_1697.dat .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/Zhang2019/amrvac .
echo ${SLURM_JOBID}

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par  
mpiexec amrvac -i zheng19.par 
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/* /data/$USER/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0


