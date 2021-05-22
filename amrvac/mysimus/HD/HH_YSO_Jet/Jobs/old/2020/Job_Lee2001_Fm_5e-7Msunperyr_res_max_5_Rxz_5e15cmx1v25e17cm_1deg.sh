#!/bin/bash
#SBATCH --job-name=Lee2001_Fm_5e-7Msunperyr_res_max_5_Rxz_5e15cmx1v25e17cm_1deg
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --time=2880
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --mem=16gb
#SBATCH --tmp=25gb

## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${HOSTNAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir Lee2001_res_max_5_R=5e15cm_1deg

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/Lee2001_res_max_5_R=5e15cm_1deg
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/tests/Lee2001_Fm_5e-7Msunperyr_res_max_5_Rxz\=5e15cmx1v25e17cm_1deg.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par  
mpiexec amrvac -i Lee2001_Fm_5e-7Msunperyr_res_max_5_Rxz\=5e15cmx1v25e17cm_1deg.par 
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/Lee2001_res_max_5_R=5e15cm_1deg/* /data/$USER/Output/Lee2001_res_max_5_R=5e15cm_1deg/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-${SLURM_JOBID}.out /data/$USER/Lee2001_res_max_5_R=5e15cm_1deg/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0

