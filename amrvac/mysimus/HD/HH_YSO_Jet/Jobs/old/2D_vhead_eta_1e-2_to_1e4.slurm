#!/bin/bash
#SBATCH --job-name=2D_vHead_eta_1e-2_to_1e4
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=4320
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=slurm-2D_vHead_eta_1e-2_to_1e4-%A_%a.out
#SBATCH --mem=32gb
#SBATCH --tmp=100gb
#SBATCH --array=0-6%3

initarray=("1" "1e-2" "1e-1" "1e1" "1e2" "1e3" "1e4")
#initarray = ("1" "1e-2" "1e-1" "1e1" "1e2" "1e3" "1e4")
## Définir le répertoire scratch et recopier les fichiers nécessaires à l'exécution
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.$SLURM_ARRAY_TASK_ID.${HOSTNAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir 2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p Output/2D_vHead
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Parfiles/Article/mi_these/2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/amrvac .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more 2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}.par
mpiexec amrvac -i 2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}.par
#srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/2D_vHead/Jet_Shang06_0100.dat /data/$USER/Output/2D_vHead/.
srun --ntasks=$SLURM_JOB_NUM_NODES mv Output/2D_vHead/* /data/$USER/Output/2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp $AMRVAC_DIR/mysimus/HD/HH_YSO_Jet/Jobs/slurm-2D_vHead_eta_1e-2_to_1e4-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /data/$USER/Output/2D_vHead_eta_${initarray[$SLURM_ARRAY_TASK_ID]}/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
