#!/bin/sh
#SBATCH --nodes=1
#SBATCH -J dup99B_spr
#SBATCH --exclusive
#SBATCH --output=%x.out
#SBATCH --error=%x.out
#SBATCH -p hpc_k80,odnl_k80
##SBATCH --array=7501-7800%4

#export XLA_FLAGS="--xla_gpu_force_compilation_parallelism=1"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/programs/x86_64-linux/cuda/10.2.89/lib64
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
ID=${SLURM_ARRAY_TASK_ID}

echo "SLURM JOB ID      : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "SLURM JOB ID      : ${SLURM_JOB_ID}"
echo "SLURM JOB NAME    : $pref"
echo "SLURM Node list   : ${SLURM_JOB_NODELIST}"

model="multimer"
seqdir="./"
outdir="${seqdir}"
#preflist='preflist01.txt'
#pref=`head -$ID preflist.txt | tail -1`
pref="${SLURM_JOB_NAME}"
out="${seqdir}/${pref}.out"
fasta="${seqdir}/${pref}.fasta"

json="${seqdir}/${pref}/ranking_debug.json"
mpkl="${seqdir}/${pref}/result_model_4_multimer_v2_pred_0.pkl"
fpkl="${seqdir}/${pref}/features.pkl"
#if [ -f $json ];then
if [ -f $mpkl ] ;then
    echo "Prediction for ${pref} finished , skipping... "
else
    echo "Prediction for ${pref} on ${SLURM_JOB_NODELIST}"
    sh alphafold_fast.sh $fasta $model $outdir
fi
