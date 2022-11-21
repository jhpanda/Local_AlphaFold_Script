#!/bin/sh
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -J mip1_Anc3
#SBATCH --output=%x.out
#SBATCH --error=%x.out
#SBATCH -p hpc_a100

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

conda activate cuda11

#export XLA_FLAGS="--xla_gpu_force_compilation_parallelism=1"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/programs/x86_64-linux/cuda/10.2.89/lib64
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
ID=${SLURM_ARRAY_TASK_ID}

echo "SLURM JOB ID      : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "SLURM JOB ID      : ${SLURM_JOB_ID}"
echo "SLURM JOB NAME    : $pref"
echo "SLURM Node list   : ${SLURM_JOB_NODELIST}"

model="multimer"
seqdir="./"
outdir="${seqdir}"
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
conda deactivate
