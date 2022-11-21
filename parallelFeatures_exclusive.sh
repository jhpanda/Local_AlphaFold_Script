#!/bin/sh
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -p hpc,odnl
#SBATCH -J preflist02
#SBATCH -o %x_%a.features.out
#SBATCH -e %x_%a.features.out
#SBATCH --array=81-88%5 ## pref02 not finished
##SBATCH --array=63-70%5 ## pref02 not finished


ID=${SLURM_ARRAY_TASK_ID}
#let "ID=ID+243684"
#pref=`head -$ID preflist02.txt | tail -1`
seqdir="./length500/"
#fasta="${seqdir}/${pref}.fasta"
outdir="${seqdir}"
model="multimer"

DB_PATH="/ru-auth/local/home/jpeng/scratch/softwares/alphafold_database"
bfd="${DB_PATH}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
uniclust30="${DB_PATH}/uniclust30/UniRef30_2022_02/UniRef30_2022_02"
wkdir=$(pwd)
echo "working dir $wkdir"

tmpdir="/tmp/jpeng/"
mkdir -p $tmpdir
rsync -aPv ${bfd}_cs219.ff{data,index} $tmpdir
cd $tmpdir
ln -s ${bfd}_a3m.ffdata
ln -s ${bfd}_a3m.ffindex
ln -s ${bfd}_hhm.ffdata
ln -s ${bfd}_hhm.ffindex
bfd_tmp="${tmpdir}/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

rsync -aPv ${uniclust30}_cs219.ff{data,index} $tmpdir
ln -s ${uniclust30}_a3m.ffdata
ln -s ${uniclust30}_a3m.ffindex
ln -s ${uniclust30}_hhm.ffdata
ln -s ${uniclust30}_hhm.ffindex
uniclust30_tmp="${tmpdir}/UniRef30_2022_02"

cd ${wkdir}

echo "SLURM JOB ID      : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "SLURM JOB ID      : ${SLURM_JOB_ID}"
echo "SLURM JOB NAME    : $pref"
echo "SLURM Node list   : ${SLURM_JOB_NODELIST}"
#echo ""
#if [ -f "${seqdir}/${pref}/features.pkl" ];then
#    echo "Finished $pref"
#else
#    echo "sh get_features.sh $fasta $model $outdir"
#fi
#preflist="preflist02_${ID}.txt"

preflist="${SLURM_JOB_NAME}_${ID}.txt"
#cat $preflist | parallel -j 6 echo {%} {}
cat $preflist | parallel -j 6 'sh get_features.sh {}'


# 1 ask for one node with many jobs (100 for example)
##SBATCH --exclusive ##
# try to maximize it to 2 days!
# 2 use gnu parallel instead of array (6 at a time for example: 6*4)
# before gnu parallel, rsync; after that, remove rsynced data
# 3 beyond that use slurm array (10 X 100 to 1000)

# 4. watch 'tail -10 output_500-1000.out ; ls -l'

rm -f $bfd_tmp $uniclust30_tmp
