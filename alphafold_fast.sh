#!/bin/sh
#SBATCH -N 1
#SBATCH --exclusive
##SBATCH -p hpc_k80,hpc,bigmem
##SBATCH -p hpc,bigmem
DB_PATH="/ru-auth/local/home/jpeng/scratch/softwares/alphafold_database"

fasta_paths="$1"
model="$2"
output_dir="$3"
msa_flag="0"
if [ -z "$output_dir" ];then
    output_dir="${wkdir}"
fi
#models="model_1,model_2,model_3,model_4,model_5"

if [ "$model" != "monomer" ] && [ "$model" != "multimer" ];then
    echo "Usage: sh alphafold.sh [sequence_name] [monomer|multimer]"
    exit "Invalid Option: $model"
fi

data_dir="$DB_PATH"
uniref90="${DB_PATH}/uniref90/uniref90.fasta"
mgnify="${DB_PATH}/mgnify/mgy_clusters.fa"
bfd="${DB_PATH}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
uniclust30="${DB_PATH}/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
pdb70="${DB_PATH}/pdb70/pdb70"
mmcif="${DB_PATH}/pdb_mmcif/mmcif_files"
pdb_seqres="${DB_PATH}/pdb_seqres/pdb_seqres.txt"
uniprot="${DB_PATH}/uniprot/uniprot_trembl.fasta"
obsolete_pdbs="${DB_PATH}/pdb_mmcif/obsolete.dat"
max_template_date="2021-11-01"
db_preset="full_dbs"


binary_path="/ru-auth/local/home/jpeng/scratch/softwares/"
hhblits="${binary_path}/hh-suite/bin/hhblits"
hhsearch="${binary_path}/hh-suite/bin/hhsearch"
hmmsearch="${binary_path}/hmmer/bin/hmmsearch"
hmmbuild="${binary_path}/hmmer/bin/hmmbuild"
jackhmmer="${binary_path}/hmmer/bin/jackhmmer"
kalign="/ru-auth/local/home/jpeng/scratch/miniconda3/bin/kalign"

## special tricks to bfd and uniclust30 ##
## https://github.com/soedinglab/hh-suite/issues/281#issuecomment-888689484 ##
echo "constructing temporary bfd and uniclust database"
bfd="${DB_PATH}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
uniclust30="${DB_PATH}/uniclust30/UniRef30_2022_02/UniRef30_2022_02"
tmpdir=/tmp/jpeng/
wkdir=`pwd`
mkdir -p $tmpdir
rsync -aPv ${bfd}_cs219.ff{data,index} $tmpdir
rsync -aPv ${uniclust30}_cs219.ff{data,index} $tmpdir
cd $tmpdir
ln -s ${bfd}_a3m.ffdata
ln -s ${bfd}_a3m.ffindex
ln -s ${bfd}_hhm.ffdata
ln -s ${bfd}_hhm.ffindex
ln -s ${uniclust30}_a3m.ffdata
ln -s ${uniclust30}_a3m.ffindex
ln -s ${uniclust30}_hhm.ffdata
ln -s ${uniclust30}_hhm.ffindex
bfd_tmp="${tmpdir}/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
uniclust30_tmp="${tmpdir}/UniRef30_2022_02"
cd $wkdir
echo "Done!"

echo "Now start prediction"
db_preset="full_dbs"

#af_script="${binary_path}/alphafold/run_alphafold_fast.py"
af_script="${binary_path}/alphafold/run_alphafold.py"
#af_script="${binary_path}/alphafold/run_from_msa.py"
if [ "$model" = "monomer" ];then
    model_preset="monomer_ptm"
    python $af_script \
        --fasta_paths=$fasta_paths \
        --output_dir=$output_dir \
        --use_precomputed_msas=$msa_flag \
        --data_dir=$data_dir \
        --uniref90_database_path=$uniref90 \
        --mgnify_database_path=$mgnify \
        --uniclust30_database_path=$uniclust30 \
        --bfd_database_path=$bfd \
        --pdb70_database_path=$pdb70 \
        --template_mmcif_dir=$mmcif \
        --max_template_date=$max_template_date \
        --obsolete_pdbs_path=$obsolete_pdbs \
        --model_preset=$model_preset \
        --db_preset=$db_preset \
        --hhblits_binary_path=$hhblits \
        --hhsearch_binary_path=$hhsearch \
        --hmmsearch_binary_path=$hmmsearch \
        --hmmbuild_binary_path=$hmmbuild \
        --jackhmmer_binary_path=$jackhmmer \
        --kalign_binary_path=$kalign \
        --use_gpu_relax=0

elif [ "$model" = "multimer" ];then
    model_preset="multimer"
    echo "$model_preset"
    python $af_script \
        --fasta_paths=$fasta_paths \
        --use_precomputed_msas=${msa_flag} \
        --output_dir=$output_dir \
        --data_dir=$data_dir \
        --uniref90_database_path=$uniref90 \
        --mgnify_database_path=$mgnify \
        --uniclust30_database_path=$uniclust30 \
        --bfd_database_path=$bfd \
        --pdb_seqres_database_path=$pdb_seqres \
        --template_mmcif_dir=$mmcif \
        --max_template_date=$max_template_date \
        --obsolete_pdbs_path=$obsolete_pdbs \
        --uniprot_database_path=$uniprot \
        --model_preset=$model_preset \
        --db_preset=$db_preset \
        --hhblits_binary_path=$hhblits \
        --hhsearch_binary_path=$hhsearch \
        --hmmsearch_binary_path=$hmmsearch \
        --hmmbuild_binary_path=$hmmbuild \
        --jackhmmer_binary_path=$jackhmmer \
        --kalign_binary_path=$kalign \
        --num_multimer_predictions_per_model=5 \
        --use_gpu_relax=0
fi
