#!/bin/bash
#PBS -S /bin/bash
#PBS -N fun
#PBS -l nodes=1:ppn=40
#PBS -e el.funann.err
#PBS -o el.funann.log
#PBS -q Mcu
#PBS -V

##### env
conda_activate="/public/zhis/miniconda3/bin/activate"
env="funannotate"

source "${conda_activate}" "${env}"

##### par setting
homedir="/public/zhis/all/026.230705_allergen_exploration/03.HDHGenome_anno.250409"
cd ${homedir}

datadir="${homedir}/data"
resdir="${homedir}/result"
#mkdir -p "${resdir}"

thread="40"

Genome_fa="${datadir}/HDH_Genome.fa"
Genome_pref="HDH_Genome"

BUSCO_DB="/public/zhis/ref/BUSCO_ref_database/mollusca_odb10"

##### 1. Genome preparation

resdir1_prep="${resdir}/05.Genome_prep"
mkdir -p "${resdir1_prep}"

### 1.1. Cleaning
cd "${resdir1_prep}"

cl_fa="${Genome_pref}.clean.fa"
cl_st_fa="${Genome_pref}.cl.sort.fa"
cl_st_mk_fa="${Genome_pref}.cl.st.mk.fa"

  funannotate clean -i "${Genome_fa}" -o "${cl_fa}" -p 95 -c 95 -m 500

### 1.2. Sorting & Rename
  funannotate sort -i "${cl_fa}" -o "${cl_st_fa}" -b "scaffold" --minlen 0

### 1.3. RepeatMasking
dbname="${Genome_pref}"

  BuildDatabase -name "${dbname}" "${cl_st_fa}"

  RepeatModeler -database "${dbname}" -threads "${thread}" -LTRStruct

  funannotate mask \
	-i "${cl_st_fa}" \
	-o "${cl_st_mk_fa}" \
	-m repeatmasker \
	-l "${dbname}-families.fa" \
	--cpus ${thread} \
	--debug


cleanfq_loc="${resdir}/01.Trimming"
cd ${cleanfq_loc}

##### 2. Funannotate
### 2.1. train
  funannotate train \
	--stranded no \
	-i "${resdir1_prep}/${cl_st_mk_fa}" \
	-o "${resdir1_prep}/fun" \
	--left mRNA1_1_1_val_1.fq.gz mRNA1_2_1_val_1.fq.gz mRNA1_3_1_val_1.fq.gz mRNA1_4_1_val_1.fq.gz mRNA2_1_1_val_1.fq.gz mRNA2_2_1_val_1.fq.gz mRNA2_3_1_val_1.fq.gz mRNA2_4_1_val_1.fq.gz mRNA3_1_1_val_1.fq.gz mRNA3_2_1_val_1.fq.gz mRNA3_3_1_val_1.fq.gz mRNA3_4_1_val_1.fq.gz \
	--right mRNA1_1_2_val_2.fq.gz mRNA1_2_2_val_2.fq.gz mRNA1_3_2_val_2.fq.gz mRNA1_4_2_val_2.fq.gz mRNA2_1_2_val_2.fq.gz mRNA2_2_2_val_2.fq.gz mRNA2_3_2_val_2.fq.gz mRNA2_4_2_val_2.fq.gz mRNA3_1_2_val_2.fq.gz mRNA3_2_2_val_2.fq.gz mRNA3_3_2_val_2.fq.gz mRNA3_4_2_val_2.fq.gz \
	--no_trimmomatic \
	--pasa_db sqlite \
	--max_intronlen 500000 \
	--species "Haliotis discus hannai" \
	--cpus ${thread}

tmpdir="${resdir1_prep}/fun/tmppp"
mkdir -p "${tmpdir}"

## Trans/Prot evidence:
evidence_dir="${datadir}/genome_pred_material"
est="${evidence_dir}/total.ESTs_hal.fa"
uni_close="${evidence_dir}/Uniprot_closely_species.fa"
prot1="${evidence_dir}/protein_HalAsinina.fa"

### 2.2. predict
  funannotate predict \
	--organism other \
	--repeats2evm \
	--busco_db "${BUSCO_DB}" \
	-i "${resdir1_prep}/${cl_st_mk_fa}" \
	-o "${resdir1_prep}/fun" \
	-s "Haliotis discus hannai" \
	--transcript_evidence "${resdir}/02.DenovoAssembly/Result.Trinity.fasta" "${est}" \
	--protein_evidence "${FUNANNOTATE_DB}/uniprot_sprot.fasta" "${uni_close}" "${prot1}" \
	--max_intronlen 500000 \
	--keep_evm \
	--genemark_mode ES \
	--tmpdir "${tmpdir}" \
	--cpus ${thread}

## some errors that required to be manally fixed..
predict_resdir="${resdir1_prep}/fun/predict_results"

  cat "${predict_resdir}/Haliotis_discus_hannai.models-need-fixing.txt" | sed 1d | cut -f 1 > ${predict_resdir}/fun.gene_model.rm.txt
  funannotate fix \
	-d ${predict_resdir}/fun.gene_model.rm.txt \
	-i ${predict_resdir}/Haliotis_discus_hannai.gbk \
	-t ${predict_resdir}/Haliotis_discus_hannai.tbl \
	-o ${predict_resdir}


### 2.3. update
  funannotate update \
	-i "${resdir1_prep}/fun" \
	--max_intronlen 500000 \
	--species "Haliotis discus hannai" \
	--cpus ${thread}

## manually fix, same as previous...
update_resdir="${resdir1_prep}/fun/update_results"

  cat "${update_resdir}/Haliotis_discus_hannai.models-need-fixing.txt" | sed 1d | cut -f 1 > ${update_resdir}/fun.gene_model.rm.txt
  funannotate fix \
        -d ${update_resdir}/fun.gene_model.rm.txt \
        -i ${update_resdir}/Haliotis_discus_hannai.gbk \
        -t ${update_resdir}/Haliotis_discus_hannai.tbl \
        -o ${update_resdir}

### 2.4. annotation
IPRSCAN="/public/zhis/software/interproscan/interproscan-5.73-104.0/interproscan.sh"

  funannotate iprscan -i "${resdir1_prep}/fun" -m local --iprscan_path ${IPRSCAN} -c ${thread}
  funannotate annotate -i "${resdir1_prep}/fun" --cpus ${thread} --busco_db "${BUSCO_DB}"
