#!/usr/bin/bash env 

################################################################################
# Annotate the piRNA reads
# 1. piRNA clusters (unique, unistrand: flam, cluster2, cluster11)
# 2. gene (unique, gene-body) 
# 3. TE (unique + multi; sense)
# 4. TE (unique + multi; anti)
#
## merge bed
#
################################################################################


################################################################################
## global variables
piRC_idx="${HOME}/data/genome/dm6/bowtie_index/dm6_piRNA_cluster"
# piRC_ex="piRC_exclude_uni.bed" # !! absolute path
piRC_ex="/data/yulab/wangming/work/devel_pipeline/pirna_seq/scripts/anno/piRC_exclude_uni.bed"
genome_idx="${HOME}/data/genome/dm6/bowtie_index/dm6"
te_idx="${HOME}/data/genome/dm6/bowtie_index/dm6_transposon"
gene_bed="${HOME}/data/genome/dm6/annotation_and_repeats/dm6.ensembl.bed"
cpu=8
merge_gap=1
################################################################################



################################################################################
## functions
function fx_name() {
    local fname=$(basename ${1/.gz})
    fname=${fname/.fast[aq]}
    fname=${fname/.f[aq]}
    fname=${fname/.not_piRC}
    fname=${fname/.not_genome}
    fname=${fname/.not_te}
    fname=${fname/.to_piRC}
    fname=${fname/.to_gene}
    fname=${fname/.to_te}
    fname=${fname/.fixed}
    fname=${fname/.bed}
    fname=${fname/_sens}
    fname=${fname/_anti}
    echo ${fname}
}
export -f fx_name


function count_bed() {
    local fname=$(fx_name $1)
    local n1=$(cat $1 | wc -l)
    local n2=$(cut -f 4 $1 | awk -F"-" '{s+=$2}END{print s}')
    echo "${fname},${n1},${n2}"
}
export -f count_bed


function count_fx() {
    local fname=$(fx_name $1)
    local n1=$(bioawk -cfastx '{print $name}' $1 | wc -l)
    local n2=$(bioawk -cfastx '{split($name, a, "-"); s+=a[2]}END{print s}' $1)
    echo "${fname},${n1},${n2}"
}
export -f count_fx


#################################################################################
## 1. to piRNA cluster
function to_piRC() {
    # remove reads on cluster2, cluster11, flam (antisense)
    # fq="input.fq.gz"
    local fq=$1
    local out_dir=$2 # "01.to_piRC"
    out_dir="${out_dir}/01.to_piRC"
    local fname=$(fx_name $fq)
    [[ -d ${out_dir} ]] || mkdir -p ${out_dir}
    # alignment
    local bam="${out_dir}/${fname}.to_piRC.bam"
    local unmap="${out_dir}/${fname}.not_piRC.fq"
    local log="${out_dir}/${fname}.to_piRC.bowtie.log"
    [[ -f ${bam} ]] ||
        (bowtie -x ${piRC_idx} -m 1 -v 2 -S -p ${cpu} --no-unal --un ${unmap} ${fq} \
        2> ${log} | \
        samtools view -bhS - | \
        samtools sort -o ${bam} - && \
        samtools index ${bam})
    # convert to bed and merge
    local bed="${out_dir}/${fname}.to_piRC.bed"
    local bed_final="${out_dir}/${fname}.to_piRC.final.bed"
    # remove cluster2, cluster11, flam (anti)
    [[ -f ${bed} ]] || \
        (bedtools intersect -v -s -a ${bam} -b ${piRC_ex} -bed > ${bed} && \
        bedtools intersect -u -s -a ${bam} -b ${piRC_ex} | \
        samtools fastq - >> ${unmap})
    # merge
    [[ -f ${bed_final} ]] || bedtools merge -i ${bed} -d ${merge_gap} -c 4,5,6 -o first -s > ${bed_final}
    # stat
    stat="${out_dir}/${fname}.to_piRC.stat"
    [[ -f ${stat} ]] || (s=$(count_bed ${bed}) && echo ${s},"piRC" > ${stat})
    # output
    echo ${unmap}
}


## 2. to gene (genome => gene)
function to_gene() {
    # fq="input.fq.gz"
    local fq=$1
    local out_dir=$2 # "02.to_gene"
    out_dir="${out_dir}/02.to_gene"
    local fname=$(fx_name $fq)
    [[ -d ${out_dir} ]] || mkdir -p ${out_dir}
    # alignment => genome
    local bam="${out_dir}/${fname}.to_genome.bam"
    local unmap="${out_dir}/${fname}.not_genome.fq"
    local log="${out_dir}/${fname}.to_genome.bowtie.log"
    [[ -f ${bam} ]] ||
        (bowtie -x ${genome_idx} -m 1 -v 2 -S -p ${cpu} --no-unal --un ${unmap} ${fq} \
        2> ${log} | \
        samtools view -bhS - | \
        samtools sort -o ${bam} - && \
        samtools index ${bam})
    # extract gene bed
    local bed="${out_dir}/${fname}.to_gene.bed"
    local bed_final="${out_dir}/${fname}.to_gene.final.bed"
    [[ -f ${bed} ]] || \
        (bedtools intersect -u -s -a ${bam} -b ${gene_bed} -bed > ${bed} && \
        bedtools intersect -v -s -a ${bam} -b ${gene_bed} | \
        samtools fastq - >> ${unmap})
    # merge
    [[ -f ${bed_final} ]] || bedtools merge -i ${bed} -d ${merge_gap} -c 4,5,6 -o first -s > ${bed_final}
    # stat
    stat="${out_dir}/${fname}.to_gene.stat"
    [[ -f ${stat} ]] || (s=$(count_bed ${bed}) && echo ${s},"gene" > ${stat})
    # output
    echo ${unmap}
}


## 3. to te_sens (unique + multi)
function to_te_sens() {
    local fq=$1
    local out_dir=$2 # "03.to_te_sens"
    out_dir="${out_dir}/03.to_te_sens"
    local fname=$(fx_name $fq)
    [[ -d ${out_dir} ]] || mkdir -p ${out_dir}
    # alignment => genome
    local bam="${out_dir}/${fname}.to_te.bam"
    local unmap="${out_dir}/${fname}.not_te.fq"
    local log="${out_dir}/${fname}.to_te.bowtie.log"
    [[ -f ${bam} ]] ||
        (bowtie -x ${te_idx} -k 1 --best -v 2 -S -p ${cpu} --no-unal --un ${unmap} ${fq} \
        2> ${log} | \
        samtools view -bhS - | \
        samtools sort -o ${bam} - && \
        samtools index ${bam})
    # extract gene bed
    local bed="${out_dir}/${fname}.to_te_sens.bed"
    local bed_final="${out_dir}/${fname}.to_te_sens.final.bed"
    [[ -f ${bed} ]] || \
        bedtools bamtobed -i <(samtools view -bhS -F 16 ${bam}) > ${bed} # forward strand
    # merge
    [[ -f ${bed_final} ]] || bedtools merge -i ${bed} -d ${merge_gap} -c 4,5,6 -o first -s > ${bed_final}
    # stat
    stat="${out_dir}/${fname}.to_te_sens.stat"
    [[ -f ${stat} ]] || (s=$(count_bed ${bed}) && echo ${s},"te_sens" > ${stat})
    # output
    echo ${unmap}
}


## 4. to te_anti (unique + multi)
function to_te_anti() {
    local fq=$1
    local out_dir=$2 # "04.to_te_anti"
    out_dir="${out_dir}/04.to_te_anti"
    local fname=$(fx_name $fq)
    [[ -d ${out_dir} ]] || mkdir -p ${out_dir}
    # alignment => genome
    local bam="${out_dir}/${fname}.to_te.bam"
    local unmap="${out_dir}/${fname}.not_te.fq"
    local log="${out_dir}/${fname}.to_te.bowtie.log"
    [[ -f ${bam} ]] ||
        (bowtie -x ${te_idx} -k 1 --best -v 2 -S -p ${cpu} --no-unal --un ${unmap} ${fq} \
        2> ${log} | \
        samtools view -bhS - | \
        samtools sort -o ${bam} - && \
        samtools index ${bam})
    # extract gene bed
    local bed="${out_dir}/${fname}.to_te_anti.bed"
    local bed_final="${out_dir}/${fname}.to_te_anti.final.bed"
    [[ -f ${bed} ]] || \
        bedtools bamtobed -i <(samtools view -bhS -f 16 ${bam}) > ${bed} # reverse strand
    # merge
    [[ -f ${bed_final} ]] || bedtools merge -i ${bed} -d ${merge_gap} -c 4,5,6 -o first -s > ${bed_final}
    # stat
    stat="${out_dir}/${fname}.to_te_anti.stat"
    [[ -f ${stat} ]] || (s=$(count_bed ${bed}) && echo ${s},"te_anti" > ${stat})
    # output
    echo ${unmap}
}


## 5. all-in-one
function main() {
    # input_fq="input.fq.gz"
    # out_dir="demo"
    input_fq=$1
    out_dir=$2
    fname=$(fx_name $input_fq)
    out_dir="${out_dir}/${fname}" # update out_dir
    
    echo "1. to piRNA clusters"
    not_piRC=$(to_piRC ${input_fq} ${out_dir})

    echo "2. to genome"
    not_gene=$(to_gene ${not_piRC} ${out_dir})

    echo "3. to te sens"
    not_te=$(to_te_sens ${not_gene} ${out_dir})
    
    echo "4. to te anti"
    not_te=$(to_te_anti ${not_gene} ${out_dir})

    ## 5. wrap stat
    s1=$(count_fx ${input_fq}) # input
    s2=$(count_fx ${not_te})   # unmap
    stat="${out_dir}/05.stat"
    echo ${s1},"input" > ${stat}
    cat ${out_dir}/*/*stat >> ${stat}
    echo ${s2},"unmap" >> ${stat}
}


[[ $# -lt 2 ]] && echo "Usage: anno.sh <in.fq> <out_dir>" && exit 1
main $1 $2

#