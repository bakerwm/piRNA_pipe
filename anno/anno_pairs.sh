#!/usr/bin/bash env 

################################################################################
# Annotate the piRNA reads
# check read counts, overlapped between groups
# query + subject
################################################################################



################################################################################
# global variables
ovl_frac=1
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


function get_groups() {
    # return the group list
    # piRC, gene, te_sens, te_anti, ...
    # format: 01.to_gene, ...
    local dir=$1
    [[ $(is_anno_dir ${dir}) -gt 0 ]] && echo "not anno_dir: ${dir}" && exit 1
    ls ${dir} | while read i
    do
        [[ ${i} == 0[0-9].to_* ]] || continue # format
        g=${i/0[1-9].to_}
        echo ${g}
    done
}
export -f get_groups


function is_anno_dir() {
    # return code
    # true: 0; false: 1
    # required files:
    # (.bed, .final.bed)
    g_list=(piRC gene te_sens te_anti)
    b1=()
    b2=()
    for g in ${g_list[@]}
    do
        b1+=($(ls ${1}/0[1-9].*/*.to_${g}.bed))
        b2+=($(ls ${1}/0[1-9].*/*.to_${g}.final.bed))
    done
    [[ ${#b1[@]} -eq ${#g_list[@]} && ${#b2[@]} -eq ${#g_list[@]} ]] && echo 0 || echo 1
}
export -f is_anno_dir



function cmp_anno() {
    # compare:
    # query:.bed vs subject:.final.bed
    local q_dir=$1
    local s_dir=$2
    # local g=$3 # piRC/gene/te_sens/te_anti
    [[ $(is_anno_dir ${q_dir}) -gt 0 ]] && echo "not anno_dir: ${q_dir}" && exit 1
    [[ $(is_anno_dir ${s_dir}) -gt 0 ]] && echo "not anno_dir: ${s_dir}" && exit 1
    # output files
    local sname=$(basename ${s_dir})
    local ovl_dir="${q_dir}/overlap/${sname}"
    # check group list piRC/gene/te_sens/te_anti
    g_list=($(get_groups ${q_dir}))
    for g in ${g_list[@]}
    do
        # list files
        local query=$(ls ${q_dir}/*/*to_${g}.bed)
        local subject=$(ls ${s_dir}/*/*to_${g}.final.bed) # final
        [[ -f ${query} && -f ${subject} ]] || (echo "file not found" && exit 1)
        # output dir
        local sub_dir="${ovl_dir}/$(basename $(dirname ${query}))"
        local bed="${sub_dir}/$(basename ${query})"
        [[ -d ${sub_dir} ]] || mkdir -p ${sub_dir}
        [[ -f ${bed} ]] || bedtools intersect -f ${ovl_frac} -s -u -a ${query} -b ${subject} > ${bed}
        local stat=${bed/.bed/.stat}
        [[ -f ${stat} ]] || (s=$(count_bed ${bed}) && echo ${s},${g} > ${stat})
    done
    # summary
    local stat="${ovl_dir}/05.stat"
    [[ -f ${stat} ]] || cat ${ovl_dir}/*/*stat > ${stat}
}
export -f cmp_anno



# q="/data/yulab/wangming/work/yu_2021/piwi_lxh/results/hiseq/smRNA/piRNA_locus_classification/results/shPiwi1_3_8h"
# q="/data/yulab/wangming/work/yu_2021/piwi_lxh/results/hiseq/smRNA/piRNA_locus_classification/results"
# s="/data/yulab/wangming/work/yu_2021/piwi_lxh/results/hiseq/smRNA/piRNA_locus_classification/results_v2/Piwi_IP_3_8h"
# q="/data/yulab/wangming/work/yu_2021/piwi_lxh/results/hiseq/smRNA/piRNA_locus_classification/results_v2/total_3_8h"
# q_list=($(ls ${q}))
# echo ${q_list[@]}
# parallel cmp_subdir "${q}"/{1} ${s} {2} ::: ${q_list[@]} ::: piRC gene te
# cmp_anno ${q} ${s}

[[ $# -lt 2 ]] && echo "Usage: bash anno_pair.sh <query> <subject>" && exit 1
cmp_anno $1 $2

#