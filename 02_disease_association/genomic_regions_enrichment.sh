#!/usr/bin/env bash

set -e -u -o pipefail


bedtools sort -i ${input.query_region} | bedtools merge | bedtools coverage -a ${input.target_region} -b /dev/stdin | sort -k 4,4 > ${output.observed_tsv}


bedtools groupby -i ${input.observed_tsv} -g 4 -c 6,7 -o sum,sum | awk -v cov=${wildcards.cov} '$2/$3>cov' | cut -f1 | sort -u | wc -l > ${output.observed_lst}


mkdir -p ${RESULTS_DIR}/temp/${wildcards.exp_id}/shuffled

for num in $(seq 0 $(({REP_TIMES}-1)))
do
    bedtools sort -i ${params.shuffled_dir}/APG_HPRC_HGSVC.${params.query_id}.shuffle_${num}.bed | bedtools merge | bedtools coverage -a ${input.target_region} -b /dev/stdin | sort -k 4,4 > ${RESULTS_DIR}/temp/${wildcards.exp_id}/shuffled/${wildcards.exp_id}.No${num}.intersect_shuffled.tsv
done


for num in $(seq 0 $(({REP_TIMES}-1)))
do
    bedtools groupby -i {RESULTS_DIR}/temp/{wildcards.exp_id}/shuffled/{wildcards.exp_id}.No${{num}}.intersect_shuffled.tsv -g 4 -c 6,7 -o sum,sum | awk -v cov={wildcards.cov} '$2/$3>cov' | cut -f1 | sort -u | wc -l
done > ${output.shuffled_lst}
