#!/bin/sh

gene=$1
Scripts_dir=$2
Res_dir=$3

wkdir=${Res_dir}/${gene}_GREX

Rscript ${Scripts_dir}/compute_grex.R ${gene} ${wkdir}
