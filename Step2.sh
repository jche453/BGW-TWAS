#!/bin/sh
gene=$1
geneFile=$2
Res_dir=$3
p_thresh=$4
max_blocks=$5


echo ${gene}
# Score statistics directory

Score_dir=${Res_dir}/${gene}_scores
cd ${Score_dir}

for fname in $(ls *.score*)
do
zcat $fname | awk -v var=$fname 'NR == 1 || $13 < min {line = $0; min = $13} END{print var, min}' >> ${gene}_ranked_segments.txt
done

# get gene info to provide to the makefile arguments

target_chr=$(grep ${gene} ${geneFile} | awk 'FS {print $1}'); echo $target_chr
start_pos=$(grep ${gene} ${geneFile} | awk 'FS {print $2}'); echo $start_pos
end_pos=$(grep ${gene} ${geneFile} | awk 'FS {print $3}'); echo $end_pos


## for gene_ranked_segs, FS=_.,
cat ${gene}_ranked_segments.txt | while read line ; do
chr=$(echo $line | awk -F[:._/] '{print $3}' )
start=$(echo $line | awk -F[:._/] '{print $4}' )
end=$(echo $line | awk -F[:._/] '{print $5}' )
study=$(echo $line | awk -F[:._/] '{print $2}' )
Name=$(echo $line | awk -F[:._/] '{print $1}' )
pval=$(echo $line | awk -F " " '{print $2}' )
if ((( $chr == ${target_chr})) && (($end > ($start_pos - 1000000))) && (($start < ($start_pos - 1000000) ))); then
printf "${Name}_${study}_${chr}_${start}_${end}\n" >> ${gene}_signif_segments.txt
elif ((( $chr == ${target_chr})) && (($start < ($end_pos + 1000000))) && (($end > ($end_pos + 1000000)))); then
printf "${Name}_${study}_${chr}_${start}_${end}\n" >> ${gene}_signif_segments.txt
elif (($(echo "$pval < $p_thresh" | bc -l ))) ; then
printf "${Name}_${study}_${chr}_${start}_${end}\n" >> ${gene}_signif_segments.txt
else
continue;
fi
done
sort -g -k 2 -u ${gene}_signif_segments.txt | head -${num_blocks} > ${gene}_signif_segmentstmp.txt && mv ${gene}_signif_segmentstmp.txt ${gene}_signif_segments.txt
filehead=${gene}_signif_segments.txt

echo Step 2 complete!

exit
