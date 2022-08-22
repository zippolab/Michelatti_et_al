#!/bin/bash

### Annotate ALL interactions per sample
BEDPE=$(ls FitHiChIP/ATAC/*.interactions_FitHiC_Q0.01.bedpe)

for f in $BEDPE; do
	echo $f
	NAME=$(basename $f _FitHiC_Q0.01.bedpe)
	echo $NAME
	/home/Programs/pgltools/sh/pgltools intersect1D -a $f -b GRCh37.p13.promoter.uniq.bed -wa | /home/Programs/pgltools/sh/pgltools merge -stdInA -o collapse -c 7 -noH > FitHiChIP_results/annotation/$NAME"_promoter.bedpe"
	/home/Programs/pgltools/sh/pgltools intersect1D -a $f -b regulatory.element.bed -wa | /home/Programs/pgltools/sh/pgltools merge -stdInA -o collapse -c 7 -noH > FitHiChIP_results/annotation/$NAME"_reg_el.bedpe"

	/home/Programs/pgltools/sh/pgltools intersect -a FitHiChIP_results/annotation/$NAME"_promoter.bedpe" -b FitHiChIP_results/annotation/$NAME"_reg_el.bedpe" -allA > FitHiChIP_results/annotation/$NAME"_annotation.bedpe"
	/home/Programs/pgltools/sh/pgltools intersect -a FitHiChIP_results/annotation/$NAME"_promoter.bedpe" -b FitHiChIP_results/annotation/$NAME"_reg_el.bedpe" -v > FitHiChIP_results/annotation/$NAME"_promoter_only.bedpe"
	/home/Programs/pgltools/sh/pgltools intersect -a FitHiChIP_results/annotation/$NAME"_reg_el.bedpe" -b FitHiChIP_results/annotation/$NAME"_promoter.bedpe" -v > FitHiChIP_results/annotation/$NAME"_reg_el_only.bedpe"
	
	cp $f FitHiChIP_results/comparison/
	mv FitHiChIP_results/comparison/$(basename $f) FitHiChIP_results/comparison/$NAME.bedpe
done


### Count different type of P-E interactions and save them into a report file
BEDPE=$(ls FitHiChIP_results/annotation/*_annotation.bedpe)

for f in $BEDPE; do
	NAME=$(basename $f _annotation.bedpe)
	echo $f >> FitHiChIP_results/annotation/report_annotation.txt
	echo "All interactions to be annotated: " $(wc -l  FitHiChIP_results/comparison/$NAME.bedpe | cut -d ' ' -f 1) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "All interactions overlapping P: " $(wc -l FitHiChIP_results/annotation/$NAME"_promoter.bedpe" | cut -d ' ' -f 1) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "All interactions overlapping E: " $(wc -l FitHiChIP_results/annotation/$NAME"_reg_el.bedpe" | cut -d ' ' -f 1) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "All annotated P-E interactions: " $(wc -l $f | cut -d ' ' -f 1) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "One to one P-E interactions: " $(awk '{if (($7=="A" && $8=="B") || ($7=="B" && $8=="A")) print $0}' $f | wc -l) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "One to many P-E interactions: " $(awk '{if (($7=="A" || $7=="B") && ($8 ~/,/)) print $0}' $f | wc -l) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "Many to one P-E interactions: " $(awk '{if (($7 ~/,/) && ($8=="A" || $8=="B")) print $0}' $f | wc -l) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "Many to Many P-E interactions: " $(awk '{if (($7 ~/,/) && ($8 ~/,/)) print $0}' $f | wc -l) >> FitHiChIP_results/annotation/report_annotation.txt
	echo "Single P-E in same bin: " $(awk '{if (($7=="A" && $8=="A") || ($7=="B" && $8=="B")) print $0}' $f | wc -l) >> FitHiChIP_results/annotation/report_annotation.txt
	echo >> FitHiChIP_results/annotation/report_annotation.txt
	echo >> FitHiChIP_results/annotation/report_annotation.txt
done
