#!/bin/sh


OUT_DIR="/xdisk/uschuch/corykeith/ViCAT_P003/out"
NCBI_REF="/xdisk/uschuch/corykeith/BLAST/NCBI_virus/viraldb_1line.fsa"
mauve="/xdisk/uschuch/corykeith/tools/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"

while IFS=$'\t' read -r -a myArray; do
     #Variables for each contig in list.
     SMPLE="${myArray[0]%_filtered_scaffolds.fasta}"
     NODE_NAME="${myArray[1]%_length*}"
     NODE="${myArray[1]}"
     ACCESSION="${myArray[2]}"

     ##Informative sample info
     suff="${SMPLE#*146201_*}"
     pre="${suff%_i5*}"

     ## Make directory for mauve
     rm -rf ${OUT_DIR}/${SMPLE}/mauve/$NODE_NAME
     mkdir -p ${OUT_DIR}/${SMPLE}/mauve/$NODE_NAME
     ## Pull accession number from viraldb_1line.fsa
     grep -A 1 $ACCESSION $NCBI_REF > ${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${ACCESSION}.fasta

     #Variables for gather_contigs.py
     FILE="${OUT_DIR}/${SMPLE}/spades_filtered/${myArray[0]}"
     echo $NODE > contig.txt
     CONTIG="contig.txt"
     OUT="${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${pre}_${NODE_NAME}.fasta"
     echo $OUT out $FILE file $CONTIG contig > ${OUT_DIR}/${SMPLE}/mauve/out_name.txt
     python gather_contigs.py $FILE $CONTIG $OUT

     ## Variables for mauve progressive alignment
     MAUVEOUT="${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${NODE_NAME}"
     REF="${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${ACCESSION}.fasta"
     CONTIGS="$OUT"
     ## Mauve whole genome alignment
     $mauve --output ${MAUVEOUT}.xmfa --backbone-output ${MAUVEOUT}.backbone $REF $CONTIGS

     ## Grabs lines of xmfa to get local collinear blocks (LCBs) coordinates
     grep ">\ 2:*" ${MAUVEOUT}.xmfa >> ${MAUVEOUT}_LCB.txt
     ## Extract LCBs
     count=0
     while IFS=$' ' read -r -a myArray; do
         if [[ ${myArray[0]} == ">" ]]; then
         count=$(( count + 1 ))
         #echo ${myArray[1]}
         STRIP="${myArray[1]#*:}"
         #echo $STRIP strip
         FIRST="${STRIP%-*}"
         SECOND="${STRIP#*-}"
         ## extract LCBs from contigs
         samtools faidx $CONTIGS $NODE:$FIRST-$SECOND >> ${MAUVEOUT}_LCB_${count}.fasta
         fi
     done < ${MAUVEOUT}_LCB.txt
     ## Combine all blocks to multiline fasta
     cat ${MAUVEOUT}_LCB_* > ${MAUVEOUT}_combined.fasta
     #minimap variables and  minimap mapping of LCBs     
     SAM="${MAUVEOUT}.sam"
     BAM="$MAUVEOUT"
     MINIMAP="/xdisk/uschuch/corykeith/ORIENTATION/minimap2/minimap2"
     $MINIMAP -ado $REF ${MAUVEOUT}_combined.fasta > $SAM
     # Convert sam to bam
     samtools view -S -b $SAM > ${BAM}.bam
     # Sort the alignment
     samtools sort ${BAM}.bam -o ${BAM}_sorted.bam
     # Get consensus fastq file
     samtools mpileup -uf $REF ${BAM}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${BAM}_tmp_consensus.fastq
     #change header to name of contig instead of reference
     cat ${BAM}_tmp_consensus.fastq | seqkit replace -p $ACCESSION -r $NODE_NAME > ${BAM}_consensus.fastq
     rm ${BAM}_tmp_consensus.fastq
     # Convert .fastq to .fasta
     seqtk seq -a ${BAM}_consensus.fastq > ${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${pre}_${NODE_NAME}_oriented.fasta
     cp ${OUT_DIR}/${SMPLE}/mauve/${NODE_NAME}/${pre}_${NODE_NAME}_oriented.fasta /xdisk/uschuch/corykeith/ViCAT_P003/out/oriented_alphas
done < sortedalphalist.txt

