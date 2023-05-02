#!/bin/sh


###Mauve progressive alignment
mauve="/xdisk/uschuch/corykeith/tools/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"
REF="/xdisk/uschuch/corykeith/ALPHA_SEARCH/refs/NC_023292.fasta"
CONTIGS="/xdisk/uschuch/corykeith/ALPHA_SEARCH/fasta/NODE_61_length_804_cov_308.707250.fasta"
MAUVEOUT="P003_WA05_NODE_61"
OUT_DIR="$MAUVEOUT"
CONTIG_FOLDER="/xdisk/uschuch/corykeith/ALPHA_SEARCH/fasta/"
CWD="/xdisk/uschuch/corykeith/ORIENTATION"

mkdir $OUT_DIR
$mauve --output ${OUT_DIR}/${MAUVEOUT}.xmfa --backbone-output ${OUT_DIR}/${MAUVEOUT}.backbone $REF $CONTIGS

## Grabs lines that have the contig coordinates for LCBs.
grep "NODE*" ${OUT_DIR}/${MAUVEOUT}.xmfa >> ${OUT_DIR}/${MAUVEOUT}_LCB.txt


count=0
## Takes lines of LCBs and strips the string to have just the numbers of the coordinates.
while IFS=$' ' read -r -a myArray; do
    if [[ ${myArray[0]} == ">" ]]; then
    count=$(( count + 1 ))
    #echo ${myArray[1]}
    STRIP="${myArray[1]#*:}"
    echo $STRIP strip
    FIRST="${STRIP%-*}"
    SECOND="${STRIP#*-}"
    CONTIG="NODE_61_length_804_cov_308.707250"
    echo $FIRST $SECOND 
    samtools faidx $CONTIGS $CONTIG:$FIRST-$SECOND >> ${OUT_DIR}/${MAUVEOUT}_LCB_${count}.fasta
    fi
done < ${OUT_DIR}/${MAUVEOUT}_LCB.txt

cat ${OUT_DIR}/*.fasta > ${OUT_DIR}/${MAUVEOUT}_combined.fasta

cd $OUT_DIR
#minimap2
SAM="${MAUVEOUT}.sam"
BAM="$MAUVEOUT"
MINIMAP="/xdisk/uschuch/corykeith/ORIENTATION/minimap2/minimap2"

$MINIMAP -ado $REF ../${OUT_DIR}/${MAUVEOUT}_combined.fasta > $SAM

# Convert sam to bam
samtools view -S -b $SAM > ${BAM}.bam

# Sort the alignment
samtools sort ${BAM}.bam -o ${BAM}_sorted.bam

# Get consensus fastq file
samtools mpileup -uf $REF ${BAM}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${BAM}_tmp_consensus.fastq


REF_NAME="NC_023292"
CONTIG_NAME="NODE_61"
#change header to name of contig instead of reference
cat ${BAM}_tmp_consensus.fastq | seqkit replace -p $REF_NAME -r $CONTIG_NAME > ${BAM}_consensus.fastq

rm ${BAM}_tmp_consensus.fastq
# Convert .fastq to .fasta
seqtk seq ${BAM}_consensus.fastq > SAMPLE1.fasta
