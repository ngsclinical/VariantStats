if [ $3 = bwa ]
then
highCutoff="50%"
lowCutoff="5%"

bwa index -p $1 -a is ../$1

bwa aln -t 2 $1 $2_R1_001.fastq.gz -f $2_1.sai && bwa aln -t 2 $1 $2_R2_001.fastq.gz -f $2_2.sai

bwa sampe -r '@RG\tID:STS131\tSM:1886\tPL:ILLUMINA\n@PG\tID:BWA\tPN:BWA\tVN:0.5.9rc' $1 $2_1.sai $2_2.sai $2_R1_001.fastq.gz $2_R2_001.fastq.gz | samtools view -b -S -h -F 4 -q 5 - | samtools sort 
- $2_$3

samtools index $2_$3.bam

samtools faidx ../$1

java -Xmx15G -jar /usr/share/java/GenomeAnalysisTK.jar -T DepthOfCoverage -R $1 -baseCounts -I $2_$3.bam -dfrac 1 -dt ALL_READS -ct 10 -o $2_$3_BamStats.txt

perl ~/bin/CumulativeVarCall_v3.1.pl $2_$3_BamStats.txt $1 0 50 10 $2_$3 y

perl ~/bin/variantstatistics_v2.1.pl ${2}_${3}_"MINOR_VARIANT_GATK" 50 1 ${2}_${3}_{highCutoff} n

perl ~/bin/variantstatistics_v2.1.pl ${2}_${3}_"MINOR_VARIANT_GATK" 5 .1 ${2}_${3}_{lowCutoff} n

elif [ $3 = novo ]
then
highCutoff="50%"
lowCutoff="5%"

#ref=`echo $1 | cut -d\. -f 1`;
tmpref=`echo $1 | cut -d\/ -f 2`;ref=`echo $tmpref | cut -d\. -f 1`;echo $ref;

/usr/local/bin/novoalign -f $2_R1_001.fastq $2_R2_001.fastq -i PE 500,100 -c 4 -oSAM $'@RG\tID:${Plasmid}\tSM:${SRR054742}\tPU:ILLUMINA\tLB:PlasmidPrep' -d ../${ref} | samtools view -S -F 4 -b -h -
q 5 - | samtools sort - $2_$3 

samtools index $2_$3.bam

samtools faidx $1

java -Xmx15G -jar /usr/share/java/GenomeAnalysisTK.jar -T DepthOfCoverage -R $1 -baseCounts -I $2_$3.bam -dfrac 1 -dt ALL_READS -ct 10 -o $2_$3_BamStats.txt

perl ~/bin/CumulativeVarCall_v3.1.pl $2_$3_BamStats.txt $1 0 50 10 $2_$3 y

perl ~/bin/variantstatistics_v2.1.pl ${2}_${3}_"MINOR_VARIANT_GATK" 50 1 ${2}_${3}_{highCutoff} n

perl ~/bin/variantstatistics_v2.1.pl ${2}_${3}_"MINOR_VARIANT_GATK" 5 .1 ${2}_${3}_{lowCutoff} n

fi