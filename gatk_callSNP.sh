cat AddSample.list| while read id ;do
spath=/gss1/home/wujun01/data/ColorGWAS
 java -Xmx4G -jar /gss1/home/wujun01/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE  -threads 10  $spath/01.fq/$id/$id'_1.fq.gz' $spath/01.fq/$id/$id'_2.fq.gz' $spath/01.fq/$id/$id'_clean_1.fastq.gz' un1.fq.gz $spath/01.fq/$id/$id'_clean_2.fastq.gz' un2.fq.gz ILLUMINACLIP:/gss1/home/wujun01/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 
-e bwa mem -t 10 -k 32 -M -R "'"@RG'\\t'ID:$id'\\t'SM:$id"'"  $spath/00.index/Pdanv2 $spath/01.fq/$id/$id'_clean_1.fastq.gz' $spath/01.fq/$id/$id'_clean_2.fastq.gz' -o $spath/03.bam/$id.sam
samtools view -F 4 -o $spath/03.bam/$id'.bam' -b -S  -q 1 $spath/03.bam/$id'.sam'
java -Xms1g -Xmx8g -XX:ParallelGCThreads=2 -jar ~/bin/picard.jar SortSam I=$spath/03.bam/$id'.bam' O=$spath/03.bam/$id'.sorted.bam' SORT_ORDER=coordinate 
java -Xms1g -Xmx8g -XX:P:warallelGCThreads=2 -jar ~/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=$spath/03.bam/$id'.sorted.bam' OUTPUT=$spath/03.bam/$id'.repeatmark.bam' METRICS_FILE=$spath/03.bam/$id'.bam.metrics'
samtools index  $spath/03.bam/$id'.repeatmark.bam'>>05.sh/$id'.sh'
rm $spath/03.bam/$id'.sorted.bam' $spath/03.bam/$id'.bam' >> 05.sh/$id'.sh'
rm $spath/03.bam/$id'.sam'
gatk --java-options '-Xms2g -Xmx2g -XX:ParallelGCThreads=1' HaplotypeCaller -I $spath/03.bam/$id'.repeatmark.bam' -ERC GVCF -R $spath/00.index/Pdan  -O $spath/04.vcf/$id.g.vcf
done