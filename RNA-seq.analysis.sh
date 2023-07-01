###STAR RNA-seq pipeline
#0. fq trim 
trimmomatic PE -threads 10 ~/02.data/02.SRA/Pyrus_HBL/01.Genome_Raw_data/X101SC20124165-Z01-F003/Illumina/clean/FDSW210004050-1r_L2_1_clean.rd.fq.gz ~/02.data/02.SRA/Pyrus_HBL/01.Genome_Raw_data/X101SC20124165-Z01-F003/Illumina/clean/FDSW210004050-1r_L2_2_clean.rd.fq.gz  HBL_NGS2_1.fq.gz unPair1.fq.gz HBL_NGS2_2.fq.gz  unPair2.fq.gz ILLUMINACLIP:/localstorage/home/manyi/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:3:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50

#1. index构建
STAR --runMode genomeGenerate --genomeDir Pdan_hap1 --runThreadN 20 --genomeFastaFiles PdanHap1.fa --sjdbGTFfile PdanHap1.gtf --sjdbOverhang 49 --genomeSAindexNbases 13

2. STAR 比对
cat sample.list|while read id \ndo \necho STAR --genomeDir 00.index/ --runThreadN 20 --readFilesIn ~/01.WorkSapce/02.YH1H/01.trim 05.sh/$id'_Pair_1.fq.gz' ~/01.WorkSapce/02.YH1H/01.trim 05.sh/$id'_Pair_2.fq.gz'   --readFilesCommand zcat --outFileNamePrefix  ~/01.WorkSapce/02.YH1H/02.bam/$id --outSAMtype BAM SortedByCoordinate --alignIntronMax 20000 --alignMatesGapMax 25000 --outFilterMultimapNmax 1 --outBAMsortingThreadN 10 >>05.sh/STAR.sh\ndone

#3. stingtie 定量
cat sample.list |while read id \ndo \nstringtie -p 10 -G ../05.Slocus/01.JZYLRNA/01.index/Pdan_hap1.gtf  -A $id'AddSlocus.fpkm -o $id'AddSlocus.stringtie.gff' -e  ../02.bam/$id'Aligned.sortedByCoord.out.bam'\ndone

#4. 表达矩阵构建
prepDE.py -i gff.list -g gene_count_matrix.csv -t transcript_count_matrix.csv