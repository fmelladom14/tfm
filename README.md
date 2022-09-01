# Transcriptomic analysis of disease-resistance response in Cucumis melo L.

<div style="text-align: justify"> In this repository, code generated to perform the analysis is included. It allows you to perform RNA-seq analysis, GO and KEGG enrichment,  and create a co-expression network using the R package igraph. </div>

## Example of analysis

<div style="text-align: justify"> This part of the analysis can be integrated into job submission scripts, where you can specify more options. Then, using the `sbatch` command, _SLURM_ runs your job when its your turn. Yet most commands are simple enough to use the command line. First, you need to download data, using _fastq-dump_. Yoo could do a simple loop to download all of the samples at once. </div>

`$ fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir path/reads/ i`

<div style="text-align: justify"> The next step is quality control._FastQC_ will accept multiple file names as input, so we can use the *.fq or *fastq.gz  wildcard in the same directory where your samples are located. </div>

`$ fastqc *.fastq.gz`

<div style="text-align: justify"> If your reads have good enough quality you can skip the next part of the process where you make sure that all the reads in the dataset have a chance to map/align to the genome, trimming unwanted information off from every read. We will not be performing this step because data did not have an appreciable amount of leftover adapter sequences or other contaminating sequences based on FastQC. Then, you have to align your reads to your reference genome. For this task we used HISAT2. Before the alignment, you have to run `hisat2-build`, which builds a HISAT2 index from the reference genome in fasta format. You can retreive this data from ENSEMBL. </div>

`$ hisat2-build -f DHL92_genome_v4.fa index`

<div style="text-align: justify"> Staying in the same directory, run the following command to align **paired reads**: </div>

`$ hisat2 -f -x $path/to/index/index -1 $path/to/reads/read_1.fq -2 $path/to/reads/read_2.fq -S eg2.sam`

<div style="text-align: justify"> You can use _SAMtools_ together with this step, which is a collection of tools for manipulating and analyzing SAM and BAM alignment files. You can simply pipe the following commands to the previus line of code `| samtools view -Sb | samtools sort -o /desired/path`. Now you have a sorted BAM file called eg2.sorted.bam for each pair of reads. You can create a new directory and move the BAM files in it, you can even include it in the SLURM script `mv *eg2.sorted.bam* ../alignment/sorted/`, for example. The last step before moving onto statistical analyisis will be counting how many reads mapped to genes. For this purpose you can use multiple softwares. In particular we used _featureCounts_, which is part of the _subreads_ software package. </div>

`$ featureCounts -p -O -t exon -g gene_id n -a /path/to/indes/DHL92_v4.gtf -o example_output_file.txt eg2.sorted.bam`

<div style="text-align: justify"> Running these commands will produce feature count tables. To make the count matrix of every bam file you can do the following:  </div>

`$ featureCounts -a /path/to/indes/DHL92_v4.gtf -o example_output_file.txt *bam \
| cut -f1,7- | sed 1d > counts.mx`

<div style="text-align: justify"> These combined feature count tables can be used for differential expression (DE) analysis. An example DE analysis script is included in this project, DESeq2. This script uses the R programming language with the DESeq2 stadistical package, the count matrix, the information about the samples (the columns of the count matrix), and the design formula. Then, GO and KEGG anotation can be performed using the example script. Finally, an example gene co-expression analysis script using R package igraph is included as well. </div>

## Credits

<div style="text-align: justify"> This analysis is a reevaluation of the original paper by Polonio et al.  [[1]](#1). It has been supervised by Manuel Franco Nicolas and Juana María Vivo Molina (Universidad de Murcia). </div>

## Bibliography

<a id="1">[1]</a> 
<div style="text-align: justify"> Alvaro Polonio et al. “RNA-seq analysis and fluorescence imaging of melon powdery mildew dis-
ease reveal an orchestrated reprogramming of host physiology”. In: Scientific Reports 9 (May 2019),
p. 7978. DOI: 10.1038/s41598-019-44443-5. </div>





