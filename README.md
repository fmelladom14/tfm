# Transcriptomic analysis of disease-resistance response in Cucumis melo L.

In this repository, code generated to perform the analysis is included. It allows you to perform RNA-seq analysis, GO and KEGG enrichment,  and create a co-expression network using the R package igraph. 

## Example of analysis

This part of the analysis can be integrated into job submission scripts, where you can specify more options. Then, using the `sbatch` command, SLURM runs your job when its your turn. Yet most commands are simple enough to use the command line. First, you need to download data, using fastq-dump. Yoo could do a simple loop to download all of the samples at once.

`$ fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir path/reads/ i`

The next step is quality control. FastQC will accept multiple file names as input, so we can use the *.fq or *fastq.gz  wildcard in the same directory where your samples are located.

`$ fastqc *.fastq.gz`

If your reads have good enough quality you can skip the next part of the process where you make sure that all the reads in the dataset have a chance to map/align to the genome, trimming unwanted information off from every read. We will not be performing this step becausee data did not have an appreciable amount of leftover adapter sequences or other contaminating sequences based on FastQC. Then, you have to align your reads to your reference genome. For this task we used HISAT2. hisat2-build builds a HISAT2 index from the reference genome in fasta format.

`$ hisat2-build -f DHL92_genome_v4.fa index`

Staying in the same directory, run the following command:

`$ hisat2 -f -x $path/to/index/index -1 $path/to/reads/read_1.fq -2 $path/to/reads/read_2.fq -S eg2.sam`

You can use SAMtools together with this step, which is a collection of tools for manipulating and analyzing SAM and BAM alignment files. You can simply pipe the following commands to the previus line of code `| samtools view -Sb | samtools sort -m 2G -o /desired/path`. Now you have a sorted BAM file called eg2.sorted.bam for each pair of reads, forward and reverse. You can create a new directory and move the BAM files in it, you can even include it in the job script `mv *eg2.sorted.bam* ../alignment/sorted/` for example.

The last step before moving onto the R notebooks contained in the repository will be counts how many reads mapped to genes. For this purpose you can use multiple softwares. In particular we used _featureCounts_, which is part of the subreads package.

`featureCounts -p -O -t exon -g gene_id n -a /path/to/indes/DHL92_v4.gtf -o example_output_file.txt eg2.sorted.bam`

To make the count matrix of every bam file you can do the following: 

`featureCounts -a /path/to/indes/DHL92_v4.gtf -o example_output_file.txt *bam \
| cut -f1,7- | sed 1d > counts.mx`

Providing the count matrix, the information about the samples (the columns of the count matrix), and the design formula, the data is now ready to analyse with the stadistical packages in R _(DESeq2, igraph)_.


## Credits

This analysis is a reevaluation of the original paper by Polonio et al.  [[1]](#1).

## Bibliography

<a id="1">[1]</a> 
Alvaro Polonio et al. “RNA-seq analysis and fluorescence imaging of melon powdery mildew dis-
ease reveal an orchestrated reprogramming of host physiology”. In: Scientific Reports 9 (May 2019),
p. 7978. DOI: 10.1038/s41598-019-44443-5.





