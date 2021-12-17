# ClonENet
ClonENet is a method that supports single-sample and multi-sample to infer subclonal population, using single nucleotide polymorphism and copy number variation information. <br>
The script uses human genes as an example to show how to obtain the input file of ClonENet, that is, from the gene sequencing file to the generation of the ClonENet input file.
## Preprocess steps
* Using BWA to align the sequencing data. and the SAM file will be generated. Converting SAM files to BAM files with samtools.
```Bash
    $ bwa mem -t 24 -M -Y -R "@RG\tID:$pfx\tPL:illumina\tLB:WGS\tSM:$pfx" $reference $fastq1 $fastq2 | samtools view -Sb > $pfx.bam

    $ samtools sort -@ 24 -o $pfx.sorted.bam $pfx.bam
```
* Using gatk to process BAM files.
```Bash
    $ gatk MarkDuplicates -I $pfx.sorted.bam -O $pfx.sorted.markdup.bam \
    -M $pfx.sorted.markdup.txt -REMOVE_DUPLICATES true

    $ gatk BuildBamIndex -I $pfx.sorted.markdup.bam -O $pfx.sorted.markdup.bai

    $ gatk BaseRecalibrator -R $reference -I $pfx.sorted.markdup.bam \
    --known-sites $indel1 --known-sites $dbsnp \
    --known-sites $indel2 -O $pfx.table

    $ gatk ApplyBQSR --bqsr-recal-file $pfx.table \
    -R $reference -I $pfx.sorted.markdup.bam \
    -O $pfx.sorted.markdup.bqsr.bam
```
* Using gatk mutect2 to detect SNP in filtered BAM files 
```Bash
    $ gatk Mutect2 -R ${ref} \
    -I $pfx.sorted.markdup.bqsr.bam -I $normal.sorted.markdup.bqsr.bam \
    -tumor ${pfx} -normal ${normal} -L ${interval_list} -O ${pfx}.mutect2.vcf

    $ gatk FilterMutectCalls -V ${pfx}.mutect2.vcf \
    -O ${pfx}.somatic.vcf -R ${ref}
```
* Using CNVkit to analyze the copy number of the BAM file
```Bash
    $ cnvkit.py batch ${pfxc}.sorted.markdup.bqsr.bam -n ${pfxn}.sorted.markdup.bqsr.bam \
    -m wgs -f ${reference} --annotate ${refFlat} \
    --output-dir ${work_dir} --diagram --scatter
```
* Convert the detection results of SNP and CN into the format of the input file
```Bash
    $ python toInput.py -S ${pfx}.somatic.vcf -C ${pfx}.sorted.call.cnr -O ${work_dir}
```
ClonENet is a method to infer subclonal populations using major and minor copy number (CN), variant allelic frequency (VAF) of somatic mutation. Therefore, the input file needs to provide relevant information. If there is no major and minor copy number, set minor with 0, and the value of major is the absolute CN.
We give some alternative software to detect CNV and SNP:
* Software for detecting SNP:
    * gatk HaplotypeCaller;
    * Varscan;
    * svmSomatic
* Software for detecting CNV:
    * CNVnator;
    * Battenberg;
    * Sclust CN;
    * IFTV-CNV;
## Run ClonENet
