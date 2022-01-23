# ClonENet
ClonENet is a method that supports single-sample to infer subclonal population, using single nucleotide polymorphism and copy number variation information. <br>
The script uses human genes as an example to show how to obtain the input file of ClonENet, that is, from the gene sequencing file to the generation of the ClonENet input file.
## Preprocess steps
* Using BWA to align the sequencing data. and the SAM file will be generated. Converting SAM files to BAM files with samtools.
```Bash
    $ bwa mem -t 24 -M -Y -R "@RG\tID:$pfx\tPL:illumina\tLB:WGS\tSM:$pfx" $reference $fastq1 $fastq2 | samtools view -Sb > $pfx.bam

    $ samtools sort -@ 24 -o ${pfx}.sorted.bam ${pfx}.bam
```
* Using gatk to process BAM files.
```Bash
    $ gatk MarkDuplicates -I ${pfx}.sorted.bam -O ${pfx}.sorted.markdup.bam \
    -M $pfx.sorted.markdup.txt -REMOVE_DUPLICATES true

    $ gatk BuildBamIndex -I ${pfx}.sorted.markdup.bam -O ${pfx}.sorted.markdup.bai

    $ gatk BaseRecalibrator -R ${reference} -I ${pfx}.sorted.markdup.bam \
    --known-sites $indel1 --known-sites $dbsnp \
    --known-sites $indel2 -O ${pfx}.table

    $ gatk ApplyBQSR --bqsr-recal-file ${pfx}.table \
    -R ${reference} -I ${pfx}.sorted.markdup.bam \
    -O ${pfx}.sorted.markdup.bqsr.bam
```
* Using gatk mutect2 to detect SNP in filtered BAM files 
```Bash
    $ gatk Mutect2 -R ${ref} \
    -I ${pfx}.sorted.markdup.bqsr.bam -I ${normal}.sorted.markdup.bqsr.bam \
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
    * [gatk HaplotypeCaller](https://github.com/broadinstitute/gatk);
    * [Varscan](https://github.com/dkoboldt/varscan);
    * [svmSomatic](https://github.com/BDanalysis/svmSomatic).
* Software for detecting CNV:
    * [CNVnator](https://github.com/abyzovlab/CNVnator);
    * [Battenberg](https://github.com/cancerit/cgpBattenberg);
    * [Sclust CN](https://www.clab-cosine.net/cun-web/);
    * IFTV-CNV.
* Tumor purity (Optional)
    * [ABSOLUTE](http://archive.broadinstitute.org/cancer/cga/absolute);
    * [Estimation](https://bioinformatics.mdanderson.org/estimate/index.html).
## Install ClonENet
   * Python >=3.7
   * numpy
   * sklearn 
   * matplotlib (Optional)
## Run ClonENet
* Subclonal population inference with single sample 
    * parameter '-c' represents CNV information file.
    * parameter '-s' represents SNV information file.
    * parameter '-p' represents tumor purity. If the purity of the sample is uncertain, it can be set to - 1.
    * parameter '-o' represents the directory of output
    * parameter '-x' represents the name or prefix of file
```Python
    python single_process.py \
        -c ./Sample/EGAF00000057355_cnv.txt \
        -s ./Sample/EGAF00000057355_snp.txt \
        -p -1 -o output_dir -x EGAF00000057355
```
* Estimation of the subclonal copy number
    * parameter '-i' represents the directory of input. Directory for saving the results of subclonal population inference
    * parameter '-c' represents CNV information file.
    * parameter '-p' represents tumor purity. If the purity of the sample is uncertain, it can be set to -1.
    * parameter '-o' represents the directory of output.
    * parameter '-x' represents the name or prefix of file.
```Python
    python subclone_copy_number.py \
        -i ./output_dir/ \
        -c ./Sample/EGAF00000057355_cnv.txt \
        -o ./output_dir/ -x EGAF00000057355 -p -1
```
* Cluster diagram (need the package named 'matplotlib')
    * parameter '-i' represents the directory of input. Directory for saving the results of subclonal population inference
    * parameter '-o' represents the directory of output.
    * parameter '-x' represents the name or prefix of file.
```Python
    python drow.py \
        -i ./output_dir/ \
        -o ./output_dir/ \
        -x EGAF00000057355
```
## File format
* [CNV file format](https://github.com/BDanalysis/ClonENet/blob/main/Sample/EGAF00000057355_cnv.txt)<br/>
chrom   start   end major_CN    minor_CN
* [SNV file format](https://github.com/BDanalysis/ClonENet/blob/main/Sample/EGAF00000057355_snp.txt)<br/>
chrom   position    allele_read_depth    total_read_depth
