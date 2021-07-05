## Edit the values of the variables below to fit your setup

# trimming scripts and primer files are from
# https://github.com/ItokawaK/Alt_nCov2019_primers
# just clone this repo to the path of this script or edit these paths to point
# to the appropriate files
trimming_script="Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py"
coverage_script="Alt_nCov2019_primers/tools/plot_depth.py"
primer_file="Alt_nCov2019_primers/Primers/ver_N1/primer.bed"
# For this we use the bed file from the Itokawa et al. 2020 paper and converted
# it to a simple annotation format file so it can be used with featureCounts
# see https://doi.org/10.1371/journal.pone.0239403.s003
fc_annotation="fc_annotation.saf"
# specify your output folder here
output_dir="results"
# bwa index needs to be build and specified
index_bwa="MN908947.3/MN908947.3.fa"
# the number of threads to be used by bwa
threads=8
# specify your R1 and R2 here.
# file names are expected to end with either _R1.fastq.gz or _R2.fastq.gz
# if the files end differently the script will fail. 
fq_read1="some_file_R1.fastq.gz"
fq_read2="some_file_R2.fastq.gz"
# this analysis workflow is inspired by the repo from this paper
# repo: https://github.com/ItokawaK/Alt_nCov2019_primers
# paper: https://doi.org/10.1371/journal.pone.0239403.s003
example_wf (){
    fq1=$1
    fq2=$2
    mkdir -p $output_dir
    # first pass alignment required for primer trimming as in the paper
    bwa mem -t $threads $index_bwa $fq1 $fq2 | \
      tee $output_dir/${fq1%_R1.fastq.gz}.sam | \
      python $trimming_script --gzip $primer_file \
        $output_dir/${fq1%.fastq.gz}_trimmed.fastq \
        $output_dir/${fq2%.fastq.gz}_trimmed.fastq
    # second pass alignment with the primer trimmed reads
    bwa mem -t $threads $index_bwa \
      $output_dir/${fq1%.fastq.gz}_trimmed.fastq.gz \
      $output_dir/${fq2%.fastq.gz}_trimmed.fastq.gz \
      > $output_dir/${fq1%_R1.fastq.gz}_trimmed.sam
    # sort and convert to bam
    samtools sort $output_dir/${fq1%_R1.fastq.gz}.sam \
      > $output_dir/${fq1%_R1.fastq.gz}.bam
    samtools index $output_dir/${fq1%_R1.fastq.gz}.bam
    rm -v $output_dir/${fq1%_R1.fastq.gz}.sam

    samtools sort $output_dir/${fq1%_R1.fastq.gz}_trimmed.sam  \
      > $output_dir/${fq1%_R1.fastq.gz}_trimmed.bam
    samtools index $output_dir/${fq1%_R1.fastq.gz}_trimmed.bam
    rm -v $output_dir/${fq1%_R1.fastq.gz}_trimmed.sam

    # only count the non primer trimmed filed as the annotation file used 
    # contains the primer binding sites
    featureCounts -a $fc_annotation -F SAF \
      -o $output_dir/${fq1%_R1.fastq.gz}.counts -p -M --fraction -s 1 \
      $output_dir/${fq1%_R1.fastq.gz}_trimmed.bam
    # here we plot the coverage. 
    # if you want the amplicon ranges and potential mutations visualized you
    # need to provide the script below with a bed and fasta file.
    # see https://github.com/ItokawaK/Alt_nCov2019_primers
    python $coverage_script -i $output_dir/${fq1%_R1.fastq.gz}_trimmed.bam \
      -o $output_dir/${fq1%_R1.fastq.gz}.pdf
}

example_wf $fq_read1 $fq_read2