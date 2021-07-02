## Edit the values of the variables to fit your setup

# trimming scripts and primer files are from https://github.com/ItokawaK/Alt_nCov2019_primers
# just colone this repo to the path of this script
# or edit these paths to point to the appropriate files
trimming_script="Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py"
coverage_script="Alt_nCov2019_primers/tools/plot_depth.py"
primer_file="Alt_nCov2019_primers/Primers/ver_N1/primer.bed"

# For this we use the bed file from the Itokawa et al. 2020 paper and 
# converted it to a saf file so it cane be used with featureCounts
# see https://doi.org/10.1371/journal.pone.0239403.s003
fc_annotation="fc_annotation.saf"

# specify your output folder here
output_dir="results"

# bwa index needs to be build and specified
index_bwa="MN908947.3/MN908947.3.fa"
threads=8

# specify your R1 and R2 here.
# file names are expected to end with either _R1.fastq.gz or _R2.fastq.gz
# if the files end differently the script will fail. 
fq_read1="some_file_R1.fastq.gz"
fq_read2="some_file_R2.fastq.gz"

# analysis enspired by the repo from this paper https://github.com/ItokawaK/Alt_nCov2019_primers
example_wf (){
    fq1=$1
    fq2=$2
    out=$output_dir
    mkdir -p $output_dir
    # first pass aligment required for primer trimming as in the paper
    bwa mem -t $threads $index_bwa $fq1 $fq2 | tee $out/${fq1%_R1.fastq.gz}.sam | python $trimming_script --gzip $primer_file ${fq1%_R1.fastq.gz}_trimmed.fastq ${fq2%.fastq.gz}_trimmed.fastq
    # second pass alignment with the primer trimmed reads
    bwa mem -t $threads $index_bwa ${fq1%_R1.fastq.gz}_trimmed.fastq.gz ${fq2%.fastq.gz}_trimmed.fastq.gz > $out/${fq1%_R1.fastq.gz}_trimmed.sam
    # sort and convert to bam
    samtools sort $out/${fq1%_R1.fastq.gz}.sam  > $out/${fq1%_R1.fastq.gz}.bam
    samtools index $out/${fq1%_R1.fastq.gz}.bam   
    rm -v $out/${fq1%_R1.fastq.gz}.sam

    samtools sort $out/${fq1%_R1.fastq.gz}_trimmed.sam  > $out/${fq1%_R1.fastq.gz}_trimmed.bam
    samtools index $out/${fq1%_R1.fastq.gz}_trimmed.bam
    rm -v $out/${fq1%_R1.fastq.gz}.sam

    # only count the non primer trimmed filed as the annotation file used 
    # contains the primer binding sites
    # TODO: this is not yet working, needs to be fixed.
    featureCounts -a $fc_annotation -F SAF -o $out/${fq1%.fastq.gz}.counts -p -M --fraction -s 1 $out/${fq1%_R1.fastq.gz}.bam

    # here we plot the coverage. 
    # if you want the amplicon ranges and potential mutations visualized you need
    # to provide it with a bed and fasta file.
    # see https://github.com/ItokawaK/Alt_nCov2019_primers
    python $coverage_script -i $out/${fq1%_R1.fastq.gz}.bam -o $out/${fq1%_R1.fastq.gz}.pdf
}

example_wf $fq_read1 $fq_read2