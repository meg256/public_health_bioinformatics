version 1.0

task primer_clip {
  input {
    File read1
    File read2
    String samplename
    File reference_genome    # single-sequence reference, e.g. NC_000962.3
    File primer_bed          # primalscheme primer.bed
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0"
  }
  command <<<
    set -euo pipefail
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version)" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # fail early if samtools is too old (written to work under pipefail)
    ampliconclip_help=$(samtools ampliconclip 2>&1 || true)
    if ! grep -q -- '--both-ends' <<< "$ampliconclip_help"; then
      echo "ERROR: samtools in this image lacks 'ampliconclip --both-ends'; use samtools >= 1.12" >&2
      exit 1
    fi

    # primalscheme BEDs use a placeholder chrom name ("reference"); rewrite it to match the reference header.
    # Only valid for a single-sequence reference (true for H37Rv). Skips header/comment lines.
    ref_chrom=$(head -n1 ~{reference_genome} | sed 's/^>//; s/[[:space:]].*//')
    awk -v c="$ref_chrom" 'BEGIN{OFS="\t"} /^#/ || NF<3 {next} {$1=c; print}' ~{primer_bed} > primers.bed
    echo "primers in BED: $(wc -l < primers.bed)"

    cp ~{reference_genome} ref.fasta
    bwa index ref.fasta

    echo "input R1 reads: $(zcat -f ~{read1} | awk 'END{print NR/4}')"

    bwa mem -t ~{cpu} ref.fasta ~{read1} ~{read2} \
      | samtools sort -@ ~{cpu} -o aligned.bam -
    samtools index aligned.bam

    # Hard-clip so primer bases are removed from SEQ (soft clips survive samtools fastq).
    # Unmapped reads pass through unchanged.
    samtools ampliconclip \
      --hard-clip \
      --both-ends \
      -b primers.bed \
      -f ~{samplename}.ampliconclip.stats.txt \
      -o clipped.bam \
      aligned.bam

    samtools sort -n -@ ~{cpu} -o clipped.namesorted.bam clipped.bam

    # samtools fastq drops secondary/supplementary by default and restores original read orientation
    samtools fastq -@ ~{cpu} -n \
      -1 ~{samplename}_1.primerclip.fastq.gz \
      -2 ~{samplename}_2.primerclip.fastq.gz \
      -s ~{samplename}.primerclip.singletons.fastq.gz \
      -0 /dev/null \
      clipped.namesorted.bam

    echo "output R1 reads: $(zcat ~{samplename}_1.primerclip.fastq.gz | awk 'END{print NR/4}')"
  >>>
  output {
    File read1_clean = "~{samplename}_1.primerclip.fastq.gz"
    File read2_clean = "~{samplename}_2.primerclip.fastq.gz"
    File ampliconclip_stats = "~{samplename}.ampliconclip.stats.txt"
    String bwa_version = read_string("BWA_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String primer_clip_docker = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}