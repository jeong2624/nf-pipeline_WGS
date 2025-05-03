nextflow.enable.dsl = 2

params.input = false
params.genome = "$baseDir/data/00-Ref/hg38_chr21.fa"
params.known_sites_vcf = "$baseDir/data/00-Ref/Known_sites/*.{vcf,vcf.gz}"

// Define the process for running Fastp
process Fastp {

    // Output directory
    publishDir "reports/Fastp", mode: "copy"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html

    script:
    def input_args = reads.size() == 1 ? "-i ${reads[0]}" : "-i ${reads[0]} -I ${reads[1]}"
    def output_args = reads.size() == 1 ? "-o ${meta.sample_id}.trimmed.fastq.gz" 
                                        : "-o ${meta.sample_id}.trimmed_R1.fastq.gz -O ${meta.sample_id}.trimmed_R2.fastq.gz"
    """
    fastp \\
    $input_args \\
    $output_args \\
    -j ${meta.sample_id}.json \\
    -h ${meta.sample_id}.html \\
    -R ${meta.sample_id}
    """
    }

// Define the process for running BWA (index)
process BWA_index {

    // Output directory
    publishDir "reports/BWA", mode: "copy"

    input:
    path fasta

    output:
    path("*.amb"), emit: amb
    path("*.ann"), emit: ann
    path("*.bwt"), emit: bwt
    path("*.pac"), emit: pac
    path("*.sa"), emit: sa

    script:
    """
    bwa index \\
    $fasta
    """
}

process BWA_mem {

    // Output directory
    publishDir "reports/BWA", mode: "symlink"

    input:
    tuple val(meta), path(reads)
    val index

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def input_read = reads.size() == 1 ? "${reads[0]}": "${reads[0]} ${reads[1]}"
    def output_bam = "${meta.sample_id}.bam"
    def read_group = "@RG\\tID:${meta.patient_id}\\tSM:${meta.sample_id}\\tPL:${meta.platform}\\tLB:${meta.lane}"

    """
    bwa mem \\
    -R '$read_group' \\
    $index \\
    $input_read | \\
    samtools view -bS - > $output_bam
    """
}

process SortSam {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam

    script:
    """
    picard \\
    SortSam \\
    -I $bam \\
    -O ${meta.sample_id}.sorted.bam \\
    --SORT_ORDER coordinate \\
    """
}

process MarkDuplicates {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("*.sorted.markdup.bam"), emit: bam
    tuple val(meta), path("*.markdup.metrics.txt"), emit: txt

    script:
    """
    picard MarkDuplicates \\
    I=$sorted_bam \\
    O=${meta.sample_id}.sorted.markdup.bam \\
    M=${meta.sample_id}.markdup.metrics.txt
    """
}

process BaseRecalibrator {
    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(markdup_bam)
    val fasta
    val known_sites_vcf

    output:
    path("recall_data.table"), emit: table

    script:
    def known_sites_args = known_sites_vcf.collect { "--known-sites ${it}" }.join(' ')
    """
    gatk \\
    BaseRecalibrator \\
    -I $markdup_bam \\
    -R $fasta \\
    $known_sites_args \\
    -O recall_data.table
    """
}

process ApplyBSQR {
    // Output directory
    publishDir "reports/GATK4", mode: "copy"

    input:
    tuple val(meta), path(markdup_bam)
    val fasta
    path table

    output:
    tuple val(meta), path("*.sorted.markdup.recal.bam"), emit: bam

    script:
    """
    gatk \\
    ApplyBQSR \\
    -R $fasta \\
    -I $markdup_bam \\
    --bqsr-recal-file $table \\
    -O ${meta.sample_id}.sorted.markdup.recal.bam
    """
}

process HaplotypeCaller {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(bam)
    val fasta

    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: vcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk \\
    HaplotypeCaller \\
    -R $fasta \\
    -I $bam \\
    -O ${meta.sample_id}.g.vcf.gz \\
    -ERC GVCF \\
    -OVI true
    """
}

process GenotypeGVCFs {

    // Output directory
    publishDir "reports/GATK4", mode: "copy"

    input:
    tuple val(meta), path(gvcf), path(tbi)
    val fasta

    output:
    tuple val(meta), path("*.rawVariants.vcf.gz"), emit: vcf
    tuple val(meta), path("*.rawVariants.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk \\
    GenotypeGVCFs \\
    -V $gvcf \\
    -R $fasta \\
    -OVI true \\
    -O ${meta.sample_id}.rawVariants.vcf.gz
    """
}

process Select_SNP {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val fasta

    output:
    tuple val(meta), path("*.rawSNPs.vcf.gz"), emit: vcf
    tuple val(meta), path("*.rawSNPs.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk \\
    SelectVariants \\
    -R $fasta \\
    -V $vcf \\
    --select-type-to-include SNP \\
    -O ${meta.sample_id}.rawSNPs.vcf.gz \\
    -OVI true
    """
}

process Select_INDEL {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val fasta

    output:
    tuple val(meta), path("*.rawINDELs.vcf.gz"), emit: vcf
    tuple val(meta), path("*.rawINDELs.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk \\
    SelectVariants \\
    -R $fasta \\
    -V $vcf \\
    --select-type-to-include INDEL \\
    -O ${meta.sample_id}.rawINDELs.vcf.gz \\
    -OVI true
    """
}

process Filter_SNP {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.snps_filtered.vcf.gz"), emit : vcf
    tuple val(meta), path("*.snps_filtered.vcf.gz.tbi"), emit : tbi

    script:
    """
    gatk \\
    VariantFiltration \\
    -V $vcf \\
    -filter "QD < 2.0" --filter-name "QD2" \\
    -filter "QUAL < 30.0" --filter-name "QUAL30" \\
    -filter "SOR > 3.0" --filter-name "SOR3" \\
    -filter "FS > 60.0" --filter-name "FS60" \\
    -filter "MQ < 40.0" --filter-name "MQ40" \\
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
    -filter "ExcessHet > 54.69" --filter-name ExcessHet \\
    -O ${meta.sample_id}.snps_filtered.vcf.gz \\
    -OVI true
    """
}

process Filter_INDEL {

    // Output directory
    publishDir "reports/GATK4", mode: "symlink"
    
    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.indels_filtered.vcf.gz"), emit : vcf
    tuple val(meta), path("*.indels_filtered.vcf.gz.tbi"), emit : tbi

    script:
    """
    gatk \\
    VariantFiltration \\
    -V $vcf \\
    -filter "QD < 2.0" --filter-name "QD2" \\
    -filter "QUAL < 30.0" --filter-name "QUAL30" \\
    -filter "FS > 200.0" --filter-name "FS200" \\
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
    -filter "ExcessHet > 54.69" --filter-name ExcessHet \\
    -O ${meta.sample_id}.indels_filtered.vcf.gz \\
    -OVI true
    """
}

process MergeVCF {

    // Output directory
    publishDir "reports/", mode: "copy"

    input:
    tuple val(meta), path(snp), path(snp_tbi)
    tuple val(meta2), path(indel), path(indel_tbi)

    output:
    tuple val(meta), path("*.output_variants.vcf.gz")
    tuple val(meta), path("*.output_variants.vcf.gz.tbi")

    script:
    """
    picard \\
    MergeVcfs \\
    -I $snp \\
    -I $indel \\
    -O ${meta.sample_id}.output_variants.vcf.gz \\
    --CREATE_INDEX true
    """
}

workflow {
    // Create a channel for the metadata and associated read files
    ch_sample_info = Channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map {row ->
                                  def meta = [patient_id: row.patient,
                                              sample_id: row.sample,
                                              lane: row.lane,
                                              platform: row.platform]
                                  def reads = [file(row.fastq_1), file(row.fastq_2)]
                                  def fasta = params.genome
                                  def known_site = file(params.known_sites_vcf)
                                  return [meta, reads, fasta, known_site]}

    // Running Fastp
    def ch_fastp_reads = ch_sample_info.map { $it -> [$it[0], $it[1]] }
    trimmed_reads = Fastp(ch_fastp_reads).reads
                                         .map{row ->
                                              def meta = row[0]
                                              def reads = row[1]
                                              return [meta, reads]}

    // Running Reference genome indexing using BWA
    def ref_fasta = ch_sample_info.map { it -> it[2] }
    ch_bwa_index = BWA_index(ref_fasta).amb.map{ it -> "${it.toString().replace(".amb", "")}" }

    // Running sequence alingment using BWA-mem
    ch_bam = BWA_mem(trimmed_reads, ch_bwa_index).bam
                                                 .map {row ->
                                                       def meta = row[0]
                                                       def bam = row[1].toString()
                                                       return [meta, bam]}
    // Sorting bam
    ch_sorted_bam = SortSam(ch_bam).bam
                                   .map {row ->
                                         def meta = row[0]
                                         def bam = row[1].toString()
                                         return [meta, bam]}
    // Mark-duplicate
    ch_markdup_bam = MarkDuplicates(ch_sorted_bam).bam
                                                  .map {row ->
                                                        def meta = row[0]
                                                        def bam = row[1].toString()
                                                        return [meta, bam]}

    // BaseRecalibrator
    def known_site_vcf = ch_sample_info.map { it -> it[3] }
    bsqr_table = BaseRecalibrator(ch_markdup_bam, ref_fasta, known_site_vcf).table
    
    // ApplyBQSR
    final_bam = ApplyBSQR(ch_markdup_bam, ref_fasta, bsqr_table).bam
                                                                .map {row ->
                                                                      def meta = row[0]
                                                                      def bam = row[1].toString()
                                                                      return [meta, bam]}

    // HaplotypeCaller
    ch_gVCF = HaplotypeCaller(final_bam, ref_fasta).vcf
                                                   .map {row ->
                                                         def meta = row[0]
                                                         def vcf = row[1].toString()
                                                         def tbi = row[1].toString() + ".tbi"
                                                         return [meta, vcf, tbi]}

    // GenotypeGVCFs
    chr_raw_variant = GenotypeGVCFs(ch_gVCF, ref_fasta).vcf
                                                   .map {row ->
                                                         def meta = row[0]
                                                         def vcf = row[1].toString()
                                                         def tbi = row[1].toString() + ".tbi"
                                                         return [meta, vcf, tbi]}

    // Hard filtering SNPs
    ch_snps = Select_SNP(chr_raw_variant, ref_fasta).vcf
                                                    .map {row ->
                                                          def meta = row[0]
                                                          def vcf = row[1].toString()
                                                          def tbi = row[1].toString() + ".tbi"
                                                          return [meta, vcf, tbi]}
    
    ch_filter_snps = Filter_SNP(ch_snps).vcf
                                        .map {row ->
                                              def meta = row[0]
                                              def vcf = row[1].toString()
                                              def tbi = row[1].toString() + ".tbi"
                                              return [meta, vcf, tbi]}

    // Hard filtering INDELs
    ch_indels = Select_INDEL(chr_raw_variant, ref_fasta).vcf
                                                        .map {row ->
                                                              def meta = row[0]
                                                              def vcf = row[1].toString()
                                                              def tbi = row[1].toString() + ".tbi"
                                                              return [meta, vcf, tbi]}
    
    ch_filter_indels = Filter_INDEL(ch_indels).vcf
                                              .map {row ->
                                                    def meta = row[0]
                                                    def vcf = row[1].toString()
                                                    def tbi = row[1].toString() + ".tbi"
                                                    return [meta, vcf, tbi]}

    // Merge filtered variants
    ch_final_vcf = MergeVCF(ch_filter_snps, ch_filter_indels)

}