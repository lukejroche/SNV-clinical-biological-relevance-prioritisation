#!/usr/bin/env nextflow

params.samples = "test_samples.txt"
params.vcf     = "/home/lroche/clinical_genomics/variant_detect/1KG.chr22.anno.vcf.gz"
params.clinvar = "/home/lroche/clinical_genomics/variant_detect/clinvar.vcf.gz"
params.gnomad  = "/home/lroche/references/gnomAD/gnomad.exomes.v4.1.1.sites.chr22.vcf.bgz"

workflow {

    vcf_tuple = tuple(
        file(params.vcf),
        file("${params.vcf}.tbi")
    )

    sample_ch = Channel.fromPath(params.samples)
        .splitText()
        .map { it.trim() }
        .filter { it }

    split_input = sample_ch.map { sample ->
        tuple(sample, vcf_tuple[0], vcf_tuple[1])
    }

    split_out = extract_sample(split_input)
    filtered_out = filter_variants(split_out)
    snpeff_out = snpeff(filtered_out)
    clinvar_out = clinvar(snpeff_out)
    gnomad_out = gnomad(clinvar_out)
    tsv_out = make_tsv(gnomad_out)
    score_script=Channel.value(file("scripts/score_variants.py"))
    score(tsv_out, score_script)
}

process extract_sample {

    input:
    tuple val(sample), path(vcf), path(vcf_index)

    output:
    tuple val(sample),
          path("${sample}.vcf.gz"),
          path("${sample}.vcf.gz.tbi")

    script:
    """
    bcftools view \
        -s $sample \
        $vcf \
        -Oz \
        -o ${sample}.vcf.gz

    tabix -f -p vcf ${sample}.vcf.gz
    """
}

process filter_variants {

    input:
         tuple val(sample), path(vcf), path(tbi)

    output:
      tuple val(sample),
      path("${sample}.filtered.vcf.gz"),
      path("${sample}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools filter \
      -i 'INFO/AF<1 && FORMAT/GQ>20' \
      $vcf \
      -Oz \
      -o ${sample}.filtered.vcf.gz

    tabix -f -p vcf ${sample}.filtered.vcf.gz
    """
}

process snpeff {

    memory '3 GB'

    input:
      tuple val(sample), path(vcf), path(tbi)

    output:
      tuple val(sample),
      path("${sample}.snpeff.vcf.gz"),
      path("${sample}.snpeff.vcf.gz.tbi")

    script:
    """
    snpEff -Xmx3g ann GRCh37.75 $vcf > ${sample}.snpeff.vcf
    bgzip -f ${sample}.snpeff.vcf
    tabix -f -p vcf ${sample}.snpeff.vcf.gz
    """
}

process clinvar {

    input:
      tuple val(sample), path(vcf), path(tbi)

    output:
      tuple val(sample),
      path("${sample}.clinvar.vcf.gz"),
      path("${sample}.clinvar.vcf.gz.tbi")

    script:
    """
    bcftools annotate \
      -a ${params.clinvar} \
      -c INFO/CLNSIG,INFO/CLNDN \
      $vcf \
      -Oz -o ${sample}.clinvar.vcf.gz

    tabix -f -p vcf ${sample}.clinvar.vcf.gz
    """
}

process gnomad {

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
      tuple val(sample),
          path("${sample}.gnomad.filtered.vcf.gz"),
          path("${sample}.gnomad.filtered.vcf.gz.tbi")

    script:
    """
    bcftools annotate \
      -a ${params.gnomad} \
      -c INFO/AF \
      $vcf \
      -Oz -o ${sample}.gnomad.vcf.gz

    tabix -f -p vcf ${sample}.gnomad.vcf.gz

    bcftools filter \
      -i 'INFO/AF<0.001 || INFO/AF="."' \
      ${sample}.gnomad.vcf.gz \
      -Oz -o ${sample}.gnomad.filtered.vcf.gz

    tabix -f -p vcf ${sample}.gnomad.filtered.vcf.gz
    """
}

process make_tsv {

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.variants.tsv")

    script:
    """
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/AF\n' \
      $vcf > ${sample}.variants.tsv
    """
}


process score {

    publishDir "results/ranked", mode: "copy", overwrite:"true"

    input:
    tuple val(sample), path(tsv)
    path script

    output:
    path "${sample}_ranked.tsv"

    script:
    """
    python3 $script \
        $tsv \
        ${sample}_ranked.tsv
    """
}
