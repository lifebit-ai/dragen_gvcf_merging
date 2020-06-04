#!/usr/bin/env nextflow
// Nextflow version of the dragen gvcfgenotyper workflow

Channel
    .fromPath(params.gvcf_list)
    .splitText(by: params.gvcf_lines, file: true)
    .into{ first_joint_aggregation_gvcfs; second_joint_aggregation_gvcfs }
    //.map{ vcf -> [vcf, file("${vcf}.csi")] }
Channel
    .fromPath(params.regions_file)
    .splitText()
    .set{ region }
Channel
    .fromPath(params.reference)
    .into{ first_joint_aggregation_reference; merge_consensus_sites_reference; second_joint_aggregation_reference }
Channel
    .value(params.max_alleles)
    .set{ max_alleles }
Channel
    .value(params.threads)
    .set{ threads }

first_joint_aggregation_channel = first_joint_aggregation_gvcfs
    .combine(region)
    .combine(first_joint_aggregation_reference)
    .combine(threads)
    .combine(max_alleles)

process first_joint_aggregation {
    tag "$gvcf_subset $fixed_region"
    input:
        tuple file(gvcf_subset), val(region), file(reference), val(num_threads), val(max_alleles) from first_joint_aggregation_channel

    output:
        tuple val(region), file("${gvcf_subset_name}_${fixed_region}_first_aggregation") into sample_consensus_sites

    // module "singularity/3.2.1"

    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    gvcf_subset_name = gvcf_subset.baseName
    // touch ${gvcf_subset_name}_${fixed_region}_first_aggregation
    // ulimit -u 65535
    // ulimit -n 65535
    """
    ulimit -a
    cat ${gvcf_subset} | awk -F'/' '{print \$NF}' > local_vcfs.txt

    dragen --sw-mode \
     --enable-gvcf-genotyper=true \
     --enable-map-align=false \
     --gg-extra-params="--drop-genotypes" \
     --gg-enable-indexing=true \
     --output-directory . \
     --output-file-prefix ${gvcf_subset_name}_${fixed_region}_first_aggregation \
     --variant-list local_vcfs.txt \
     --ht-reference ${reference} \
     --gg-max-alternate-alleles ${max_alleles} \
     --gg-regions ${region} \
     --num-threads ${num_threads}
    """
}

sample_consensus_sites
    .groupTuple()
    .set{ grouped_sample_consensus_sites }

process merge_consensus_sites {

    input:
        tuple val(region), file(gvcf_subset) from grouped_sample_consensus_sites
        each file(reference) from merge_consensus_sites_reference

    output:
        tuple val(region), file("sites_with_consensus_samples_${fixed_region}.bcf") into consensus_sites_for_second_aggregation
    
    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")

    // touch ${fixed_region}_consensus_sites
    """
    echo "${gvcf_subset.join("\n")}" >> ${fixed_region}_consensus_site_file_list

    bcftools concat -a \
     -f ${fixed_region}_consensus_site_file_list \
     -Ov \
    | bcftools norm -m+any \
     -f ${reference} \
     -d none \
     -Ob \
     -o sites_with_consensus_samples_${fixed_region}.bcf
    """   
}

second_joint_aggregation_channel = consensus_sites_for_second_aggregation
    .combine(second_joint_aggregation_gvcfs)
    .combine(second_joint_aggregation_reference)
    .combine(threads)
    .combine(max_alleles)


process second_joint_aggregation {

    input:
        tuple val(region), file(merged_consensus_sites), file(gvcf_subset), file(gvcf_subset_index), file(reference), val(num_threads), val(max_alleles) from second_joint_aggregation_channel

    output:
        tuple val(region), file("${gvcf_subset_name}_second_aggregation_output") into merge_samples

    // module "singularity/3.2.1"

    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    gvcf_subset_name = "gvcf_subset_${fixed_region}"

    // echo ${consensus_sites.baseName} >> ${gvcf_subset_name}__second_aggregation_input
    // cat "${gvcf_subset}" >> ${gvcf_subset_name}__second_aggregation_input
    // ulimit -u 65535
    // ulimit -n 65535
    """
    ulimit -a

    dragen --sw-mode \
     --enable-gvcf-genotyper=true \
     --enable-map-align=false \
     --gg-enable-indexing=true \
     --variant-list ${gvcf_subset} \
     --gg-extra-params="--allele-list ${merged_consensus_sites}" \
     --gg-regions ${region} \
     --ht-reference ${reference} \
     --gg-max-alternate-alleles ${max_alleles} \
     --num-threads ${threads} \
     --output-directory . \
     --output-file-prefix ${gvcf_subset_name}_${fixed_region}_second_aggregation \
     --gg-enable-concat=false
    """
}

process merge_chunks {
    
    publishDir "./results", mode: "copy"

    input:
        tuple val(region), file("samples") from merge_samples.groupTuple()
    output:
        tuple file("all_chunks_merged_${fixed_region}.vcf.gz"), file("all_chunks_merged_${fixed_region}.vcf.gz.")
    
    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    
    """
    ls ${samples} >> "${fixed_region}_second_joint_aggregation_file_list"

    bcftools merge -m+all \
     -l ${fixed_region}_second_joint_aggregation_file_list
     -Oz \
     -o all_chunks_merged_${fixed_region}.vcf.gz

    bcftools index -t all_chunks_merged_${fixed_region}.vcf.gz
    """
}