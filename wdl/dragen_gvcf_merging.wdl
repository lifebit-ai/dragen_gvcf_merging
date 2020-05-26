# Dragen gVCF merging, with optional joint calling
# This workflow is for testing on Helix, dragen running in Singularity container
# Cloud workflow will have dragen running in docker container, will need minor tweaks

import "dragen_scatter_over_gvcfs.wdl" as gvcf_scatter

workflow dragen_aggregate_gvcfs {
    # Static inputs
    File dragen_container
    File reference
    File dragen_hash_table
    Int gvcf_lines
    Int max_alleles
    Int threads
    String singularity_module

    # Variable Inputs
    File gvcf_list
    File regions_file
    Array[String] regions = read_lines(regions_file)
    Boolean joint_call

    call split {
        input:  gvcf_list=gvcf_list,
                gvcf_lines=gvcf_lines,
    }

    scatter (region in regions) {
        call gvcf_scatter.scatter_over_gvcfs{
            input:  gvcf_chunk_list=split.gvcf_chunk_list,
                    reference=reference,
                    max_alleles=max_alleles,
                    region=region,
                    singularity_module=singularity_module,
                    dragen_container=dragen_container,
                    dragen_hash_table=dragen_hash_table,
                    threads=threads
        }
    }
    output {
        Array[File] merged_regions = scatter_over_gvcfs.merged_region
        Array[File] merged_regions_index = scatter_over_gvcfs.merged_region_index
    }
}

task split {
    File gvcf_list
    Int gvcf_lines
    command {
        split -d \
        -a 6 \
        --lines ${gvcf_lines} \
        ${gvcf_list} \
        gvcf_chunk_
    }
    runtime {
        cpus: 1
        memory: "1 GB"
        lsf_queue: "high"
        lsf_project: "bio"
    }
    output {
        Array[File] gvcf_chunk_list = glob("gvcf_chunk_*")
    }
}

