workflow scatter_over_gvcfs {
    # Static inputs
    File reference
    File dragen_hash_table
    File dragen_container
    Int max_alleles
    Int threads
    String singularity_module
    
    # Variable inputs
    Array[File] gvcf_chunk_list
    String region
    String fixed_region = sub(region, "[:-]", "_")

    scatter (gvcf_subset in gvcf_chunk_list) {
        call first_joint_aggregation {
            input:  gvcf_subset=gvcf_subset,
                    reference=reference,
                    region=region,
                    fixed_region=fixed_region,
                    max_alleles=max_alleles,
                    singularity_module=singularity_module,
                    dragen_container=dragen_container,
                    threads=threads
        }
    }
    call merge_consensus_sites {
        input:  first_joint_aggregation_output=first_joint_aggregation.first_joint_aggregation_output,
                reference=reference,
                dragen_container=dragen_container,
                singularity_module=singularity_module,
                fixed_region=fixed_region
    }
    scatter (gvcf_subset in gvcf_chunk_list) {
        call second_joint_aggregation {
            input:  gvcf_subset=gvcf_subset,
                    fixed_region=fixed_region,
                    region=region,
                    max_alleles=max_alleles,
                    reference=reference,
                    singularity_module=singularity_module,
                    dragen_container=dragen_container,
                    merged_consensus_sites=merge_consensus_sites.merged_consensus_sites,
                    threads=threads
        }
    }
    call merge_chunks {
        input:  second_joint_aggregation_output_files=second_joint_aggregation.second_joint_aggregation_output,
                fixed_region=fixed_region,
                singularity_module=singularity_module,
                dragen_container=dragen_container
    }
    # call normalise {
    #     input:  merged_region=merge_chunks.merged_region,
    #             vt=vt,
    #             reference=reference,
    #             bcftools=bcftools,
    #             fixed_region=fixed_region
    # }
    output {
        File merged_region = merge_chunks.merged_region
        File merged_region_index = merge_chunks.merged_region_index
    }
}

# 1st joint-merging pass for each chunk, to build consensus sites
task first_joint_aggregation {
    File gvcf_subset
    File reference
    File dragen_container
    String gvcf_subset_name = basename(gvcf_subset)
    String region
    String fixed_region
    String singularity_module
    Int max_alleles
    Int threads
    runtime {
        cpus: "${threads}"
        memory: "22 GB"
        lsf_queue: "high"
        lsf_project: "bio"
    }
    command {
        set -eou pipefail;

        ${singularity_module}

        singularity run ${dragen_container} dragen --sw-mode \
         --enable-gvcf-genotyper=true \
         --enable-map-align=false \
         --gg-extra-params="--drop-genotypes" \
         --gg-enable-indexing=true \
         --output-directory . \
         --output-file-prefix ${gvcf_subset_name}_${fixed_region}_first_aggregation \
         --variant-list ${gvcf_subset} \
         --ht-reference ${reference} \
         --gg-max-alternate-alleles ${max_alleles} \
         --gg-regions ${region} \
         --num-threads ${threads}
    }
    output {
        File first_joint_aggregation_output = "${gvcf_subset_name}_${fixed_region}_first_aggregation.vcf.gz"
        File first_joint_aggregation_output_index = "${gvcf_subset_name}_${fixed_region}_first_aggregation.vcf.gz.tbi"
    }
}

# merge ~concatenate~ consensus sites from all chunks
task merge_consensus_sites {
    Array[File] first_joint_aggregation_output
    File create_consensus_file_list = write_lines(first_joint_aggregation_output)
    File reference
    File dragen_container
    String fixed_region
    String singularity_module
    runtime {
        cpus: 1
        memory: "1 GB"
        lsf_queue: "high"
        lsf_project: "bio"
    }
    command {
        ${singularity_module}

        singularity run ${dragen_container} bcftools concat -a \
         -f ${create_consensus_file_list} \
         -Ov \
        | bcftools norm -m+ any \
         -f ${reference} \
         -d none \
         -Ob \
         -o sites_with_consensus_samples_${fixed_region}.bcf

        singularity run ${dragen_container} bcftools index sites_with_consensus_samples_${fixed_region}.bcf
    }
    output {
        File merged_consensus_sites = "sites_with_consensus_samples_${fixed_region}.bcf"
    }
}

# Second gvcf-merging pass for each chunk, forcing genotypes for each site
# in consensus allele list from the previous step
task second_joint_aggregation {
    File gvcf_subset
    File reference
    File merged_consensus_sites
    File dragen_container
    String gvcf_subset_name = basename(gvcf_subset)
    String region
    String fixed_region
    String singularity_module
    Int max_alleles
    Int threads
    runtime {
        cpus: "${threads}"
        memory: "22 GB"
        lsf_queue: "high"
        lsf_project: "bio"
    }
    command {
        set -eou pipefail;
        ${singularity_module}

        singularity run ${dragen_container} dragen --sw-mode \
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
    }
    output {
        File second_joint_aggregation_output = "${gvcf_subset_name}_${fixed_region}_second_aggregation.vcf.gz"
        File second_joint_aggregation_output_index = "${gvcf_subset_name}_${fixed_region}_second_aggregation.vcf.gz.tbi"
    }
}

# Merge chunks per region
task merge_chunks {
    Array[File] second_joint_aggregation_output_files
    File second_joint_aggregation_file_list = write_lines(second_joint_aggregation_output_files)
    File dragen_container
    String fixed_region
    String singularity_module
    runtime {
        cpus: 1
        memory: "30 GB"
        lsf_queue: "high"
        lsf_project: "bio"
    }
    command {
        ${singularity_module}

        singularity run ${dragen_container} bcftools merge -m all \
         -l ${second_joint_aggregation_file_list}
         -Oz \
         -o all_chunks_merged_${fixed_region}.vcf.gz \
        
        singularity run ${dragen_container} bcftools index all_chunks_merged_${fixed_region}.vcf.gz
    }
    output {
        File merged_region = "all_chunks_merged_${fixed_region}.vcf.gz"
        File merged_region_index = "all_chunks_merged_${fixed_region}.vcf.gz.csi"
    }
}

# Joint calling - optional
# task joint_calling {
#     File dragen_container
#     File merged_chunks
#     File dragen_hash_table
#     String singularity_module
#     String fixed_region
#     runtime {
#         cpus: "${threads}"
#         memory: "22 GB"
#         lsf_queue: "short"
#         lsf_project: "bio"
#     }
#     command {
#         ${singularity_module}

#         singularity run ${dragen_container} dragen --sw-mode \
#          --enable-joint-genotyping=true \
#          --num-threads \
#          --variant ${merged_chunks} \
#          --ref-dir ${dragen_hash_table} \
#          --output-directory .
#          --output_file_prefix all_chunks_joint_called_${fixed_region}

#         singularity run ${dragen_container} bcftools index all_chunks_joint_called_${fixed_region}.vcf.gz
#     }
#     output {

#     }
# }





# Normalise regions
# task normalise {
#     File merged_region
#     File reference
#     String vt
#     String bcftools
#     String fixed_region
#     runtime {
#         cpus: 1
#         memory: "30 GB"
#         lsf_queue: "medium"
#         lsf_project: "bio"
#     }
#     command {
#         ${vt}
#         ${bcftools}

#         vt decompose -s ${merged_region} | \
#         vt normalize -n \
#         -w 10000 \
#         -r ${reference} - | \
#         vt decompose_blocksub - | \
#         bcftools view -Ob -o all_chunks_merged_norm_${fixed_region}.bcf ;
#         bcftools index all_chunks_merged_norm_${fixed_region}.bcf
#     }
#     output {
#         File norm_region = "all_chunks_merged_norm_${fixed_region}.bcf"
#         File norm_region_index = "all_chunks_merged_norm_${fixed_region}.bcf.csi"
#     }
# }