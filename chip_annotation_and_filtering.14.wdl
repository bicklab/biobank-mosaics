version 1.0

## This WDL filters mutect2 wgs vcfs to exons and annotates them with Annovar.
## ** ANNOVAR **
## Annovar functionally annotates genetic variants detected from diverse genomes.
## Given a list of variants with chromosome, start position, end position, reference nucleotide
## and observed nucleotides, Annovar can perform gene-based annotation, region-based annotation,
## filter-based annotation, and more.
##
## See ANNOVAR documentation to fully understand functionality:
## https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
##
## annovar_zip: the zipped folder with all of the needed files to run Annovar
##              NOTE: This file path is set on Terra - The file must be in the Workspace's bucket
## annovar_vcf_input: the Tables/sample column containing the vcf output files from a run of Mutect2
##                    NOTE: This is set on Terra (ex. this.filtered_vcf)
## annovar_protocols: the specificed protocols needed to run annovar (default = refGene,cosmic70)
##                    NOTE: You must add the needed file paths to annovar_data_sources
## annovar_operation: the specified operations needed to run annovar (default = g,f)
##                    NOTE: They must match up with annovar_protocols
## ref_name: the reference name needed for annovar to run (default = hg38)
## annovar_docker: the docker image to be used in the Annovar task
##
##
##
## Distributed under terms of the MIT License
## Copyright (c) 2025 Yash Pershad 

workflow passenger_mutations_exon_filtering_and_annotation {
    input {
      String vcf_directory
      File exon_bed
      File annovar_zip
      String docker = "gcr.io/bick-aps2/chip_annotation_and_filtering:1.0"
    }

    call find_files {
        input:
            directory = vcf_directory
    }

    scatter(idx in range(length(find_files.vcf_files))) {
      call filter_exons {
        input:
          vcf = find_files.vcf_files[idx],
          exon_bed = exon_bed
      }
    }

    Array[File] singleton_exons = select_all(filter_exons.singleton_exon)

    call Annovar {
        input:
            docker = docker,
            vcf_inputs = singleton_exons,
            annovar_zip = annovar_zip
    }
    
    output {
      Array[File] annovar_annotated_file_vcf = Annovar.annovar_output_files_vcf
      Array[File] annovar_annotated_file_table = Annovar.annovar_output_files_table
    }
}

# Find all VCF files in a directory
task find_files {
    input {
        String directory
        Int cpu = 1
        Float memory = 2
    }

    command {
        gsutil ls ${directory}*.vcf.gz > vcf_files.txt
    }
    
    output {
        Array[File] vcf_files = read_lines("vcf_files.txt")
    }

    runtime {
        docker: "google/cloud-sdk:latest"
        memory: memory + " GB"
        cpu: cpu
    }
}

task filter_exons {
    input {
        File vcf 
        File exon_bed
        Float memory = 8
        Int cpu = 4
        Int preemptible = 2
        Int maxRetries = 3
    }   
    
    Int disk_size = ceil(20 + 3 * size(vcf, "GiB"))

    String base = basename(vcf, ".vcf.gz")
    
    command <<<
        set -eu -o pipefail

        # Step 0: Index
        tabix -p vcf ~{vcf}

        #Step 1: Filter exons
        bcftools view -R ~{exon_bed} ~{vcf} -o ~{base}_exons.vcf
    >>>

    output {
        File singleton_exon = "~{base}_exons.vcf"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}


task Annovar {
    input {
      # Annovar inputs
      Array[File] vcf_inputs
      File annovar_zip
      String ref_name = "hg38"
      String annovar_protocols = "refGene,cosmic70"
      String annovar_operation = "g,f"

      # Runtime
      String docker
      Float memory = 4.0
      Int annovar_disk_space = 300
      Int cpu = 1
      Int preemptible = 1
      Int maxRetries = 0
    }

    command <<<
        set -euxo pipefail

        # Setup
        cp ~{annovar_zip} .
        unzip annovar_files.zip

        chmod +x annovar_files/convert2annovar.pl
        chmod +x annovar_files/table_annovar.pl
        chmod +x annovar_files/annotate_variation.pl
        chmod +x annovar_files/coding_change.pl
        chmod +x annovar_files/retrieve_seq_from_fasta.pl
        chmod +x annovar_files/variants_reduction.pl

        # Process each VCF file
        for vcf_input in ~{sep=' ' vcf_inputs}; do
            # Extract sample_id without extension
            sample_id=$(basename "${vcf_input}" .vcf | cut -d'.' -f1)
            file_prefix=${sample_id}.annovar_out

            # Annovar
            if ! perl annovar_files/table_annovar.pl ${vcf_input} annovar_files \
              -buildver ~{ref_name} \
              -out ${file_prefix} \
              -remove \
              -protocol ~{annovar_protocols} \
              -operation ~{annovar_operation} \
              -nastring . -vcfinput 2>/dev/null; then
                echo "ANNOVAR failed for ${vcf_input}"
                continue
            fi
        done
    >>>

    runtime {
      docker: docker
      memory: memory + " GiB"
      disk: "local-disk " + annovar_disk_space + " HDD"
      cpu: cpu
      preemptible: preemptible
      maxRetries: maxRetries
    }

    output {
      Array[File] annovar_output_files_vcf = glob("*.hg38_multianno.vcf")
      Array[File] annovar_output_files_table = glob("*.hg38_multianno.txt")
    }
}