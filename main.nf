#!/usr/bin/env nextflow

reads = Channel.fromFilePairs(params.input_reads)

process extract_PETs {
    input:
        tuple val(sample), file(reads) from reads

    output:
        tuple val(sample), file("*.pet{1,2}.fq") into pets

    "expet $reads --linker $params.linker --enzyme $params.enzyme --output_prefix $sample -t $params.cpus -b"
}


process build_bedpe {
    input:
        tuple val(sample), file(pets) from pets

    output:
        tuple val(sample), file("*.uniq.bedpe") into bedpes

    "dlohic build_bedpe $pets $sample --bwa-index $params.bwa_index_prefix -p $params.cpus"
}

rest_sites_file = Channel.create()
rest_file_exists = params.rest_file != "" && file(params.rest_file).exists()
if (rest_file_exists) {
    rest_sites_file << params.rest_file
}
if (params.enzyme_name == "") {
    params.enzyme_name = "UNKNOW"
}


process extract_fragments {
    input:
        val x from Channel.of(1)
    
    output:
        file "${params.enzyme_name}.hdf5" into rest_sites_file
    
    when:
        !rest_file_exists

    "dlohic extract_fragments ${params.fasta_path} ${params.enzyme_name}.hdf5 -f hdf5 -p $params.cpus -r $params.enzyme"
}

process noise_reduce {
    input:
        tuple val(sample), file(bedpe) from bedpes
        file rest_file from rest_sites_file
    
    output:
        tuple val(sample), file("*.nr.bedpe") into nr_bedpes

    "dlohic noise_reduce $bedpe ${sample}.nr.bedpe -r $rest_file -p $params.cpus"
}

process build_pairs {
    input:
        tuple val(sample), file(nr_bedpe) from nr_bedpes
    
    output:
        tuple val(sample), file("*.pairs") into pairs_files

    "dlohic bedpe2pairs $nr_bedpe ${sample}.pairs --keep --remove-redundancy --ncpu $params.cpus"
}

juicer_tools_jar = Channel.create()
jar_exists = params.juicer_tools_jar != "" && file(params.juicer_tools_jar).exists()
if (jar_exists) {
    juicer_tools_jar << params.juicer_tools_jar
}

process download_juicer_tools_jar {
    input:
        val x from Channel.of(1)

    output:
        file "juicer_tools.jar" into juicer_tools_jar

    when:
        !jar_exists

    "wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar -O juicer_tools.jar"
}

process build_dot_hic {
    input:
        tuple val(sample), file(pairs) from pairs_files
        val jar from juicer_tools_jar

    output:
        tuple val(sample), file("*.hic") into dot_hics

    "java -jar $jar pre $pairs ${sample}.hic ${params.chrom_file} -r $params.resolutions"
}

dot_hics.println()