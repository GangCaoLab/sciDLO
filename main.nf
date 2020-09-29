#!/usr/bin/env nextflow

Channel.fromFilePairs(params.input_reads)
    .set{reads}

process unzip {
    input:
        tuple val(lib_id), file(reads) from reads
    output:
        tuple val(lib_id), file("*_{1,2}.fq") into unziped
    
    publishDir 'result/tmp/unzip'

    """
    if [[ ${reads[0]} == *.gz ]]; then
        unpigz ${reads[0]} -p $task.cpus -c > ${lib_id}_1.fq
        unpigz ${reads[1]} -p $task.cpus -c > ${lib_id}_2.fq
    else
        ln -s ${reads[0]} ${lib_id}_1.fq
        ln -s ${reads[1]} ${lib_id}_2.fq
    fi
    """
}

process extract_PETs {
    input:
        tuple val(lib_id), file(reads) from unziped

    output:
        tuple val(lib_id), file("*.pet{1,2}.fq") into pets

    publishDir 'result/tmp/extract_PETs'

    script:
    if (params.SE_mode)
        """
        expet --fq1 ${reads[0]} --linker $params.linker --enzyme $params.enzyme --output_prefix $lib_id -t $task.cpus -b --adapter $params.adapter
        """
    else
        """
        expet --fq1 ${reads[0]} --fq2 ${reads[1]} --linker $params.linker --enzyme $params.enzyme --output_prefix $lib_id -t $task.cpus -b
        """
}

process build_bedpe {
    input:
        tuple val(lib_id), file(pets) from pets

    output:
        tuple val(lib_id), file("*.uniq.bedpe") into bedpes

    publishDir 'result/tmp/build_bedpe'

    "dlohic build_bedpe $pets $lib_id --bwa-index $params.bwa_index_prefix -p $task.cpus"
}

rest_sites_file = Channel.value()
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

    "dlohic extract_fragments ${params.fasta_path} ${params.enzyme_name}.hdf5 -f hdf5 -p $task.cpus -r $params.enzyme"
}

process noise_reduce {
    input:
        tuple val(lib_id), file(bedpe) from bedpes
        val rest_file from rest_sites_file
    
    output:
        tuple val(lib_id), file("*.nr.bedpe") into nr_bedpes

    publishDir 'result/tmp/noise_reduce'

    "dlohic noise_reduce $bedpe ${lib_id}.nr.bedpe -r $rest_file -p $task.cpus"
}

process build_pairs {
    input:
        tuple val(lib_id), file(nr_bedpe) from nr_bedpes
    
    output:
        tuple val(lib_id), file("*.pairs") into pairs_files_1, pairs_files_2

    publishDir 'result/tmp/build_pairs'

    "dlohic bedpe2pairs $nr_bedpe ${lib_id}.pairs --keep --remove-redundancy --ncpu $task.cpus"
}

juicer_tools_jar = Channel.value()
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

process split_cells {
    input:
        tuple val(lib_id), file(pairs) from pairs_files_1
    output:
        tuple val(lib_id), file("${lib_id}.cell.*.pairs") into _pairs_per_cell 

    publishDir 'result/spcell'
    
    "spcell $pairs $params.barcodes_file -o ${lib_id}.cell -t $task.cpus"
}

_pairs_per_cell
    .flatMap { lib_id, cell_pairs ->
        def cells = []
        for (p in cell_pairs) {
            String fname = p.name
            _parts = fname.split("\\.")
            cell_id = _parts[0..(_parts.length-2)].join(".")
            cells << [cell_id, p]
        }
        return cells
    }
    .into{pairs_per_cell; _pairs_per_cell_2}

process sort_pairs_per_cell {
    input:
        tuple val(cell_id), file(pairs) from pairs_per_cell
    output:
        tuple val(cell_id), file("*.cell.*.sorted.pairs") into sorted_pairs_per_cell

    publishDir 'result/pairs_cell'
    
    """
    echo "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2" > ${cell_id}.sorted.pairs
    sort --parallel=${task.cpus} -k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7 $pairs >> ${cell_id}.sorted.pairs
    """
}

process build_dot_hic_per_cell {
    input:
        tuple val(cell_id), file(pairs) from sorted_pairs_per_cell
        val jar from juicer_tools_jar
    output:
        tuple val(cell_id), file("*.cell.*.hic") into dot_hics_per_cell

    publishDir 'result/dot_hic_cell'
    
    "java -jar $jar pre $pairs ${cell_id}.hic ${params.chrom_file} -r $params.resolutions"
}

_pairs_per_cell_2
    .map { cell_id, cell_pairs ->
        _parts = cell_id.split("\\.")
        String lib_id = _parts[0]
        return [lib_id, cell_pairs]
    }
    .groupTuple()
    .set{grouped_cells_one_lib}

process merge_pairs_per_library {
    input:
        tuple val(lib_id), file(pairs_files) from grouped_cells_one_lib
    output:
        tuple val(lib_id), file("*.merge.pairs") into merged_pairs, merged_pairs_2

    publishDir 'result/pairs_lib'
    """
    echo "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2" > ${lib_id}.merge.pairs
    cat $pairs_files | sort --parallel=${task.cpus} -k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7 >> ${lib_id}.merge.pairs
    """
}

process build_dot_hic_per_library {
    input:
        tuple val(lib_id), file(pairs) from merged_pairs
        val jar from juicer_tools_jar

    output:
        tuple val(lib_id), file("*.hic") into dot_hics

    publishDir 'result/dot_hic_lib'

    "java -jar $jar pre $pairs ${lib_id}.hic ${params.chrom_file} -r $params.resolutions"
}

merged_pairs_2
    .map{lib_id, f -> return f }
    .collect()
    .set{for_merge_all}

process merge_all_pairs {
    input:
        file(pairs_files) from for_merge_all
    output:
        file("merge_all.pairs") into merge_all_pairs
    
    publishDir 'result/pairs_all'
    """
    echo "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2" > merge_all.pairs
    cat $pairs_files | sort --parallel=${task.cpus} -k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7 >> merge_all.pairs
    """
}

process build_dot_hic_all {
    input:
        file(pairs) from merge_all_pairs
        val jar from juicer_tools_jar
    output:
        file("merge_all.hic") into merge_all_hic
    
    publishDir 'result/dot_hic_all'

    "java -jar $jar pre $pairs merge_all.hic ${params.chrom_file} -r $params.resolutions"
}
