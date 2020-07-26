#!/usr/bin/env nextflow

reads = Channel.fromFilePairs(params.input_reads)

process extract_PETs {
    input:
        tuple val(lib_id), file(reads) from reads

    output:
        tuple val(lib_id), file("*.pet{1,2}.fq") into pets

    "expet $reads --linker $params.linker --enzyme $params.enzyme --output_prefix $lib_id -t $params.cpus -b"
}


process build_bedpe {
    input:
        tuple val(lib_id), file(pets) from pets

    output:
        tuple val(lib_id), file("*.uniq.bedpe") into bedpes

    "dlohic build_bedpe $pets $lib_id --bwa-index $params.bwa_index_prefix -p $params.cpus"
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
        tuple val(lib_id), file(bedpe) from bedpes
        file rest_file from rest_sites_file
    
    output:
        tuple val(lib_id), file("*.nr.bedpe") into nr_bedpes

    "dlohic noise_reduce $bedpe ${lib_id}.nr.bedpe -r $rest_file -p $params.cpus"
}

process build_pairs {
    input:
        tuple val(lib_id), file(nr_bedpe) from nr_bedpes
    
    output:
        tuple val(lib_id), file("*.pairs") into pairs_files_1, pairs_files_2

    "dlohic bedpe2pairs $nr_bedpe ${lib_id}.pairs --keep --remove-redundancy --ncpu $params.cpus"
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
        tuple val(lib_id), file("*.cell.*.pairs") into _pairs_per_cell 
    
    "spcell $pairs $params.barcodes_file -o ${lib_id}.cell"
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
    
    """
    echo "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2" > ${cell_id}.sorted.pairs
    sort --parallel=${params.cpus} -k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7 $pairs >> ${cell_id}.sorted.pairs
    """
}

process build_dot_hic_per_cell {
    input:
        tuple val(cell_id), file(pairs) from sorted_pairs_per_cell
        val jar from juicer_tools_jar
    output:
        tuple val(cell_id), file("*.cell.*.hic") into dot_hics_per_cell
    
    "java -jar $jar pre $pairs ${cell_id}.hic ${params.chrom_file} -r $params.resolutions"
}

dot_hics_per_cell.println()

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
        tuple val(lib_id), file("*.merge.pairs") into merged_pairs
    """
    echo "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2" > ${lib_id}.merge.pairs
    cat $pairs_files | sort --parallel=${params.cpus} -k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7 >> ${lib_id}.merge.pairs
    """
}

process build_dot_hic_per_library {
    input:
        tuple val(lib_id), file(pairs) from merged_pairs
        val jar from juicer_tools_jar

    output:
        tuple val(lib_id), file("*.hic") into dot_hics

    "java -jar $jar pre $pairs ${lib_id}.hic ${params.chrom_file} -r $params.resolutions"
}


dot_hics.println()
