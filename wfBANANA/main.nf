nextflow.enable.dsl = 2

/*
 * BaNaNA Nextflow/Singularity implementation.
 *
 * This workflow mirrors the original Snakemake pipeline while keeping the
 * existing Python helper scripts unchanged.
 */

params.samplesheet = params.samplesheet ?: null
params.barcode_dir = params.barcode_dir ?: null
params.outdir = params.outdir ?: 'results'
params.threads = params.threads ?: 8

params.min_len_filtering = params.min_len_filtering ?: 2000
params.max_len_filtering = params.max_len_filtering ?: 6000
params.min_mean_quality_filering = params.min_mean_quality_filering ?: 90
params.min_mean_quality_polishing = params.min_mean_quality_polishing ?: 20

params.rrnas = params.rrnas ?: '18S_rRNA:1000,28S_rRNA:400,5_8S_rRNA:90'
params.chosen_rrna = params.chosen_rrna ?: '18S_rRNA'

params.db_location = params.db_location ?: null
params.db_id = params.db_id ?: 0.7
params.db_query_cov = params.db_query_cov ?: 0.9

def truthy(value) {
    if (value == null) {
        return true
    }
    return value.toString().toLowerCase() in ['true', '1', 'yes', 'y']
}

def requireParam(name, value) {
    if (value == null || value.toString().trim() == '') {
        error "Missing required parameter: --${name}"
    }
}

def samplesheetToChannel(samplesheet) {
    Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample = row.sample?.toString()
            def fastq = row.fastq?.toString()

            if (!sample) {
                error "Samplesheet row is missing the 'sample' value: ${row}"
            }
            if (!fastq) {
                error "Samplesheet row is missing the 'fastq' value: ${row}"
            }
            if (sample.contains('_')) {
                error "Sample '${sample}' contains '_'. BaNaNA abundance parsing assumes sample names do not contain underscores."
            }

            tuple(sample, file(fastq, checkIfExists: true))
        }
}

def inputDirToChannel(inputDir) {
    Channel
        .fromPath("${inputDir}/*", checkIfExists: true)
        .filter { inputPath ->
            inputPath.isDirectory() || inputPath.getName().toString() ==~ /.*\.(fastq|fq)(\.gz)?$/
        }
        .map { inputPath ->
            def sample = inputPath.isDirectory()
                ? inputPath.getName().toString()
                : inputPath.getName().toString().replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')

            if (sample.contains('_')) {
                error "Sample '${sample}' contains '_'. BaNaNA abundance parsing assumes sample names do not contain underscores."
            }

            tuple(sample, inputPath)
        }
}

process FILTLONG {
    tag "$sample"
    label 'filtlong'
    publishDir "${params.outdir}/intermediate/01_filtlong", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("filtlong_${sample}.fastq"), emit: reads

    script:
    """
    filtlong \\
        --min_length ${params.min_len_filtering} \\
        --max_length ${params.max_len_filtering} \\
        --min_mean_q ${params.min_mean_quality_filering} \\
        ${reads} > filtlong_${sample}.fastq
    """
}

process PREPARE_READS {
    tag "$sample"
    label 'base'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fastq"), emit: reads

    script:
    """
    if [ -d "${reads}" ]; then
        find -L "${reads}" -maxdepth 1 -type f \\( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \\) | sort > ${sample}.fastq.list
        if [ ! -s ${sample}.fastq.list ]; then
            echo "No FASTQ files found for ${sample} in ${reads}" >&2
            exit 1
        fi
        while read fastq; do
            case "\$fastq" in
                *.gz) gzip -cd "\$fastq" ;;
                *) cat "\$fastq" ;;
            esac
        done < ${sample}.fastq.list > ${sample}.fastq
    else
        case "${reads}" in
            *.gz) gzip -cd "${reads}" > ${sample}.fastq ;;
            *) cat "${reads}" > ${sample}.fastq ;;
        esac
    fi
    """
}

process FASTQ_TO_FASTA {
    tag "$sample"
    label 'base'
    publishDir "${params.outdir}/intermediate/02_fastq_to_fasta", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("filtlong_${sample}.fasta"), emit: fasta

    script:
    """
    sed -n '1~4s/^@/>/p;2~4p' ${reads} > filtlong_${sample}.fasta
    """
}

process BARRNAP {
    tag "$sample"
    label 'barrnap'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/03_barrnap", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("barrnap_${sample}.fasta"), emit: rrna_calls

    script:
    """
    if grep -q '^>' ${fasta}; then
        barrnap \\
            --kingdom euk \\
            --reject 0.1 \\
            --outseq barrnap_${sample}.fasta \\
            ${fasta} \\
            --threads ${task.cpus}
    else
        touch barrnap_${sample}.fasta
    fi
    """
}

process EXTRACT_RRNA {
    tag "$sample"
    label 'python_tools'
    publishDir "${params.outdir}/intermediate/04_rrna_extracted", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(rrna_calls)
    path scripts_dir

    output:
    tuple val(sample), path("rrna_extracted_${sample}.fasta"), emit: fasta
    path "rrna_extraction_stats_${sample}.txt", emit: stats

    script:
    """
    python ${scripts_dir}/extracting_rrna.py \\
        -i ${rrna_calls} \\
        -r '${params.rrnas}' \\
        -cr '${params.chosen_rrna}' \\
        -o rrna_extracted_${sample}.fasta \\
        > rrna_extraction_stats_${sample}.txt
    """
}

process NANOPLOT {
    tag "$sample"
    label 'nanoplot'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/05_nanoplot", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("nanoplot_${sample}"), emit: stats_dir

    script:
    """
    mkdir -p nanoplot_${sample}
    if [ -s "${reads}" ]; then
        NanoPlot \\
            --fastq ${reads} \\
            --tsv_stats \\
            -t ${task.cpus} \\
            --info_in_report \\
            -o nanoplot_${sample}
    else
        cat > nanoplot_${sample}/NanoStats.txt <<'EOF'
Metrics	dataset
number_of_reads	0
number_of_bases	0
median_read_length	0
mean_read_length	0
read_length_stdev	0
n50	0
mean_qual	0
median_qual	0
EOF
    fi
    """
}

process CLUSTER_THRESHOLD {
    tag "$sample"
    label 'python_tools'
    publishDir "${params.outdir}/intermediate/06_cluster_threshold", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(nanoplot_dir)
    path scripts_dir
    path error_table

    output:
    tuple val(sample), path("clust_file_${sample}.txt"), emit: threshold

    script:
    """
    python ${scripts_dir}/calculate_clustering_threshold.py \\
        -s ${nanoplot_dir}/NanoStats.txt \\
        -e ${error_table} \\
        -o clust_file_${sample}.txt
    """
}

process ERROR_CLUSTERING {
    tag "$sample"
    label 'python_tools'
    cpus { params.threads as int }
    
    //publishDir "${params.outdir}/intermediate/07_error_clustering", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(rrna_fasta), path(threshold)

    output:
    tuple val(sample), path("samples/clusters_error_${sample}"), path("rrna_extracted_${sample}.fasta"), emit: clusters

    script:
    """
    mkdir -p samples/clusters_error_${sample}
    if grep -q '^>' ${rrna_fasta}; then
        vsearch \\
            --cluster_fast ${rrna_fasta} \\
            -id \$(cat ${threshold}) \\
            --threads ${task.cpus} \\
            --clusters samples/clusters_error_${sample}/cluster
    fi
    """
}

process CONSENSUS {
    tag "$sample"
    label 'python_tools'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/08_consensus", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(error_clusters), path(rrna_fasta)
    path scripts_dir

    output:
    tuple val(sample), path("consensus_${sample}.fasta"), path("samples/clusters_error_${sample}"), path("rrna_extracted_${sample}.fasta"), emit: consensus

    script:
    """
    mkdir -p samples
    if [ "\$(readlink -f ${error_clusters})" != "\$(readlink -f samples/clusters_error_${sample} 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_error_${sample}
        cp -r ${error_clusters} samples/clusters_error_${sample}
    fi
    mkdir -p samples/clusters_error_${sample}/alignments
    python ${scripts_dir}/mafft_consensus.py \\
        -i samples/clusters_error_${sample} \\
        -a samples/clusters_error_${sample}/alignments/ \\
        -t ${task.cpus} \\
        -o consensus_${sample}.fasta
    """
}

process MINIMAP {
    tag "$sample"
    label 'python_tools'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/09_minimap", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(consensus), path(error_clusters), path(rrna_fasta)
    path scripts_dir

    output:
    tuple val(sample), path("minimap_out_all_${sample}.paf"), path("samples/clusters_error_${sample}"), path("rrna_extracted_${sample}.fasta"), path("consensus_${sample}.fasta"), emit: paf

    script:
    """
    mkdir -p samples
    if [ "\$(readlink -f ${error_clusters})" != "\$(readlink -f samples/clusters_error_${sample} 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_error_${sample}
        cp -r ${error_clusters} samples/clusters_error_${sample}
    fi
    mkdir -p minimap_out_${sample}
    python ${scripts_dir}/minimap.py \\
        -c ${consensus} \\
        -cl samples/clusters_error_${sample} \\
        -t ${task.cpus} \\
        -o minimap_out_${sample}
    if find minimap_out_${sample} -type f -name '*.paf' | grep -q .; then
        cat minimap_out_${sample}/*.paf > minimap_out_all_${sample}.paf
    else
        touch minimap_out_all_${sample}.paf
    fi
    """
}

process RACON {
    tag "$sample"
    label 'racon'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/10_racon", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(paf), path(error_clusters), path(rrna_fasta), path(consensus)

    output:
    tuple val(sample), path("racon_${sample}.fasta"), path("samples/clusters_error_${sample}"), emit: polished

    script:
    """
    mkdir -p samples
    if [ "\$(readlink -f ${error_clusters})" != "\$(readlink -f samples/clusters_error_${sample} 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_error_${sample}
        cp -r ${error_clusters} samples/clusters_error_${sample}
    fi
    if [ ! -s ${consensus} ] || [ ! -s ${paf} ]; then
        touch racon_${sample}.fasta
    else
        racon \\
            ${rrna_fasta} \\
            -q ${params.min_mean_quality_polishing} \\
            -w 500 \\
            -t ${task.cpus} \\
            ${paf} \\
            ${consensus} \\
            > racon_${sample}.fasta
    fi
    """
}

process ADD_SAMPLE_NAMES {
    tag "$sample"
    label 'python_tools'
    publishDir "${params.outdir}/intermediate/11_named_polished", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(polished), path(error_clusters)
    path scripts_dir

    output:
    tuple val(sample), path("racon_name_${sample}.fasta"), path("samples/clusters_error_${sample}"), emit: named

    script:
    """
    mkdir -p samples
    if [ "\$(readlink -f ${error_clusters})" != "\$(readlink -f samples/clusters_error_${sample} 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_error_${sample}
        cp -r ${error_clusters} samples/clusters_error_${sample}
    fi
    python ${scripts_dir}/add_sample_id.py \\
        -i ${polished} \\
        -sn ${sample} \\
        -o racon_name_${sample}.fasta
    """
}

process MERGE_SAMPLES {
    label 'base'
    publishDir "${params.outdir}/intermediate/12_merged", mode: 'copy', overwrite: true

    input:
    path named_fastas

    output:
    path "merged.fasta", emit: merged

    script:
    """
    cat ${named_fastas} > merged.fasta
    """
}

process CHIMERA_REMOVAL {
    label 'python_tools'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/13_chimeras", mode: 'copy', overwrite: true

    input:
    path merged
    path database

    output:
    path "nonchim_db.fasta", emit: ref_filtered
    path "nonchim_db_dn.fasta", emit: nonchimeras

    script:
    """
    if grep -q '^>' ${merged}; then
        vsearch \\
            --uchime_ref ${merged} \\
            --db ${database} \\
            --threads ${task.cpus} \\
            --nonchimeras nonchim_db.fasta

        if grep -q '^>' nonchim_db.fasta; then
            vsearch \\
                --uchime2_denovo nonchim_db.fasta \\
                --threads ${task.cpus} \\
                --nonchimeras nonchim_db_dn.fasta
        else
            touch nonchim_db_dn.fasta
        fi
    else
        touch nonchim_db.fasta nonchim_db_dn.fasta
    fi
    """
}

process FINAL_CLUSTERING {
    label 'python_tools'
    cpus { params.threads as int }
    publishDir "${params.outdir}/intermediate/14_final_clustering", mode: 'copy', overwrite: true

    input:
    path nonchimeras

    output:
    path "pre_otus.fasta", emit: pre_otus
    path "samples/clusters_final", emit: clusters

    script:
    """
    mkdir -p samples/clusters_final
    if grep -q '^>' ${nonchimeras}; then
        vsearch \\
            --cluster_fast ${nonchimeras} \\
            -id 0.99 \\
            --threads ${task.cpus} \\
            --clusters samples/clusters_final/cluster \\
            --centroids pre_otus.fasta
    else
        touch pre_otus.fasta
    fi
    """
}

process REMOVE_N_SEQS {
    label 'python_tools'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path pre_otus
    path scripts_dir

    output:
    path "otus.fasta", emit: otus

    script:
    """
    python ${scripts_dir}/remove_Nseqs.py \\
        -i ${pre_otus} \\
        -o otus.fasta
    """
}

process TAXONOMY {
    label 'python_tools'
    cpus { params.threads as int }
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path otus
    path database

    output:
    path "taxonomy.tsv", emit: taxonomy

    script:
    """
    if grep -q '^>' ${otus}; then
        vsearch \\
            --usearch_global ${otus} \\
            --db ${database} \\
            --id ${params.db_id} \\
            --threads ${task.cpus} \\
            --blast6out taxonomy.tsv \\
            --query_cov ${params.db_query_cov}
    else
        touch taxonomy.tsv
    fi
    """
}

process TAXONOMY_TABLE {
    label 'python_tools'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path taxonomy
    path scripts_dir

    output:
    path "taxonomy_table.tsv", emit: table

    script:
    """
    python ${scripts_dir}/get_taxonomy_table.py \\
        -i ${taxonomy} \\
        -o taxonomy_table.tsv
    """
}

process ABUNDANCE {
    tag "$sample"
    label 'python_tools'
    publishDir "${params.outdir}/intermediate/15_abundance", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(error_clusters)
    path otus
    path final_clusters
    path scripts_dir

    output:
    path "abundance_${sample}.tsv", emit: abundance

    script:
    """
    mkdir -p samples
    if [ "\$(readlink -f ${error_clusters})" != "\$(readlink -f samples/clusters_error_${sample} 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_error_${sample}
        cp -r ${error_clusters} samples/clusters_error_${sample}
    fi
    if [ "\$(readlink -f ${final_clusters})" != "\$(readlink -f samples/clusters_final 2>/dev/null || true)" ]; then
        rm -rf samples/clusters_final
        cp -r ${final_clusters} samples/clusters_final
    fi
    python ${scripts_dir}/abundance.py \\
        -otu ${otus} \\
        -fclu samples/clusters_final \\
        -eclu samples/clusters_error_${sample} \\
        -sn ${sample} \\
        -o abundance_${sample}.tsv
    """
}

process OTU_TABLE {
    label 'python_tools'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path taxonomy
    path abundance_files
    path scripts_dir

    output:
    path "otu_table.tsv", emit: table

    script:
    """
    mkdir -p abundance_inputs
    cp abundance_*.tsv abundance_inputs/
    python ${scripts_dir}/get_otu_table.py \\
        -t ${taxonomy} \\
        -i abundance_inputs \\
        -o otu_table.tsv
    """
}

process BIOLOGICAL_REPORT {
    label 'python_tools'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path taxonomy
    path otu_table
    path nanoplot_dirs
    path rrna_stats
    path scripts_dir

    output:
    path "biological_report.html", emit: report
    path "banana_reports", emit: report_dir

    script:
    """
    python ${scripts_dir}/biological_report.py \\
        --taxonomy ${taxonomy} \\
        --otu-table ${otu_table} \\
        --output biological_report.html \\
        --output-dir banana_reports \\
        --nanoplot-dirs ${nanoplot_dirs} \\
        --rrna-stats ${rrna_stats}
    """
}

workflow {
    def hasSamplesheet = params.samplesheet != null && params.samplesheet.toString().trim() != ''
    def hasBarcodeDir = params.barcode_dir != null && params.barcode_dir.toString().trim() != ''

    if (hasSamplesheet && hasBarcodeDir) {
        error "Use either --samplesheet or --barcode_dir, not both"
    }
    if (!hasSamplesheet && !hasBarcodeDir) {
        error "Missing input: provide --barcode_dir or --samplesheet"
    }
    if (hasBarcodeDir && !file(params.barcode_dir).isDirectory()) {
        error "--barcode_dir must be a directory: ${params.barcode_dir}"
    }

    requireParam('db_location', params.db_location)
    def enable_taxonomy_table = truthy(params.enable_optional_taxonomy_format)

    samples_ch = hasBarcodeDir ? inputDirToChannel(params.barcode_dir) : (file(params.samplesheet).isDirectory() ? inputDirToChannel(params.samplesheet) : samplesheetToChannel(params.samplesheet))
    scripts_ch = Channel.value(file("${projectDir}/scripts", checkIfExists: true))
    error_table_ch = Channel.value(file("${projectDir}/files/P_error_table.tsv", checkIfExists: true))
    database_ch = Channel.value(file(params.db_location, checkIfExists: true))

    PREPARE_READS(samples_ch)
    FILTLONG(PREPARE_READS.out.reads)
    FASTQ_TO_FASTA(FILTLONG.out.reads)
    BARRNAP(FASTQ_TO_FASTA.out.fasta)
    EXTRACT_RRNA(BARRNAP.out.rrna_calls, scripts_ch)

    NANOPLOT(FILTLONG.out.reads)
    CLUSTER_THRESHOLD(NANOPLOT.out.stats_dir, scripts_ch, error_table_ch)

    rrna_with_threshold_ch = EXTRACT_RRNA.out.fasta.join(CLUSTER_THRESHOLD.out.threshold)
    ERROR_CLUSTERING(rrna_with_threshold_ch)
    CONSENSUS(ERROR_CLUSTERING.out.clusters, scripts_ch)
    MINIMAP(CONSENSUS.out.consensus, scripts_ch)
    RACON(MINIMAP.out.paf)
    ADD_SAMPLE_NAMES(RACON.out.polished, scripts_ch)

    named_fastas_ch = ADD_SAMPLE_NAMES.out.named.map { sample, fasta, clusters -> fasta }.collect()
    MERGE_SAMPLES(named_fastas_ch)

    CHIMERA_REMOVAL(MERGE_SAMPLES.out.merged, database_ch)
    FINAL_CLUSTERING(CHIMERA_REMOVAL.out.nonchimeras)
    REMOVE_N_SEQS(FINAL_CLUSTERING.out.pre_otus, scripts_ch)
    TAXONOMY(REMOVE_N_SEQS.out.otus, database_ch)
    if (enable_taxonomy_table) {
        TAXONOMY_TABLE(TAXONOMY.out.taxonomy, scripts_ch)
    }

    error_clusters_for_abundance_ch = ADD_SAMPLE_NAMES.out.named.map { sample, fasta, clusters -> tuple(sample, clusters) }
    ABUNDANCE(error_clusters_for_abundance_ch, REMOVE_N_SEQS.out.otus, FINAL_CLUSTERING.out.clusters, scripts_ch)
    abundance_files_ch = ABUNDANCE.out.abundance.collect()
    OTU_TABLE(TAXONOMY.out.taxonomy, abundance_files_ch, scripts_ch)
    nanoplot_dirs_ch = NANOPLOT.out.stats_dir.map { sample, stats_dir -> stats_dir }.collect()
    rrna_stats_ch = EXTRACT_RRNA.out.stats.collect()
    BIOLOGICAL_REPORT(TAXONOMY.out.taxonomy, OTU_TABLE.out.table, nanoplot_dirs_ch, rrna_stats_ch, scripts_ch)
}
