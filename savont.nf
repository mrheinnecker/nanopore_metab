nextflow.enable.dsl = 2

/*
 * Standalone Savont workflow.
 *
 * Input FASTQ files are normalized per sample, trimmed to 18S while preserving
 * FASTQ qualities, and then processed by savont asv. Optional classification
 * can be run against a pre-downloaded Savont database directory.
 */

params.samplesheet = params.samplesheet ?: null
params.outdir = params.outdir ?: 'results/savont_run'
params.threads = params.threads ?: 8
params.trim_threads = params.trim_threads ?: params.threads
params.savont_threads = params.savont_threads ?: params.threads
params.taxonomy_threads = params.taxonomy_threads ?: params.threads

params.kingdom = params.kingdom ?: 'euk'
params.pooled_samples = params.pooled_samples ?: false
params.savont_db = params.savont_db ?: null
params.classifier = params.classifier ?: 'none'
params.export_qiime = params.export_qiime ?: false
params.pr2_db = params.pr2_db ?: null
params.db_id = params.db_id ?: 0.7
params.db_query_cov = params.db_query_cov ?: 0.9
params.enable_taxonomy_table = params.enable_taxonomy_table != null ? params.enable_taxonomy_table : true
params.enable_biological_report = params.enable_biological_report != null ? params.enable_biological_report : true

params.min_read_length = params.min_read_length ?: 1100
params.max_read_length = params.max_read_length ?: 2000
params.quality_value_cutoff = params.quality_value_cutoff ?: 98
params.minimum_base_quality = params.minimum_base_quality ?: 25
params.min_cluster_size = params.min_cluster_size ?: 12

params.rrna_operon = params.rrna_operon ?: false
params.hifi = params.hifi ?: false
params.single_strand = params.single_strand != null ? params.single_strand : true
params.use_hpc = params.use_hpc ?: false
params.low_polymorphism = params.low_polymorphism ?: false
params.mask_low_quality = params.mask_low_quality ?: false

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

            tuple(sample, file(fastq, checkIfExists: true))
        }
}

def savontPresetFlags() {
    def flags = []
    if (truthy(params.rrna_operon)) {
        flags << '--rrna-operon'
    }
    if (truthy(params.hifi)) {
        flags << '--hifi'
    }
    if (truthy(params.single_strand)) {
        flags << '--single-strand'
    }
    if (truthy(params.use_hpc)) {
        flags << '--use-hpc'
    }
    if (truthy(params.low_polymorphism)) {
        flags << '--low-polymorphism'
    }
    if (truthy(params.mask_low_quality)) {
        flags << '--mask-low-quality'
    }
    return flags.join(' ')
}

process PREPARE_READS {
    tag "$sample"
    label 'savont_base'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("prepared_${sample}.fastq.gz"), emit: reads

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
        done < ${sample}.fastq.list | gzip -c > prepared_${sample}.fastq.gz
    else
        case "${reads}" in
            *.gz) cp "${reads}" prepared_${sample}.fastq.gz ;;
            *) gzip -c "${reads}" > prepared_${sample}.fastq.gz ;;
        esac
    fi
    """
}

process TRIM_18S {
    tag "$sample"
    label 'savont_trim18s'
    cpus { params.trim_threads as int }
    publishDir "${params.outdir}/intermediate/01_trim18s", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)
    path trim18s_dir

    output:
    tuple val(sample), path("${sample}.fastq.gz"), emit: trimmed
    tuple val(sample), path("trim18s.summary_${sample}.tsv"), emit: stats

    script:
    """
    python ${trim18s_dir}/trim18s.py \\
        ${reads} \\
        --outdir trim18s \\
        --threads ${task.cpus} \\
        --kingdom ${params.kingdom} \\
        --overwrite

    cp trim18s/prepared_${sample}.18S.fastq.gz ${sample}.fastq.gz
    cp trim18s/trim18s.summary.tsv trim18s.summary_${sample}.tsv
    """
}

process SAVONT_ASV_SAMPLE {
    tag "$sample"
    label 'savont_asv'
    cpus { params.savont_threads as int }
    publishDir "${params.outdir}/asv", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("savont_${sample}"), path("savont_${sample}/final_asvs.fasta"), path("savont_${sample}/feature-table.tsv"), path("savont_${sample}/final_clusters.tsv"), emit: asv_dir

    script:
    def flags = savontPresetFlags()
    """
    mkdir -p inputs
    ln -s ../${reads} inputs/${sample}.fastq.gz

    savont asv \\
        inputs/${sample}.fastq.gz \\
        --output-dir savont_${sample} \\
        --threads ${task.cpus} \\
        --min-read-length ${params.min_read_length} \\
        --max-read-length ${params.max_read_length} \\
        --quality-value-cutoff ${params.quality_value_cutoff} \\
        --minimum-base-quality ${params.minimum_base_quality} \\
        --min-cluster-size ${params.min_cluster_size} \\
        ${flags}
    """
}

process SAVONT_ASV_POOLED {
    label 'savont_asv'
    cpus { params.savont_threads as int }
    publishDir "${params.outdir}/asv", mode: 'copy', overwrite: true

    input:
    path reads

    output:
    tuple path("savont_pooled"), path("savont_pooled/final_asvs.fasta"), path("savont_pooled/feature-table.tsv"), path("savont_pooled/final_clusters.tsv"), emit: asv_dir

    script:
    def flags = savontPresetFlags()
    """
    mkdir -p inputs
    for fastq in ${reads}; do
        ln -s ../"\$fastq" inputs/"\$(basename "\$fastq")"
    done

    savont asv \\
        --pooled-samples inputs/*.fastq.gz \\
        --output-dir savont_pooled \\
        --threads ${task.cpus} \\
        --min-read-length ${params.min_read_length} \\
        --max-read-length ${params.max_read_length} \\
        --quality-value-cutoff ${params.quality_value_cutoff} \\
        --minimum-base-quality ${params.minimum_base_quality} \\
        --min-cluster-size ${params.min_cluster_size} \\
        ${flags}
    """
}

process SAVONT_CLASSIFY {
    tag "$sample"
    label 'savont_classify'
    cpus { params.taxonomy_threads as int }
    publishDir "${params.outdir}/classification", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(asv_dir)
    path database

    output:
    tuple val(sample), path("classification_${sample}"), emit: classification

    script:
    """
    savont ${params.classifier} \\
        --input-dir ${asv_dir} \\
        --db ${database} \\
        --output-dir classification_${sample} \\
        --threads ${task.cpus}
    """
}

process PR2_TAXONOMY {
    tag "$sample"
    label 'savont_taxonomy'
    cpus { params.taxonomy_threads as int }
    publishDir "${params.outdir}/taxonomy", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(asv_dir), path(final_asvs)
    path database

    output:
    tuple val(sample), path("pr2_${sample}/taxonomy.tsv"), path("pr2_${sample}/final_asvs.fasta"), emit: taxonomy

    script:
    """
    mkdir -p pr2_${sample}
    cp ${final_asvs} pr2_${sample}/final_asvs.fasta
    if grep -q '^>' ${final_asvs}; then
        vsearch \\
            --usearch_global ${final_asvs} \\
            --db ${database} \\
            --id ${params.db_id} \\
            --threads ${task.cpus} \\
            --blast6out pr2_${sample}/taxonomy.tsv \\
            --query_cov ${params.db_query_cov}
    else
        touch pr2_${sample}/taxonomy.tsv
    fi
    """
}

process PR2_TAXONOMY_TABLE {
    tag "$sample"
    label 'savont_taxonomy'
    publishDir "${params.outdir}/taxonomy", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(taxonomy), path(final_asvs)
    path scripts_dir

    output:
    tuple val(sample), path("pr2_${sample}"), emit: table

    script:
    """
    mkdir -p pr2_${sample}
    if [ -s ${taxonomy} ]; then
        python ${scripts_dir}/get_taxonomy_table.py \\
            -i ${taxonomy} \\
            -o pr2_${sample}/taxonomy_table.tsv
    else
        printf 'OTU\\tPident\\tAccession\\tDomain\\tSupergroup\\tDivision\\tSubdivision\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\tLength\\n' > pr2_${sample}/taxonomy_table.tsv
    fi
    """
}

process SAVONT_BIOLOGICAL_REPORT {
    label 'savont_report'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path trim_summaries
    path asv_dirs
    path taxonomy_tables
    path scripts_dir

    output:
    path "savont_biological_report.html", emit: report
    path "savont_reports", emit: report_dir

    script:
    """
    python ${scripts_dir}/savont_biological_report.py \\
        --trim-summaries ${trim_summaries} \\
        --asv-dirs ${asv_dirs} \\
        --taxonomy-tables ${taxonomy_tables} \\
        --output savont_biological_report.html \\
        --output-dir savont_reports
    """
}

process SAVONT_EXPORT {
    label 'savont_base'
    publishDir "${params.outdir}/export", mode: 'copy', overwrite: true

    input:
    path asv_dirs

    output:
    path "savont_export", emit: export_dir

    script:
    """
    savont export \\
        --input-dirs ${asv_dirs} \\
        --output-dir savont_export
    """
}

workflow {
    requireParam('samplesheet', params.samplesheet)

    def classifier = params.classifier.toString()
    if (!(classifier in ['none', 'classify', 'sintax'])) {
        error "--classifier must be one of: none, classify, sintax"
    }
    if (classifier != 'none') {
        requireParam('savont_db', params.savont_db)
    }

    samples_ch = samplesheetToChannel(params.samplesheet)
    trim18s_ch = Channel.value(file("${projectDir}/savont_niko", checkIfExists: true))

    PREPARE_READS(samples_ch)
    TRIM_18S(PREPARE_READS.out.reads, trim18s_ch)

    if (truthy(params.pooled_samples)) {
        pooled_reads_ch = TRIM_18S.out.trimmed.map { sample, reads -> reads }.collect()
        SAVONT_ASV_POOLED(pooled_reads_ch)
        asv_dirs_ch = SAVONT_ASV_POOLED.out.asv_dir.map { dir, final_asvs, feature_table, final_clusters -> tuple('pooled', dir) }
    } else {
        SAVONT_ASV_SAMPLE(TRIM_18S.out.trimmed)
        asv_dirs_ch = SAVONT_ASV_SAMPLE.out.asv_dir.map { sample, dir, final_asvs, feature_table, final_clusters -> tuple(sample, dir) }
    }

    if (classifier != 'none') {
        database_ch = Channel.value(file(params.savont_db, checkIfExists: true))
        SAVONT_CLASSIFY(asv_dirs_ch, database_ch)
    }

    if (params.pr2_db != null && params.pr2_db.toString().trim() != '') {
        pr2_database_ch = Channel.value(file(params.pr2_db, checkIfExists: true))
        scripts_ch = Channel.value(file("${projectDir}/scripts", checkIfExists: true))
        if (truthy(params.pooled_samples)) {
            pr2_input_ch = SAVONT_ASV_POOLED.out.asv_dir.map { dir, final_asvs, feature_table, final_clusters -> tuple('pooled', dir, final_asvs) }
        } else {
            pr2_input_ch = SAVONT_ASV_SAMPLE.out.asv_dir.map { sample, dir, final_asvs, feature_table, final_clusters -> tuple(sample, dir, final_asvs) }
        }
        PR2_TAXONOMY(pr2_input_ch, pr2_database_ch)
        if (truthy(params.enable_taxonomy_table)) {
            PR2_TAXONOMY_TABLE(PR2_TAXONOMY.out.taxonomy, scripts_ch)
            if (truthy(params.enable_biological_report)) {
                trim_summaries_ch = TRIM_18S.out.stats.map { sample, summary -> summary }.collect()
                report_asv_dirs_ch = asv_dirs_ch.map { sample, dir -> dir }.collect()
                taxonomy_tables_ch = PR2_TAXONOMY_TABLE.out.table.map { sample, table -> table }.collect()
                SAVONT_BIOLOGICAL_REPORT(trim_summaries_ch, report_asv_dirs_ch, taxonomy_tables_ch, scripts_ch)
            }
        }
    }

    if (truthy(params.export_qiime)) {
        export_input_ch = asv_dirs_ch.map { sample, dir -> dir }.collect()
        SAVONT_EXPORT(export_input_ch)
    }
}
