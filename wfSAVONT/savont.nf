nextflow.enable.dsl = 2

/*
 * Standalone Savont workflow.
 *
 * Input FASTQ files are normalized per sample, trimmed to 18S while preserving
 * FASTQ qualities, and then processed by savont asv. Optional classification
 * can be run against a pre-downloaded Savont database directory.
 */

params.samplesheet = params.samplesheet ?: null
params.barcode_dir = params.barcode_dir ?: null
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
params.blast_queries = params.blast_queries ?: 0
params.blast_db = params.blast_db ?: 'core_nt'
params.blast_task = params.blast_task ?: 'megablast'
params.blast_evalue = params.blast_evalue ?: '1e-20'
params.blast_max_target_seqs = params.blast_max_target_seqs ?: 20
params.blast_per_query_retries = params.blast_per_query_retries ?: 5
params.blast_retry_delay_seconds = params.blast_retry_delay_seconds ?: 20
params.blast_poll_interval_seconds = params.blast_poll_interval_seconds ?: 20
params.blast_max_wait_seconds = params.blast_max_wait_seconds ?: 3600

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

def inputDirToChannel(inputDir) {
    Channel
        .fromPath("${inputDir}/*", checkIfExists: true, type: 'any')
        .filter { inputPath ->
            inputPath.isDirectory() || inputPath.getName().toString() ==~ /.*\.(fastq|fq)(\.gz)?$/
        }
        .map { inputPath ->
            def sample = inputPath.isDirectory()
                ? inputPath.getName().toString()
                : inputPath.getName().toString().replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample, inputPath)
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

process NCBI_REMOTE_BLAST {
    tag "$sample"
    label 'savont_blast'
    publishDir "${params.outdir}/blast", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(query)
    path plot_script
    path api_script

    output:
    tuple val(sample), path("${sample}.blast_results.json"), emit: json
    tuple val(sample), path("${sample}.blast_requests.tsv"), emit: requests
    tuple val(sample), path("${sample}.blast_hits.tsv"), emit: hits
    tuple val(sample), path("${sample}.blast_report.txt"), emit: report
    tuple val(sample), path("${sample}.blast_dashboard.html"), emit: dashboard

    script:
    """
    printf 'query_id\\tquery_length\\tsubject_accession\\ttaxonomy_ids\\tscientific_name\\tcommon_name\\tpercent_identity\\talignment_length\\tmismatches\\tgaps\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tquery_coverage\\tevalue\\tbit_score\\tsubject_title\\n' > ${sample}.blast_hits.tsv
    printf 'query_file\\trid\\n' > ${sample}.blast_requests.tsv
    printf '{"BlastOutput2":[' > ${sample}.blast_results.json
    json_count=0
    : > ${sample}.blast_report.txt

    if grep -q '^>' ${query}; then
        mkdir split_queries
        awk '
            /^>/ {
                if (output) close(output)
                count++
                output = sprintf("split_queries/query_%06d.fasta", count)
            }
            { print > output }
        ' ${query}

        processed_queries=0
        for query_file in split_queries/query_*.fasta; do
            if (( processed_queries >= ${params.blast_queries} )); then
                break
            fi
            processed_queries=\$((processed_queries + 1))

            query_name=\$(basename "\${query_file}" .fasta)
            query_json="\${query_name}.json"
            query_report="\${query_name}.txt"
            query_ok=false

            for attempt in \$(seq 1 ${params.blast_per_query_retries}); do
                rm -f "\${query_json}" "\${query_report}"
                echo "Submitting \${query_name} to NCBI (attempt \${attempt}/${params.blast_per_query_retries})" >&2

                if rid=\$(bash ${api_script} \
                    "\${query_file}" \
                    "${params.blast_db}" \
                    "${params.blast_task}" \
                    "${params.blast_evalue}" \
                    "${params.blast_max_target_seqs}" \
                    "${params.blast_poll_interval_seconds}" \
                    "${params.blast_max_wait_seconds}" \
                    "\${query_json}" \
                    "\${query_report}"); then
                    query_ok=true
                    break
                fi

                echo "NCBI returned an incomplete/error response for \${query_name}" >&2
                if (( attempt < ${params.blast_per_query_retries} )); then
                    sleep \$(( ${params.blast_retry_delay_seconds} * attempt ))
                fi
            done

            if [[ "\${query_ok}" != true ]]; then
                echo "No valid NCBI response for \${query_name} after ${params.blast_per_query_retries} attempts" >&2
                exit 75
            fi

            printf '%s\\t%s\\n' "\${query_name}" "\${rid}" >> ${sample}.blast_requests.tsv
            if (( json_count > 0 )); then printf ',' >> ${sample}.blast_results.json; fi
            jq -c '.BlastOutput2[]' "\${query_json}" >> ${sample}.blast_results.json
            json_count=\$((json_count + 1))

            jq -r '
                .BlastOutput2[].report.results.search as \$search
                | \$search.hits[]?
                | .description[0] as \$description
                | .hsps[0] as \$hsp
                | [
                    (\$search.query_title | split(" ")[0]),
                    \$search.query_len,
                    \$description.accession,
                    (\$description.taxid // ""),
                    (\$description.sciname // ""),
                    "",
                    (100 * \$hsp.identity / \$hsp.align_len),
                    \$hsp.align_len,
                    (\$hsp.align_len - \$hsp.identity - \$hsp.gaps),
                    \$hsp.gaps,
                    \$hsp.query_from,
                    \$hsp.query_to,
                    \$hsp.hit_from,
                    \$hsp.hit_to,
                    (100 * ((\$hsp.query_to - \$hsp.query_from) | fabs + 1) / \$search.query_len),
                    \$hsp.evalue,
                    \$hsp.bit_score,
                    \$description.title
                  ] | @tsv
            ' "\${query_json}" >> ${sample}.blast_hits.tsv
            cat "\${query_report}" >> ${sample}.blast_report.txt
            printf '\\n\\n' >> ${sample}.blast_report.txt
        done
    else
        echo "No ASV records detected in ${query}; skipping NCBI remote BLAST" >&2
        printf 'No ASV records detected; NCBI remote BLAST was not run.\\n' > ${sample}.blast_report.txt
    fi
    printf ']}' >> ${sample}.blast_results.json

    awk -f ${plot_script} -v sample="${sample}" ${sample}.blast_hits.tsv > ${sample}.blast_dashboard.html
    """

    stub:
    """
    printf '{"BlastOutput2":[]}' > ${sample}.blast_results.json
    printf 'query_file\\trid\\n' > ${sample}.blast_requests.tsv
    printf 'query_id\\tquery_length\\tsubject_accession\\ttaxonomy_ids\\tscientific_name\\tcommon_name\\tpercent_identity\\talignment_length\\tmismatches\\tgaps\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tquery_coverage\\tevalue\\tbit_score\\tsubject_title\\n' > ${sample}.blast_hits.tsv
    printf 'BLASTN stub output\\n' > ${sample}.blast_report.txt
    printf '<!doctype html><title>Stub BLAST dashboard</title>\\n' > ${sample}.blast_dashboard.html
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

    def classifier = params.classifier.toString()
    if (!(classifier in ['none', 'classify', 'sintax'])) {
        error "--classifier must be one of: none, classify, sintax"
    }
    if (classifier != 'none') {
        requireParam('savont_db', params.savont_db)
    }
    if ((params.blast_queries as int) < 0) {
        error "--blast_queries cannot be negative"
    }
    if (!(params.blast_task.toString() in ['blastn', 'blastn-short', 'megablast'])) {
        error "--blast_task must be blastn, blastn-short, or megablast"
    }
    if ((params.blast_max_target_seqs as int) < 1) {
        error "--blast_max_target_seqs must be at least 1"
    }
    if ((params.blast_per_query_retries as int) < 1) {
        error "--blast_per_query_retries must be at least 1"
    }
    if ((params.blast_retry_delay_seconds as int) < 0) {
        error "--blast_retry_delay_seconds cannot be negative"
    }
    if ((params.blast_poll_interval_seconds as int) < 10) {
        error "--blast_poll_interval_seconds must be at least 10"
    }
    if ((params.blast_max_wait_seconds as int) < 1) {
        error "--blast_max_wait_seconds must be at least 1"
    }

    samples_ch = hasBarcodeDir ? inputDirToChannel(params.barcode_dir) : (file(params.samplesheet).isDirectory() ? inputDirToChannel(params.samplesheet) : samplesheetToChannel(params.samplesheet))
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

    if ((params.blast_queries as int) > 0) {
        if (truthy(params.pooled_samples)) {
            blast_input_ch = SAVONT_ASV_POOLED.out.asv_dir.map { dir, final_asvs, feature_table, final_clusters -> tuple('pooled', final_asvs) }
        } else {
            blast_input_ch = SAVONT_ASV_SAMPLE.out.asv_dir.map { sample, dir, final_asvs, feature_table, final_clusters -> tuple(sample, final_asvs) }
        }
        blast_plot_ch = Channel.value(file("${projectDir}/blast_remote/scripts/blast_interactive_report.awk", checkIfExists: true))
        blast_api_ch = Channel.value(file("${projectDir}/blast_remote/scripts/ncbi_blast_url.sh", checkIfExists: true))
        NCBI_REMOTE_BLAST(blast_input_ch, blast_plot_ch, blast_api_ch)
    }
}
