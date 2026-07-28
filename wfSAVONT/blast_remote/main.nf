nextflow.enable.dsl = 2

/*
 * Standalone NCBI remote BLAST workflow for an existing ASV FASTA.
 *
 * The search is submitted once in BLAST archive format. blast_formatter then
 * creates both a machine-readable table and a human-readable alignment report
 * without submitting the sequences to NCBI a second time.
 */

params.query = params.query ?: null
params.asv_dir = params.asv_dir ?: null
params.outdir = params.outdir ?: 'results/remote_blast'
params.db = params.db ?: 'core_nt'
params.task = params.task ?: 'megablast'
params.evalue = params.evalue ?: '1e-20'
params.max_target_seqs = params.max_target_seqs ?: 20
params.per_query_retries = params.per_query_retries ?: 5
params.retry_delay_seconds = params.retry_delay_seconds ?: 20
params.dry_run_queries = params.dry_run_queries ?: 0
params.poll_interval_seconds = params.poll_interval_seconds ?: 20
params.max_wait_seconds = params.max_wait_seconds ?: 3600
params.container = params.container ?: "${projectDir}/containers/blast_remote.sif"

process NCBI_REMOTE_BLAST {
    tag "$sample"
    label 'remote_blast'

    publishDir params.outdir, mode: 'copy', overwrite: true

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
    if ! grep -q '^>' ${query}; then
        echo "Input does not contain any FASTA records: ${query}" >&2
        exit 1
    fi

    printf 'query_id\\tquery_length\\tsubject_accession\\ttaxonomy_ids\\tscientific_name\\tcommon_name\\tpercent_identity\\talignment_length\\tmismatches\\tgaps\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tquery_coverage\\tevalue\\tbit_score\\tsubject_title\\n' \
        > ${sample}.blast_hits.tsv
    printf 'query_file\\trid\\n' > ${sample}.blast_requests.tsv
    printf '{"BlastOutput2":[' > ${sample}.blast_results.json
    json_count=0
    : > ${sample}.blast_report.txt

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
        if (( ${params.dry_run_queries} > 0 && processed_queries >= ${params.dry_run_queries} )); then
            break
        fi
        processed_queries=\$((processed_queries + 1))

        query_name=\$(basename "\${query_file}" .fasta)
        query_json="\${query_name}.json"
        query_report="\${query_name}.txt"
        query_ok=false

        for attempt in \$(seq 1 ${params.per_query_retries}); do
            rm -f "\${query_json}" "\${query_report}"
            echo "Submitting \${query_name} to NCBI (attempt \${attempt}/${params.per_query_retries})" >&2

            if rid=\$(bash ${api_script} \
                "\${query_file}" \
                "${params.db}" \
                "${params.task}" \
                "${params.evalue}" \
                "${params.max_target_seqs}" \
                "${params.poll_interval_seconds}" \
                "${params.max_wait_seconds}" \
                "\${query_json}" \
                "\${query_report}"); then
                query_ok=true
                break
            fi

            echo "NCBI returned an incomplete/error response for \${query_name}" >&2
            if (( attempt < ${params.per_query_retries} )); then
                sleep \$(( ${params.retry_delay_seconds} * attempt ))
            fi
        done

        if [[ "\${query_ok}" != true ]]; then
            echo "No valid NCBI response for \${query_name} after ${params.per_query_retries} attempts" >&2
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
    printf ']}' >> ${sample}.blast_results.json

    awk -f ${plot_script} \
        -v sample="${sample}" \
        ${sample}.blast_hits.tsv \
        > ${sample}.blast_dashboard.html
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
    def hasQuery = params.query != null && params.query.toString().trim() != ''
    def hasAsvDir = params.asv_dir != null && params.asv_dir.toString().trim() != ''

    if (hasQuery && hasAsvDir) {
        error "Use either --query or --asv_dir, not both"
    }
    if (!hasQuery && !hasAsvDir) {
        error "Missing input: provide --query FASTA or --asv_dir SAVONT_RESULTS_DIR"
    }
    if (!(params.task.toString() in ['blastn', 'blastn-short', 'megablast'])) {
        error "--task must be blastn, blastn-short, or megablast for the NCBI URL API"
    }
    if ((params.max_target_seqs as int) < 1) {
        error "--max_target_seqs must be at least 1"
    }
    if ((params.per_query_retries as int) < 1) {
        error "--per_query_retries must be at least 1"
    }
    if ((params.retry_delay_seconds as int) < 0) {
        error "--retry_delay_seconds cannot be negative"
    }
    if ((params.dry_run_queries as int) < 0) {
        error "--dry_run_queries cannot be negative"
    }
    if ((params.poll_interval_seconds as int) < 10) {
        error "--poll_interval_seconds must be at least 10"
    }
    if ((params.max_wait_seconds as int) < 1) {
        error "--max_wait_seconds must be at least 1"
    }

    if (hasQuery) {
        query_file = file(params.query, checkIfExists: true)
        sample_name = query_file.parent.name.toString().replaceFirst(/^savont_/, '').replaceAll(/[^A-Za-z0-9_.-]/, '_')
        query_ch = Channel.value(tuple(sample_name, query_file))
    } else {
        if (!file(params.asv_dir, checkIfExists: true).isDirectory()) {
            error "--asv_dir must be a directory: ${params.asv_dir}"
        }
        query_ch = Channel
            .fromPath("${params.asv_dir}/**/final_asvs.fasta", checkIfExists: true)
            .map { fasta ->
                def sample = fasta.parent.name.toString().replaceFirst(/^savont_/, '').replaceAll(/[^A-Za-z0-9_.-]/, '_')
                tuple(sample, fasta)
            }
    }
    plot_script_ch = Channel.value(file("${projectDir}/scripts/blast_interactive_report.awk", checkIfExists: true))
    api_script_ch = Channel.value(file("${projectDir}/scripts/ncbi_blast_url.sh", checkIfExists: true))
    NCBI_REMOTE_BLAST(query_ch, plot_script_ch, api_script_ch)
}
