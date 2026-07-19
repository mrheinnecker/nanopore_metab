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
params.db = params.db ?: 'nt'
params.task = params.task ?: 'blastn'
params.evalue = params.evalue ?: '1e-20'
params.max_target_seqs = params.max_target_seqs ?: 20
params.container = params.container ?: "${projectDir}/containers/blast_remote.sif"

process NCBI_REMOTE_BLAST {
    tag sample
    label 'remote_blast'

    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(query)
    path plot_script

    output:
    tuple val(sample), path("${sample}.blast_archive.asn1"), emit: archive
    tuple val(sample), path("${sample}.blast_hits.tsv"), emit: hits
    tuple val(sample), path("${sample}.blast_report.txt"), emit: report
    tuple val(sample), path("${sample}.blast_dashboard.html"), emit: dashboard

    script:
    """
    if ! grep -q '^>' ${query}; then
        echo "Input does not contain any FASTA records: ${query}" >&2
        exit 1
    fi

    blastn \
        -query ${query} \
        -db ${params.db} \
        -remote \
        -task ${params.task} \
        -evalue ${params.evalue} \
        -max_target_seqs ${params.max_target_seqs} \
        -outfmt 11 \
        -out ${sample}.blast_archive.asn1

    printf 'query_id\\tquery_length\\tsubject_accession\\ttaxonomy_ids\\tscientific_name\\tcommon_name\\tpercent_identity\\talignment_length\\tmismatches\\tgaps\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tquery_coverage\\tevalue\\tbit_score\\tsubject_title\\n' \
        > ${sample}.blast_hits.tsv

    blast_formatter \
        -archive ${sample}.blast_archive.asn1 \
        -outfmt "6 qseqid qlen saccver staxids sscinames scomnames pident length mismatch gaps qstart qend sstart send qcovs evalue bitscore stitle" \
        >> ${sample}.blast_hits.tsv

    blast_formatter \
        -archive ${sample}.blast_archive.asn1 \
        -outfmt 0 \
        -out ${sample}.blast_report.txt

    awk -f ${plot_script} \
        -v sample="${sample}" \
        ${sample}.blast_hits.tsv \
        > ${sample}.blast_dashboard.html
    """

    stub:
    """
    touch ${sample}.blast_archive.asn1
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
    if (!(params.task.toString() in ['blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn'])) {
        error "--task must be a task supported by blastn"
    }
    if ((params.max_target_seqs as int) < 1) {
        error "--max_target_seqs must be at least 1"
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
    NCBI_REMOTE_BLAST(query_ch, plot_script_ch)
}
