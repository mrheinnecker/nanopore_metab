  nextflow -C savont.config run savont.nf \
    -profile singularity \
    -work-dir /tmp/savont_300_work_v2 \
    --samplesheet /mnt/c/Users/rheinnec/Documents/taxseq/test_data/ont_300_reads_three_barcodes/samplesheet.tsv \
    --outdir /tmp/savont_300_results_v2 \
    --threads 2 \
    --min_read_length 1000 \
    --max_read_length 3300