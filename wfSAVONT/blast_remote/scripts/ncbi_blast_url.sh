#!/usr/bin/env bash
set -euo pipefail

if (( $# != 9 )); then
    echo "Usage: ncbi_blast_url.sh QUERY DB TASK EVALUE MAX_HITS POLL_SECONDS MAX_WAIT_SECONDS JSON_OUT REPORT_OUT" >&2
    exit 2
fi

query=$1
database=$2
task=$3
evalue=$4
max_hits=$5
poll_seconds=$6
max_wait_seconds=$7
json_out=$8
report_out=$9
endpoint=https://blast.ncbi.nlm.nih.gov/Blast.cgi

submit_response=$(mktemp)
status_response=$(mktemp)
trap 'rm -f "$submit_response" "$status_response"' EXIT

megablast=off
short_query=false
if [[ "$task" == megablast ]]; then
    megablast=on
elif [[ "$task" == blastn-short ]]; then
    short_query=true
fi

curl --fail --silent --show-error --max-time 120 \
    --retry 3 --retry-delay 5 \
    "$endpoint" \
    --data-urlencode "CMD=Put" \
    --data-urlencode "PROGRAM=blastn" \
    --data-urlencode "DATABASE=$database" \
    --data-urlencode "MEGABLAST=$megablast" \
    --data-urlencode "SHORT_QUERY_ADJUST=$short_query" \
    --data-urlencode "EXPECT=$evalue" \
    --data-urlencode "HITLIST_SIZE=$max_hits" \
    --data-urlencode "QUERY@$query" \
    > "$submit_response"

rid=$(sed -n 's/^[[:space:]]*RID = \([^[:space:]]*\).*/\1/p' "$submit_response" | head -n 1)
if [[ -z "$rid" ]]; then
    echo "NCBI did not return a request ID (RID)" >&2
    grep -E 'Message ID#|Error:|Message:' "$submit_response" >&2 || true
    exit 75
fi
echo "NCBI RID: $rid" >&2

elapsed=0
while true; do
    sleep "$poll_seconds"
    elapsed=$((elapsed + poll_seconds))

    curl --fail --silent --show-error --max-time 60 \
        --retry 3 --retry-delay 5 \
        "$endpoint?CMD=Get&RID=$rid&FORMAT_OBJECT=SearchInfo" \
        > "$status_response"

    if grep -q 'Status=READY' "$status_response"; then
        break
    elif grep -q 'Status=FAILED' "$status_response"; then
        echo "NCBI search failed for RID $rid" >&2
        exit 75
    elif grep -q 'Status=UNKNOWN' "$status_response"; then
        echo "NCBI no longer recognizes RID $rid" >&2
        exit 75
    fi

    if (( elapsed >= max_wait_seconds )); then
        echo "Timed out after ${max_wait_seconds}s waiting for RID $rid" >&2
        exit 75
    fi
    echo "RID $rid is still queued/running (${elapsed}s)" >&2
done

curl --fail --silent --show-error --max-time 120 \
    --retry 3 --retry-delay 5 \
    "$endpoint?CMD=Get&RID=$rid&FORMAT_TYPE=JSON2_S" \
    > "$json_out"

curl --fail --silent --show-error --max-time 120 \
    --retry 3 --retry-delay 5 \
    "$endpoint?CMD=Get&RID=$rid&FORMAT_TYPE=Text" \
    > "$report_out"

printf '%s\n' "$rid"
