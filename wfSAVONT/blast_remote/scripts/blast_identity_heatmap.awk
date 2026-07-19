BEGIN {
    FS = "\t"
    cell_w = 168
    cell_h = 72
    left = 260
    top = 80
}

NR == 1 { next }

{
    query = $1
    if (!(query in seen_query)) {
        seen_query[query] = 1
        queries[++query_count] = query
    }

    rank = ++hit_count[query]
    if (rank > max_rank) {
        max_rank = rank
    }

    key = query SUBSEP rank
    identity[key] = $7 + 0
    accession[key] = $3
    organism[key] = $5
    description[key] = $18
    alignment_length[key] = $8
    mismatches[key] = $9
    gaps[key] = $10
    coverage[key] = $15
}

function xml_escape(value, escaped) {
    escaped = value
    gsub(/&/, "\\&amp;", escaped)
    gsub(/</, "\\&lt;", escaped)
    gsub(/>/, "\\&gt;", escaped)
    gsub(/"/, "\\&quot;", escaped)
    return escaped
}

function heat_color(value, fraction, red, green) {
    fraction = (value - 70) / 30
    if (fraction < 0) fraction = 0
    if (fraction > 1) fraction = 1
    red = int(245 - 190 * fraction)
    green = int(75 + 110 * fraction)
    return sprintf("rgb(%d,%d,70)", red, green)
}

function shorten(value, limit) {
    if (length(value) <= limit) return value
    return substr(value, 1, limit - 3) "..."
}

END {
    if (max_rank < 1) max_rank = 1
    if (query_count < 1) query_count = 1

    width = left + max_rank * cell_w + 80
    height = top + query_count * cell_h + 100

    print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" width "\" height=\"" height "\" viewBox=\"0 0 " width " " height "\">"
    print "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
    print "<style>text{font-family:Arial,sans-serif}.label{font-size:12px}.cell{font-size:10px;fill:white}.title{font-size:18px;font-weight:bold}.hit{font-weight:bold}</style>"
    print "<text class=\"title\" x=\"20\" y=\"28\">" xml_escape(title) "</text>"
    print "<text class=\"label\" x=\"20\" y=\"50\">Rows: ASVs; columns: BLAST hit rank; color: percent identity; cov: query coverage; len: alignment length; mm: mismatches</text>"

    for (rank = 1; rank <= max_rank; rank++) {
        x = left + (rank - 1) * cell_w + cell_w / 2
        print "<text class=\"label\" text-anchor=\"middle\" x=\"" x "\" y=\"" top - 12 "\">" rank "</text>"
    }

    if (length(queries[1]) == 0) {
        print "<text class=\"label\" x=\"20\" y=\"" top + 20 "\">No BLAST hits</text>"
    } else {
        for (row = 1; row <= query_count; row++) {
            query = queries[row]
            y = top + (row - 1) * cell_h
            print "<text class=\"label\" text-anchor=\"end\" x=\"" left - 8 "\" y=\"" y + 16 "\">" xml_escape(query) "</text>"

            for (rank = 1; rank <= hit_count[query]; rank++) {
                key = query SUBSEP rank
                value = identity[key]
                x = left + (rank - 1) * cell_w
                display_description = description[key]
                if (display_description == "") display_description = organism[key]
                tooltip = query " | hit " rank " | accession " accession[key] " | identity " value "% | query coverage " coverage[key] "% | alignment length " alignment_length[key] " | mismatches " mismatches[key] " | gaps " gaps[key] " | " display_description
                print "<g><title>" xml_escape(tooltip) "</title>"
                print "<rect x=\"" x "\" y=\"" y "\" width=\"" cell_w - 4 "\" height=\"" cell_h - 4 "\" rx=\"3\" fill=\"" heat_color(value) "\"/>"
                print "<text class=\"cell hit\" x=\"" x + 6 "\" y=\"" y + 14 "\">" xml_escape(accession[key]) "</text>"
                print "<text class=\"cell\" x=\"" x + 6 "\" y=\"" y + 29 "\">identity " sprintf("%.2f", value) "%</text>"
                print "<text class=\"cell\" x=\"" x + 6 "\" y=\"" y + 43 "\">cov " coverage[key] "% | len " alignment_length[key] "</text>"
                print "<text class=\"cell\" x=\"" x + 6 "\" y=\"" y + 57 "\">mm " mismatches[key] " | gaps " gaps[key] " | " xml_escape(shorten(display_description, 18)) "</text></g>"
            }
        }
    }

    legend_y = top + query_count * cell_h + 38
    print "<text class=\"label\" x=\"20\" y=\"" legend_y "\">Identity:</text>"
    for (step = 0; step <= 30; step++) {
        value = 70 + step
        print "<rect x=\"" 80 + step * 8 "\" y=\"" legend_y - 13 "\" width=\"8\" height=\"12\" fill=\"" heat_color(value) "\"/>"
    }
    print "<text class=\"label\" x=\"80\" y=\"" legend_y + 16 "\">≤70%</text>"
    print "<text class=\"label\" text-anchor=\"end\" x=\"" 80 + 31 * 8 "\" y=\"" legend_y + 16 "\">100%</text>"
    print "</svg>"
}
