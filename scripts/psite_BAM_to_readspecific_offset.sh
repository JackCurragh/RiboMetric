samtools view $1 | awk '{
    read_name = $1
    ps_value = ""

    for (i = 12; i <= NF; i++) {
        if ($i ~ /^PS:i:/) {
            ps_value = substr($i, 6)
            break
        }
    }

    if (ps_value != "") {
        print read_name "\t" ps_value
    }
}'> $2
