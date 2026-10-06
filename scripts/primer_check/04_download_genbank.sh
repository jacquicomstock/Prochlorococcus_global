for f in strain_genbank/*.gb
do
    echo "===== $(basename "$f") ====="
    grep -B 2 -A 6 '16S ribosomal RNA' "$f"
done
