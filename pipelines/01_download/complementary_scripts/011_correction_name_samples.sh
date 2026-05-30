for file in *_1.fastq.gz; do
    mv -i "$file" "${file/_1.fastq.gz/.fastq.gz}"
done
