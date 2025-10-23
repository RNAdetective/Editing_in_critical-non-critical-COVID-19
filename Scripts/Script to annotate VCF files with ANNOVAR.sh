#!/bin/bash

# Directories
INPUT_DIR="/home/aiswarya/Share/annovar/jacusa_post_female"   # Directory with the VCF files
OUTPUT_DIR="/home/aiswarya/Share/annovar/annovar_gene_annotation"   # Store results inside Shared directory
DB_DIR="/home/aiswarya/Share/annovar/annovar/humandb"   # ANNOVAR database directory
ANNOVAR="/home/aiswarya/Share/annovar/annovar/table_annovar.pl"   # Path to table_annovar.pl

# Output directory 
mkdir -p "$OUTPUT_DIR"

# Loop through all VCF files in the input directory
for vcf_file in "$INPUT_DIR"/*.vcf; do
    # Extract filename without extension
    filename=$(basename "$vcf_file" .vcf)

    # It will perform gene-based annotation
    perl "$ANNOVAR" "$vcf_file" "$DB_DIR" \
        -buildver hg19 \
        -out "$OUTPUT_DIR/$filename" \
        -protocol refGene \
        -operation g \
        -remove -polish -nastring . \
        -vcfinput \
        --thread 4
done

echo "Gene annotation completed yey!!!"

