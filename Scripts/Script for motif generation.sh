## convert the vcf file to bed file
awk -F"," 'BEGIN{OFS="\t"} {if(NR>1) printf "%s\t%d\t%d\t%s\t.\t.\n", $2, $3-1, $3, $12}' Present_in_REDI_AG.csv > HC_Present_in_REDI_AG.bed

## Extract the flanking 5nt region and save as fasta 
bedtools slop -b 5 -i HC_Present_in_REDI_AG.bed -g Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.fai | bedtools getfasta -fi Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa -bed - -fo HC_Present_in_REDI_AG.fasta

## Use weblogo to make the motif logo
weblogo -f HC_Present_in_REDI_AG.fasta -D fasta -o HC_Present_in_REDI_AG_logo.png -F png -U probability -s large --fineprint ""
