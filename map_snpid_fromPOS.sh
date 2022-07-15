

wget ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
mkdir /medpop/esp2/btruong/Tools/hg19_dbSNP151
gunzip -c All_20180423.vcf.gz | vcf2bed --sort-tmpdir=/medpop/esp2/btruong/Tools/hg19_dbSNP151 --max-mem=80G - > hg19.dbSNP151.bed


zcat /medpop/esp2/honigberg/HDP_GWAS/MGI/X642.1.output.gz | awk -v OFS="\t" '(FNR!=1){ print $1, ($2 - 1), $2; }' | sort-bed - > positions.bed

bedmap --echo --echo-map-id --delim '\t' positions.bed /medpop/esp2/btruong/Tools/hg38.dbSNP151.bed > MGI.bed





wget ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
gunzip -c All_20180418.vcf.gz | vcf2bed --sort-tmpdir=/medpop/esp2/btruong/Tools/hg38_dbSNP151 --max-mem=80G - > hg38.dbSNP151.bed


