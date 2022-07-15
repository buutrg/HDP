
library(data.table)
library(stringr)
library(dplyr)

options(datatable.fread.datatable=FALSE)


genelist = read.table("/broad/hptmp/btruong/cS2G/data/critical_gene_sets/list_genes_qc.txt.gz", header=T)

# gene_trait = read.table("/broad/hptmp/btruong/cS2G/data/critical_gene_sets/critical_gene_sets_geshtn.txt")

genelist = genelist %>% 
	# filter(name %in% gene_trait[,1]) %>%
	select(name, chr, ensembl.start, ensembl.end)
colnames(genelist) = c("GENE", "CHR", "START", "END")

write.table(genelist, "genecord.txt", sep="\t", col.names=T, row.names=F, quote=F)