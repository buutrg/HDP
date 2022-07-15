
library(dplyr)

dd = read.table("/broad/hptmp/btruong/HDP/pops/preec/preec_0kb_pops_step2_GeneSymbol.preds", header=T)

all_genes = read.table("/broad/hptmp/btruong/cS2G/data/critical_gene_sets/list_genes_qc.txt.gz", header=T)

inter = intersect(dd$NAME, all_genes$name)

dd = dd %>% filter(NAME %in% inter)
dd = dd[1:floor(nrow(dd)*0.1),]
write.table(data.frame(dd$NAME,1), "/broad/hptmp/btruong/cS2G/data/critical_gene_sets/critical_gene_sets_preec.txt", row.names=F, col.names=F, sep="\t", quote=F)


# all_genes1 = all_genes
# all_genes1 = all_genes1 %>% select()

