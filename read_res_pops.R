
library(data.table)
library(dplyr)

options(datatable.fread.datatable=FALSE)


setwd("/broad/hptmp/btruong/HDP/pops/geshtn")

predres = read.delim("geshtn_0kb_pops_step2.preds")
gene = read.table("/broad/hptmp/btruong/HDP/pops/gene_annot_jun10.txt", header=T)
gene = gene %>% select(ENSGID, NAME)

predres = merge(predres, gene)
predres = predres[order(predres$PoPS_Score, decreasing=T),]
head(predres)


write.table(predres, "geshtn_0kb_pops_step2_GeneSymbol.preds", row.names=F, quote=F, sep="\t")


marginalsres = read.delim("geshtn_0kb_pops_step2.marginals")
idx = which(startsWith(marginalsres[,1], "mouse"))
if (length(idx)>0) 
	marginalsres = marginalsres[-idx,]
marginalsres = marginalsres[order(marginalsres$r2, decreasing=T),]
head(marginalsres)

write.table(marginalsres, "geshtn_0kb_pops_step2.marginals", row.names=F, quote=F, sep="\t")




# snploc = fread("geshtn.clumped")
# snploc_clumped = snploc %>% select("CHR", "SNP", "BP", "P")
# snploc_sig = snploc_clumped

# geshtn
snploc = fread("geshtn.snploc_hg19")
snploc_sig = snploc %>% filter(SNP %in% c("rs16998073", "rs13154066", "rs7139122", "rs260017"))


# geshtn
snploc = fread("geshtn.snploc_hg19")
snploc_sig = snploc %>% filter(SNP %in% c("rs17367504", "rs3754357", "rs1918969", "rs16998073", "rs2442752", "rs2508367", "rs10774624", "rs7318880", "rs1421085", "rs167479", "rs259983"))



# snploc_sig = snploc %>% filter(P<=5e-9)
# colnames(snploc_sig)[3] = "BP"

gene = read.table("/broad/hptmp/btruong/HDP/pops/gene_annot_jun10.txt", header=T)

priorgenes = unique(predres$NAME)[1:500]
predres$out = paste(predres$NAME, " (", round(predres$PoPS_Score,2), ")", sep="")

snpgene = NULL
allgenes = c()
for (i in 1:nrow(snploc_sig)) {
	if (i %% 50 == 0)
		print(i)
	snp = snploc_sig[i,]
	gene1 = gene %>% filter(CHR==snp$CHR & START >= snp$BP - 500000 & END <= snp$BP + 500000)
	gene1 = paste0(intersect(gene1$NAME, priorgenes), collapse=", ")
	snpgene = rbind(snpgene, cbind(snp, gene1))
	allgenes = unique(c(allgenes, gene1))
}
colnames(snpgene)[5] = "PoPSgene_within_500kb"
snpgene = snpgene[order(snpgene$P),]
head(snpgene)


write.table(snpgene, "geshtn_snpgene.500kb_clumped.txt", row.names=F, quote=F, sep="\t")

###################################

predres_gene = predres[match(priorgenes, predres$NAME),]

snpgene = NULL
allgenes = c()
for (i in 1:nrow(snploc_sig)) {
	# i=1
	if (i %% 50 == 0)
		print(i)
	snp = snploc_sig[i,]
	gene1 = gene %>%
		filter(CHR == snp$CHR) %>%
		mutate(dis = ifelse(END < snp$POS, snp$POS-END, START-snp$POS)) %>%
		mutate(out2 = paste0(NAME, " (", dis, ")"))
		
	gene1 = merge(gene1, predres %>% select(NAME, PoPS_Score, out))
	gene1 = gene1[order(gene1$PoPS_Score, decreasing=T),]
	
	common_gene = intersect(predres_gene$NAME, gene1$NAME)
	gene1 = gene1[match(common_gene, gene1$NAME),]
	
	gene1_prioritized = paste0(gene1$out, collapse=", ")
	gene1_prioritized_loc = paste0(gene1$out2, collapse=", ")
	snpgene = rbind(snpgene, cbind(snp, gene1_prioritized, gene1_prioritized_loc))
	allgenes = unique(c(allgenes, gene1_prioritized))
}
colnames(snpgene)[5:6] = c("PoPSgene_top", "PoPSgene_distance")
snpgene = snpgene[order(snpgene$P),]
head(snpgene)

write.table(snpgene, "geshtn_snpgene.added.txt", row.names=F, quote=F, sep="\t")



###################################


snpgene = NULL
allgenes = c()
for (i in 1:nrow(snploc_sig)) {
	# i=2
	if (i %% 50 == 0)
		print(i)
	snp = snploc_sig[i,]
	gene1 = gene %>% 
		filter(CHR==snp$CHR & 
			START >= snp$POS - 500000 & 
			END <= snp$POS + 500000) %>%
		mutate(dis = ifelse(END < snp$POS, snp$POS-END, START-snp$POS)) %>%
		mutate(out2 = paste0(NAME, " (", dis, ")"))
		
	gene1 = merge(gene1, predres %>% select(NAME, PoPS_Score, out))
	gene1 = gene1[order(gene1$PoPS_Score, decreasing=T),]
	idx = which(duplicated(gene1$NAME))
	if (length(idx) > 0)
		gene1 = gene1[-idx,]
	gene1 = gene1[1:min(5,nrow(gene1)),]
	
	gene1_prioritized = paste0(gene1$out, collapse=", ")
	gene1_prioritized_loc = paste0(gene1$out2, collapse=", ")
	snpgene = rbind(snpgene, cbind(snp, gene1_prioritized, gene1_prioritized_loc))
	allgenes = unique(c(allgenes, gene1_prioritized))
}
colnames(snpgene)[5:6] = c("PoPSgene_top", "PoPSgene_distance")
snpgene = snpgene[order(snpgene$P),]
head(snpgene)

write.table(snpgene, "geshtn_snpgene_500kbAroundLeadSNP_NotTop500.txt", row.names=F, quote=F, sep="\t")





