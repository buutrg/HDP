
library(dplyr)
library(data.table)
library(stringr)
options(datatable.fread.datatable=FALSE)


# setwd("/broad/hptmp/btruong/HDP/data")
setwd("/broad/hptmp/btruong/HDP/coloc")

filelist = list.files(".")
# filelist = filelist[which(endsWith(filelist, "peqtl.1e-4_p1.1e-04_p2.1e-04_p12.1e-05.txt"))]
# filelist = filelist[which(endsWith(filelist, "_pgwas.5e-9_peqtl.1e-3_p1.1e-04_p2.1e-04_p12.1e-05.txt"))]
# filelist = filelist[which(endsWith(filelist, "_pgwas.5e-9_peqtl.1e-3_p1.1e-04_p2.1e-04_p12.1e-05.txt"))]
# filelist = filelist[which(startsWith(filelist, "coloc_fullsumstat02_preec"))]
# filelist = filelist[which(startsWith(filelist, "coloc_varNearLeadVar_preec"))]
# filelist = filelist[which(startsWith(filelist, "coloc_AlleqtlNearLeadVar_preec"))]
# filelist = filelist[which(startsWith(filelist, "coloc_allnearbygenes_alleqtl_preec"))]
filelist = filelist[which(startsWith(filelist, "coloc_allnearbygenes_alleqtl_pgwas5e-9_preec"))]
# filelist = filelist[which(startsWith(filelist, "coloc_allnearbygenes_alleqtl_preec"))]

res_all = NULL
for (file in filelist) {
	tmp = read.delim(file, header=T)
	res_all = rbind(res_all, tmp)
}
res_all = res_all %>%
	rowwise() %>%
	mutate(CHR = unlist(strsplit(GWAS.Lead.SNP, split=":"))[1]) %>%
	mutate(ENSG_gene = unlist(strsplit(ENSG_gene, split="[.]"))[1])

res_all = data.frame(res_all)

res_all$Region = str_replace(res_all$Region, ":", "-")
res_all$Region = paste(res_all$CHR, res_all$Region, sep=":")

res_all = res_all[, c("Region", "Gene", "TopSNP_eQTL", "GWAS.Lead.SNP", "Tissue", "Coloc.H0", "Coloc.H1", "Coloc.H2", "Coloc.H3", "Coloc.H4", "eQTL.P", "GWAS.P")]

colnames(res_all) = c("Region", "Gene", "qtl_lead", "gwas_lead", "Tissue", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "qtl_pval", "gwas_pval")


res_all = res_all[order(res_all$PP.H4, decreasing=T),]
dim(res_all)
head(res_all, 20)

res_all = res_all %>%
	filter(!Tissue %in% c("Testis", "Prostate"))
res_all$Coloc_res = "No coloc"
res_all$Coloc_res = ifelse(res_all$PP.H4 > 0.5 & res_all$PP.H4 < 0.7, "Weak coloc", res_all$Coloc_res)
res_all$Coloc_res = ifelse(res_all$PP.H4 >= 0.7, "Strong coloc", res_all$Coloc_res)
head(res_all, 20)


# head(res_all)
# leadsnps = unique(res_all$gwas_lead)
# neargene = data.frame(
# 	"6:31383987" = "HLA-S",
# 	"4:80263187" = "FGF5",
# 	"3:169422102" = "MECOM",
# 	"19:11416089" = "RGL3",
# 	"20:59160402" = "ZNF831",
# 	"12:111395984" = "SH2B3",
# 	"16:53767042" = "FTO",
# 	"13:28564148" = "FLT1",
# 	"11:101399851" = "TRPC6",
# 	"5:32831564" = "NPR3",
# 	"20:59143668" = "ZNF831",
# 	"12:122861523" = "HIP1R",
# 	"12:53063801" = "TNS2",
# 	"1:41430071" = "MTHFR",
# 	check.names=F
# 	)

# res_all$Locus = as.character(neargene[, res_all$gwas_lead])
# res_all = res_all %>%
# 	relocate(Locus, .before="Gene")
# head(res_all)


# write.csv(res_all, "coloc_preec.csv", row.names=F)
write.csv(res_all, "coloc_preec_p5e-9.csv", row.names=F)



