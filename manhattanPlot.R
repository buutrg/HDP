




library(ggplot2)
library(dplyr)
library(qqman)
library(data.table)




pos_list = as.data.frame(fread("/medpop/esp2/btruong/Projects/HDP/summary_all_annot.txt"))
pos_list_save = pos_list


pos_list = pos_list_save
colnames(pos_list) = c("CHR", "POS", "SNP")

pos_list[,1] = as.numeric(pos_list[,1])
pos_list[,2] = as.numeric(pos_list[,2])
pos_list[,3] = as.character(pos_list[,3])




unique(pos_list[,1])


trait = "geshtn"
trait = "comp_htn_preg"
trait = "preeclampsia"


# gwasResults1 = as.data.frame(fread(paste0("/medpop/esp2/btruong/Projects/HDP/sumstat_meta_analysis/meta_", trait,".txt")))
# gwasResults = as.data.frame(fread(paste0("/medpop/esp2/btruong/Projects/HDP/meta_", trait,"1.txt")))
gwasResults = as.data.frame(fread(paste0("/medpop/esp2/btruong/Projects/HDP/sumstat_meta_analysis/meta_", trait,"_finngenf6.txt")))
dim(gwasResults)


idx = match(gwasResults[,1], pos_list[,3])
gwasResults["CHR"] = pos_list[idx,1]
gwasResults["POS"] = pos_list[idx,2]

gwasResults_missing = gwasResults[!complete.cases(gwasResults),]


gwasResults = gwasResults[complete.cases(gwasResults),]
gwasResults = gwasResults[gwasResults$CHR <= 22,]
gwasResults_p = gwasResults[gwasResults["P-value"] <= 0.05, ]
gwasResults_p = gwasResults_p[gwasResults_p["Freq1"] >= 0.05 | gwasResults_p["P-value"] <= 5e-8, ]
gwasResults_p = gwasResults_p %>%
	rename(any_of(c(P='P-value')))

dim(gwasResults_p)

gwasResults_sig = gwasResults_p[gwasResults_p$P <= 5e-8,]
gwasResults_sig = gwasResults_sig[order(gwasResults_sig$P), ]
dim(gwasResults_sig)



write.table(gwasResults_sig[,1], paste0("gwasResults_sig_snplist", trait, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(gwasResults_sig, paste0("gwasResults_sig", trait, ".txt"), row.names=F, quote=F, sep="\t")


png(paste0("manhat_", trait,"_finngenf6.png"))
manhattan(gwasResults_p, chr="CHR", bp="POS", snp="MarkerName", p="P", annotatePval = 0.01)
dev.off()

pdf(paste0("manhat_", trait,"_finngenf6.pdf"))
manhattan(gwasResults_p, chr="CHR", bp="POS", snp="MarkerName", p="P", annotatePval = 0.01)
dev.off()



##########################

snp_annot = as.data.frame(fread("significant_snp_preeclampsia_haploreg.csv"))









