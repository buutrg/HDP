
library(data.table)
library(stringr)
library(dplyr)

options(datatable.fread.datatable=FALSE)


sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/InterPregGen/summary_stats_maternal_preeclampsia_with_rsids.txt")

newpos_lifted = fread("summary_stats_maternal_preeclampsia_with_rsids.txt_lifted.bed")

sumstat$newID = paste(sumstat$CHR, newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2], sep=":")
sumstat$newPOS = newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2]

fwrite(sumstat, file="summary_stats_maternal_preeclampsia_with_rsids_lifted.txt", sep="\t", row.names=F)


################################


sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/Estonia/O13_O14_O15_EstBB_200421_TL")
sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/Estonia/O14_O15_EstBB_080421_TL")
sumstat$locus = paste(sumstat$CHR, sumstat$POS, sep=":")

newpos_lifted = fread("O13_O14_O15_EstBB_200421_TL_lifted.bed")



sumstat$newID = paste(sumstat$CHR, newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2], sep=":")
sumstat$newPOS = newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2]

# fwrite(sumstat, file="O13_O14_O15_EstBB_200421_TL_lifted.txt", sep="\t", row.names=F)
fwrite(sumstat, file="O14_O15_EstBB_080421_TL_lifted.txt", sep="\t", row.names=F)


################################


sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/Estonia/O13_EstBB_080421_TL.gz")
sumstat = fread("/medpop/esp2/aniruddh/GH/ghtn/2022_03_01_GNH_binary_GWAS_firth_GNH0246_Preeclampsia_Eclampsia.regenie.gz")
sumstat$locus = paste(sumstat$CHROM, sumstat$GENPOS, sep=":")

newpos_lifted = fread("O13_EstBB_080421_TL.gz_lifted.bed")



sumstat$newID = paste(sumstat$CHR, newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2], sep=":")
sumstat$newPOS = newpos_lifted[match(sumstat$locus, newpos_lifted[,4]), 2]

# fwrite(sumstat, file="O13_O14_O15_EstBB_200421_TL_lifted.txt", sep="\t", row.names=F)
fwrite(sumstat, file="O13_EstBB_080421_TL_lifted.txt", sep="\t", row.names=F)


#############################

ref = fread("/medpop/esp2/btruong/Tools/hg38_common_chrpos.txt")
colnames(ref)[4]="newID"

# file = "HA_composite.composite.glm.logistic"
# sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/data/mssm/",file))


file = "BioMe_AA_Ec_saige_output.txt"
sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/HDP-mssm/",file))
# sumstat = sumstat %>% filter(!is.na(OR))

sumstat$newID = paste(sumstat$CHR, sumstat$POS, sep=":")

sumstat$newRSID = ref[match(sumstat$newID, ref$newID),3]

sumstat = sumstat %>% select(-newID)

fwrite(sumstat, paste0("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/",file ,"_newRSID"), row.names=F, sep="\t")










#################################################

ref = fread("/medpop/esp2/btruong/Tools/hg37_common_chrpos.txt")
colnames(ref)[4]="newID"

file = "O13_O14_O15_EstBB_200421_TL_lifted"
sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/data/",file, ".txt"))
sumstat$newID = paste(sumstat$CHR, sumstat$POS, sep=":")

sumstat$newRSID = ref[match(sumstat$newID, ref$newID),3]

sumstat = sumstat %>% select(-newID)

fwrite(sumstat, paste0("/medpop/esp2/btruong/Projects/HDP/data/",file ,"_newRSID.txt"), row.names=F, sep="\t")







file="2022_03_01_GNH_binary_GWAS_firth_GNH0246_Preeclampsia_Eclampsia.regenie"
file="2022_03_01_GNH_binary_GWAS_firth_GNH0248_Preeclampsia_Eclampsia_GestationalHypertension.regenie"
file="2022_03_01_GNH_binary_GWAS_firth_GNH0251_GestationalHypertension_noPreeclampsia_Eclampsia.regenie"
sumstat = fread(paste0("/medpop/esp2/aniruddh/GH/ghtn/", file, ".gz"))
sumstat$newID = paste(sumstat$CHROM, sumstat$GENPOS, sep=":")

sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.txt.gz"))
sumstat$newID = paste(sumstat$CHROM, sumstat$GENPOS, sep=":")


ref = fread("/medpop/esp2/btruong/Tools/hg38_common_chrpos.txt")
colnames(ref)[4]="newID"

sumstat$newRSID = ref[match(sumstat$newID, ref$newID),3]


sumstat = sumstat %>% select(-newID)

fwrite(sumstat, paste0("/medpop/esp2/btruong/Projects/HDP/data/",file ,"_newRSID.txt"), row.names=F, sep="\t")





sumstat  = fread("preec_metal_1.txt")
sumstat  = fread("geshtn_metal_1.txt")
sumstat  = fread("metal_preec_Asian_1.txt")
sumstat  = fread("/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_preec_European_1.txt")
sumstat  = fread("/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_preec_European.Hispanic.Asian_1.txt")
sumstat$chrpos = paste(sumstat[,"Chromosome"], sumstat[,"Position"], sep=":")

ref = fread("/medpop/esp2/btruong/Tools/hg38_common_chrpos.txt")
colnames(ref)[4]="chrpos"


sumstat$newRSID = ref$V3[match(sumstat$chrpos, ref$chrpos)]
sumstat$AF = sumstat$MinFreq
idx = which(sumstat$MinFreq == 0.5)
sumstat$AF[idx] = sumstat$MaxFreq[idx]
sumstat = sumstat %>% select(-chrpos)

sumstat$Allele1 = toupper(sumstat$Allele1)
sumstat$Allele2 = toupper(sumstat$Allele2)
sumstat = sumstat[which(!is.na(sumstat$newRSID)),]

fwrite(sumstat, "/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_preec_European.Hispanic.Asian_1_withRSID.txt", row.names=F, sep="\t", quote=F)

fwrite(sumstat, "/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_preec_European_1_newRSID.txt", row.names=F, sep="\t", quote=F)
fwrite(sumstat, "metal_preec_Asian_1_newRSID.txt", row.names=F, sep="\t", quote=F)



idx = match(sumstat$chrpos, ref$chrpos)
sumstat$newrsid = ref[idx,3]
sumstat$newpos = ref[idx,2]
sumstat$newchr = ref[idx,1]

sumstat1 = sumstat %>% filter(Chromosome != newchr)








options(scipen = 999)

sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_saved.txt.gz"))

sumstat$newID = paste(sumstat$CHR, sumstat$POS, sep=":")

newpos_lifted = fread("/medpop/esp2/btruong/Projects/HDP/data/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_lifted.bed")


sumstat$newPOS = newpos_lifted[match(sumstat$newID, newpos_lifted[,4]), 2]
sumstat$chrpos = paste(sumstat$CHR, sumstat$newPOS, sep=":")
# sumstat$POS = sumstat$newPOS
sumstat = sumstat %>% select(-newID)


options(scipen = 999)
fwrite(sumstat, "/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.txt.gz", row.names=F, sep="\t")




tmp = sumstat %>% select(CHR, POS)
tmp$newID = paste(tmp$CHR, tmp$POS, sep=":")
tmp$CHR = paste("chr", tmp$CHR, sep="")
tmp$end = tmp$POS + 1
tmp = tmp %>% select(CHR, POS, end, newID)
 
ss = format(tmp, scientific=F)
fwrite(ss, "/medpop/esp2/btruong/Projects/HDP/data/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed", row.names=F, col.names=F, sep="\t")


sumstat$chrpos = paste(sumstat[,"Chromosome"], sumstat[,"Position"], sep=":")






sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/InterPregGen/European_only/mat_all_chrALL_STERR_EU.1tbl.gz")
sumstat = fread("/medpop/esp2/honigberg/HDP_GWAS/InterPregGen/centAsian/mat_all_chrALL_STERR_ASIA.1tbl.gz")
sumstat = sumstat %>%
	rowwise() %>%
	mutate(rsid=unlist(strsplit(MarkerName, split=":"))[1])
sumstat = as.data.frame(sumstat)

ref = fread("/medpop/esp2/btruong/Tools/hg38_common_chrpos.txt")
colnames(ref)[4]="newID"

sumstat$newPOS = ref[match(sumstat$rsid, ref[,3]),2]
sumstat$CHR = ref[match(sumstat$rsid, ref[,3]),1]
sumstat$chrpos = ref[match(sumstat$rsid, ref[,3]),4]

head(sumstat)

fwrite(sumstat, "/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/mat_all_chrALL_STERR_EU.1tbl_lifted.gz", row.names=F, sep="\t")
fwrite(sumstat, "/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/mat_all_chrALL_STERR_ASIA.1tbl_lifted.gz", row.names=F, sep="\t")




sumstat = fread("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/mat_all_chrALL_STERR_EU.1tbl_lifted.gz")








