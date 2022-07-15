
library(data.table)
library(dplyr)

options(datatable.fread.datatable=FALSE)

ref = fread("/medpop/esp2/btruong/Tools/hg38_common_chrpos.txt")
colnames(ref) = c("chr", "pos37", "rsid", "newID")colnames(ref_hg37) = c("chr", "pos37", "rsid", "newID")



ref_hg37 = fread("/medpop/esp2/btruong/Tools/hg37_common_chrpos.txt")
colnames(ref_hg37) = c("chr", "pos37", "rsid", "newID")


trait = "geshtn"
# trait = "composite"
# trait = "geshtn"
# ancestry = "European.Hispanic.Asian"
ancestry = "European.Hispanic.Asian.African"
# ancestry = "European.Hispanic"
ancestry = "European"

sumstat = fread(paste0("/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_", trait, "_", ancestry, "_1.txt"))

sumstat = fread("/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/output2/preeclampsia_chr1-23.results.txt")





sumstat$AF = sumstat$MinFreq
idx = which(sumstat$MinFreq == 0.5)
sumstat$AF[idx] = sumstat$MaxFreq[idx]

sumstat$rsid = ref$V3[match(sumstat$MarkerName, ref$newID)]
sumstat = sumstat[-which(is.na(sumstat$rsid)),]

head(sumstat)
# sumstat = sumstat %>% select(-chrpos)


sumstat$pos37 = ref_hg37$pos37[match(sumstat$rsid, ref_hg37$rsid)]

head(sumstat)

fwrite(sumstat, paste0("/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_", trait, "_", ancestry, "_1_withRSID.txt"), sep="\t", row.names=F)











# tmp = unlist(sumstat["P-value"] / 2)
# z = qnorm(tmp)
# lambda = round(median(z^2) / 0.454, 4)
# lambda

sigsnps = sumstat %>% filter(get("P-value") <= 5e-9)
sigsnps = sigsnps[order(sigsnps$Chromosome),]



outdf = sigsnps  %>% select(Chromosome, Position, Allele2, Allele1)
outdf$newID = paste(outdf$Chromosome, outdf$Position, sep=":")



outdf$rsid = ref[match(outdf$newID, ref$newID),3]
outdf = outdf[complete.cases(outdf),]

fwrite(data.frame(outdf$rsid), file=paste0("sigMETAL_", trait, "_", ancestry, ".txt"), row.names=F, col.names=F)
fwrite(sigsnps %>% filter(MarkerName %in% outdf$newID), file=paste0("sigMETAL_", trait, "_", ancestry, "_full.txt"), row.names=F, sep="\t")

