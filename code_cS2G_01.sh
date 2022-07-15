conda activate ldsc

ldsc_folder=/medpop/esp2/btruong/Tools/ldsc

use Python-2.7

trait=preec
SUMSTATS=/broad/hptmp/btruong/cS2G/data/sumstat/preec_European.Hispanic.Asian.txt
allgene_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/allgenes/allgenes.
training_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/training.
training_cs2g_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/cS2G.
allgene_cs2g_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/cS2G.

valid_pref=/broad/hptmp/btruong/cS2G/data/trait_annot

FRQFILE=/broad/hptmp/btruong/cS2G/data/1000G_Phase3_frq/1000G.EUR.QC. 
WGTFILE=/broad/hptmp/btruong/cS2G/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
BASELINELD=/broad/hptmp/btruong/cS2G/data/baselineLD_v2.2/baselineLD.
# training

python ${ldsc_folder}/ldsc.py \
	--h2 $SUMSTATS \
	--ref-ld-chr $BASELINELD,${training_pref},${allgene_pref} \
	--frqfile-chr $FRQFILE \
	--w-ld-chr $WGTFILE \
	--overlap-annot --print-coefficients --print-delete-vals \
	--out sldsc_baselineS2G_outputs/$trait.baselineS2G.training

echo "Finished step 1"

python ${ldsc_folder}/ldsc.py \
	--h2 $SUMSTATS \
	--ref-ld-chr $BASELINELD,${training_pref},${allgene_pref},${cs2g_pref},${allgene_cs2g_pref} \
	--frqfile-chr $FRQFILE \
	--w-ld-chr $WGTFILE \
	--overlap-annot --print-coefficients --print-delete-vals \
	--out sldsc_baselineS2G_outputs/$trait.baselineS2G_cS2G.training

echo "Finished step 2"

# validation

python ${ldsc_folder}/ldsc.py \
	--h2 $SUMSTATS \
	--ref-ld-chr $BASELINELD,${allgene_pref},${valid_pref}/${trait}. \
	--frqfile-chr $FRQFILE \
	--w-ld-chr $WGTFILE \
	--overlap-annot --print-coefficients --print-delete-vals \
	--out sldsc_baselineS2G_outputs/$trait.baselineS2G.validation

echo "Finished step 3"

python ${ldsc_folder}/ldsc.py \
	--h2 $SUMSTATS \
	--ref-ld-chr $BASELINELD,${allgene_pref},${cs2g_pref},${allgene_cs2g_pref},${valid_pref}/${trait}. \
	--frqfile-chr $FRQFILE \
	--w-ld-chr $WGTFILE \
	--overlap-annot --print-coefficients --print-delete-vals \
	--out sldsc_baselineS2G_outputs/$trait.baselineS2G_cS2G.validation

echo "Finished step 4"




Rscript /broad/hptmp/btruong/cS2G/data/scripts/sumannots_baselineS2G.r $BASELINELD,validation_${trait}/validation_${trait}.,allgenes/allgenes. $FRQFILE validation_${trait}/validation_${trait}.sumannots

Rscript /broad/hptmp/btruong/cS2G/data/scripts/sumannots_baselineS2G.r $BASELINELD,${allgene_pref},${valid_pref} $FRQFILE ${valid_pref}/validation_${trait}.sumannots

Rscript /broad/hptmp/btruong/cS2G/data/scripts/sumannots_baselineS2G.r $BASELINELD,validation_${trait}/validation_${trait}.,validation_${trait}/cS2G.,allgenes/allgenes.,allgenes/cS2G. $FRQFILE validation_${trait}/cS2G.sumannots

Rscript /broad/hptmp/btruong/cS2G/data/scripts/sumannots_baselineS2G.r $BASELINELD,${allgene_pref},${valid_pref},${allgene_cs2g_pref} $FRQFILE ${valid_pref}/validation_${trait}.sumannots



##########################################################

gene_set_file=/broad/hptmp/btruong/cS2G/data/critical_gene_sets/critical_gene_sets_geshtn_geneset.txt
# annot=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/allgenes/allgenes.22.annot.gz
annot=/broad/hptmp/btruong/cS2G/data/trait_annot/geshtn.22.annot.gz
bimfile=/broad/hptmp/btruong/cS2G/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22.bim
gene_coord="/broad/hptmp/btruong/cS2G/data/genecord.txt"


ldsc_folder=/medpop/esp2/btruong/Tools/ldsc


python3 -m pdb ${ldsc_folder}/make_annot.py \
	--gene-set-file ${gene_set_file} \
	--gene-coord-file ${gene_coord} \
	--windowsize 100000 \
	--bimfile ${bimfile} \
	--annot-file ${annot}



with open("bin.dat", "wb") as f:
	pickle.dump(iter_bim, "bin.dat")



./ldsc.py \
--print-snps /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/list.txt \
--ld-wind-cm 1.0 \
--out ../cS2G/baselineS2Gannots/allgenes/allgenes.7 \
--bfile /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.7 \
--thin-annot  \
--annot ../cS2G/baselineS2Gannots/allgenes/allgenes.7.annot.gz \
--l2








##########################################################


#!/bin/sh
#$ -pe smp 2 -R y -binding linear:2
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=20:00:00
#$ -N make_annot

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3
conda init bash

conda activate ldsc

cd /broad/hptmp/btruong/cS2G/data
use BEDTools
use Python-2.7

trait=preec

# chr=$SGE_TASK_ID
# chr=22

gene_set_file=/broad/hptmp/btruong/cS2G/data/critical_gene_sets/critical_gene_sets_${trait}_geneset.txt
# annot=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/allgenes/allgenes.22.annot.gz
annot=/broad/hptmp/btruong/cS2G/data/trait_annot/${trait}.${chr}.annot.gz
# bimprefix=/broad/hptmp/btruong/cS2G/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}
bimprefix=/broad/hptmp/btruong/cS2G/data/hm3snp/hm3_1kg_chr${chr}
gene_coord="/broad/hptmp/btruong/cS2G/data/genecord.txt"


ldsc_folder=/medpop/esp2/btruong/Tools/ldsc


# python3 ${ldsc_folder}/make_annot.py \
# 	--gene-set-file ${gene_set_file} \
# 	--gene-coord-file ${gene_coord} \
# 	--windowsize 100000 \
# 	--bimfile ${bimprefix}.bim \
# 	--annot-file ${annot}

# out=/broad/hptmp/btruong/cS2G/data/trait_annot/${trait}.${chr}
# python ${ldsc_folder}/ldsc.py \
# 	--l2 \
# 	--bfile ${bimprefix} \
# 	--ld-wind-kb 500 \
# 	--chunk-size 500000 \
# 	--annot ${annot} \
# 	--thin-ann \
# 	--out ${out}


SUMSTATS=/broad/hptmp/btruong/cS2G/data/sumstat/preec_European.Hispanic.Asian.txt
allgene_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/allgenes/allgenes.
training_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/training.
training_cs2g_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/cS2G.
allgene_cs2g_pref=/broad/hptmp/btruong/cS2G/data/baselineS2Gannots/training/cS2G.


valid_pref=/broad/hptmp/btruong/cS2G/data/trait_annot

FRQFILE=/broad/hptmp/btruong/cS2G/data/1000G_Phase3_frq/1000G.EUR.QC. 
WGTFILE=/broad/hptmp/btruong/cS2G/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
BASELINELD=/broad/hptmp/btruong/cS2G/data/baselineLD_v2.2/baselineLD.


out=/broad/hptmp/btruong/cS2G/data/sldsc_baselineS2G_outputs/${trait}.baselineS2G.validation
python ${ldsc_folder}/ldsc.py \
	--h2 $SUMSTATS \
	--ref-ld-chr $BASELINELD,${allgene_pref},${valid_pref}/${trait}. \
	--frqfile-chr $FRQFILE \
	--w-ld-chr $WGTFILE \
	--overlap-annot --print-coefficients --print-delete-vals \
	--out ${out}

##########################################################
