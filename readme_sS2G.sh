#Last edit: 4/5/22 by Steven Gazal

#This file and directory contain the lines of code to replicate the results of Gazal et al. medrxiv

#Note that the name format of the S2G-derived SNP annotations should be
#- S2G.geneset for the critical gene set S2G-derived SNP annotation
#- S2G.allgenes for the S2G-derived SNP annotation for all 19,995 genes
#- for S2G annotations linking all SNPs to at least one gene (such as Closest TSS), no .allgenes annotation is required; instead, the base annotation will be used as an equivalent annotation

##########################################################################################
#1) s-ldsc outputs for all the annotations of the baseline-S2G and the baseline-S2G model with the cS2G-derived SNP annotation (available in the baselineS2Gannots/sldsc_baselineS2G_outputs directory of https://alkesgroup.broadinstitute.org/cS2G/) were obtained using the following code:
SUMSTAT_PATH=/PATH/sumstats_63 # this variable is the path to the sumstats_63 directory available at https://alkesgroup.broadinstitute.org/cS2G/
BASELINELD=/PATH/baselineLD_v2.2/baselineLD. # this variable is the path to the baseline-LD v2.2 files available at https://alkesgroup.broadinstitute.org/LDSCORE/
FRQFILE=/PATH/1000G.EUR.QC. # this variable is the path to the reference files available at https://alkesgroup.broadinstitute.org/LDSCORE/
WGTFILE=/PATH/weights.hm3_noMHC. # this variable is the path to the weight files available at https://alkesgroup.broadinstitute.org/LDSCORE/
while read LINE; do
    ARR=( $LINE ); TRAIT=${ARR[0]}
    SUMSTATS=SUMSTAT_PATH/$TRAIT.sumstats.gz
    python ldsc.py \
        --h2 $SUMSTATS \
        --ref-ld-chr $BASELINELD,training/training.,allgenes/allgenes. \
        --frqfile-chr $FRQFILE \
        --w-ld-chr $WGTFILE \
        --overlap-annot --print-coefficients --print-delete-vals \
        --out sldsc_baselineS2G_outputs/$TRAIT.baselineS2G.training
    python ldsc.py \
        --h2 $SUMSTATS \
        --ref-ld-chr $BASELINELD,validation_${TRAIT}/validation_${TRAIT}.,allgenes/allgenes. \
        --frqfile-chr $FRQFILE \
        --w-ld-chr $WGTFILE \
        --overlap-annot --print-coefficients --print-delete-vals \
        --out sldsc_baselineS2G_outputs/$TRAIT.baselineS2G.validation
    python ldsc.py \
        --h2 $SUMSTATS \
        --ref-ld-chr $BASELINELD,training/training.,training/cS2G.,allgenes/allgenes.,allgenes/cS2G. \
        --frqfile-chr $FRQFILE \
        --w-ld-chr $WGTFILE \
        --overlap-annot --print-coefficients --print-delete-vals \
        --out sldsc_baselineS2G_outputs/$TRAIT.baselineS2G_cS2G.training
    python ldsc.py \
        --h2 $SUMSTATS \
        --ref-ld-chr $BASELINELD,validation_${TRAIT}/validation_${TRAIT}.,validation_${TRAIT}/cS2G.,allgenes/allgenes.,allgenes/cS2G. \
        --frqfile-chr $FRQFILE \
        --w-ld-chr $WGTFILE \
        --overlap-annot --print-coefficients --print-delete-vals \
        --out sldsc_baselineS2G_outputs/$TRAIT.baselineS2G_cS2G.validation
done < SUMSTAT_PATH/traits63.txt
##########################################################################################


##########################################################################################
#2) sumannots file (used in next sections) were generated using the following code
while read LINE; do
    ARR=( $LINE ); TRAIT=${ARR[0]}
    Rscript sumannots_baselineS2G.r $BASELINELD,validation_${TRAIT}/validation_${TRAIT}.,allgenes/allgenes. $FRQFILE validation_${TRAIT}/validation_${TRAIT}.sumannots
    Rscript sumannots_baselineS2G.r $BASELINELD,validation_${TRAIT}/validation_${TRAIT}.,validation_${TRAIT}/cS2G.,allgenes/allgenes.,allgenes/cS2G. $FRQFILE validation_${TRAIT}/cS2G.sumannots
done < SUMSTAT_PATH/traits63.txt
Rscript sumannots_baselineS2G.r $BASELINELD,training/training.,allgenes/allgenes. $FRQFILE training/training.sumannots
Rscript sumannots_baselineS2G.r $BASELINELD,training/training.,training/cS2G.,allgenes/allgenes.,allgenes/cS2G.$FRQFILE training/cS2G.sumannots
#.sumannots files contain for all S2G-derived SNP annotations the number of SNPs in all the annotations of the baseline-S2G model
#The script sumannots_baselineS2G.r takes 3 arguments:
#argument 1: The list of annotation files (same argument provided to --ref-ld-chr ldsc option,the first one should be the baseline-LD model)
#argument 2: The frequency file of the reference panel (same argument provided to --frqfile-chr ldsc option)
#argument 3: The output filename
##########################################################################################


##########################################################################################
#3) h2 coverage, precision and recall for the 50 S2G-derived SNP annotation were computed using the following code
#3a) create validation.traits file containing the 63 trait names, corresponding ldsc baseline-S2G output, and sumannots file
S2G_PATH=/PATH/cS2G/baselineS2Gannots # this variable is the path to the baselineS2Gannots directory available at https://alkesgroup.broadinstitute.org/cS2G/
head -63 SUMSTAT_PATH/traits63.txt > tmp
while read LINE; do
    ARR=( $LINE ); TRAIT=${ARR[0]}
    echo "$TRAIT $S2G_PATH/sldsc_baselineS2G_outputs/$TRAIT.baselineS2G.validation $S2G_PATH/validation_$TRAIT/validation_$TRAIT.sumannots" >> validation.traits
done < tmp
rm tmp
#3b) for each of the 50 S2G-derived SNP annotation we compute h2 coverage, precision and recall 
for S2G in Exon Promoter Genebody Genebody100kb ClosestTSS ClosestTSS_1kb ClosestTSS_1kb_5kb ClosestTSS_5kb_10kb ClosestTSS_10kb_50kb ClosestTSS_50kb_100kb ClosestTSS_100kb_500kb ClosestTSS_500kb_1000kb Closest2ndTSS Closest3rdTSS Closest4thTSS Closest5thTSS Closest6thTSS Closest7thTSS Closest8thTSS Closest9thTSS Closest10thTSS Closest11thTSS Closest12thTSS Closest13thTSS Closest14thTSS Closest15thTSS Closest16thTSS Closest17thTSS Closest18thTSS Closest19thTSS Closest20thTSS ciseQTLs_GTeX ciseQTLs_GTeXblood ciseQTLs_eQTLGen finemappedciseQTLs_GTeX finemappedciseQTLs_GTeXblood finemappedciseQTLs_eQTLGen Roadmap Roadmapblood EpiMap EpiMapblood ABC ABCblood HiC_ClosestTSS HiC_distance PCHiC_Jung PCHiC_Javierreblood Ciceroblood GeneHancer OpenTargets; do
    Rscript meta_analysis_baselineS2G.r validation.traits $S2G baselineS2G.$S2G.txt
done
#The script meta_analysis_baselineS2G.r takes 3 arguments:
#argument 1: The file containing the traits to meta-analyze (column 1: trait ID; column 2: path to ldsc outputs; column 3: path to sumannot file)  
#argument 2: The name of the S2G strategy to analyze
#argument 3: The output filename
##########################################################################################


##########################################################################################
#4) h2 coverage, precision and recall for the cS2G-derived SNP annotation was computed using the following code
head -63 SUMSTAT_PATH/traits63.txt > tmp
while read LINE; do
    ARR=( $LINE ); TRAIT=${ARR[0]}
    echo "$TRAIT $S2G_PATH/sldsc_baselineS2G_outputs/$TRAIT.baselineS2G_cS2G.validation $S2G_PATH/validation_$TRAIT/cS2G.sumannots" >> validation.cS2G.traits
done < tmp
rm tmp
Rscript meta_analysis_baselineS2G.r validation.cS2G.traits cS2G baselineS2G_cS2G.cS2G.txt
##########################################################################################


##########################################################################################
#5) The weights of the S2G strategies of the cS2G strategy was computed using the following code
#5a) first we computed the per-snp heritability of each SNP by meta-analyzing regression coefficients of the baseline-S2G model
head -63 SUMSTAT_PATH/traits63.txt > tmp
while read LINE; do
    ARR=( $LINE ); TRAIT=${ARR[0]}
    echo "$TRAIT $S2G_PATH/sldsc_baselineS2G_outputs/$TRAIT.baselineS2G.training" >> training.traits
done < tmp
rm tmp
Rscript optimize_weights_persnph2.r training.traits $BASELINELD,$S2G_PATH/training/training.,$S2G_PATH/allgenes/allgenes. $FRQ baselineS2G.training.persnph2.txt
#5b) the weights were optimized using the following code; it uses 10 optimization algorithms with 10 different seeds
echo -e "Exon\nPromoter\nfinemappedciseQTLs_GTeX\nfinemappedciseQTLs_eQTLGen\nRoadmap\nEpiMap\nABC\nPCHiC_Jung\nPCHiC_Javierreblood\nCiceroblood" > listS2G.txt
for SEED in {1..10}; do
    Rscript optimize_weights.r baselineS2G.training.persnph2.txt listS2G.txt $S2G_PATH/training/training.,$S2G_PATH/allgenes/allgenes. $FRQ $SEED cS2G.seed$SEED
done
#The script optimize_weights_persnph2.r takes 4 arguments:
#argument 1: The file containing the traits to meta-analyze (column 1: trait ID; column 2: path to ldsc outputs)  
#argument 2: The list of annotation files (same argument provided to --ref-ld-chr ldsc option)
#argument 3: The frequency file of the reference panel (same argument provided to --frqfile-chr ldsc option)
#argument 4: The output filename
#The script optimize_weights.r takes 6 arguments:
#argument 1: The per-SNP heritability file from 5a) 
#argument 2: The file containing list of S2G startegies to consider; the first strategy needs to be Exon
#argument 3: The list of annotation files containining the S2G-derived SNP annotations (same argument provided to --ref-ld-chr ldsc option, but the baseline-LD model is not necessary)
#argument 4: The frequency file of the reference panel (same argument provided to --frqfile-chr ldsc option)
#argument 5: The seed number for the algorithm
#argument 6: The output filename header

##########################################################################################
#6) The omnigencity analyses were performed using the following code:
BETA2_PATH=/PATH/beta2 # this variable is the path to the beta2 estimates available at https://alkesgroup.broadinstitute.org/cS2G/pergeneh2_cS2G_UKBB/beta2
TRAIT=bp_SYSTOLICadjMEDz    #the trait to analyze; for example here bp_SYSTOLICadjMEDz
S2G_UKBB=/PATH/cS2G_UKBB/cS2G # this variable is the path to the cS2G links for UKBB SNPs available at https://alkesgroup.broadinstitute.org/cS2G/cS2G_UKBB
S2G_1000GEUR=/PATH/cS2G_1000GEUR/cS2G # this variable is the path to the cS2G links for 1000G EUR SNPs available at https://alkesgroup.broadinstitute.org/cS2G/cS2G_1000GEUR
GENEFILE=/PATH/list_genes_qc.txt.gz available at https://alkesgroup.broadinstitute.org/cS2G/critical_gene_sets
SUMSTATUKBB_PATH=/PATH/sumstats_UKBB122K # this variable is the path to the sumstats_UKBB122K directory available at https://alkesgroup.broadinstitute.org/cS2G/
BASELINELD=/PATH/baselineLD_v2.2/baselineLD. # this variable is the path to the baseline-LD v2.2 files available at https://alkesgroup.broadinstitute.org/LDSCORE/
FRQFILE=/PATH/1000G.EUR.QC. # this variable is the path to the reference files available at https://alkesgroup.broadinstitute.org/LDSCORE/
WGTFILE=/PATH/weights.hm3_noMHC. # this variable is the path to the weight files available at https://alkesgroup.broadinstitute.org/LDSCORE/
PRINTFILE=/PATH/hapmap3_snps.tgz # this variable is the path to the HapMap3 SNP file available at https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz

#6a) First we computed unadjusted per-gene heritabilities using cS2G 
perl get_h2_per_gene_init.pl $BETA2_PATH $TRAIT $S2G_UKBB $GENEFILE $TRAIT.cS2G.unadjusted.txt

#6b) Second we computed adjusted per-gene heritabilities using cS2G 
perl get_h2_per_gene.pl      $BETA2_PATH $TRAIT $S2G_UKBB $GENEFILE $TRAIT.cS2G.unadjusted.txt $TRAIT.cS2G.txt
rm $TRAIT.cS2G.init.txt

#6c) create annotations for S-LDSC
mkdir annotations
for SIZE in 100 200 500 1000 2000 5000 10000 19995; do
    #create annotations
    head -n $(($SIZE + 1)) $TRAIT.cS2G.txt | tail -n $SIZE | awk '{print $1, 1}' > annotations/$TRAIT.cS2G.$SIZE.txt
    echo $SIZE
    perl create_S2G_geneset.pl annotations/$TRAIT.cS2G.$SIZE $S2G_1000GEUR $FRQFILE
    #compute LD scores
    for CHR in {1..22}; do
        python ldsc.py \
          --l2 \
          --bfile $FRQFILE$CHR \
          --ld-wind-cm 1 \
          --annot annotations/$TRAIT.cS2G.$SIZE.$CHR.annot.gz \
          --thin-annot \
          --out annotations/$TRAIT.cS2G.$SIZE \
          --print-snps $PRINTFILE
    done
done

#6d) run S-LDSC
SUMSTATS=SUMSTATUKBB_PATH/UKB_122K.$TRAIT.sumstats.gz
python ldsc.py \
    --h2 $SUMSTATS \
    --ref-ld-chr $BASELINELD,annotations/$TRAIT.cS2G.100.,annotations/$TRAIT.cS2G.200.,annotations/$TRAIT.cS2G.500.,annotations/$TRAIT.cS2G.1000.,annotations/$TRAIT.cS2G.2000.,annotations/$TRAIT.cS2G.5000.,annotations/$TRAIT.cS2G.10000.,annotations/$TRAIT.cS2G.19995., \
    --frqfile-chr $FRQFILE \
    --w-ld-chr $WGTFILE \
    --overlap-annot --print-coefficients --print-delete-vals \
    --out sldsc_omnigenic_outputs/$TRAIT.cS2G.omnigenic

#6e) Finally compute Ge, Ge_c, Ge_lf in R
R
TRAIT="" #Fill with the trait
data = read.table(paste0(TRAIT,".cS2G.txt"),h=T)
Ge <- function (per_gene_beta2) { 3*19995/(mean(per_gene_beta2**2)/(mean(per_gene_beta2)**2)) }
Ge_all    = Ge(data$h2gene)
Ge_common = Ge(data$h2gene_common)
Ge_lowfrq = Ge(data$h2gene_lowfreq)
c(Ge_all,Ge_common,Ge_lowfrq)