use Regenie


ln -s /medpop/esp2/projects/MGB_Biobank/imputation/30K_GSA/result/bit8_merged .
ln -s /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged .
ln -s /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38 .

#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N regenie_assoc
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=4G
#$ -o regenie_step1.log
#$ -e regenie_step1.log

use .regenie-2.0.2

plink2 \
  --bfile /medpop/esp2/btruong/Projects/HDP/MGB/merged/GSA_30K.tag \
  --keep /medpop/esp2/btruong/Projects/HDP/MGB/pheno_preeclampsia.tsv \
  --maf 0.01 \
  --make-bed \
  --out /medpop/esp2/btruong/Projects/HDP/MGB/preeclampsia


plink2 \
  --bfile /medpop/esp2/btruong/Projects/HDP/MGB/merged/GSA_30K.tag \
  --keep /medpop/esp2/btruong/Projects/HDP/MGB/pheno_geshtn.tsv \
  --maf 0.01 \
  --make-bed \
  --out /medpop/esp2/btruong/Projects/HDP/MGB/geshtn


plink2 \
  --bfile /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38 \
  --maf 0.01 \
  --write-snplist \
  --out snps_pass


ls /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release


plink2 \
  --bfile /broad/ukbb/genotype/ukb_s_chr1_v2 \
  --maf 0.01 \
  --write-snplist \
  --out snps_pass


rm mergelist.txt
for ((chr=1; chr<=22; chr++)); do
  echo "/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam" >> mergelist.txt
done

qsub -l h_rt=5:00:00 -l h_vmem=32G -wd /medpop/esp2/btruong/Projects/logjobs -V tmp.sh


#######################################################


#/medpop/esp2/btruong/Projects/HDP/scripts/regenie_step1.sh

echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N rgn_${trait}
#$ -l h_rt=24:00:00
#$ -l s_rt=24:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=2G
#$ -j y

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2
use .regenie-2.0.2


mkdir -p ${wdir}
cd ${wdir}

regenie \
  --step 1 \
  --bed ${bed} \
  --phenoFile ${phenofile} \
  --covarFile ${covarfile} \
  --bsize 1000 \
  --extract ${snpextract} \
  --bt \
  --covarColList ${covarCol} \
  --catCovarList ${catcovar} \
  --out ${out}
" > tmp.sh
cat tmp.sh


qsub tmp.sh

trait=composite
trait=preeclampsia
trait=composite

wdir=/medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}
phenofile=/medpop/esp2/btruong/Projects/HDP/MGB/pheno_${trait}.tsv
covarfile=/medpop/esp2/btruong/Projects/HDP/MGB/basic_covariates.txt
covarCol=age,PC{1:10}
catcovar=batch,sex
out=${trait}_regenie_step1
snpextract=/medpop/esp2/btruong/Projects/HDP/MGB/snps_pass.snplist
bed=/medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38
trait=${trait} wdir=${wdir} bed=${bed} phenofile=${phenofile} covarfile=${covarfile} snpextract=${snpextract} covarCol=${covarCol} catcovar=${catcovar} $WDIR/scripts/regenie_step1.sh







###############################
for chr in {1..22}; do 
  ln -s /broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed ukb_chr${chr}_v2.bed
  ln -s /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim ukb_chr${chr}_v2.bim
  ln -s /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam ukb_chr${chr}_v2.fam
done





trait=composite
trait=preec
trait=geshtn

wdir=/broad/hptmp/btruong/HDP/regenie/${trait}
phenofile=/medpop/esp2/btruong/Projects/HDP/data/preec_UKB_2.txt
covarfile=/medpop/esp2/btruong/Projects/HDP/data/preec_UKB_2.txt
covarCol=age,PC{1:10}
catcovar=sex
out=${trait}_regenie_step1
snpextract=/medpop/esp2/btruong/Projects/HDP/MGB/snps_pass.snplist
bed=/medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38

mkdir -p ${wdir}
trait=${trait} wdir=${wdir} bed=${bed} phenofile=${phenofile} covarfile=${covarfile} snpextract=${snpextract} covarCol=${covarCol} catcovar=${catcovar} $WDIR/scripts/regenie_step1.sh





#######################################################

awk '{print $1,$2,$2,$2}' /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38.fam > GSA_30K_hg38_sameFID.IID.txt
plink --bfile /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38 --update-ids GSA_30K_hg38_sameFID.IID.txt --make-bed --out GSA_30K_hg38_sameFID.IID

#######################################################

awk '{if(FNR!=1) $1=$2; print $0}' /medpop/esp2/btruong/Projects/HDP/MGB/pheno_${trait}.tsv > /medpop/esp2/btruong/Projects/HDP/MGB/pheno_${trait}_1.tsv
awk '{if(FNR!=1) $1=$2; print $0}' /medpop/esp2/btruong/Projects/HDP/MGB/basic_covariates.txt > /medpop/esp2/btruong/Projects/HDP/MGB/basic_covariates_1.txt

#######################################################

echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N rS2_${trait}
#$ -l h_rt=10:00:00
#$ -l s_rt=10:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=4G
#$ -j y
#$ -t 1-22


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2
use .regenie-2.0.2

chr=\$SGE_TASK_ID
echo \${chr}

cd /medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}
pwd

regenie \
  --step 2 \
  --bgen /medpop/esp2/projects/MGB_Biobank/imputation/30K_GSA/result/bit8_merged/30K_GSA.chr\${chr}.bit8.dose.bgen \
  --phenoFile /medpop/esp2/btruong/Projects/HDP/MGB/pheno_${trait}.tsv \
  --covarFile /medpop/esp2/btruong/Projects/HDP/MGB/basic_covariates.txt \
  --sample /medpop/esp2/btruong/Projects/HDP/MGB/MGB_FIDeq0.sample \
  --bsize 200 \
  --bt \
  --firth \
  --approx \
  --pThresh 0.01 \
  --maxiter-null 10000 \
  --maxstep-null 2 \
  --pred ${trait}_regenie_step1_pred.list \
  --catCovarList batch,sex \
  --out /medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}_regenie_step2.firth.chr\${chr}

" > tmp.sh
cat tmp.sh
qsub tmp.sh



#######################


trait=composite
trait=preeclampsia
trait=composite

trait=${trait} /medpop/esp2/btruong/Projects/HDP/scripts/regenie_step2.sh


#######################

gunzip -c /medpop/esp2/btruong/Tools/All_20180418.vcf.gz | vcf2bed - > hg38.dbSNP151.bed1

awk -vOFS="\t" '{if(FNR!=1) print $1, ($2 - 1), $2; }' composite_regenie_step2.firth.chr22_Y1.regenie | sort-bed - > composite_regenie_step2.firth.chr22_Y1.regenie.bed

use BEDTools
bedtools intersect -a composite_regenie_step2.firth.chr22_Y1.regenie_new -b /medpop/esp2/btruong/Tools/All_20180418.vcf.gz -wa -wb > tmp



################################

#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N regenie_preeclampsia
#$ -l h_rt=2:00:00
#$ -l s_rt=2:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=4G
#$ -j y
#$ -t 1-22


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

trait=preeclampsia
chr=$SGE_TASK_ID
echo $chr

use BEDTools

cd /medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}

awk -vOFS="\t" '{if(FNR!=1) print $1, ($2 - 1), $2; }' ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie | sort-bed - > ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie.bed


bedmap --echo --echo-map-id --delim '\t' ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie.bed /medpop/esp2/btruong/Tools/hg38.dbSNP151.common.bed > ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie_withRSID.bed

  

################################

#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N regenie_composite
#$ -l h_rt=2:00:00
#$ -l s_rt=2:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=4G
#$ -j y
#$ -t 1-22


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

trait=composite
chr=$SGE_TASK_ID
echo $chr

use BEDTools

cd /medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}

Rscript /medpop/esp2/btruong/Projects/HDP/scripts/posTOrsid_sumstat_regenie.R \
  --sumstat ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie \
  --annot ${trait}_regenie_step2.firth.chr${chr}_Y1.regenie_withRSID.bed \
  --out ${trait}_regenie_step2.firth.chr${chr}_Y1_withRSID.regenie



######################################

trait=preeclampsia
# trait=geshtn
# trait=composite

cd /medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}

cat ${trait}_regenie_step2.firth.chr*_Y1_withRSID.regenie > ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie

awk '{if (FNR==1 || $1!="CHROM") print $0}' ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie > ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie_1
sort -k1n -k2n ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie_1 > ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie

rm ${trait}_regenie_step2.firth.allchr_Y1_withRSID.regenie_1




############################################
#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -l h_rt=1:00:00
#$ -l s_rt=1:00:00
#$ -l h_vmem=32G


awk '(FNR!=1){
    split($1,a,"|")
    split(a[1],b,":")
    print "chr"b[1],b[2],b[2]+1,b[1]":"b[2]
}' OFS="\t" /medpop/esp2/SarahUrbut/for_sarah/ALL_lipid_CAD_associations_merged.txt > ALL_lipid_CAD_associations_merged_rsid4lift.bed

head ALL_lipid_CAD_associations_merged_rsid4lift.bed


liftOver ALL_lipid_CAD_associations_merged_rsid4lift.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz ALL_lipid_lifted.bed ALL_lipid_unlifted.bed 


###########################################

awk '(FNR!=1){
  $1=gsub ("^0*", "", $1)
    print "chr"$1,$2,$2+1,$1":"$2
}' /medpop/esp2/honigberg/HDP_GWAS/ukb_HDP_full_results.tsv > ukb_HDP_full_results.tsv.bed

head ukb_HDP_full_results.tsv.bed

liftOver ukb_HDP_full_results.tsv.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz ukb_HDP_full_results.tsv_lifted.bed ukb_HDP_full_results.tsv_unlifted.bed




awk '(FNR!=1){
    print "chr"$16,$17,$17+1,$16":"$17
}' /medpop/esp2/honigberg/HDP_GWAS/InterPregGen/summary_stats_maternal_preeclampsia_with_rsids.txt > summary_stats_maternal_preeclampsia_with_rsids.txt.bed

head summary_stats_maternal_preeclampsia_with_rsids.txt.bed

liftOver summary_stats_maternal_preeclampsia_with_rsids.txt.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz summary_stats_maternal_preeclampsia_with_rsids.txt_lifted.bed summary_stats_maternal_preeclampsia_with_rsids.txt_unlifted.bed





awk '(FNR!=1){
    print "chr"$1,$2,$2+1,$1":"$2
}' /medpop/esp2/honigberg/HDP_GWAS/Estonia/O13_O14_O15_EstBB_200421_TL > O13_O14_O15_EstBB_200421_TL.bed

head O13_O14_O15_EstBB_200421_TL.bed

liftOver O13_O14_O15_EstBB_200421_TL.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz O13_O14_O15_EstBB_200421_TL_lifted.bed O13_O14_O15_EstBB_200421_TL_unlifted.bed




zcat /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.txt.gz | awk '(FNR!=1){
    sprintf("chr%d\t%d\t%d\t%d:%d",$2,$3,$3+1,$2,$3)
    # a=sprintf("chr%d",$2)
    # print a
    # print "chr"$2,$3,$3+1,$2":"$3
}' > GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed

head GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed

liftOver GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_lifted.bed GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_unlifted.bed



zcat /medpop/esp2/honigberg/HDP_GWAS/InterPregGen/European_only/mat_all_chrALL_STERR_EU.1tbl.gz | awk '(FNR!=1){
    sprintf("chr%d\t%d\t%d\t%d:%d",$2,$3,$3+1,$2,$3)
    # a=sprintf("chr%d",$2)
    # print a
    # print "chr"$2,$3,$3+1,$2":"$3
}' > GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed

head GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed

liftOver GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_lifted.bed GWASsummary_PreEclampsia_Japanese_SakaueKanai2020.auto_unlifted.bed





zcat /medpop/esp2/honigberg/HDP_GWAS/Estonia/O13_EstBB_080421_TL.gz | awk '(FNR!=1){
    print "chr"$1,$2,$2+1,$1":"$2
}' > O13_EstBB_080421_TL.gz.bed

head O13_EstBB_080421_TL.gz.bed

liftOver O13_EstBB_080421_TL.gz.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz O13_EstBB_080421_TL.gz_lifted.bed O13_EstBB_080421_TL.gz_unlifted.bed



use Bcftools

wget -O hg37.db151.vcf.gz ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget -O hg37.db151.vcf.gz.tbi ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi


wget -O hg38.db151.vcf.gz ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
wget -O hg38.db151.vcf.gz.tbi ftp://ftp.ncbi.nih.gov/snp/organisms//human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.tbi



bcftools query -f '%CHROM\t%POS\t%ID\n' hg38.db151.vcf.gz > hg37_common.txt

awk '{print $0,$1":"$2}' OFS="\t" /medpop/esp2/btruong/Tools/hg37_common.txt > /medpop/esp2/btruong/Tools/hg37_common_chrpos.txt

hg=hg38
echo "bcftools query -f '%CHROM\t%POS\t%ID\n' ${hg}.db151.vcf.gz > ${hg}_allsnp.txt" > tmp.sh
qsub -l h_vmem=16G -l h_rt=3:00:00 -wd /medpop/esp2/btruong/Projects/logjobs -V tmp.sh

hg=hg37
hg=hg38
echo "awk '{print \$0,\$1\":\"\$2}' OFS=\"\t\" /medpop/esp2/btruong/Tools/${hg}_allsnp.txt > /medpop/esp2/btruong/Tools/${hg}_allsnp_chrpos.txt" > tmp.sh
qsub -l h_vmem=16G -l h_rt=3:00:00 -wd /medpop/esp2/btruong/Projects/logjobs -V tmp.sh




#########################################


echo "#!/bin/bash -l
#$ -N ${trait}_${ancestry}
#$ -pe smp 4 -R y -binding linear:4
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -R y
#$ -l h_rt=3:30:00
#$ -l s_rt=3:30:00
#$ -l h_vmem=8G

cd /medpop/esp2/btruong/Projects/HDP/data

Rscript /medpop/esp2/btruong/Projects/HDP/scripts/write_for_metal.R --trait=${trait} --ancestry=${ancestry} --out=${out}

metal ${out}running.sh
" > tmp.sh



cat ${out}running.sh
cat tmp.sh


# qsub tmp.sh


##################################

trait=composite
trait=preec
# trait=geshtn

# ancestry=European
# ancestry=European.Hispanic
# ancestry=Asian

trait_list="preec composite geshtn"
trait_list="preec"
# trait_list="composite"
# trait_list="geshtn"
# ancestry_list="European European.Hispanic European.Hispanic.Asian"
ancestry_list="European.Hispanic.Asian.African"
ancestry_list="European.African"
# ancestry_list="European"


cd /medpop/esp2/btruong/Projects/HDP/jobs

for trait in $trait_list; do
  echo ${trait}
  for ancestry in $ancestry_list; do
    echo "---${ancestry}"
    out=metal_ukb.pmbb.hunt_${trait}_${ancestry}_
    trait=${trait} ancestry=${ancestry} out=${out} /medpop/esp2/btruong/Projects/HDP/scripts/make_metal.sh
    cat tmp.sh
  done
done




filename=/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/output2/preeclampsia_chr1-23.results_hg38.txt
filename=/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/preec_allchr_maf0.001_hg38.regenie
filename=/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/geshtn_allchr_maf0.001_hg38.regenie
awk '{$1=$1;print $0}' OFS="\t" ${filename} > ${filename}1
mv ${filename}1 ${filename}



for jobs in {32159152..32159477}; do
  # qalter ${jobs} -l h_vmem=2G
  qalter ${jobs} -l h_rt=5:00:00:00
done








out=metal_${trait}_${ancestry}_

trait=${trait} ancestry=${ancestry} out=${out} /medpop/esp2/btruong/Projects/HDP/scripts/make_metal.sh

snp="16:53763996"
awk -v snp=${snp} '{
  if (FNR==1) {
    cc=0
    for(i=1;i<=NF;i++) { cc+=1; if($i=="chrpos") col=cc; }
    print col;
  } else {
    print $2
    # if ($col==snp) print $0
  }
}' 2022_03_01_GNH_binary_GWAS_firth_GNH0248_Preeclampsia_Eclampsia_GestationalHypertension.regenie_newRSID.txt | head


grep -o -a -m 1 -h -r "16:53763996" 2022_03_01_GNH_binary_GWAS_firth_GNH0248_Preeclampsia_Eclampsia_GestationalHypertension.regenie_newRSID.txt | head -1

grep -rIl "16:53763996" 2022_03_01_GNH_binary_GWAS_firth_GNH0248_Preeclampsia_Eclampsia_GestationalHypertension.regenie_newRSID.txt

cat ${out}running.sh

/medpop/esp2/honigberg/HDP_GWAS/finngen/summary_stats_finngen_R6_O15_PRE_OR_ECLAMPSIA.gz
41836590

/medpop/esp2/btruong/Projects/HDP/data/O14_O15_EstBB_080421_TL_lifted.txt
42403619

rs370119665

rs113211757



10      42403619        rs370119665     T       A       3214.98974609375        0.0501981340348721      0.99923 32023   0.0736346751490625      0.0887032291474243      9.35844117668625        0.406468706078803       0.406468706078803       1       127.092856154271     128.946401342343 0.0541285714433928      0.0500184501885114      1400    30623   10:42403619     10:41836588     41836588




#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N qcsumstat
#$ -l h_rt=1:00:00
#$ -l s_rt=1:00:00
#$ -l h_vmem=8G
#$ -j y
#$ -t 9-17


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

row=$SGE_TASK_ID

Rscript /medpop/esp2/btruong/Projects/HDP/scripts/preprocess_sumstat.R --row=${row}




java -jar /medpop/esp2/btruong/Tools/Metasoft.jar 

python plink2metasoft.py outputfile /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/O13_EstBB_080421_TL_lifted.txt /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/EA_GH.glm.logistic_newRSID


while read -d, i;
do
  echo ${i[1]}
  # echo ${p[1]}
done < /medpop/esp2/btruong/Projects/HDP/data/sumstat_list.csv


while IFS=, read p; do
  p=($p)
  # echo ${p[1]}
done < /medpop/esp2/btruong/Projects/HDP/data/sumstat_list.csv


ancestry=European
snp="16:53763996"
# snp="53763996"
trait="preec"
awk -F ',' -v trait=${trait} -v ancestry=${ancestry} '{
  # snp="53763996"
  # cmd="grep "snp" "$16
  if ($2==trait & $3==ancestry) {
    print $16
  #   system(cmd)
  }
}' /medpop/esp2/btruong/Projects/HDP/data/sumstat_list.csv


rs146780608:9816170:C:T

1:9756112



while read files; do
    echo $files
    grep $snp $files
    
    
done




#######################################

echo "#!/bin/bash -l
#$ -pe smp 2 -R y -binding linear:2
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N ${jobtitle}
#$ -l h_rt=1:00:00
#$ -l s_rt=1:00:00
#$ -l h_vmem=8G
#$ -j y


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2


Rscript /medpop/esp2/btruong/Projects/HDP/scripts/manhattanPlot.R --file=${file} --title=$title --out=$out
" > tmp.sh
cat tmp.sh
qsub tmp.sh



for trait in $trait_list; do
  for ancestry in $ancestry_list; do
    # echo $ancestry
    file=/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_${trait}_${ancestry}_1.txt
    title=${trait}":"${ancestry}
    out=/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_${trait}_${ancestry}_1.pdf
    jobtitle=${trait}.${ancestry}
    file=$file out=$out title=$title jobtitle=$jobtitle /medpop/esp2/btruong/Projects/HDP/scripts/plotmanhat.sh
  done
done


trait=geshtn
trait=preec
trait=composite
ancestry=European.Hispanic.Asian

use .perl-5.28.0

vep \
  -i sigsnp_${trait}_${ancestry}.txt \
  -o res_sigsnp_${trait}_${ancestry}.txt \
  --max_af --everything \
  --force_overwrite --cache \
  --pick --tab --plugin Phenotypes \
  --fields Uploaded_variation,SYMBOL,CANONICAL,Consequence,BIOTYPE,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,MAX_AF,PHENOTYPES \
  --dir_cache $veppath/cache \
  --dir_plugin /medpop/esp2/btruong/Tools/VEPplugins 


filter_vep -i res_sigsnp_${trait}_${ancestry}.txt --force_overwrite -filter "MAX_AF > 0.01 or not MAX_AF" -o res_sigsnp_${trait}_${ancestry}_filtered.txt






#####################



use .perl-5.28.0

veppath=/medpop/esp2/projects/software/vep105/ensembl-vep-release-105

export PERL5LIB=$veppath/cpanm/lib/perl5:$veppath/loftee

#tabix should match with vep version
PATH=$veppath/htslib:$PATH

#loftee resource
loftee_resource=/medpop/esp2/projects/software/loftee_resource
gerp=$loftee_resource/gerp_conservation_scores.homo_sapiens.GRCh38.bw
anc=$loftee_resource/human_ancestor.fa.gz




#/medpop/esp2/btruong/Projects/HDP/scripts/regenie_step1.sh

echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N rgn_${trait}
#$ -l h_rt=24:00:00
#$ -l s_rt=24:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=2G
#$ -j y

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2
use .regenie-2.0.2


mkdir -p ${wdir}
cd ${wdir}

regenie \
  --step 1 \
  --bed ${bed} \
  --phenoFile ${phenofile} \
  --phenoCol ${phenoCol} \
  --covarFile ${covarfile} \
  --bsize 1000 \
  --extract ${snpextract} \
  --bt \
  --maxCatLevels 50 \
  --covarColList ${covarCol} \
  --catCovarList ${catcovar} \
  --out ${out}
" > tmp.sh
cat tmp.sh


qsub tmp.sh

trait=composite
trait=preeclampsia
trait=composite

wdir=/medpop/esp2/btruong/Projects/HDP/MGB/regenie/${trait}
phenofile=/medpop/esp2/btruong/Projects/HDP/MGB/pheno_${trait}.tsv
covarfile=/medpop/esp2/btruong/Projects/HDP/MGB/basic_covariates.txt
covarCol=age,PC{1:10}
catcovar=batch,sex
out=${trait}_regenie_step1
snpextract=/medpop/esp2/btruong/Projects/HDP/MGB/snps_pass.snplist
bed=/medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38

trait=${trait} wdir=${wdir} bed=${bed} phenofile=${phenofile} phenoCol=${phenoCol} covarfile=${covarfile} snpextract=${snpextract} covarCol=${covarCol} catcovar=${catcovar} $WDIR/scripts/regenie_step1.sh

plink2 \
  --bfile /medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38 \
  --maf 0.01 \
  --write-snplist \
  --out snps_pass

rm mergelist.txt
for ((chr=1; chr<=22; chr++)); do
  echo "ukb_hm3_chr${chr}" >> mergelist.txt
done




plink --merge-list mergelist.txt --memory 16000 --make-bed --out ukb_hm3




qsub -pe smp 4 -R y -binding linear:4 -l h_rt=24:00:00 -l h_vmem=8G -V -wd /medpop/esp2/btruong/Projects/logjobs tmp.sh



phenofile=/medpop/esp2/btruong/Projects/HDP/data/geshtn_UKB_1.txt
phenoCol=gh1
covarfile=/medpop/esp2/btruong/Projects/HDP/data/geshtn_UKB_1.txt




/broad/hptmp/btruong/ukb_hm3



echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N PRS_meta
#$ -l h_rt=10:00:00
#$ -l s_rt=10:00:00
#$ -l h_vmem=16G

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

chr=\$SGE_TASK_ID
# cd /medpop/esp2/btruong/Projects/HDP/validation_ukb
cd /broad/hptmp/btruong/validation_ukb


# /medpop/esp2/btruong/Tools/qctool/qctool -g /broad/ukbb/imputed_v3/ukb_imp_chr\${chr}_v3.bgen -og ukb_${trait}_${ancestry}_chr\${chr} -s /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr15_v3_s487395.sample -ofiletype binary_ped -incl-rsids /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_${ancestry}_chr\${chr}_snplist.txt

plink2 --bfile ukb_${trait}_${ancestry}_chr\${chr} --indep-pairwise 500 50 0.8 --rm-dup force-first --score /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_${ancestry}_chr\${chr}_w.txt 1 2 3 cols=+scoresums --out ukb_${trait}_${ancestry}_chr\${chr}


plink2 --bfile ukb_${trait}_${ancestry}_chr\${chr} --indep-pairwise 500 50 0.8 --rm-dup force-first --score /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_${ancestry}_chr\${chr}_w_sig.txt 1 2 3 cols=+scoresums --out ukb_${trait}_${ancestry}_chr\${chr}_sig
" > tmp.sh
cat tmp.sh

qsub -t 1-1 tmp.sh



trait="geshtn"
ancestry="European.Hispanic.Asian"

trait="preec"
ancestry="European.Hispanic.Asian"



plink2 --bfile ukb_${trait}_${ancestry}_chr${chr} --indep-pairwise 500 50 0.1 --rm-dup force-first --score /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_${ancestry}_chr${chr}_w.txt 1 2 3 cols=+scoresums --out ukb_${trait}_${ancestry}_chr${chr}_w500.s50.r20.8


trait=${trait} ancestry=${ancestry} /medpop/esp2/btruong/Projects/HDP/scripts/job_PRS_from_meta.sh
qsub -t 1-22 tmp.sh

chr=22

plink2 --bfile ukb_${trait}_${ancestry}_chr${chr} --rm-dup force-first --out ukb_${trait}_${ancestry}_chr${chr}_nodup


plink2 --bfile ukb_${trait}_${ancestry}_chr${chr} --rm-dup force-first --score /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_${ancestry}_chr${chr}_w.txt 1 2 3 cols=+scoresums --out ukb_${trait}_${ancestry}_chr${chr}








### QC ancestry MGBB
/medpop/esp2/skoyama/ukbb_sampleqc/src/qc.R





vdb-decrypt --ngc /medpop/esp2/btruong/Tools/prj_6213_D34507.ngc

vdb-decrypt --ngc /medpop/esp2/btruong/Tools/prj_6213_D34507.ngc phg001583.v1.TOPMed_WGS_HCHS_SOL_v2_frz9.genotype-qc.MULTI.tar.ncbi_enc
vdb-decrypt --ngc /medpop/esp2/btruong/Tools/prj_6213_D34507.ngc phg001583.v1.TOPMed_WGS_HCHS_SOL_v2_frz9.genotype-calls-vcf.c1.HMB-NPU.tar.ncbi_enc


vdb-decrypt --ngc /medpop/esp2/btruong/Tools/prj_6213_D34507.ngc Release_Notes.phs001395.TOPMed_WGS_HCHS_SOL.v2.p1.MULTI.pdf.ncbi_enc.aspx




awk '(FNR!=1){print "chr"$2,$3,$3+1,$2":"$3}' /medpop/esp2/honigberg/GWAS/MVP/sub20190109/sbp_MVP-only_transethnic_12052018.txt > sbp_MVP-only_transethnic_12052018_forlift.bed

liftOver sbp_MVP-only_transethnic_12052018_forlift.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz sbp_MVP-only_transethnic_12052018_forlift_lifted.bed sbp_MVP-only_transethnic_12052018_forlift_unlifted.bed 


awk '(NR==FNR){newpos[$4]=$2;next}{
  oldid=$2":"$3
  if (FNR==1) {
    print $0"\tnewpos"
    } else {
      print $0"\t"newpos[oldid]
    }
}' sbp_MVP-only_transethnic_12052018_forlift_lifted.bed /medpop/esp2/honigberg/GWAS/MVP/sub20190109/sbp_MVP-only_transethnic_12052018.txt > sbp_MVP-only_transethnic_12052018_newposHg38.txt



zcat /medpop/esp2/btruong/Projects/HDP/data/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | awk '(FNR!=1){print "chr"$1,$2,$2+1,$1":"$2}' > Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_forlift.bed

liftOver Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_forlift.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_forlift_lifted.bed Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_forlift_unlifted.bed 


zcat /medpop/esp2/btruong/Projects/HDP/data/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | \
  awk '(NR==FNR){newpos[$4]=$2;next}{
    oldid=$1":"$2
    if (FNR==1) {
      print $0"\tnewpos"
      } else {
        print $0"\t"newpos[oldid]
      }
  }' Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_forlift_lifted.bed - > meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_newposHg38.txt




python /PRScs/PRScs.py  \
  --ref_dir=/data/alh-admmt/Projects/Collider_bias/PRSCS/data/ldblk_ukbb_eur \
  --bim_prefix=/ukbb_hrc_chrautosom \
  --sst_file=/data/alh-admmt/Projects/Collider_bias/PRSCS/BMI/BMI_PRSCSsumstat.txt \
  --n_gwas=700000 \
  --out_dir=/data/alh-admmt/Projects/Collider_bias/PRSCS/BMI/BMI \
  --chrom=${chr} \
  --seed=${chr}


echo "#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -l h_vmem=20G
#$ -l h_rt=20:00:00
#$ -t 1-22
#$ -N prscs


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3

base_pop=eur
chr=\${SGE_TASK_ID}

cd /broad/hptmp/btruong/validation_ukb/PRSCS/posterior/


# testplink=/broad/hptmp/btruong/validation_ukb/ukb_${trait}_European.Hispanic.Asian_chr\${chr}_newids_nodup


PRScs=/medpop/esp2/yruan/tools/PRScs/PRScs.py
ref_dir=/medpop/esp2/yruan/raw.data/ld.ref.prs-csx.ukbb/ldblk_ukbb_\${base_pop}


Rscript \$WDIR/codes/munge_sumstat_PRSCS.R --sumstat ${sumstat} --SNP ${SNP} --A1 ${A1} --A2 ${A2} --BETA ${BETA} --P ${P} --out ${out}_chr\${chr}_w.txt

python \$PRScs \
 --ref_dir=\${ref_dir} \
 --bim_prefix=${bim_prefix} \
 --sst_file=${out}_chr\${chr}_w.txt \
 --n_gwas=${n_trn} \
 --chrom=\$chr \
 --phi=$phi \
 --out_dir=${out}_chr\${chr}_posterior

plink --bfile ${testplink} --score ${out}_chr\${chr}_posterior_pst_eff_a1_b0.5_phi${phi}_chr\${chr}.txt 2 4 6 sum --out ${out}_chr\${chr}

" > tmp.sh
cat tmp.sh

# Rscript \$WDIR/codes/munge_sumstat_PRSCS.R --sumstat ${sumstat} --SNP ${rsid} --A1 ${A1} --A2 ${A2} --BETA ${beta} --P ${P} --out ${trait}_European.Hispanic.Asian_chr\${chr}_w_prscs.txt






trait=preec
n_trn=400000

trait=geshtn
n_trn=200000



philist=(1e-06 1e-04 1e-02 1e+00)



for ((chr=1; chr<=22; chr++)); do
  # mv posterior/PRS-CS_preec_phi.1_pst_eff_a1_b0.5_phi1e+00_chr${chr}.txt posterior/PRS-CS_preec_phi.1e+00_pst_eff_a1_b0.5_phi1e+00_chr${chr}.txt

done





zcat /medpop/esp2/smjcho/outgoing/SBP_MVP_transethnic.results.gz | \
  awk '(FNR==NR){snp[$2]=1;next}{
    if (FNR==1 || snp[$1]) {print $0; next}
  }' /broad/hptmp/btruong/ukb_hm3/ukb_hm3.bim - > SBP_MVP_transethnic.results.txt



A1=EA
A2=Allele2
SNP=SNP_ID
BETA=EffectEstimate
P=Pvalue
n_trn=250000
sumstat=/medpop/esp2/btruong/Projects/HDP/data/SBP_MVP_transethnic.results.txt
sumstat=/medpop/esp2/btruong/Projects/HDP/data/DBP_MVP_transethnic.results.txt
# sumstat=/medpop/esp2/smjcho/outgoing/DBP_MVP_transethnic.results.gz


bim_prefix=/broad/hptmp/btruong/ukb_hm3/ukb_hm3

testplink=/broad/hptmp/btruong/ukb_hm3/ukb_hm3_chr\${chr}



philist=(1e-06 1e-04 1e-02 1e+00)

for ((phi_i=3; phi_i<=3; phi_i++)); do
  # phi_i=0
  phi=${philist[phi_i]}
  # Rscript /medpop/esp2/btruong/scripts/concatenatePRS.R --file PRSCS_${trait}_phi.${phi}_chr{chr}.profile --out PRSCS_${trait}_phi.${phi}_allchr.profile
  
  # sumstat=/medpop/esp2/smjcho/outgoing/DBP_MVP_transethnic.results.gz
  # out=PRSCS_SBP.MVP_phi.${phi}
  out=PRSCS_DBP.MVP_phi.${phi}
  
  sumstat=${sumstat} testplink=${testplink} bim_prefix=${bim_prefix} trait=${trait} A1=${A1} A2=${A2} SNP=${SNP} BETA=${BETA} P=${P} n_trn=${n_trn} phi=${phi} out=${out} $WDIR/Projects/HDP/scripts/PRSCSx.sh
  qsub tmp.sh
  
  # cd /broad/hptmp/btruong/validation_ukb/PRSCS/posterior/
  # Rscript /medpop/esp2/btruong/scripts/concatenatePRS.R --file ${out}_chr{chr}.profile --out ${out}_allchr.profile
  
done





cp PRSCS_DBP.MVP*allchr.profile $WDIR/Projects/HDP/results/


Rscript $WDIR/codes/munge_sumstat_PRSCS.R --sumstat /medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_European.Hispanic.Asian_chr${chr}_w.txt --SNP rsid --A1 Allele1 --A2 Allele2 --BETA Effect --P "P-value" --out ${trait}_European.Hispanic.Asian_chr${chr}_w_prscs.txt


/medpop/esp2/btruong/Projects/HDP/validation_ukb/${trait}_European.Hispanic.Asian_chr${chr}_w.txt


################################################ 

PRS-CSx


echo "#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -l h_vmem=20G
#$ -l h_rt=20:00:00
#$ -t 1-22
#$ -N prscsx


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3

chr=\${SGE_TASK_ID}

cd /broad/hptmp/btruong/validation_ukb/PRSCS/posterior/

PRScsx=/medpop/esp2/yruan/tools/PRScsx/PRScsx.py
ref_dir=/medpop/esp2/yruan/raw.data/ld.ref.prs-csx.ukbb
bim_prefix=/broad/hptmp/btruong/ukb_hm3/ukb_hm3

# Rscript \$WDIR/codes/munge_sumstat_PRSCS.R --sumstat ${sumstat1} --SNP ${SNP} --A1 ${A1} --A2 ${A2} --BETA ${BETA} --P ${P} --out ${sumstat1}_preproc1
# Rscript \$WDIR/codes/munge_sumstat_PRSCS.R --sumstat ${sumstat2} --SNP ${SNP} --A1 ${A1} --A2 ${A2} --BETA ${BETA} --P ${P} --out ${sumstat2}_preproc1

python \$PRScsx \
 --ref_dir=\${ref_dir} \
 --bim_prefix=\${bim_prefix} \
 --sst_file=${sumstat1}_preproc1,${sumstat2}_preproc1\
 --pop=${base_pop1},${base_pop2} \
 --n_gwas=${n_trn1},${n_trn2} \
 --chrom=\${chr} \
 --phi=${phi} \
 --out_dir=/broad/hptmp/btruong/validation_ukb/PRSCS/posterior/ \
 --out_name=${out}_chr\${chr}_${base_pop1}.${base_pop2}



plink --bfile ${testplink} --score ${out}_chr\${chr}_${base_pop1}.${base_pop2}_${base_pop1}_pst_eff_a1_b0.5_phi${phi}_chr\${chr}.txt 2 4 6 sum --out ${out}_${base_pop1}_chr\${chr}

plink --bfile ${testplink} --score ${out}_chr\${chr}_${base_pop1}.${base_pop2}_${base_pop2}_pst_eff_a1_b0.5_phi${phi}_chr\${chr}.txt 2 4 6 sum --out ${out}_${base_pop2}_chr\${chr}



" > tmp.sh
cat tmp.sh





SNP=newRSID
A1=Allele1
A2=Allele2
BETA=Effect
P="P-value"

sumstat1=/medpop/esp2/btruong/Projects/HDP/data/metal_preec_European_1_withRSID.txt_preproc1
sumstat2=/medpop/esp2/btruong/Projects/HDP/data/metal_preec_Asian_1_withRSID.txt_preproc1

sumstat=/medpop/esp2/btruong/Projects/HDP/results/MEGA_META/metal_preec_European_1_withRSID.txt

SNP=rsid
sumstat=$WDIR/Projects/HDP/results/MEGA_META/metal_preec_European.Hispanic.Asian_1_withRSID.txt

Rscript $WDIR/codes/munge_sumstat_PRSCS.R --sumstat ${sumstat} --SNP ${SNP} --A1 ${A1} --A2 ${A2} --BETA ${BETA} --P ${P} --out metal_preec_European.Hispanic.Asian_1_withRSID.txt_preproc1




# must have file with correpsonding *preproc1
sumstat1=/medpop/esp2/btruong/Projects/HDP/data/metal_preec_European_1_withRSID.txt
sumstat1=metal_preec_European.Hispanic.Asian_1_withRSID.txt
base_pop1=eur
n_trn1=400000
# n_trn1=350000

sumstat2=/medpop/esp2/btruong/Projects/HDP/data/metal_preec_Asian_1_withRSID.txt
sumstat2=metal_preec_European.Hispanic.Asian_1_withRSID.txt
base_pop2=eas
# n_trn2=100000
n_trn2=400000

# sumstat=/medpop/esp2/smjcho/outgoing/DBP_MVP_transethnic.results.gz


bim_prefix=/broad/hptmp/btruong/ukb_hm3/ukb_hm3

testplink=/broad/hptmp/btruong/ukb_hm3/ukb_hm3_chr\${chr}


base_pop1=eur
base_pop2=eas

philist=(1e-06 1e-04 1e-02 1e+00)

for ((phi_i=0; phi_i<=3; phi_i++)); do
  # phi_i=0
  phi=${philist[phi_i]}
  # Rscript /medpop/esp2/btruong/scripts/concatenatePRS.R --file PRSCS_${trait}_phi.${phi}_chr{chr}.profile --out PRSCS_${trait}_phi.${phi}_allchr.profile
  
  # sumstat=/medpop/esp2/smjcho/outgoing/DBP_MVP_transethnic.results.gz
  # out=PRSCS_SBP.MVP_phi.${phi}
  out=PRSCSx_metaanalysis_preec_phi.${phi}
  # chr=22
  
  # sumstat2=${sumstat2} base_pop2=${base_pop2} n_trn2=${n_trn2} sumstat1=${sumstat1} base_pop1=${base_pop1} n_trn1=${n_trn1} testplink=${testplink} bim_prefix=${bim_prefix} trait=${trait} A1=${A1} A2=${A2} SNP=${SNP} BETA=${BETA} P=${P} n_trn=${n_trn} phi=${phi} base_pop1=${base_pop1} base_pop2=${base_pop2} out=${out} $WDIR/Projects/HDP/scripts/PRSCSx.sh
  # qsub tmp.sh
  
  cd /broad/hptmp/btruong/validation_ukb/PRSCS/posterior/
  Rscript /medpop/esp2/btruong/scripts/concatenatePRS.R --file ${out}_${base_pop1}_chr{chr}.profile --out ${out}_${base_pop1}_allchr.profile
  Rscript /medpop/esp2/btruong/scripts/concatenatePRS.R --file ${out}_${base_pop2}_chr{chr}.profile --out ${out}_${base_pop2}_allchr.profile
  
done



cp PRSCSx*allchr.profile /medpop/esp2/btruong/Projects/HDP/results/UKB_PRS/PRSCSx

# preecl

trait=preec
sumstat=/broad/hptmp/btruong/HDP/data/metal_preec_European.Hispanic.Asian_1_leadvar500kb.txt

# geshtn
trait=geshtn
sumstat=/broad/hptmp/btruong/HDP/data/metal_geshtn_European.Hispanic.Asian_1_leadvar500kb.txt

# gene_name_list="ABHD16A ACAD10 ACP5 ACTRT3 AIF1 AKTIP ALDH2 ANGPTL8 ANTXR2 APOM ATP5F1E ATP6V1G2 ATXN2 BAG6 BRAP C6orf15 C6orf47 CCDC159 CCHCR1 CDSN CLIC1 CNN1 CSNK2B CTSZ CUX2 DDAH2 DDR1 DDX39B DOCK6 ECSIT EDN3 ELAVL3 ELOF1 EPOR FGF5 FLT1 FTO GNAS GPANK1 GTF2H4 HLA-B HLA-C HSPA1A HSPA1B HSPA1L KANK2 LDLR LRRC31 LRRC34 LRRIQ4 LSM2 LST1 LTA LTB LY6G5B LY6G5C LY6G6C LY6G6D LY6G6F MCCD1 MECOM MICA MICB MPIG6B MSH5 MUC21 MUC22 MUCL3 MYL2 MYNN NCR3 NELFCD NEU1 NFKBIL1 NPEPL1 ODAD3 PAN3 PGR PHETA1 PLPPR2 POMP POU5F1 PRDM8 PRELID3B PRKCSH PRRC2A PSORS1C1 PSORS1C2 RAB3D RBL2 RGL3 RPGRIP1L SAPCD1 SFTA2 SH2B3 SLC44A4 SLC46A3 SMARCA4 SPC24 SWSAP1 TCF19 TIMM29 TMEM205 TNF TRPC6 TSPAN16 TUBB1 VARS1 VARS2 VWA7 YIPF2 ZNF439 ZNF440 ZNF441 ZNF491 ZNF627 ZNF653 ZNF69 ZNF823 ZNF831"


# gene_name_list="AAAS ABCB9 AMHR2 ANTXR2 ARL6IP4 ATF7 ATP5F1E CCDC62 CDK2AP1 CSAD CTPS1 CTSZ DENR EDN2 EDN3 EIF4B ESPL1 FGF5 FOXO6 GNAS HCAR1 HCAR2 HCAR3 HIP1R HIVEP3 IGFBP6 ITGB7 KNTC1 KRT1 KRT18 KRT3 KRT4 KRT76 KRT77 KRT78 KRT79 KRT8 MAP3K12 MFSD5 MPHOSPH9 MTRFR MYG1 NELFCD NPEPL1 NPFF NPR3 OGFOD2 PCBP2 PFDN5 PITPNM2 PRDM8 PRELID3B PRR13 RARG RSRC2 SCMH1 SLFNL1 SOAT2 SP1 SP7 SPRYD3 STX16 SUB1 TARBP2 TNS2 TUBB1 VPS37B ZCCHC8 ZFR ZNF740 ZNF831"
# gene_name_list="FGF5"

# gene_name=ZNF831
# gene_name=SH2B3
# gene_name=RGL3
# gene_name=MECOM
# gene_name=CSAD

# tissue_file=/medpop/esp2/projects/GTEx/v8/eQTL_all_associations/Artery_Coronary.allpairs.txt.gz
# N_eqtl_tissue=213
# tissue=Artery_Coronary


# tissue_file=/medpop/esp2/projects/GTEx/v8/eQTL_all_associations/Adipose_Subcutaneous.allpairs.txt.gz
# N_eqtl_tissue=581
# tissue=Adipose_Subcutaneous


# tissue=Artery_Coronary
# tissue=Adipose_Subcutaneous
# tissue_list="Artery_Coronary Adipose_Subcutaneous Kidney_Cortex"
# tissue_list="Artery_Coronary,Adipose_Subcutaneous,Kidney_Cortex,Adipose_Visceral,Artery_Tibial,Uterus,Thyroid,Vagina"
# tissue_list="Kidney_Cortex"

# tissue_list="Whole_Blood,Vagina,Uterus,Thyroid,Testis,Stomach,Spleen,Small_Intestine_Terminal_Ileum,Skin_Sun_Exposed_Lower_leg,Skin_Not_Sun_Exposed_Suprapubic,Prostate,Pituitary,Pancreas,Ovary,Nerve_Tibial,Muscle_Skeletal,Minor_Salivary_Gland,Lung,Liver,Kidney_Cortex,Heart_Left_Ventricle,Heart_Atrial_Appendage,Esophagus_Muscularis,Esophagus_Mucosa,Esophagus_Gastroesophageal_Junction,Colon_Transverse,Colon_Sigmoid,Cells_EBV-transformed_lymphocytes,Cells_Cultured_fibroblasts,Breast_Mammary_Tissue,Brain_Substantia_nigra,Brain_Spinal_cord_cervical_c-1,Brain_Putamen_basal_ganglia,Brain_Nucleus_accumbens_basal_ganglia,Brain_Hypothalamus,Brain_Hippocampus,Brain_Frontal_Cortex_BA9,Brain_Cortex,Brain_Cerebellum,Brain_Cerebellar_Hemisphere,Brain_Caudate_basal_ganglia,Brain_Anterior_cingulate_cortex_BA24,Brain_Amygdala,Artery_Tibial,Artery_Coronary,Artery_Aorta,Adrenal_Gland,Adipose_Visceral_Omentum,Adipose_Subcutaneous"

# tissue_list_list="Whole_Blood Vagina Uterus Thyroid Testis Stomach Spleen Small_Intestine_Terminal_Ileum Skin_Sun_Exposed_Lower_leg Skin_Not_Sun_Exposed_Suprapubic Prostate Pituitary Pancreas Ovary Nerve_Tibial Muscle_Skeletal Minor_Salivary_Gland Lung Liver Kidney_Cortex Heart_Left_Ventricle Heart_Atrial_Appendage Esophagus_Muscularis Esophagus_Mucosa Esophagus_Gastroesophageal_Junction Colon_Transverse Colon_Sigmoid Cells_EBV-transformed_lymphocytes Cells_Cultured_fibroblasts Breast_Mammary_Tissue Brain_Substantia_nigra Brain_Spinal_cord_cervical_c-1 Brain_Putamen_basal_ganglia Brain_Nucleus_accumbens_basal_ganglia Brain_Hypothalamus Brain_Hippocampus Brain_Frontal_Cortex_BA9 Brain_Cortex Brain_Cerebellum Brain_Cerebellar_Hemisphere Brain_Caudate_basal_ganglia Brain_Anterior_cingulate_cortex_BA24 Brain_Amygdala Artery_Tibial Artery_Coronary Artery_Aorta Adrenal_Gland Adipose_Visceral_Omentum Adipose_Subcutaneous"



# tissue_list_list="a Whole_Blood Vagina Uterus Thyroid Testis Stomach Spleen Small_Intestine_Terminal_Ileum Skin_Sun_Exposed_Lower_leg Skin_Not_Sun_Exposed_Suprapubic Prostate Pituitary Pancreas Ovary Nerve_Tibial Muscle_Skeletal Minor_Salivary_Gland Lung Liver Kidney_Cortex Heart_Left_Ventricle Heart_Atrial_Appendage Esophagus_Muscularis Esophagus_Mucosa Esophagus_Gastroesophageal_Junction Colon_Transverse Colon_Sigmoid Cells_EBV-transformed_lymphocytes Cells_Cultured_fibroblasts Breast_Mammary_Tissue Brain_Substantia_nigra Brain_Spinal_cord_cervical_c-1 Brain_Putamen_basal_ganglia Brain_Nucleus_accumbens_basal_ganglia Brain_Hypothalamus Brain_Hippocampus Brain_Frontal_Cortex_BA9 Brain_Cortex Brain_Cerebellum Brain_Cerebellar_Hemisphere Brain_Caudate_basal_ganglia Brain_Anterior_cingulate_cortex_BA24 Brain_Amygdala Artery_Tibial Artery_Coronary Artery_Aorta Adrenal_Gland Adipose_Visceral_Omentum Adipose_Subcutaneous"


eqtl_p_threshold=1
gwas_p_threshold=5e-9
# gwas_p_threshold=5e-8

coloc_p1=1e-04
coloc_p2=1e-04
coloc_p12=1e-05


cd /broad/hptmp/btruong/HDP/coloc/jobs            

# for gene_name in ${gene_name_list}; do
#   echo $gene_name

#   out=coloc_allnearbygenes_alleqtl_${trait}_${gene_name}_pgwas.${gwas_p_threshold}_peqtl.${eqtl_p_threshold}_p1.${coloc_p1}_p2.${coloc_p2}_p12.${coloc_p12}
#   out=coloc_allnearbygenes_alleqtl_pgwas5e-8_${trait}_${gene_name}_pgwas.${gwas_p_threshold}_peqtl.${eqtl_p_threshold}_p1.${coloc_p1}_p2.${coloc_p2}_p12.${coloc_p12}
#   # out=coloc_AlleqtlNearLeadVar_${trait}_${gene_name}_${tissue}_pgwas.${gwas_p_threshold}_peqtl.${eqtl_p_threshold}_p1.${coloc_p1}_p2.${coloc_p2}_p12.${coloc_p12}.txt
  
#   trait=${trait} sumstat=${sumstat} coloc_p1=${coloc_p1} coloc_p2=${coloc_p2} coloc_p12=${coloc_p12} gwas_p_threshold=${gwas_p_threshold} eqtl_p_threshold=${eqtl_p_threshold} tissue_list1=${tissue_list1} gene_name=${gene_name} out=${out} $WDIR/codes/job_coloc.sh
#   qsub -t 1-49 tmp.sh
  
# done


gene_nearby_file=/broad/hptmp/btruong/HDP/coloc/gene_nearby_p${gwas_p_threshold}_${trait}.tsv
out=coloc_allnearbygenes_alleqtl_pgwas${gwas_p_threshold}_${trait}


gene_nearby_file=${gene_nearby_file} trait=${trait} sumstat=${sumstat} coloc_p1=${coloc_p1} coloc_p2=${coloc_p2} coloc_p12=${coloc_p12} gwas_p_threshold=${gwas_p_threshold} eqtl_p_threshold=${eqtl_p_threshold} tissue_list1=${tissue_list1} gene_list=${gene_list} out=${out} $WDIR/codes/job_coloc.sh
cat tmp.sh



qsub -t 1-49 tmp.sh





cd /broad/hptmp/btruong/HDP/coloc

SGE_TASK_ID=1


gene_list=""
while read gene;
do
    gene_list=${gene_list}" ${gene}"
done < /broad/hptmp/btruong/HDP/coloc/gene_nearby_p5e-9_geshtn.tsv


lineno=0
while read line;
do
    params=($line)
    let lineno+=1
    if [[ ${lineno} == ${SGE_TASK_ID} ]];
    then
        tissue=${line}
    fi
done < tissue_gtex_list.tsv

echo ${tissue}
tissue=Minor_Salivary_Gland

Rscript $WDIR/codes/coloc_GTEx_allnearbyvariants.R   --trait geshtn   --sumstat /broad/hptmp/btruong/HDP/data/metal_geshtn_European.Hispanic.Asian_1_leadvar500kb.txt   --coloc_p1 1e-04   --coloc_p2 1e-04   --coloc_p12 1e-05   --gwas_p_threshold 5e-9   --eqtl_p_threshold 1   --tissue_list ${tissue}   --tissue_file ""   --gene_list "${gene_list}"   --out coloc_allnearbygenes_alleqtl_pgwas5e-9_geshtn










tissue=Whole_Blood
gene_list=" ABHD16A ACAD10 ACP5 ACTRT3 AIF1 AKTIP ALDH2 ANGPTL8 ANTXR2 APOM ATP5F1E ATP6V1G2 ATXN2 BAG6 BRAP C6orf15 C6orf47 CCDC159 CCHCR1 CDSN CLIC1 CNN1 CSNK2B CTSZ CUX2 DDAH2 DDR1 DDX39B DOCK6 ECSIT EDN3 ELAVL3 ELOF1 EPOR FGF5 FLT1 FTO GNAS GPANK1 GTF2H4 HLA-B HLA-C HSPA1A HSPA1B HSPA1L KANK2 LDLR LRRC31 LRRC34 LRRIQ4 LSM2 LST1 LTA LTB LY6G5B LY6G5C LY6G6C LY6G6D LY6G6F MCCD1 MECOM MICA MICB MPIG6B MSH5 MUC21 MUC22 MUCL3 MYL2 MYNN NCR3 NELFCD NEU1 NFKBIL1 NPEPL1 ODAD3 PAN3 PGR PHETA1 PLPPR2 POMP POU5F1 PRDM8 PRELID3B PRKCSH PRRC2A PSORS1C1 PSORS1C2 RAB3D RBL2 RGL3 RPGRIP1L SAPCD1 SFTA2 SH2B3 SLC44A4 SLC46A3 SMARCA4 SPC24 SWSAP1 TCF19 TIMM29 TMEM205 TNF TRPC6 TSPAN16 TUBB1 VARS1 VARS2 VWA7 YIPF2 ZNF439 ZNF440 ZNF441 ZNF491 ZNF627 ZNF653 ZNF69 ZNF823 ZNF831"

Rscript $WDIR/codes/coloc_GTEx_allnearbyvariants.R   --trait geshtn   --sumstat /broad/hptmp/btruong/HDP/data/metal_geshtn_European.Hispanic.Asian_1.txt   --coloc_p1 1e-04   --coloc_p2 1e-04   --coloc_p12 1e-05   --gwas_p_threshold 5e-9   --eqtl_p_threshold 1   --tissue_list ${tissue}   --tissue_file ""   --gene_list "${gene_list}"   --out coloc_allnearbygenes_alleqtl_pgwas5e-9_geshtn










while read tissue;
do
    
  gene_nearby_file=${gene_nearby_file} trait=${trait} sumstat=${sumstat} coloc_p1=${coloc_p1} coloc_p2=${coloc_p2} coloc_p12=${coloc_p12} gwas_p_threshold=${gwas_p_threshold} eqtl_p_threshold=${eqtl_p_threshold} tissue_list1=${tissue_list1} gene_list=${gene_list} out=${out} $WDIR/codes/job_coloc.sh
  
done < /broad/hptmp/btruong/HDP/data/tissue_gtex_list.tsv




trait=geshtn

SGE_TASK_ID=3

lineno=0
while read line;
do
    params=(\$line)
    let lineno+=1
    if [[ \${lineno} == \${SGE_TASK_ID} ]];
    then
        tissue=${line}
    fi
done < tissue_gtex_list.tsv



for job in {31510598..31510695}; do
  qdel $job
done

######## $WDIR/codes/job_coloc.sh

#$ -pe smp 4 -R y -binding linear:4

echo "#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -l h_vmem=16G
#$ -l h_rt=20:00:00
#$ -N ${trait}


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3


cd /broad/hptmp/btruong/HDP/coloc



gene_list=\"\"
while read gene;
do
    gene_list=\${gene_list}\" \${gene}\"
done < ${gene_nearby_file}


lineno=0
while read line;
do
    params=(\$line)
    let lineno+=1
    if [[ \${lineno} == \${SGE_TASK_ID} ]];
    then
        tissue=\${line}
    fi
done < tissue_gtex_list.tsv

echo \${tissue}

Rscript \$WDIR/codes/coloc_GTEx_allnearbyvariants.R \
  --trait ${trait} \
  --sumstat ${sumstat} \
  --coloc_p1 ${coloc_p1} \
  --coloc_p2 ${coloc_p2} \
  --coloc_p12 ${coloc_p12} \
  --gwas_p_threshold ${gwas_p_threshold} \
  --eqtl_p_threshold ${eqtl_p_threshold} \
  --tissue_list \${tissue} \
  --tissue_file \"\" \
  --gene_list \"\${gene_list}\" \
  --out ${out}
" > tmp.sh
chmod +x tmp.sh
cat tmp.sh


# --gene_name ${gene_name} \
# --tissue_file /broad/hptmp/btruong/HDP/data/pgen.1007799.s007/eQTL_placenta_hg38.txt \



##############################

$WDIR/codes/coloc_GTEx_allnearbyvariants.R \
  --trait ${trait} \
  --sumstat ${sumstat} \
  --coloc_p1 ${coloc_p1} \
  --coloc_p2 ${coloc_p2} \
  --coloc_p12 ${coloc_p12} \
  --gwas_p_threshold ${gwas_p_threshold} \
  --eqtl_p_threshold ${eqtl_p_threshold} \
  --tissue_list ${tissue} \
  --tissue_file "" \
  --gene_list "${gene_list}" \
  --out ${out}



cd /broad/hptmp/btruong/HDP/data


tissue_list_list=(a Whole_Blood Vagina Uterus Thyroid Testis Stomach Spleen Small_Intestine_Terminal_Ileum Skin_Sun_Exposed_Lower_leg Skin_Not_Sun_Exposed_Suprapubic Prostate Pituitary Pancreas Ovary Nerve_Tibial Muscle_Skeletal Minor_Salivary_Gland Lung Liver Kidney_Cortex Heart_Left_Ventricle Heart_Atrial_Appendage Esophagus_Muscularis Esophagus_Mucosa Esophagus_Gastroesophageal_Junction Colon_Transverse Colon_Sigmoid Cells_EBV-transformed_lymphocytes Cells_Cultured_fibroblasts Breast_Mammary_Tissue Brain_Substantia_nigra Brain_Spinal_cord_cervical_c-1 Brain_Putamen_basal_ganglia Brain_Nucleus_accumbens_basal_ganglia Brain_Hypothalamus Brain_Hippocampus Brain_Frontal_Cortex_BA9 Brain_Cortex Brain_Cerebellum Brain_Cerebellar_Hemisphere Brain_Caudate_basal_ganglia Brain_Anterior_cingulate_cortex_BA24 Brain_Amygdala Artery_Tibial Artery_Coronary Artery_Aorta Adrenal_Gland Adipose_Visceral_Omentum Adipose_Subcutaneous)

tissue_list=${tissue_list_list[${SGE_TASK_ID}]}
# tissue_list=placenta


Rscript $WDIR/codes/coloc_GTEx_allnearbyvariants.R   --trait preec   --sumstat /broad/hptmp/btruong/HDP/data/metal_preec_European.Hispanic.Asian_1.txt   --coloc_p1 1e-04   --coloc_p2 1e-04   --coloc_p12 1e-05   --gwas_p_threshold 5e-8   --eqtl_p_threshold 1   --tissue_list ${tissue_list}   --tissue_file    --gene_name ZNF653   --out coloc_allnearbygenes_alleqtl_pgwas5e-8_preec_ZNF653_pgwas.5e-8_peqtl.1_p1.1e-04_p2.1e-04_p12.1e-05




freq_diff

--maf 0.001 \


gcta \
--diff-freq 1 \
--bfile /broad/hptmp/btruong/HDP/pops/1kg_chr/g1000_eur.hg38.6 \
--chr 6 \
--cojo-file intermediate/eQTL_Whole_Blood_chr6_pos_31383987_ABHD16A.txt \
--cojo-cond intermediate/preec_chr6_pos_31383987_ABHD16A_Whole_Blood_topsnp.txt \
--out intermediate/Whole_Blood_chr_6_pos_31383987_ABHD16A_conditionalP






##############################


tissue_list_list=(Whole_Blood Vagina Uterus Thyroid Testis Stomach Spleen Small_Intestine_Terminal_Ileum Skin_Sun_Exposed_Lower_leg Skin_Not_Sun_Exposed_Suprapubic Prostate Pituitary Pancreas Ovary Nerve_Tibial Muscle_Skeletal Minor_Salivary_Gland Lung Liver Kidney_Cortex Heart_Left_Ventricle Heart_Atrial_Appendage Esophagus_Muscularis Esophagus_Mucosa Esophagus_Gastroesophageal_Junction Colon_Transverse Colon_Sigmoid Cells_EBV-transformed_lymphocytes Cells_Cultured_fibroblasts Breast_Mammary_Tissue Brain_Substantia_nigra Brain_Spinal_cord_cervical_c-1 Brain_Putamen_basal_ganglia Brain_Nucleus_accumbens_basal_ganglia Brain_Hypothalamus Brain_Hippocampus Brain_Frontal_Cortex_BA9 Brain_Cortex Brain_Cerebellum Brain_Cerebellar_Hemisphere Brain_Caudate_basal_ganglia Brain_Anterior_cingulate_cortex_BA24 Brain_Amygdala Artery_Tibial Artery_Coronary Artery_Aorta Adrenal_Gland Adipose_Visceral_Omentum Adipose_Subcutaneous)


tissue_list=${tissue_list_list[${SGE_TASK_ID}]}




use .perl-5.28.0


vep \
  -i sigsnp_preec_European.Hispanic.Asian.txt \
  -o resvep_sigsnp_preec.txt \
  --max_af --everything \
  --force_overwrite --cache \
  --pick --tab --plugin Phenotypes \
  --fields Uploaded_variation,SYMBOL,CANONICAL,Consequence,BIOTYPE,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,MAX_AF,PHENOTYPES \
  --dir_cache $veppath/cache \
  --dir_plugin /medpop/esp2/btruong/Tools/VEPplugins 


filter_vep -i resvep_sigsnp_preec.txt --force_overwrite -filter "MAX_AF > 0.01 or not MAX_AF" -o resvep_sigsnp_preec_filtered.txt

  
echo -e "CHR\tSNP\tPOS\teffect_allele\tnoneffect_allele\tBETA" > PRSCS_preec_phi.1e-04_allchr.txt
echo -e "CHR\tSNP\tPOS\teffect_allele\tnoneffect_allele\tBETA" >  PRSCS_geshtn_phi.1e-06_allchr.txt
echo -e "CHR\tSNP\tPOS\teffect_allele\tnoneffect_allele\tBETA" >  PRSCS_SBP.MVP_phi.1e-04_allchr.txt

for chr in {1..22};
do
  cat PRS-CS_preec_phi.1e-04_pst_eff_a1_b0.5_phi1e-04_chr${chr}.txt >>PRSCS_preec_phi.1e-04_allchr.txt
  cat PRS-CS_geshtn_phi.1e-06_pst_eff_a1_b0.5_phi1e-06_chr${chr}.txt >>PRSCS_geshtn_phi.1e-06_allchr.txt
  cat PRSCS_SBP.MVP_phi.1e-04_chr${chr}_posterior_pst_eff_a1_b0.5_phi1e-04_chr${chr}.txt >> PRSCS_SBP.MVP_phi.1e-04_allchr.txt
done

cp PRSCS_SBP.MVP_phi.1e-04_allchr.txt $WDIR/Projects/HDP/results/PRSCS_weights
cp PRSCS_geshtn_phi.1e-06_allchr.txt $WDIR/Projects/HDP/results/PRSCS_weights
cp PRSCS_preec_phi.1e-04_allchr.txt $WDIR/Projects/HDP/results/PRSCS_weights







awk '(FNR!=1){
  split($5,a,":")
  print "chr"a[1]"\t"a[2]"\t"a[2]+1"\t"$5
}' eQTL_placenta_final.txt > eQTL_placenta_final_forlift.txt

liftOver eQTL_placenta_final_forlift.txt /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz eQTL_placenta_final_lifted.txt eQTL_placenta_final_unlifted.txt




awk '(NR==FNR){newpos[$4]=$2;next}{
  if (FNR==1) print "Beta\tSE\tpval_nominal\tgene_id\tvariant_id\tmaf"
  if (newpos[$5]) {
    split($5,a,":")
    newid = "chr"a[1]"_"newpos[$5]"_"a[4]"_"a[3]
    $5 = newid
    print $0
  }
}' OFS="\t" eQTL_placenta_final_lifted.txt eQTL_placenta_final.txt > eQTL_placenta_hg38.txt





/magma\
  --bfile /medpop/esp2/SarahUrbut/1kgEUR/1000G.EUR \
  --gene-annot magma_0kb.genes.annot\
  --pval AFib.sumstats ncol=N\
  --gene-model snp-wise=mean\
  --out AFib



SNP ID, chromosome, and base pair position

3 1 2


######################################################


awk '(NR==FNR){ensg[$2]=$1;next}{
  if (ensg[$6]) {
    $1=ensg[$6]
    print $0
  }
}' OFS="\t" /broad/hptmp/btruong/HDP/pops/gene_annot_jun10.txt /broad/hptmp/btruong/HDP/pops/NCBI37.3.gene.loc > /broad/hptmp/btruong/HDP/pops/NCBI37.3.gene.loc_1




trait=geshtn
# trait=preec
cd /broad/hptmp/btruong/HDP/pops/${trait}

# $(NF-1)>0.01 && $(NF-1)<0.99 && 

awk -v FS='\t' '{
  if (FNR==1) {print "SNP\tCHR\tPOS\tP"; next}
  if ($NF!="NA" && $2!="NA" && $1!="NA" && $NF!="" && $2!="" && $1!="" && $10!="" && $10!="NA") {
    $1=sprintf("%d",$1)
    $2=sprintf("%d",$2)
    print $NF,$1,$2,$10
  }
  }' OFS="\t" /broad/hptmp/btruong/HDP/data/metal_${trait}_European.Hispanic.Asian_1_withRSID.txt > ${trait}.snploc


trait=geshtn
awk 'BEGIN{minp=1}($1==5 && minp>$10){minp=$10;rsid=$3}END{print rsid," ",minp}' /broad/hptmp/btruong/HDP/data/metal_${trait}_European.Hispanic.Asian_1_withRSID.txt 

awk ''

awk '(FNR!=1){print "chr"$2,$3,$3+1,$1}' OFS="\t" ${trait}.snploc > ${trait}.snploc_4lift
liftOver ${trait}.snploc_4lift /medpop/esp2/btruong/Tools/hg38ToHg19.over.chain.gz ${trait}.snploc_4lift_lifted ${trait}.snploc_4lift_unlifted


awk '(NR==FNR){pos[$4]=$2;next}{
  if (FNR==1) {print $0;next}
  if (pos[$1]) {
    $3=pos[$1]
    print $0
  }
}' OFS="\t" ${trait}.snploc_4lift_lifted ${trait}.snploc > ${trait}.snploc_hg19



geneloc=/broad/hptmp/btruong/HDP/pops/NCBI37.3.gene.loc


/medpop/esp2/btruong/Tools/magma --annotate window=500,500 --snp-loc ${trait}.snploc_hg19 --gene-loc ${geneloc}_1 --out magma_annot_${trait}
/medpop/esp2/btruong/Tools/magma --annotate --snp-loc ${trait}.snploc_hg19 --gene-loc ${geneloc}_1 --out magma_annot_${trait}_0kb


cp */*.preds $WDIR/Projects/HDP/results/pops
cp */*.coefs $WDIR/Projects/HDP/results/pops
cp */*.marginals $WDIR/Projects/HDP/results/pops


# python /medpop/esp2/btruong/Tools/pops/munge_feature_directory.py \
#    --gene_annot_path gene_annot_jun10.txt \
#    --nan_policy zero \
#    --feature_dir /broad/hptmp/btruong/HDP/pops/gtex_expression \
#    --save_prefix ./features_munged/pops_features \
#    --max_cols 5000

        
# python /medpop/esp2/btruong/Tools/pops/munge_feature_directory.py \
#    --gene_annot_path gene_annot_jun10.txt \
#    --nan_policy zero \
#    --feature_dir /broad/hptmp/btruong/HDP/pops/features_files/ \
#    --save_prefix ./features_munged_1/pops_features \
#    --max_cols 5000



# /medpop/esp2/btruong/Tools/magma \
#   --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur \
#   --gene-annot magma_annot_${trait}.genes.annot \
#   --pval ${trait}.snploc N=450000 \
#   --annotate window=500,500 \
#   --gene-model snp-wise=mean\
#   --out ${trait}


        
# python -m pdb /medpop/esp2/btruong/Tools/pops/pops.py \
# python /medpop/esp2/btruong/Tools/pops/pops.py \
#  --gene_annot_path ../gene_annot_jun10.txt \
#  --feature_mat_prefix /broad/hptmp/btruong/HDP/pops/features_munged/pops_features \
#  --num_feature_chunks 4 \
#  --magma_prefix ${trait} \
#  --out_prefix ${trait}_pops_step2


####################################

# $ -pe smp 4 -R y -binding linear:4


####################################

#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=16G
#$ -l h_rt=20:00:00
#$ -N preec


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3


trait=preec

cd /broad/hptmp/btruong/HDP/pops/${trait}

geneloc=/broad/hptmp/btruong/HDP/pops/NCBI37.3.gene.loc


# /medpop/esp2/btruong/Tools/magma \
#   --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur \
#   --gene-annot magma_annot_${trait}_0kb.genes.annot \
#   --pval ${trait}.snploc_hg19 N=450000 \
#   --batch ${SGE_TASK_ID} 27 \
#   --gene-model snp-wise=mean\
#   --out ${trait}_0kb

# echo "finished magma"

python /medpop/esp2/btruong/Tools/pops/pops.py \
 --gene_annot_path /broad/hptmp/btruong/HDP/pops/gene_annot_jun10.txt \
 --feature_mat_prefix /broad/hptmp/btruong/HDP/pops/features_munged_1/pops_features \
 --num_feature_chunks 12 \
 --magma_prefix ${trait}_0kb_combine \
 --feature_selection_p_cutoff 5e-9 \
 --control_features_path /medpop/esp2/projects/software/PoPS/data/control.features \
 --out_prefix ${trait}_0kb_p5e-9_pops_step2

echo "finished pops"

###########################################################

qsub -t 9,24 pops.sh



trait=preec
trait=geshtn
cd /broad/hptmp/btruong/HDP/pops/${trait}
/medpop/esp2/btruong/Tools/magma --merge ${trait} --out ${trait}_combine











# python /medpop/esp2/btruong/Tools/pops/munge_feature_directory.py \
#    --gene_annot_path gene_annot_jun10.txt \
#    --nan_policy zero \
#    --feature_dir /broad/hptmp/btruong/HDP/pops/features_files/ \
#    --save_prefix ./features_munged_1/pops_features \
#    --max_cols 5000


chr=22

python /medpop/esp2/SarahUrbut/pops/pops.predict_scores.py \
  --gene_loc /medpop/esp2/projects/software/PoPS/data/gene_loc.txt \
  --features /broad/hptmp/btruong/HDP/pops/PoPS.features.parquet \
  --gene_results ${trait} \
  --chromosome ${chr} \
  --out ${trait}


###########################


qsub -t 1-22 pops.sh
qsub -t 1-49 tmp.sh



dd = read.delim("preec_pops_step2.preds")
dd = dd[order(dd$PoPS_Score, decreasing=T),]
head(dd)

head -n 1000 geshtn.snploc_hg19 | awk 'BEGIN{OFS=" "}{$1=$1}1' > geshtn.snploc_hg19_sub

/medpop/esp2/btruong/Tools/magma \
  --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur \
  --gene-annot magma_annot_geshtn.genes.annot \
  --pval geshtn.snploc_hg19_sub1 N=450000 \
  --batch 24 27 \
  --gene-model snp-wise=mean \
  --out geshtn1



awk 'BEGIN{count=0}(NR==FNR){snp[$2]=1;next}{if (snp[$1]) print $0}' OFS="\t" /broad/hptmp/btruong/HDP/pops/g1000_eur.bim geshtn.snploc_hg19_sub > geshtn.snploc_hg19_sub1





plink2 --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur --freq --out /broad/hptmp/btruong/HDP/pops/g1000_eur


awk 'BEGIN{count=0}(NR==FNR){snp[$1]=1;next}{if(snp[$2])count++}END{print count}' ${trait}.snploc_hg19 /broad/hptmp/btruong/HDP/pops/g1000_eur.bim 






trait=geshtn
cd /broad/hptmp/btruong/HDP/pops/${trait}
plink --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur --clump ${trait}.snploc_hg19 --clump-p1 5e-9 --clump-r2 0.1 --clump-kb 500 --out ${trait}_p.5e-9_r2.0.1_kb.500



plink --bfile g1000_eur --recode tab --out g1000_eur.tab

python /medpop/esp2/btruong/Tools/liftOverPlink/liftOverPlink.py --map g1000_eur.tab.map --out g1000_eur.hg38 --chain /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz
python /medpop/esp2/btruong/Tools/liftOverPlink/rmBadLifts.py --map g1000_eur.hg38.map --out good_g1000_eur.hg38.map --log bad_g1000_eur.hg38.dat

cut -f 2 bad_g1000_eur.hg38.dat > to_exclude.dat
cut -f 4 g1000_eur.hg38.bed.unlifted | sed "/^#/d" >> to_exclude.dat 

# Note: this will clobber the g1000_eur.hg38 MAP file generated by `liftOverPlink`:
plink --file g1000_eur.tab --recode --out g1000_eur.hg37 --exclude to_exclude.dat

plink --ped g1000_eur.hg37.ped --map g1000_eur.hg38.map --recode --allow-extra-chr --make-bed --out g1000_eur.hg38


sed 's/^chrM\s/25\t/g; s/^chrX\s/23\t/g; s/^chrY\s/24\t/g; s/^chr//g' g1000_eur.hg38.bim > fixed.bim

mv g1000_eur.hg38-temporary.bed g1000_eur.hg38.bed
mv g1000_eur.hg38-temporary.bim g1000_eur.hg38.bim
mv g1000_eur.hg38-temporary.fam g1000_eur.hg38.fam



plink --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur.hg38 --allow-extra-chr --exclude remove_snp.txt --make-bed --out g1000_eur.hg38.new


mv g1000_eur.hg38.new.bed g1000_eur.hg38.bed
mv g1000_eur.hg38.new.bim g1000_eur.hg38.bim
mv g1000_eur.hg38.new.fam g1000_eur.hg38.fam


cd /broad/hptmp/btruong/HDP/pops/1kg_chr
for chr in {1..22};
do
  plink --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur.hg38 --chr ${chr} --make-bed --out g1000_eur.hg38.${chr}
done




########################
mqtl
pqtl /medpop/esp2/projects/pQTL/Summary_sig/TOPMed/JHS/JHS_1E-5_pQTL.txt





#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=16G
#$ -l h_rt=20:00:00
#$ -N greml


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3

cd /broad/hptmp/btruong/HDP/greml

bfile=/medpop/esp2/projects/MGB_Biobank/genotype/30K_GSA/result/merged_hg38/GSA_30K_hg38

chr=$SGE_TASK_ID

gcta64 --bfile $bfile --autosome --maf 0.01 --make-grm --out grm_mgbb --thread-num 10 --chr $chr




ls -lt /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/HDP-mssm/*_output.txt | awk '{print $NF}' > all_saige.txt



ls -lt /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/*logistic_newRSID | awk '{print $NF}' > all_logistic.txt


rm pmin_biome.txt
rm pmin_biome_logistic.txt

while read pp; do
  # pp=/medpop/esp2/btruong/Projects/HDP/corrected_sumstat/HDP-mssm/BioMe_HA_composite_saige_output.txt
  echo $pp
  file=$(basename $pp .txt)
  # awk '{
  #   if (FNR==1) { print $0,"chrpos" }
  #   else {
  #     chr=substr($1,4)
  #     print $0,chr":"$2
  #   }
  # }' $pp > /medpop/esp2/btruong/Projects/HDP/corrected_sumstat/${file}_withCHRPOS
  awk -v file=$(basename $pp) 'BEGIN{minp=1}{
    if (FNR==1) { 
      # for (i=1;i<=NF;i++) if ($i=="p.value") pcol=i
      for (i=1;i<=NF;i++) if ($i=="P") pcol=i
    }
    else {
      if (minp>$pcol) {
        minp=$pcol
        snp=$0
      }
    }
  }END{print file,snp}' OFS=" " $pp >> pmin_biome_logistic.txt
# done < all_saige.txt
done < all_logistic.txt


ls *logistic_newRSID





########### HUNT


awk '(FNR!=1){print "chr"$1,$2,$2+1,$1":"$2}' preeclampsia_chr1-23.results.txt > preeclampsia_chr1-23.results.txt_forlift

liftOver preeclampsia_chr1-23.results.txt_forlift /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz preeclampsia_chr1-23.results.txt_lifted preeclampsia_chr1-23.results.txt_unlifted

awk '(NR==FNR){snp[$4]=$2;next}{
  if (FNR==1) {print $0,"chrpos","newpos"}
  else {
    chrpos=$1":"$2
    if (snp[chrpos]) print $0,$1":"snp[chrpos],snp[chrpos]
  }
}' preeclampsia_chr1-23.results.txt_lifted preeclampsia_chr1-23.results.txt > preeclampsia_chr1-23.results_hg38.txt











################ REGENIE UKB


echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N rgn_${trait}
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=8G
#$ -j y

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2


mkdir -p ${wdir}
cd ${wdir}


sif=/medpop/esp2/projects/software/singularity/regenie/v3.0.3/regenie.v3.0.3.sif

singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp \$sif \
regenie \
  --step 1 \
  --bed ${bed} \
  --phenoFile ${phenofile} \
  --phenoCol ${trait} \
  --covarFile ${covarfile} \
  --bsize 1000 \
  --extract ${snplist} \
  --keep ${samplist} \
  --bt \
  --maxCatLevels 1000 \
  --covarColList ${covarCol} \
  --catCovarList ${catcovar} \
  --out ${out}
" > tmp.sh
cat tmp.sh


qsub tmp.sh

###############################
rm mergelist.txt
for chr in {1..22}; do 
  ln -s /broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed ukb_chr${chr}_v2.bed
  ln -s /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim ukb_chr${chr}_v2.bim
  ln -s /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam ukb_chr${chr}_v2.fam
  echo "ukb_chr${chr}_v2" >> mergelist.txt
done


plink --merge-list mergelist.txt --make-bed --memory 64000 --out ukb_cal



plink2 \
  --bfile ukb_cal \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out qc_pass




trait=composite
trait=geshtn
trait=preec


wdir=/broad/hptmp/btruong/HDP/regenie/${trait}
phenofile=/medpop/esp2/btruong/Projects/HDP/data/${trait}_UKB_2.txt
covarfile=/medpop/esp2/btruong/Projects/HDP/data/${trait}_UKB_2.txt
covarCol=age,PC{1:10}
catcovar=batch,array
snplist=/broad/hptmp/btruong/HDP/regenie/plink_ukb/qc_pass.snplist
samplist=/broad/hptmp/btruong/HDP/regenie/plink_ukb/qc_pass.id
out=${trait}
bed=/broad/hptmp/btruong/HDP/regenie/plink_ukb/ukb_cal

mkdir -p ${wdir}
cd ${wdir}

trait=${trait} wdir=${wdir} bed=${bed} phenofile=${phenofile} covarfile=${covarfile} snpextract=${snpextract} covarCol=${covarCol} catcovar=${catcovar} snplist=${snplist} samplist=${samplist} out=${out} $WDIR/scripts/regenie_step1.sh

qsub tmp.sh




echo "#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N rgn_${trait}
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=8G
#$ -j y
#$ -t 1-22

source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2


mkdir -p ${wdir}
cd ${wdir}

chr=\${SGE_TASK_ID}


sif=/medpop/esp2/projects/software/singularity/regenie/v3.0.3/regenie.v3.0.3.sif

singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp,/broad/ukbb/:/broad/ukbb/ \$sif \
regenie \
  --step 2 \
  --bgen /broad/ukbb/imputed_v3/ukb_imp_chr\${chr}_v3.bgen \
  --phenoFile ${phenofile} \
  --phenoCol ${trait} \
  --covarFile ${covarfile} \
  --covarColList ${covarCol} \
  --catCovarList ${catcovar} \
  --maxCatLevels 10000 \
  --sample /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr15_v3_s487395.sample \
  --bsize 200 \
  --bt \
  --firth \
  --approx \
  --pThresh 0.01 \
  --maxiter-null 10000 \
  --maxstep-null 2 \
  --pred ${trait}_pred.list \
  --out ${trait}_chr\${chr}


" > tmp.sh
cat tmp.sh




trait=composite
trait=geshtn
trait=preec


wdir=/broad/hptmp/btruong/HDP/regenie/${trait}
phenofile=/medpop/esp2/btruong/Projects/HDP/data/${trait}_UKB_2.txt
covarfile=/medpop/esp2/btruong/Projects/HDP/data/${trait}_UKB_2.txt
covarCol=age,PC{1:10}
catcovar=batch,array
snplist=/broad/hptmp/btruong/HDP/regenie/plink_ukb/qc_pass.snplist
samplist=/broad/hptmp/btruong/HDP/regenie/plink_ukb/qc_pass.id
out=${trait}
bed=/broad/hptmp/btruong/HDP/regenie/plink_ukb/ukb_cal

mkdir -p ${wdir}
cd ${wdir}

trait=${trait} wdir=${wdir} bed=${bed} phenofile=${phenofile} covarfile=${covarfile} snpextract=${snpextract} covarCol=${covarCol} catcovar=${catcovar} snplist=${snplist} samplist=${samplist} out=${out} $WDIR/scripts/regenie_step2.sh

qsub -hold_jid 32091422 tmp.sh





chr=22
sif=/medpop/esp2/projects/software/singularity/regenie/v3.0.3/regenie.v3.0.3.sif

singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp,/broad/ukbb/:/broad/ukbb/ $sif regenie   --step 2   --bgen /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen   --phenoFile /medpop/esp2/btruong/Projects/HDP/data/geshtn_UKB_2.txt   --phenoCol geshtn   --covarFile /medpop/esp2/btruong/Projects/HDP/data/geshtn_UKB_2.txt   --covarColList age,PC{1:10}   --catCovarList sex,batch   --sample /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr15_v3_s487395.sample   --bsize 200   --bt   --firth   --approx   --pThresh 0.01   --maxiter-null 10000   --maxstep-null 2   --pred geshtn_pred.list   --out geshtn_chr${chr}


CXX=.g++-7.3.0 CC=.gcc-7.3.0 cmake ..

use .eigen-3.3.7
use .spectra-0.8.1
use .boost-build-2.0
use .mkl-2019.3.199
export CC=`which gcc`
export CXX=`which g++`

export EIGEN3_INCLUDE_DIR=/broad/software/free/Linux/redhat_7_x86_64/pkgs/eigen_3.3.7/
export SPECTRA_LIB=/broad/software/free/Linux/redhat_7_x86_64/pkgs/spectra-0.8.1/
export BOOST_LIB=/broad/software/free/Linux/redhat_7_x86_64/pkgs/boost_build-2.0/
export MKLROOT=/broad/software/free/Linux/redhat_7_x86_64/pkgs/mkl_2019.3.199/
export EXP=.spectra-0.8.1



rm -r *
cmake ..


cmd="gcta64 --bfile /broad/hptmp/btruong/HDP/pops/1kg_chr/g1000_eur.hg38.12 --chr 12 --maf 0.001 --diff-freq 1 --cojo-file intermediate/eQTL_Minor_Salivary_Gland_chr12_pos_53063801_AAAS.txt --cojo-cond intermediate/geshtn_chr12_pos_53063801_AAAS_Minor_Salivary_Gland_topsnp.txt --out intermediate/Minor_Salivary_Gland_chr_12_pos_53063801_AAAS_conditionalP"








head /broad/hptmp/btruong/HDP/data/metal_geshtn_European.Hispanic.Asian_1_leadvar500kb.txt
head /broad/hptmp/btruong/HDP/data/metal_preec_European.Hispanic.Asian_1_leadvar500kb.txt

tail -n +2 /broad/hptmp/btruong/HDP/data/metal_preec_European.Hispanic.Asian_1_leadvar500kb.txt | cut -f3 > preec_lead1.txt
tail -n +2 /broad/hptmp/btruong/HDP/data/metal_geshtn_European.Hispanic.Asian_1_leadvar500kb.txt | cut -f3  > geshtn_lead1.txt

cat preec_lead1.txt geshtn_lead1.txt | sort -u > both_lead.txt


awk '{
  split($1,a,":")
  print "chr"a[1],a[2],a[2]+1,$0
}' OFS="\t" both_lead.txt > extract_both_lead_hg38.txt


liftOver extract_both_lead_hg38.txt /medpop/esp2/btruong/Tools/hg38ToHg19.over.chain.gz extract_both_lead_hg38_lifted.txt extract_both_lead_hg38_unlifted.txt







# awk '(NR==FNR){snp[$1]=1;next}{
#   if (snp[$4]) print $3
# }' OFS="\t" both_lead.txt /medpop/esp2/btruong/Tools/hg38_common_chrpos.txt > extract_both_lead_hg38.txt

# grep rs10747667 /medpop/esp2/btruong/Tools/hg38_common_chrpos.txt


awk '{
  split($4,a,":")
  print a[1]":"$2 > "/broad/hptmp/btruong/HDP/extract_ukb_10k/leadsnp500kb_hg37_chr"a[1]".txt"
}' extract_both_lead_hg38_lifted.txt








#!/bin/bash -l
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -N extract_lead500kb
#$ -l h_rt=5:00:00
#$ -l s_rt=5:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=8G
#$ -j y
#$ -t 1-22


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2


cd /broad/hptmp/btruong/HDP/extract_ukb_10k/

chr=${SGE_TASK_ID}
/medpop/esp2/btruong/Tools/qctool/qctool -g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen -og extract_lead500kb_chr${chr} -s /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr15_v3_s487395.sample -ofiletype binary_ped -incl-positions leadsnp500kb_hg37_chr${chr}.txt 



-incl-samples UKB_white_10k.txt



########### R
library(dplyr)
library(data.table)


related = read.table("/medpop/esp2/wallace/projects/ukbb/relateds/data/removed.kin.3rd.gz")


file <- fread('/medpop/esp2/aniruddh/Lpa/UKBB_TG_BackgroundVariables_22020.txt') # background variables; not updated with Lpa (check notebook for that updated info)
file <- as.data.frame(file)

white <- subset(file, file$ethnicity==1 | file$ethnicity==1001 | file$ethnicity==1002 | file$ethnicity==1003)
white$race_specific <- "White"

fam = read.table("/medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr15_v3_s487395.sample")
fam = fam %>% filter(!V2 %in% related[,1])
fam = fam %>% filter(V2 %in% white$eid)
set.seed(1)
fam_sel = fam[sample(1:nrow(fam), 10000),]



write.table(fam_sel, "UKB_white_10k.txt", row.names=F, sep="\t", col.names=F, quote=F)





trait=preec
trait=geshtn

head -n1 ${trait}_chr1_${trait}.regenie > ${trait}_allchr.regenie

for chr in {1..22}; do
  echo $chr
  # cat geshtn_chr${chr}_geshtn.regenie | sort -nk13,13 | head -1 >> summary_5e-8.txt
  # awk '(-$13<=-7.30103){print $0}' geshtn_chr${chr}_geshtn.regenie >> summary_5e-8.txt
  awk '(FNR!=1){if ($NF!="TEST_FAIL") print $0}' ${trait}_chr${chr}_${trait}.regenie >> ${trait}_allchr.regenie
done


awk '{
  if (FNR==1) {print $0,"P"; next}
  if ($6>0.001 && $6<1-0.001) 
    print $0,10**(-$13)
  }' ${trait}_allchr.regenie > ${trait}_allchr_maf0.001.regenie


awk '(FNR!=1){print "chr"$1,$2,$2+1,$3}' ${trait}_allchr_maf0.001.regenie > ${trait}_forlift.bed


liftOver ${trait}_forlift.bed /medpop/esp2/btruong/Tools/hg19ToHg38.over.chain.gz ${trait}_lifted.bed ${trait}_unlifted.bed




awk '(NR==FNR){newpos[$NF]=$2;next}{
  if (FNR==1) print $0,"newpos","chrpos"
  if (newpos[$3])
    print $0,newpos[$3],$1":"newpos[$3]
}' ${trait}_lifted.bed ${trait}_allchr_maf0.001.regenie > ${trait}_allchr_maf0.001_hg38.regenie


cp ${trait}_allchr_maf0.001_hg38.regenie $WDIR/Projects/HDP/corrected_sumstat








