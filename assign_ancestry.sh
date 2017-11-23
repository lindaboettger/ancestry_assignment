#!/bin/bash

reuse PLINK2
reuse Tabix
reuse VCFtools
reuse GridEngine8
reuse R-3.3
reuse Python-2.7


### input file is binary plink, snps names are bp locations, check strandedness with ref dataset

#THIS IS THE RIGHT FILE

######### VARIABLES FED INTO THIS SCRIPT ###########
## example variables
# query_dogs_plink=/seq/vgb/linda/mutts/all_ancestry_dogs/DD_MM_10xg_refsnps_filt
# query_dog_names=/seq/vgb/linda/mutts/all_ancestry_dogs/DD_MM_10x_DogIDs_Fixed.txt
# outfile=/seq/vgb/linda/mutts/breed_calling_output/testn
# win=80
# fix_phase=T
# min_seg=1
# min_pct=0.01
# min_hap_len=1e6
# min_prob=0.5

#old variables not used in this version:
# rm_bad_post=T
#posterior_cutoff=0.7

#script_dir=/seq/vgb/linda/mutts/scripts/ancestry_assignment_pipeline_Eagle10it_v2

#set vcftools
vcftools=/home/unix/boettger/bin/vcftools_0.1.12b/bin/vcftools

## set reference datasets and paths
#plink_ref_forPhase=/seq/vgb/linda/mutts/breed_reference/allOstrander_Boyko_merged_my12perbreed_6forPHASE

## this didnt work because it sees it as in job script dir
# script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#error opening /var/spool/sge/vert03/job_scripts/query_phase_bychr_eagle.sh
outdir=$(dirname "${outfile}")

#make output dir
mkdir -p $outdir

####################### Filter SNPS #############################
# keep only snps in the reference dataset
# get rid of known regions with bad posteriors, if desired
# posteriors from this specific query dataset are not used because
# they will vary widely with sample size 

# #get snps in the reference set we will use
ref=( $(ls $refpath/*.tped) ) # var is an array witch contain the result of ls 
awk '{print $2}' ${ref[0]} > $outfile'_ref_snps.txt'


# keep only snps that are in our reference file
# get rid of missing data >0.1 (single bad snps can cause phase switching)
# had to add maf filter.... otherwize eagle "ERROR: Samples are not diploid"
plink --bfile $query_dogs_plink \
  --extract $outfile'_ref_snps.txt' \
  --bmerge $plink_ref_forPhase'.bed' $plink_ref_forPhase'.bim' $plink_ref_forPhase'.fam' \
  --dog \
  --make-bed \
  --recode vcf-iid \
  --geno 0.10 \
  --out $outfile'_filt'

# # which chromosomes do we want to run?
cut -f1 $outfile'_filt'.bim | sort -n | uniq > $outfile'_mychrs.txt'
cut -f1 $outfile'_filt'.bim | sort -n | uniq | awk '{ printf "%02i\n", $1}' > $outfile'_mychrs_alphanum.txt'

# set number of chromosomes
nchr=`cat $outfile'_mychrs.txt' | wc -l`

#get query snps
awk '{print $2}' $outfile'_filt'.bim > $outfile'_filt'_snps.txt

#get map for query snps
plink --bfile $plink_ref_forPhase \
  --extract $outfile'_filt'_snps.txt \
  --dog \
  --make-bed \
  --out $outfile'_genmap'


#make indexed bgzipped file 
bgzip $outfile'_filt'.vcf
#make index file
tabix -p vcf $outfile'_filt'.vcf.gz


#phase all chromosomes
while read c; do
qsub -S /bin/bash -v chr=$c,targ=$outfile'_filt' \
   -o /seq/vgb/linda/mutts/out/mutt_out \
   -e /seq/vgb/linda/mutts/err/mutt_err \
   -N 'chr_'${c} \
   -p -10 \
   -wd /seq/vgb/linda/mutts/err \
   $script_dir/query_phase_bychr_eagle.sh
done < $outfile'_mychrs.txt'

## check that length of dog chrs == len of outfiles before moving on

# check that there are still jobs in qstat and total out files < nchrs
while ! [ $nchr == `ls $outfile'_filt_'*.phased.vcf.gz -1 | wc -l` ] && ! [ `qstat | wc -l` == '0' ]
do 
    sleep 30
done

# if nchrs not equal out files now, then die with error
nfiles=`ls $outfile'_filt_'*.phased.vcf.gz -1 | wc -l`
if ! [ $nchr == $nfiles ]; then 
    echo "Number of chroms ($nchr) does not match number of phased outfiles ($nfiles)"
    exit 1
else
    echo "Phasing complete, all $nchr outfiles present."
fi

#get list of files to concatinate
ls -1v $outfile'_filt_'*.phased.vcf.gz > $outfile'_filt_phased_chrs.txt'

#### Combine all chromosomes
vcf-concat \
--files $outfile'_filt_phased_chrs.txt' \
 | bgzip > $outfile'_filt_allchr.phased.vcf.gz'

## remove breed phasing reference dogs
awk '{print $2}' $plink_ref_forPhase'.fam' > $outdir/phase_purebreds.txt



$vcftools --gzvcf $outfile'_filt_allchr.phased.vcf.gz' \
    --plink-tped \
    --remove $outdir/phase_purebreds.txt \
    --out $outfile'_filt_allchr_mydogs.phased'



## remove individual chromosomes
rm $outfile'_filt_'[0-9]*
rm $outfile'_filt_phased_chrs.txt'
rm $outfile'_filt_allchr.phased.vcf.gz'


#################  Convert to SupportMix Input #####################

# add genetic map (col 3 of $outfile'_genmap'.bim)
echo "$(awk 'FNR==NR{a[NR]=$3;next}{$3=a[FNR]}1' $outfile'_genmap.bim' $outfile'_filt_allchr_mydogs.phased.tped')" \
  > $outfile'_filt_allchr_mydogs.phased.tped'

# supportmix uses alphanumeric chr sorting, so append 0s
echo "$(awk '{$1 = sprintf("%02d", $1); print}' $outfile'_filt_allchr_mydogs.phased.tped')" > $outfile'_filt_allchr_mydogs.phased.tped'


# remove unneeded files
rm $outfile'_filt_snps.txt'
rm $outfile'_filt.bed'
rm $outfile'_filt.bim'
rm $outfile'_filt.fam'
rm $outfile'_filt.nosex'
rm $outfile'_filt.log'
rm $outfile'_genmap.bed'
rm $outfile'_genmap.fam'
rm $outfile'_genmap.nosex'
rm $outfile'_genmap.log'
rm $outfile'_filt.vcf.gz.tbi'
rm $outfile*.vcf.gz

# Make Supportmix input file for our reference panel #
echo 'supportmix=/seq/vgb/linda/bin/SupportMixDistribution/Application/SupportMix' > $outfile'_run_supportmix_bychr.sh'
printf 'spmix_command=' >> $outfile'_run_supportmix_bychr.sh'
printf '"' >> $outfile'_run_supportmix_bychr.sh'
echo '$supportmix \' >> $outfile'_run_supportmix_bychr.sh'
echo '    --save ${out}_${chr} \' >> $outfile'_run_supportmix_bychr.sh'
echo '    --chromosome=$chr \' >> $outfile'_run_supportmix_bychr.sh'
echo '    --window=$win \' >> $outfile'_run_supportmix_bychr.sh'
printf '    ' >> $outfile'_run_supportmix_bychr.sh'
# find all tpeds and add them as supportmix arguments (removing newlines)
ls -1 $refpath/*.tped | tr '\n' ' ' >> $outfile'_run_supportmix_bychr.sh'
# adding query dataset ($query.tped) to end of file
printf '$query.tped' >> $outfile'_run_supportmix_bychr.sh'
printf '"' >> $outfile'_run_supportmix_bychr.sh'
printf "\n" >> $outfile'_run_supportmix_bychr.sh'
# print the supportmix command and run it.
printf 'printf "$spmix_command"' >> $outfile'_run_supportmix_bychr.sh'
printf "\n" >> $outfile'_run_supportmix_bychr.sh'
printf '`$spmix_command`' >> $outfile'_run_supportmix_bychr.sh'
printf "\n" >> $outfile'_run_supportmix_bychr.sh'


### Start Supportmix run by chromosome ####
#script=$script_dir/run_supportmix_bychr_may21ref.sh
script=$outfile'_run_supportmix_bychr.sh'
query=$outfile'_filt_allchr_mydogs.phased'
out=$outfile'_SM'

### Start Supportmix run by chromosome ####
while read c; do
    echo "I got to chr $c"
    python $script_dir/qsub_watch.py $win $query $out $c $script &      
done < $outfile'_mychrs_alphanum.txt'
#/seq/vgb/linda/mutts/scripts/dog_chrs_alphanum.txt


# check that there are still python jobs with this outfile name is `ps -ef`; there will always
# be at least 1 because the original grep command is returned
while ! [ `ps -ef | grep $out | wc -l` == '1' ]
do 
    sleep 60
done

# if nchrs not equal out files now, then die with error
nfiles=`ls $outfile'_SM_'*.cfg -1 | wc -l`
if ! [ $nchr == $nfiles ]; then 
    echo "Number of chroms ($nchr) does not match number of Supportmix outfiles ($nfiles)"
    exit 1
else
    echo "Supportmix run complete, all $nchr outfiles present."
fi

###### Combine SM output across chromosomes ######

#### breed tpeds
cat \
 $outfile'_SM_01.tped' \
 $outfile'_SM_02.tped' \
 $outfile'_SM_03.tped' \
 $outfile'_SM_04.tped' \
 $outfile'_SM_05.tped' \
 $outfile'_SM_06.tped' \
 $outfile'_SM_07.tped' \
 $outfile'_SM_08.tped' \
 $outfile'_SM_09.tped' \
 $outfile'_SM_10.tped' \
 $outfile'_SM_11.tped' \
 $outfile'_SM_12.tped' \
 $outfile'_SM_13.tped' \
 $outfile'_SM_14.tped' \
 $outfile'_SM_15.tped' \
 $outfile'_SM_16.tped' \
 $outfile'_SM_17.tped' \
 $outfile'_SM_18.tped' \
 $outfile'_SM_19.tped' \
 $outfile'_SM_20.tped' \
 $outfile'_SM_21.tped' \
 $outfile'_SM_22.tped' \
 $outfile'_SM_23.tped' \
 $outfile'_SM_24.tped' \
 $outfile'_SM_25.tped' \
 $outfile'_SM_26.tped' \
 $outfile'_SM_27.tped' \
 $outfile'_SM_28.tped' \
 $outfile'_SM_29.tped' \
 $outfile'_SM_30.tped' \
 $outfile'_SM_31.tped' \
 $outfile'_SM_32.tped' \
 $outfile'_SM_33.tped' \
 $outfile'_SM_34.tped' \
 $outfile'_SM_35.tped' \
 $outfile'_SM_36.tped' \
 $outfile'_SM_37.tped' \
 $outfile'_SM_38.tped' \
 $outfile'_SM_41.tped' \
 > $outfile'_SM_allchr.tped'


 #### Posterior Probabilities
cat \
 $outfile'_SM_01.Probs.tped' \
 $outfile'_SM_02.Probs.tped' \
 $outfile'_SM_03.Probs.tped' \
 $outfile'_SM_04.Probs.tped' \
 $outfile'_SM_05.Probs.tped' \
 $outfile'_SM_06.Probs.tped' \
 $outfile'_SM_07.Probs.tped' \
 $outfile'_SM_08.Probs.tped' \
 $outfile'_SM_09.Probs.tped' \
 $outfile'_SM_10.Probs.tped' \
 $outfile'_SM_11.Probs.tped' \
 $outfile'_SM_12.Probs.tped' \
 $outfile'_SM_13.Probs.tped' \
 $outfile'_SM_14.Probs.tped' \
 $outfile'_SM_15.Probs.tped' \
 $outfile'_SM_16.Probs.tped' \
 $outfile'_SM_17.Probs.tped' \
 $outfile'_SM_18.Probs.tped' \
 $outfile'_SM_19.Probs.tped' \
 $outfile'_SM_20.Probs.tped' \
 $outfile'_SM_21.Probs.tped' \
 $outfile'_SM_22.Probs.tped' \
 $outfile'_SM_23.Probs.tped' \
 $outfile'_SM_24.Probs.tped' \
 $outfile'_SM_25.Probs.tped' \
 $outfile'_SM_26.Probs.tped' \
 $outfile'_SM_27.Probs.tped' \
 $outfile'_SM_28.Probs.tped' \
 $outfile'_SM_29.Probs.tped' \
 $outfile'_SM_30.Probs.tped' \
 $outfile'_SM_31.Probs.tped' \
 $outfile'_SM_32.Probs.tped' \
 $outfile'_SM_33.Probs.tped' \
 $outfile'_SM_34.Probs.tped' \
 $outfile'_SM_35.Probs.tped' \
 $outfile'_SM_36.Probs.tped' \
 $outfile'_SM_37.Probs.tped' \
 $outfile'_SM_38.Probs.tped' \
 $outfile'_SM_41.Probs.tped' \
 > $outfile'_SM_allchr.Probs.tped'

 # make new tfam 
 cp $outfile'_SM_01.tfam' $outfile'_SM_allchr.tfam'


# Fix population legend #
cp $outfile'_SM_01.population_legend.txt' $outfile'_SM_allchr.population_legend.txt'
popheader="popnum\tbreed_guess"
sed -i "s|$refpath||g" $outfile'_SM_allchr.population_legend.txt'
sed -i "s|.tped||g" $outfile'_SM_allchr.population_legend.txt'
sed -i "1s/.*/$popheader/" $outfile'_SM_allchr.population_legend.txt'
sed -i 's|/||g' $outfile'_SM_allchr.population_legend.txt'



## remove individual chromosome files
rm $outfile'_SM_'[0-9]*

### Make folders for plots ###
mkdir -p $outfile'_pie'
mkdir -p $outfile'_bar'
mkdir -p $outfile'_paint'
mkdir -p $outfile'_stats'

########## Call R to make plots ############
# Variables to pass:
# 1. breed legend
# 2. file prefix for Supportmix output files
# 3. markers file
# 4. the names file (id in tfam and actual name separated by tab)
# 5. minimum number of segments for breed to be included in plots
# 6. minimum percent of genome for breed to be included in plots
# 7. minimum length of haplotype to be included in plot
# 8. minimum posterior probability of haplotype breed call to be included in plot
# 9. outpath for bar charts
# 10. outpath for pie charts
# 11. outpath for chr paintings
# 12. mean posterior probability cutoff


Rscript $script_dir/make_dog_plots.R \
    $outfile'_SM_allchr.population_legend.txt' \
    $outfile'_SM_allchr' \
    $outfile'_genmap.bim' \
    $query_dog_names \
    $min_seg \
    $min_pct \
    $min_hap_len \
    $min_prob \
    $outfile'_bar/' \
    $outfile'_pie/' \
    $outfile'_paint/' \
    $outfile'_stats/' \
    $script_dir 
    #$posterior_cutoff

    ### it starts over???
#     Making plots for DX11 : Vera_Lynn_MM
# Making plots for DX02 : Hunter_MM
# Making plots for D157 : Beskow_MM
# Making plots for DX03 : Reily_MM
# Empty data.table (0 rows) of 1 col: ID
# Making plots for a550771-4306235-112317-918_A01.CEL : Hunter
# Making plots for a550771-4306235-112317-918_A02.CEL : Daytona
# Making plots for a550771-4306235-112317-918_A03.CEL : Lily2
# Making plots for a550771-4306235-112317-918_A04.CEL : Scarlet
# Making plots for a550771-4306235-112317-918_A05.CEL : Peso
# Making plots for a550771-4306235-112317-918_A06.CEL : Fendi
# Making plots for a550771-4306235-112317-918_A07.CEL : Bella
# Making plots for a550771-4306235-112317-918_A08.CEL : Andy
# Making plots for a550771-4306235-112317-918_A09.CEL : Lucy2
# Making plots for a550771-4306235-112317-918_A10.CEL : Barney
# Making plots for a550771-4306235-112317-918_A11.CEL : Moose
# Making plots for a550771-4306235-112317-918_A12.CEL : Tesla


# #concatinate the mean posterior files
# for filename in $outfile'_stats/'*_meanPosterior.txt
# do
#     dogname=$(echo $filename | perl -pe 's/.*\/([\w]+)_meanPosterior.txt/$1/')
#     perl -pe "s/^1/$dogname/" $filename | grep -Pv '^x'
# done > cmb_meanPost.txt

#delete the intermediate files


# ## re-add underscores in breed pct file
# pct_files=$outfile'_stats/'*_breedPcts.txt
# python $script_dir/fix_joint_stats.py $pct_files


# Adjust the phase to account for breeds if fix_phase==T
if [ $fix_phase == T ]; then
    # RUN THE R SCRIPT HERE THAT FIXES THE PHASING 
    echo "adjusting the phasing to maximize size of breed haplotypes"


    ########## Call R to fix phase ############
    # Variables to pass:
    # 1. breed legend
    # 2. file prefix for old tped from SM
    # 3. markers file 
    # 4. tped with snp calls (pre-SM)
    # 5. new tped with phase adjusted for breed calls


    Rscript $script_dir/fix_phase.R \
        $outfile'_SM_allchr.population_legend.txt' \
        $outfile'_SM_allchr' \
        $outfile'_genmap.bim' \
        $outfile'_allchr.phased' \
        $outfile'_phasefix'
    


    # make new tfam 
     cp $outfile'_allchr.phased.tfam' $outfile'_phasefix.tfam'

    # make chr alphanumeric
    echo "$(awk '{$1 = sprintf("%02d", $1); print}' $outfile'_phasefix.tped')" > $outfile'_phasefix.tped'


    script=$outfile'_run_supportmix_bychr.sh'
    query=$outfile'_phasefix'
    out=$outfile'_SM_phasefix'

    ### Start Supportmix run by chromosome ####
    while read c; do
    python $script_dir/qsub_watch.py $win $query $out $c $script &      
    done < /seq/vgb/linda/mutts/scripts/dog_chrs_alphanum.txt


    # check that there are still python jobs with this outfile name is `ps -ef`; there will always
    # be at least 1 because the original grep command is returned
    while ! [ `ps -ef | grep $out | wc -l` == '1' ]
    do 
        sleep 60
    done

    # if nchrs not equal out files now, then die with error
    nfiles=`ls $outfile'_SM_phasefix'*.cfg -1 | wc -l`
    if ! [ $nchr == $nfiles ]; then 
        echo "Number of chroms ($nchr) does not match number of Supportmix outfiles ($nfiles)"
        exit 1
    else
        echo "Supportmix run complete, all $nchr outfiles present."
    fi



    #### breed tpeds
    cat \
     $outfile'_SM_phasefix_01.tped' \
     $outfile'_SM_phasefix_02.tped' \
     $outfile'_SM_phasefix_03.tped' \
     $outfile'_SM_phasefix_04.tped' \
     $outfile'_SM_phasefix_05.tped' \
     $outfile'_SM_phasefix_06.tped' \
     $outfile'_SM_phasefix_07.tped' \
     $outfile'_SM_phasefix_08.tped' \
     $outfile'_SM_phasefix_09.tped' \
     $outfile'_SM_phasefix_10.tped' \
     $outfile'_SM_phasefix_11.tped' \
     $outfile'_SM_phasefix_12.tped' \
     $outfile'_SM_phasefix_13.tped' \
     $outfile'_SM_phasefix_14.tped' \
     $outfile'_SM_phasefix_15.tped' \
     $outfile'_SM_phasefix_16.tped' \
     $outfile'_SM_phasefix_17.tped' \
     $outfile'_SM_phasefix_18.tped' \
     $outfile'_SM_phasefix_19.tped' \
     $outfile'_SM_phasefix_20.tped' \
     $outfile'_SM_phasefix_21.tped' \
     $outfile'_SM_phasefix_22.tped' \
     $outfile'_SM_phasefix_23.tped' \
     $outfile'_SM_phasefix_24.tped' \
     $outfile'_SM_phasefix_25.tped' \
     $outfile'_SM_phasefix_26.tped' \
     $outfile'_SM_phasefix_27.tped' \
     $outfile'_SM_phasefix_28.tped' \
     $outfile'_SM_phasefix_29.tped' \
     $outfile'_SM_phasefix_30.tped' \
     $outfile'_SM_phasefix_31.tped' \
     $outfile'_SM_phasefix_32.tped' \
     $outfile'_SM_phasefix_33.tped' \
     $outfile'_SM_phasefix_34.tped' \
     $outfile'_SM_phasefix_35.tped' \
     $outfile'_SM_phasefix_36.tped' \
     $outfile'_SM_phasefix_37.tped' \
     $outfile'_SM_phasefix_38.tped' \
     $outfile'_SM_phasefix_41.tped' \
     > $outfile'_SM_phasefix_allchr.tped'


    #### Posterior Probabilities
    cat \
     $outfile'_SM_phasefix_01.Probs.tped' \
     $outfile'_SM_phasefix_02.Probs.tped' \
     $outfile'_SM_phasefix_03.Probs.tped' \
     $outfile'_SM_phasefix_04.Probs.tped' \
     $outfile'_SM_phasefix_05.Probs.tped' \
     $outfile'_SM_phasefix_06.Probs.tped' \
     $outfile'_SM_phasefix_07.Probs.tped' \
     $outfile'_SM_phasefix_08.Probs.tped' \
     $outfile'_SM_phasefix_09.Probs.tped' \
     $outfile'_SM_phasefix_10.Probs.tped' \
     $outfile'_SM_phasefix_11.Probs.tped' \
     $outfile'_SM_phasefix_12.Probs.tped' \
     $outfile'_SM_phasefix_13.Probs.tped' \
     $outfile'_SM_phasefix_14.Probs.tped' \
     $outfile'_SM_phasefix_15.Probs.tped' \
     $outfile'_SM_phasefix_16.Probs.tped' \
     $outfile'_SM_phasefix_17.Probs.tped' \
     $outfile'_SM_phasefix_18.Probs.tped' \
     $outfile'_SM_phasefix_19.Probs.tped' \
     $outfile'_SM_phasefix_20.Probs.tped' \
     $outfile'_SM_phasefix_21.Probs.tped' \
     $outfile'_SM_phasefix_22.Probs.tped' \
     $outfile'_SM_phasefix_23.Probs.tped' \
     $outfile'_SM_phasefix_24.Probs.tped' \
     $outfile'_SM_phasefix_25.Probs.tped' \
     $outfile'_SM_phasefix_26.Probs.tped' \
     $outfile'_SM_phasefix_27.Probs.tped' \
     $outfile'_SM_phasefix_28.Probs.tped' \
     $outfile'_SM_phasefix_29.Probs.tped' \
     $outfile'_SM_phasefix_30.Probs.tped' \
     $outfile'_SM_phasefix_31.Probs.tped' \
     $outfile'_SM_phasefix_32.Probs.tped' \
     $outfile'_SM_phasefix_33.Probs.tped' \
     $outfile'_SM_phasefix_34.Probs.tped' \
     $outfile'_SM_phasefix_35.Probs.tped' \
     $outfile'_SM_phasefix_36.Probs.tped' \
     $outfile'_SM_phasefix_37.Probs.tped' \
     $outfile'_SM_phasefix_38.Probs.tped' \
     $outfile'_SM_phasefix_41.Probs.tped' \
     > $outfile'_SM_phasefix_allchr.Probs.tped'

     # make new tfam 
     cp $outfile'_SM_phasefix_01.tfam' $outfile'_SM_phasefix_allchr.tfam'

     # Fix population legend #
     cp $outfile'_SM_phasefix_01.population_legend.txt' $outfile'_SM_phasefix_allchr.population_legend.txt'
     popheader="popnum\tbreed_guess"
     sed -i "s|$refpath/||g" $outfile'_SM_phasefix_allchr.population_legend.txt'
     sed -i "s|.tped||g" $outfile'_SM_phasefix_allchr.population_legend.txt'
     sed -i "1s/.*/$popheader/" $outfile'_SM_phasefix_allchr.population_legend.txt'
     sed -i 's|/||g' $outfile'_SM_phasefix_allchr.population_legend.txt'


    ### Make folders for the fixphase plots ###
     mkdir $outfile'_pie_fixphase'
     mkdir $outfile'_bar_fixphase'
     mkdir $outfile'_paint_fixphase'
     mkdir $outfile'_stats_fixphase'

    ## remove individual chromosome files
    rm $outfile'_SM_phasefix_'[0-9]*

    ########## Call R to make plots ############
    # Variables to pass:
    # 1. breed legend
    # 2. file prefix for Supportmix output files
    # 3. markers file
    # 4. the names file (id in tfam and actual name separated by tab)
    # 5. minimum number of segments for breed to be included in plots
    # 6. minimum percent of genome for breed to be included in plots
    # 7. minimum length of haplotype to be included in plot
    # 8. minimum posterior probability of haplotype breed call to be included in plot
    # 9. outpath for bar charts
    # 10. outpath for pie charts
    # 11. outpath for chr paintings
    # 12. dog stats
    # 13. script directory


    Rscript $script_dir/make_dog_plots.R \
        $outfile'_SM_phasefix_allchr.population_legend.txt' \
        $outfile'_SM_phasefix_allchr' \
        $outfile'_genmap.bim' \
        $query_dog_names \
        $min_seg \
        $min_pct \
        $min_hap_len \
        $min_prob \
        $outfile'_bar_fixphase/' \
        $outfile'_pie_fixphase/' \
        $outfile'_paint_fixphase/' \
        $outfile'_stats_fixphase/' \
        $script_dir

fi














