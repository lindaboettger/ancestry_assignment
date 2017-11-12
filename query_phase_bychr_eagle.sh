#!/bin/bash

## submit this chr to eagle
/seq/vgb/linda/bin/Eagle_v2.3.2/eagle --vcf=$targ'.vcf.gz' \
 --geneticMapFile=/seq/vgb/linda/mutts/genetic_map_wg/gen_map_canfam3.1 \
 --chrom=$chr \
 --chromX=42 \
 --pbwtIters 10 \
 --expectIBDcM 20 \
 --outPrefix=$targ'_'$chr'.phased'

#could add --genoErrProb arg (=0.003)
#--expectIBDcM arg (=2.0)         expected length of haplotype copying (cM)

#2 mb is the expected length

#move thi s up to 20?