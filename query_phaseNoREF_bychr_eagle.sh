#!/bin/bash

## submit this chr to eagle
/seq/vgb/linda/bin/Eagle_v2.3.2/eagle --vcf=$targ'.vcf' \
 --geneticMapFile=/seq/vgb/linda/mutts/genetic_map_wg/gen_map_canfam3.1 \
 --chrom=$chr \
 --chromX=42 \
 --pbwtIters \
 --outPrefix=$ref'_'$chr'.phased'

