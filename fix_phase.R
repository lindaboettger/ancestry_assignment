### 0. get arguments #######################################

args = commandArgs(trailingOnly=TRUE)

### 1. Load packages and functions  ########################

require(data.table)
require(ggplot2)
require(ggrepel)
require(base)
require(scales)
require(base)

# get directory name
# Old way, doesnt always work: script.dir <- dirname(sys.frame(1)$ofile)
# New way: 
#  https://stackoverflow.com/questions/1815606/
#   rscript-determine-path-of-the-executing-script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
script.dir <- dirname(thisFile())

# load recombination functions from scripts dir
source(file.path(script.dir, 'recomb_fxn.R'))

# load dog phase functions from scripts dir
source(file.path(script.dir, 'fix_phase_fxn.R'))

# load dog plot functions from scripts dir
source(file.path(script.dir, 'dog_plot_fxn.R'))


### 2. Load Darwin's Dogs Data ########################

# load or re-run fix phasing
load_old <- F #<- set to T if we want to use old fixphase data

breed.leg <- data.table(read.table(args[1],header=T))
        
# get tped, tfam, probs, and marker files
out.tped.file <- paste(args[2],'.tped',sep='')
out.tfam.file <- paste(args[2],'.tfam',sep='')
out.probs.file <- paste(args[2],'.Probs.tped',sep='')

# get markers file (same snps as before)
markers.file <- args[3]
markers <- data.table(read.table(markers.file, header=F)[,c(2,1,3:6)])

# add marker numbers
markers[, marker := 1:nrow(markers)]

# add column names
names(markers) <- c('snp', 'chr','morgans','pos', 'allele1', 'allele2', 'marker')

# convert cM to M
markers[, morgans := morgans/100]
setkey(markers, marker)

# get chr ends from markers
chr_ends <- markers[, list(marker=max(marker)), by=chr]
setkey(chr_ends, marker)

# if load_old is false, perform phase flipping on all dogs, else load existing
if (!load_old) {
  
  # Load the tped.probs for this dataset.
  tped.probs.dt <- load_supportmix_output(
    out.tped.file, out.tfam.file, out.probs.file, markers, breed.leg)
  tped.probs.dt[, hap.orig := hap]

  # Make into breed haplotype blocks
  dog.block <- make.dog.block(tped.probs.dt)
  
  # helper table of percentages of breeds in this indiv
  bp.per.dog <- dog.block[indiv==indiv[1], sum(as.numeric(pos.win.len))]
  breed.pcts <- dog.block[, list(
    pct=sum(as.numeric(pos.win.len))/bp.per.dog),
    by=c('indiv','breed')]
  names(breed.pcts) <- c('indiv','chr.breed','chr.pct')

  dog.block.phasefixed <- dog.block[order(chr),
    flip.chr.haps(
      .SD,
      breed.pcts[indiv==.BY[[1]]],
      make.debug.plots=F,
      debug.plot.prefix=paste0(
        '/Users/boettger/postdoc/mutts/figs/fixphase_may_debugplots/hubble.chr',
        sprintf('%02d',as.numeric(.BY[2])),'.'),
      lookahead.n=12,
      max.flips=4),
    by=c('indiv','chr')]
  write.table(dog.block.phasefixed,
    paste(args[5],'_oldhap.txt',sep=''),
      quote=F, row.names=F)
} else {
  dog.block.phasefixed <- data.table(
    read.table(dog.block.phasefixed,
    paste(args[5],'_oldhap.txt',sep=''), header=T))
}

######## load actual tped with bases #########
tped.orig <- data.table(read.table(paste(args[4],'.tped', sep='')))
tfam.orig <-  data.table(read.table(paste(args[4],'.tfam', sep='')))

# merge dog.block.phasefixed with tfam row #s to get tped col numbers
# so that we can flip columns corresponding to flipped haplotypes
dog.block.phasefixed <- dog.block.phasefixed[
  tfam.orig[, list(indiv=V1, tped.col=4+(.I*2)-1)],
  on='indiv']

# copy of original tped to flip, used by flip.tped.win fxn below
tped.phasefixed <- copy(tped.orig)

# list of window sets to be flipped in the new tped
markers.to.flip <- unique(dog.block.phasefixed[
  hap != hap.orig,
  list(mkr.win.start, mkr.win.end, tped.col), 
  by=c('indiv','marker')][order(marker)])

# perform the flips (ignore output)
markers.to.flip[, flip.tped.win(
 .BY[[2]], mkr.win.end, tped.col),
  by = c('indiv','mkr.win.start')]

write.table(tped.phasefixed, paste(args[5], '.tped', sep=''),
    quote=F, row.names=F, col.names=F)

