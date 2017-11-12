library(RCurl)
library(data.table)
library(zoo)
library(plyr)
library(reshape2)

############################## FUNCTIONS #################################

# load and melt tped/tfam ##############
# Inputs:
#   tped file & tfam file paths
# 
# Outputs:
#   a melted tped/tfam representation with columns for individual name,
#   snp/chr location, haploid genome (pops are phased), and allele
load_melt_tped <- function(tped.file, tfam.file, markers) {

  # load tped & tfam
  tped.dt <- data.table(read.table(tped.file))
  tfam.dt <- data.table(read.table(tfam.file))
  
  # remove cols 2 & 3 and melt by indivs
  # wtf this removes all but the SNP col and would melt on the first base???
  tped.dt <- tped.dt[,-c(1,3,4), with=F]
  tped.dt <- melt(tped.dt, id.vars=c('V2'))
  setnames(tped.dt, names(tped.dt), c('snp','indiv','allele'))
  
  # get number of snps
  n.snps <- length(tped.dt[, levels(snp)])
  
  # get indiv names from tfam and number of indivs
  indiv.names <- tfam.dt$V2
  n.indivs <- length(indiv.names)
  
  #check that tfam and tped lengths match
  stopifnot(length(indiv.names)*2*n.snps == nrow(tped.dt))
  
  # replace the blank indiv names Vx - Vy with the individuals
  tped.dt[, indiv := rep(indiv.names, each=2*n.snps, times=1)]
  
  # add haplotype genome number (1 or 2)
  tped.dt[, hap := rep(c(1,2), each=n.snps, times=n.indivs)]
  
  # add marker numbers
  tped.dt[, marker := rep(1:n.snps, times=n.indivs)]
  
  # add keys for indexing
  setkey(tped.dt, 'marker')
  
  # merge with marker data
  tped.dt <- tped.dt[markers[, list(chr,marker,morgans,pos)]]
  
  setkey(tped.dt, 'indiv','marker','hap')
  
  return(tped.dt)
}

# load and melt output tped/tfam/probs ##############
# Inputs:
#   tped/probs file & tfam file paths that were outputs from supportmix
# 
# Outputs:
#   a melted tped/tfam representation with columns for individual name,
#   snp/chr location, haploid genome (pops are phased), and allele
load_melt_output_tped <- function(out.tped.file, out.tfam.file, markers) {

  # load tped & tfam
  tped.out.dt <- data.table(read.table(out.tped.file))
  tfam.out.dt <- data.table(read.table(out.tfam.file))
  
  # remove cols 2 & 3 and melt by indivs
  # wtf this removes all but the SNP col and would melt on the first base???
  tped.out.dt <- tped.out.dt[,-c(1,3,4), with=F]
  tped.out.dt <- melt(tped.out.dt, id.vars=c('V2'))
  setnames(tped.out.dt, names(tped.out.dt), c('snp','indiv','allele'))
  
  # get number of snps
  n.snps <- length(tped.out.dt[, levels(snp)])
  
  # get indiv names from tfam and number of indivs
  indiv.names <- tfam.out.dt$V2
  n.indivs <- length(indiv.names)
  
  #check that tfam and tped lengths match
  stopifnot(length(indiv.names)*2*n.snps == nrow(tped.out.dt))
  
  # replace the blank indiv names Vx - Vy with the individuals
  tped.out.dt[, indiv := rep(indiv.names, each=2*n.snps, times=1)]
  
  # add haplotype genome number (1 or 2)
  tped.out.dt[, hap := rep(c(1,2), each=n.snps, times=n.indivs)]
  
  # add marker numbers
  #tped.out.dt[, marker := rep(1:n.snps, times=n.indivs)]
  
  # add keys for merging
  setkey(tped.out.dt, 'snp')
  setkey(markers, 'snp')
  
  # merge with marker data to get numbers
  tped.out.dt <- na.omit(tped.out.dt[markers,list(chr,snp,indiv,allele,hap,marker)])
  
  setkey(tped.out.dt, 'indiv','marker','hap')
  
  return(tped.out.dt)
}

# load and melt output tped/tfam/probs ##############
# Inputs:
#   tped/probs file & tfam file paths that were outputs from supportmix
#   markers file
# 
# Outputs:
#   a melted tped/tfam/probs representation with columns for individual name,
#   snp/chr location, haploid genome (pops are phased), and allele
load_supportmix_output <- function(
  out.tped.file, out.tfam.file, out.probs.file, markers, breed.leg) {

  #load output tped.file, tfam.file and probabilities file
  out.tped.dt <- load_melt_output_tped(out.tped.file, out.tfam.file, markers)
  #change prob colname
  setnames(out.tped.dt, 'allele', 'popnum')

  #load probs file
  out.probs.dt <- load_melt_output_tped(out.probs.file, out.tfam.file, markers)
  #change prob colname
  setnames(out.probs.dt, 'allele', 'prob')

  #merge output tped and probs
  tped.probs.dt <- out.tped.dt[
    out.probs.dt,
    list(chr,snp,indiv,popnum,hap,marker,prob)]

  #replace numbers with populations and replace indiv factors with characters
  tped.probs.dt <- merge(tped.probs.dt, breed.leg, by='popnum')[,-'popnum', with=F]
  tped.probs.dt[, indiv := as.character(indiv)]
  setkey(tped.probs.dt, 'indiv','marker','hap')
  
  # split test individual pairs I001_I002 to IN-1_INNN into haploid individuals
  # with hap == NA.
  tped.probs.dt[
    J(grep('^I\\d+_I\\d+',unique(tped.probs.dt$indiv), value=T)),
    c('indiv','hap') := list(
      rep(unlist(strsplit(indiv[1],'_')), .N/2),
      0),
    by= indiv]
  
  # make the window calls be labelled by their LAST marker instead of their
  # first marker, to match the haplotype truth calls. The last window
  # is labelled by the end of the chr.
  # marker[2:.N]-1 : this grabs the marker which starts the next window
  # and subtracts one from it. 
  # .BY[3] : .BY is a vector of the grouped-by values, the third is chr.
  # this lets us grab the current chr we are on and get its last marker.
  tped.probs.dt[, 
    marker := as.integer(
      c(marker[2:.N]-1, chr_ends[chr==.BY[3], marker])),
    by=c('indiv','hap','chr')]
  
  setkey(tped.probs.dt, 'indiv','marker')
  
  # do this because assignment fuction (:=) hides the return() otherwise
  tped.probs.dt[]

  return(tped.probs.dt)
}






# get breakpoint positions (recomb + chr ends) ##############
# Inputs:
#   i: generation number, for bookkeeping purposes, does not change recomb pos
#   pop: two populations as 'strings' which we are recombining.
# 
# Outputs:
#   a data.frame with markers that mark breakpoint positions (either
#   recombination sites or chromosome ends). Breakpoint positions are
#   associated with populations (pop) which correspond to the identity
#   of the preceeding haplotype. The gen column is arbitrary and used
#   to keep track of bookkeeping elsewhere.
do_rcmb <- function (pop) {
  
  # ensure markers are ordered
  stopifnot(markers[,!is.unsorted(marker)])
  
  # no more than two pops
  stopifnot(length(pop) <= 2)
  
  # if only one population, just return chr ends.
  if(length(pop) == 1) {
    return(cbind(chr_ends, parent=pop))
  } else {
    # if two populations, return recombination sites and chr ends.
    return(rbind(
      markers[, 
        list(
          'marker'=marker[which(
            (runif(.N-1) < (1-exp(-diff(morgans)))))]),
        by=chr],
        chr_ends)[order(chr, marker), 
      list(
        'marker'=marker, 
        'parent'=rep_len(sample(pop), .N)),
        by=chr])
  }
}


# Recombine the markers for two individuals, a and b. ############
# Inputs:
#   a: recombination marker set (e.g. output from do_rcmb()) for 
#      population/individual a
#   b: same for population b
# 
# Outputs:
#   a recombined market set from a and b
rcmb_markers <- function (a, b) {
  
  # combine markers for both parents
  ab <- rbind(a,b)
  ab$id <- c(rep('a',nrow(a)), rep('b',nrow(b)))
  names(ab) <- paste0('ab.',names(ab))
  
  # generate recombination between individuals a & b
  # mark these recombination breakpoints as being from 
  # generation i
  rcmb <- do_rcmb(pop=c('a','b'))
  # keep track of the previous markers
  rcmb$prev_marker <- c(0,rcmb$marker[1:nrow(rcmb)-1])
  
  rcmb <- rcmb[, get_ab_breaks(
      ab=ab, 
      chr=chr, 
      this_marker=marker,
      parent=parent,
      prev_marker=prev_marker), by=c('chr','marker','parent')][,
        list(
          'chr'=ab.chr, 
          'marker'=ab.marker,
          'parent'=ab.parent)]
  
  #simplify adjacent breaks from same pop and return
  return(simplify_breaks(rcmb))
}

# Simplify breaks for a recombination ############
# removes breaks if adjacent to another break from 
# same population/indiv
simplify_breaks <- function(rcmb) {
  rcmb[
    parent != c(NA,parent[1:nrow(rcmb)-1]) | 
    marker %in% chr_ends$marker,][order(chr, marker)] 
}

# Get breaks for a recombination ############
# Inputs:
#   ab: joined a and b recomb markers
#   chr: chrom 
#   this_marker: marker to recombine at
#   parent: population to recombine from (either a or b)
#   prev_marker: marker before this_marker
# 
# Outputs:
#   markers that go between this marker and the 
#   last marker from the parent
get_ab_breaks <- function(
  ab, chr, this_marker, parent, prev_marker) {
  
  parent_markers <- ab[
    ab.chr == chr & 
      ab.id == parent & 
      ab.marker < this_marker &
      ab.marker > prev_marker,
    list(ab.chr, ab.marker, ab.parent)]
  
  # get genotype of this ab break
  final_marker_parent <- ab[
      ab.chr == chr & 
      ab.id == parent & 
      ab.marker >= this_marker]$ab.parent[1]
  
  # return old breaks and this new break,
  return(rbind(
    parent_markers,
    list(chr, this_marker, final_marker_parent)))
}



# old: crosses.dt <- sample_for_crosses(minor.breed, indiv.list, num.generations)
# now, instead of sample_for_crosses, we want: sample one dog from every breed (57)
# NEW: sampled.indivs.dt <- sample_for_indivs(indiv.list)
#      (sampled.indivs.dt will contain 57 individuals, one per breed)
#      sample_for_indivs() <- could even be a call to sample_for_crosses() with 
#      num.generations == 57. Both 
# Sample one individual from each breed ############
#
# Inputs:
# list of all individuals and breeds
# 
# Outputs:
#   one individual per breed 
#   last marker from the parent
sample_for_indivs <- function(indiv.list) {
  
  # get 1 individual per breed
  sampled.indivs.dt <- indiv.list[, sample(indiv, size=1), by=breed]
  colnames(sampled.indivs.dt)[2] <- "indiv"

  return(sampled.indivs.dt)
}

# Sample individuals for one admixed individual ############
#
# Inputs:
#   sampled.indivs.dt : list of individuals sampled in sample_for_indivs()
#   num.generations : number of successive admixtures for this individual
#
# Outputs:
#
#   crosses.dt, which has an ordered list of admixture individuals,
#   of length num.generations + 1, each with cross # and expected
#   percentage ancestry of the final admixed haploidgenome
sample_for_crosses <- function(sampled.indivs.dt, num.generations) {

  crosses.dt <- sampled.indivs.dt[sample(.N,size=num.generations+1)]
  crosses.dt[, cross := c(-(num.generations), (-(num.generations):-1))]
  crosses.dt[, exp.pct := 2^cross]
  
  # needed due to bug: http://stackoverflow.com/questions/32988099
  crosses.dt[] 
  return(crosses.dt)
}

# Do recombination ############
# Inputs:
#   crosses.dt : a data.table from sample_for_crosses()
# Output:
#    a marker set for a recombined individual after the
#    specified crosses. The last individual in crosses.dt
#    is not included in the cross, she gets her own haploid
#    genome adjacent to this one.
do_n_rcmb <- function(crosses.dt) {
  
  first.indiv.markers <- do_rcmb(crosses.dt[1, indiv])
  
  # don't use last row of crosses, she has her own
  # haploid genome
  remaining.indiv.markers.list <- lapply(
    crosses.dt[2:(nrow(crosses.dt)-1), indiv],
    do_rcmb)
  
  # how reduce works:
  # c <- rcmb_markers(init,x1)
  # c <- rcmb_markers(c,x2)
  # c <- rcmb_markers(c,x3)
  # c <- rcmb_markers(c,x4)
  # c <- rcmb_markers(c,x5)
  recomb.indiv <- Reduce(
    f= rcmb_markers, 
    x= remaining.indiv.markers.list,
    init= first.indiv.markers)
  
  return(recomb.indiv)
}

# Map alleles to recombined individual from a tped ############
# Inputs:
#   recombination marker sets across many indiv.i individuals
# Outputs:
#  (none, updates input dt)
assign_random_haps <- function(all.haps.dt) {
  all.haps.dt[ , hap := sample(1:2, 1), keyby=c('indiv.i', 'chr', 'parent')]
}

# Map alleles to recombined individual from a tped ############
# Inputs:
#   a recombined marker set from an individual, via do_n_rcmb() or similar
# Outputs:
#  an haploid allele/marker table, with columns:
#  chr, marker, allele, morgans, indiv, haploid genome #, and snp pos
#
# (note) recomb.indiv includes a col called indiv.i, so that
# this function can take in multiple individuals simultaneously
get_hap_markers <- function(recomb.indiv, tped.dt) {
  setkey(recomb.indiv, indiv.i, marker)
  # make full size genome using original markers file (fill in ancestry of empty markers)
  # merge against all markers and do NOCB (next obvervation carried backward)
  recomb.indiv.sample.hap <- recomb.indiv[, 
    # individually do a rolling join the .SD subsets of the indiv.i 'by' grouping
    .SD[markers[, list(marker)], on='marker', roll=-Inf],
    by=indiv.i]
  
  # remove chr column, they will be added back from the tped
  recomb.indiv.sample.hap[, 'chr' := NULL, with=F]
  # rename parent to indiv, so we can join with tped.dt 
  setnames(recomb.indiv.sample.hap, 'parent','indiv')
  setkey(recomb.indiv.sample.hap, indiv, marker, hap)
  
  return(tped.dt[
    # pull out the alleles for the indiv/marker/hap combo in tped using
    # supa-fast key, and ditch all cols except marker, allele, and indiv.i
    recomb.indiv.sample.hap,
    c('marker','allele','indiv.i'),
    with=F][
      # reorder by indiv.i and marker, remove other columns
      order(indiv.i, marker), c('allele','indiv.i'), with=F]
  )
}

# Map alleles to recombined individual from a tped ############
# Inputs:
#   all.haps.dt: recombined marker sets from many individuals,
#    grouped by indiv.i, the individual ID
#   all.crosses.dt: matching list of crosses from many individuals
#   tped.dt: input tped file, melted
# Outputs:
#  write a an haploid allele/marker tped table, with columns:
#    chr, marker, allele, morgans, (... indiv allele columsns ...)
alleles_to_tped <- function(all.haps.dt, all.crosses.dt, tped.dt, outfile.prefix) {
  
  # get the haploid markers for the recombined chromosome
  indiv.allele.dt <- get_hap_markers(all.haps.dt, tped.dt)  
  
  # number of admixed indivs in all.crosses.dt
  num.admixed.indiv <- all.crosses.dt[, max(indiv.i)]
  
  # check that each individual has the same number of markers
  # so we can do the matrix() resize safely below
  stopifnot(length(unique(indiv.allele.dt[, .N, by=indiv.i]$N)) == 1)
  # also check that indiv.i's go from 1 to max
  stopifnot(length(unique(all.haps.dt$indiv.i))==num.admixed.indiv)
  
  # grab sampled purebred individuals from all.crosses.dt
  sampled.indiv.list <- unique(all.crosses.dt$indiv)
  
  # make maker columns and allele columns for the purebreds sampled from
  # (to serve as our controls)
  tped.out.dt <- dcast.data.table(tped.dt[J(sampled.indiv.list)],
        formula=chr+pos+morgans+snp~indiv+hap,
        value.var='allele')
  
  # swap column order of snp and pos to appease tped gods
  setcolorder(tped.out.dt, 
    names(tped.out.dt)[c(1,4,3,2,5:ncol(tped.out.dt))])
  
  purebred.indiv.names <- unique(
    gsub('_[12]$','',names(tped.out.dt)[5:ncol(tped.out.dt)]))
  
  # add new columns for each admixed indiv named I1 to I'n'
  tped.out.dt[, 
    paste0('I',sprintf('%03d',1:num.admixed.indiv)) := 
      # do a fake 'dcast' using a matrix() resize on the admixed allele list
      data.table(indiv.allele.dt[, matrix(allele,
        nrow=.N/num.admixed.indiv,
        ncol=num.admixed.indiv)])]
  
  admixed.indiv.names <- paste0(
    'I',
    sprintf('%03d',seq(1,num.admixed.indiv,2)),
    '_I',
    sprintf('%03d',seq(2,num.admixed.indiv,2)))
  
  # supportmix uses alphanumeric chr sorting, so append 0s
  tped.out.dt[, chr.i := sprintf('chr%02d',chr)]
  tped.out.dt[, chr := NULL]
  setnames(tped.out.dt, 'chr.i','chr')
  setcolorder(tped.out.dt, 
    names(tped.out.dt)[
      c(ncol(tped.out.dt),1:(ncol(tped.out.dt)-1))])

  # write to a tped file (eventually in tmp dir)
  fwrite(
    tped.out.dt, 
    paste0(outfile.prefix, '.tped'),
    col.names=F, quote=F, sep='\t')
  
  # write to a tfam file (eventually in tmp dir)
  write_tfam(
    indiv.names=c(purebred.indiv.names, admixed.indiv.names),
    outfile=paste0(outfile.prefix, ".tfam"))

}

# Create reference pop files for breeds used (without recomb indvs) ############
# Inputs:
#   corresponding breeds for each indiv (indiv.list)
#   data for all individuals (tped.dt)
# Outputs:
#  reference population files for each breed we used 
`%ni%` <- Negate(`%in%`) 
## note this hasn't been tested... it takes too long!!!
write_new_ref_tped_tfams <- function(
    tped.dt, 
    sampled.indivs.dt, 
    indiv.list, 
    outfile.prefix) {

  # get all indivs from breeds we will output
  #mybreeds.indivs <- indiv.list[breed %in% crosses.dt$breed]
  
  # remove individuals that we used in the crosses
  setkey(indiv.list, indiv)
  ref.indiv.dt <- indiv.list[indiv %ni% sampled.indivs.dt$indiv]
  setkey(indiv.list, 'indiv','breed')
    
  # helper function to print a tped for one breed
  write_tped_breed <- function(breed.i, tped.subset.dt) {
    print(paste("Writing ref tped for breed",breed.i))
      
    # must cast with chr,pos rather than chr,snp b/c snps will be 
    # sorted alphanumerically. Then change order back to
    # chr,snp with setcolorder() function below
    tped.subset.cast.dt <- dcast.data.table(
        tped.subset.dt,
        formula=chr+pos+morgans+snp~indiv+hap,
        value.var='allele')

    # swap column order of snp and pos for tped format
    # (see comment above)
    setcolorder(tped.subset.cast.dt, 
      names(tped.subset.cast.dt)[
        c(1,4,3,2,5:ncol(tped.subset.cast.dt))])
    
    # supportmix uses alphanumeric chr sorting, so append 0s
    tped.subset.cast.dt[, chr.i := sprintf('chr%02d',chr)]
    tped.subset.cast.dt[, chr := NULL]
    setnames(tped.subset.cast.dt, 'chr.i','chr')
    setcolorder(tped.subset.cast.dt, 
      names(tped.subset.cast.dt)[
        c(ncol(tped.subset.cast.dt),1:(ncol(tped.subset.cast.dt)-1))])
    
    fwrite(
      x= tped.subset.cast.dt,
      file= paste0(outfile.prefix, '.', breed.i,".tped"), 
      col.names = F, sep='\t')

    write_tfam(
      indiv.names=tped.subset.dt[, unique(indiv)],
      outfile=paste0(outfile.prefix, '.', breed.i,".tfam"))
  }
  
  # run the helper function by each breed in for loop
  for (breed.i in unique(sampled.indivs.dt$breed)) {
    write_tped_breed(
      breed= breed.i, 
      tped.subset.dt= tped.dt[ref.indiv.dt[breed==breed.i]]
    )
  }
}

# get real majority indiv for each window and % of window that is that indiv
# compare to supportmix call
merge_truth_supportmix <- function(all.haps.dt, tped.probs.dt, indiv.list) {
  # NOTE: we ran this to test: 
  # all.haps.dt.orig <- all.haps.dt
  # all.haps.dt <- all.haps.dt[chr==9]
  
  # this will be a list of haplotypes, labelled by their last marker,
  # with the breed/truth column being the true assignment for that
  # haplotype
  truth.mutts.dt <- indiv.list[all.haps.dt, on=c(indiv="parent")][,
    indiv := sprintf('I%03d',indiv.i)][,
    c('indiv.i','hap') := list(NULL,NULL)][]
  setnames(truth.mutts.dt,'breed','truth')
  setkey(truth.mutts.dt, 'indiv','marker')
  
  #merge truth haplotypes with supportmix prediction windows
  truth.windows.dt <- truth.mutts.dt[,-'chr', with=F][
    tped.probs.dt, roll=-Inf][, snp := NULL]
    
  # mark the end of true haplotypes of mutts(for bookkeeping)
  truth.windows.dt[, hap_end := F]
  
  # add haplotype rows back to truth.window rows
  truth.windows.dt <- rbind(
    truth.windows.dt, 
    # join marker chr:pos back to truth.mutts.dt,
    # (but remove chromosome ends from haplotype list)
    truth.mutts.dt[marker %ni% chr_ends$marker,
      list(
        indiv,
        truth,
        marker, 
        chr,
        hap=0,
        prob=NA,
        breed_guess=NA,
        hap_end=T)])
  
  # sort the whole table again to mix in the recombination breakpoints
  setkey(truth.windows.dt, indiv, hap, chr, marker)
  
  # fill in the truth column for non-admixed control indivs by
  # mapping their breed from the indiv.list
  truth.windows.dt[hap > 0, truth := indiv.list[indiv, breed]]
  
  # map the probs and breed guesses to the recombination breakpoints
  # using na.locf
  truth.windows.dt[,
    c('breed_guess','prob') := list(
      na.locf(breed_guess, fromLast=T),
      na.locf(prob, fromLast=T)),
    by=c('indiv','chr','hap')]
  
  # similarly, add chr ends to hap_end==T for all
  truth.windows.dt[chr_ends, hap_end := T, on=c('chr','marker'), nomatch=0L]
  
  # flag windows which cross a recombination event
  truth.windows.dt[, trans_window := F]
  truth.windows.dt[, 
    trans_window := c(F, hap_end[1:(.N-1)]), 
    by=c('indiv','chr','hap')]
  
  truth.windows.dt

  return(truth.windows.dt)
}

  
# write tfam is a tiny helper function used for making
# any tfam we need, it's called by write_tped_breed,
# (fxn to write recomb tped), etc
write_tfam <- function (indiv.names, outfile) {
  fwrite(
    data.table(
      indiv.names,
      indiv.names,
      0,
      0,
      0,
      0),
    file=outfile,
    col.names = F, sep='\t')
}


# determine breed for each individual 
load_indiv_to_breed <- function(breed_file, indiv_file) {

  # load files
  breed.dt <- data.table(read.table(breed.file, header=T))
  indiv.dt <- data.table(read.table(indiv.file, header=F))
  colnames(indiv.dt) <- 'IID'
  
  # merge breeds and dogs
  return(breed.dt[indiv.dt, on='IID'][, list(indiv=FID, breed=BREED)])
}

