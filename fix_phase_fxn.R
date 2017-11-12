library(cowplot)
library(gridExtra)

#' This function takes in one chromosome for one individual of dog.block and
#' returns a modified set of rows for that chromosome with haplotypes fipped so
#' that each breed is haploid on the same chromosome.
flip.chr.haps <- function(
    dog.block.chr.orig, breed.pcts, make.debug.plots=F, debug.plot.prefix='',
    lookahead.n=12, max.flips=4) {
  
  # max possible flips is number of unique block-ending markers - 1
  max.flips <- min(
    max.flips, 
    length(unique(dog.block.chr.orig[, mkr.end]))-1)
  if (max.flips == 0) return(dog.block.chr.orig)
  
  # (Each incoming row is a breed haplotype on this chr for this dog.)
  # make a copy that we can edit and add cols to.
  dog.block.chr <- copy(dog.block.chr.orig)
  setkey(dog.block.chr, 'marker')
  
  # 1. Add total called length for each breed
  dog.block.chr[, 
    tot.len := sum(as.numeric(pos.len)), by=c('breed')]
  
  # 2. Add homozygosity flag for each window
  dog.block.chr[, homozyg := breed[hap==1] == breed[hap==2], by=marker]

  # 3. Calc majority hap for each breed
  dog.block.chr[, hap.maj := which.max(c(
    sum(as.numeric(pos.len[hap==1])), 
    sum(as.numeric(pos.len[hap==2])))), by=c('breed')]
  
  ##########
  # perform iterative greedy chromosome flips to improve global chrom phasing:
  
  # This iteratively picks the best marker location at which to flip the phasing
  # across the chromosome according the sum of the scores calculated in 
  # delta.scores(). At every flip it prints a debug plot of why it chose the 
  # flip it did and also prints a final plot. 
  
  # add scores across the markers for this chromsome
  dog.block.chr <- dog.block.chr[
    dog.block.chr[, delta.scores(dog.block.chr,.BY[1], breed.pcts), by=mkr.end],
    on='mkr.end']
    
  # helper function tests all single additional flips, starting with 
  #an existing flip set, which may be empty.
  test.next.flip.set <- function(starting.flip.set) {
    dog.block.chr[mkr.end != max(mkr.end), 
      delta.scores(
        dog.block.chr,
        unlist(c(starting.flip.set, .BY[1])),
        breed.pcts), 
        by=mkr.end][, list(
          mkr.end, 
          delta.total, 
          n.flips=length(unlist(starting.flip.set))+1)]
  }
  
  curr.flips <- 1
  
  # initialize best first flips and make the flip marker set
  best.flip.sets <- test.next.flip.set(c())[
    order(-delta.total, n.flips)]
  best.flip.sets <- best.flip.sets[1:min(nrow(best.flip.sets),lookahead.n)]
  flip.mkr.set.list <- sapply(best.flip.sets$mkr.end, list)
  
  # this datatable will keep a running list of the best n flip sets
  best.flip.sets <- data.table(
    mkr.end = best.flip.sets$mkr.end, 
    delta.total = best.flip.sets$delta.total,
    n.flips = best.flip.sets$n.flips,
    flip.mkr.set.list,
    flip.set.id = 1:nrow(best.flip.sets))
  
  # two stop conditions:
  # 1. we have made the maximum amount of flips
  # 2. none of the new flips we performed last round made the top N flips sets
  while (curr.flips < max.flips & curr.flips == best.flip.sets[, max(n.flips)]) {
    
    curr.flips <- curr.flips + 1
  
    # check for additional flips on top of these best n that improve the
    # score  
    next.flip.sets <- best.flip.sets[n.flips == curr.flips - 1, 
      test.next.flip.set(flip.mkr.set.list), by=flip.set.id][
         best.flip.sets[, list(flip.mkr.set.list, flip.set.id)], 
         on='flip.set.id', nomatch=0]
    
    # update the flip sets with the new flips from this round
    new.lists <- sapply(1:nrow(next.flip.sets),
      function (i) list(append(
        next.flip.sets[i, flip.mkr.set.list][[1]],
        next.flip.sets[i, mkr.end])))
    if (length(new.lists) == 1) {
      next.flip.sets$flip.mkr.set.list <- list(new.lists)
    } else {
      next.flip.sets$flip.mkr.set.list <- new.lists 
    }
    
    # merge next.flip.sets with the previous best.flip.sets and take 
    # top n again, update ids for merging next round
    best.flip.sets <- rbind(best.flip.sets, next.flip.sets)[
      order(-delta.total, n.flips), .SD[1], by=delta.total]
    best.flip.sets <- best.flip.sets[1:min(nrow(best.flip.sets),lookahead.n)][,
      flip.set.id := 1:.N]
  }
  
  # identify the best flip set from lookaheads
  best.flip.set <- best.flip.sets[1]$flip.mkr.set.list[[1]]
  
  # if this best.flip.set has a negative score, then don't do it,
  # just return.
  if (best.flip.sets[1]$delta.total < 0) {
    best.flip.set <- c()
  }
  
  # perform the flips in order, optionally printing out debugs at 
  # each step.
  flip.count <- 1
  for(flip.mkr in best.flip.set) {
    
    # save a plot with the stats
    if (make.debug.plots) {
      fn <- sprintf('%sflip.%02d.pdf', debug.plot.prefix, flip.count)
      message(fn)
      ggsave(fn, plot.chr.phase.stats(dog.block.chr, flip.mkr),
             width=8, height=6)
    }
    
    # flip the haplotypes at all markers beyond this flip.mkr
    dog.block.chr[marker > flip.mkr, hap := ifelse(hap==1,2,1)]
    dog.block.chr[, hap.maj := which.max(c(
      sum(as.numeric(pos.win.len[hap==1])), 
      sum(as.numeric(pos.win.len[hap==2])))), by=c('breed')]
    
    # update the blocks
    dog.block.chr[order(hap, marker),
      block := rleid(as.vector(breed)), by=c('hap')]
    dog.block.chr[, 
      c('mkr.start', 'mkr.end', 'pos.start', 'pos.end') := list(
          min(mkr.win.start), 
          max(mkr.win.end), 
          min(pos.win.start),
          max(pos.win.end)), by=c('hap','block')]
    
    # remove old scores
    remove.delta.cols <- names(dog.block.chr)[
      grep('delta',names(dog.block.chr))]
    set(dog.block.chr, j=remove.delta.cols,value=NULL)

    # update scores across the markers for this chromsome
    dog.block.chr <- dog.block.chr[
      dog.block.chr[, delta.scores(dog.block.chr,.BY[1], breed.pcts), 
        by=mkr.end], on='mkr.end']
    
    flip.count <- flip.count + 1
  }
  
  # save a final plot with the stats
  if (make.debug.plots) ggsave(
    sprintf('%sflip.final.pdf', debug.plot.prefix),
    plot.chr.phase.stats(dog.block.chr, best.flip.set),
    width=8, height=6)
  
  # return copy with original columns
  return(dog.block.chr[, names(dog.block.chr.orig), with=F])
}

# relative weightings of the 3 scores (black, red, green) described below.
subscore_weights <- c(1.5*1e-7,1,0.1)

# this helper function performs the haplopid chr flips at one marker
# (saltimbanqui = acrobat)
saltimbanqui <- function(dog.block.test, flip.mkr) {
  dog.block.test[marker > flip.mkr, hap.test := ifelse(hap.test==1,2,1)]
  dog.block.test[, hap.test.maj := which.max(c(
    sum(as.numeric(pos.win.len[hap.test==1])), 
    sum(as.numeric(pos.win.len[hap.test==2])))), by=c('breed')]
}

# This function calculates several distinct scores for a current phasing
# and how those scores would change if the phasing was flipped at a set
# marker.
# 1. Change in bp that are on the 'major' haploid chr for each breed
# (black in debug plot below)
# 2. Change in number of contiguous breed haplotypes created or destroyed
# by a flip
# (red in debug plot below)
# 3. Change in probability of windows that are on the major chromosome for
# a breed. This only applies to homozygous regions. (green)
# In the debug plot blue is the sum of all 3.
delta.scores <- function(dog.block.test, flip.mkr.set, breed.pcts) {
  
  # if duplicate flips in set ('flip-backs') return an empty DT.
  if (length(unique(flip.mkr.set)) < length(flip.mkr.set)) {
    return(list())
  }
  
  # remove last marker if it exists in the list of flips
  flip.mkr.set <- setdiff(
    flip.mkr.set, dog.block.test[, max(mkr.end)])
  
  # if no flips in set after removal, return an empty DT
  if (length(flip.mkr.set) == 0) {
    return(list(
      delta.maj=0,
      delta.haps.wt=0,
      delta.prob=0,
      delta.total=0))
  }
  
  if ('hap.test' %in% names(dog.block.test)) {
    dog.block.test$hap.test <- NULL
    dog.block.test$hap.test.maj <- NULL
    dog.block.test$block.test <- NULL
  }
  
  dog.block.test[, hap.test := hap]
  
  # make every flip in the flip set
  for (flip.mkr in flip.mkr.set) {
    dog.block.test <- saltimbanqui(dog.block.test, flip.mkr)
  }
  
  # calculate change in bps that are on the major haploid chrom
  orig_score <- dog.block.test[hap != hap.maj][
    breed.pcts, on=c(breed='chr.breed'), nomatch=0,
    sum(pos.win.len*chr.pct)]
  new_score <- dog.block.test[hap.test != hap.test.maj][
    breed.pcts, on=c(breed='chr.breed'), nomatch=0,
    sum(pos.win.len*chr.pct)]
  # change in score, positive is better:
  delta.maj <- subscore_weights[1] * (orig_score - new_score)
  
  # calculate change in number of distinct breed haplotypes
  orig_haps <- sum(dog.block.test[, length(unique(block)), by='hap']$V1)
  
  dog.block.test[order(hap.test, marker),
    block.test := rleid(as.vector(breed)), by=c('hap.test')]
  new_haps <- sum(
    dog.block.test[, length(unique(block.test)), by='hap.test']$V1)
  
  # weight new haplotype breakpoints by the breed it breaks, bigger pct of
  # chrom means worse to create new breakpoint
  delta.haps.wt.neg <- unique(dog.block.test[, list(
    orig.hap=paste(hap, block, sep='.'),
    new.hap=paste(hap.test, block.test,sep='.'), 
    breed)])[, list(breed=breed[1], freq=.N-1), by=orig.hap][freq > 0][
      breed.pcts, on=c(breed='chr.breed'), 
      breed.hap.score := chr.pct * freq][, sum(breed.hap.score)]
  
  delta.haps.wt.pos <- unique(dog.block.test[, list(
    orig.hap=paste(hap, block, sep='.'),
    new.hap=paste(hap.test, block.test,sep='.'), 
    breed)])[, list(breed=breed[1], freq=.N-1), by=new.hap][freq > 0][
      breed.pcts, on=c(breed='chr.breed'), 
      breed.hap.score := chr.pct * freq][, sum(breed.hap.score)]

  # if homozygous, weight high-posterior windows to go on the major chromosome
  orig.prob.major.hap <- sum(
    dog.block.test[homozyg==T, 
      prob[hap.maj==hap]-prob[hap.maj!=hap], 
      by=marker]$V1)
  new.prob.major.hap <- sum(
    dog.block.test[homozyg==T, 
      prob[hap.test.maj==hap.test]-prob[hap.test.maj!=hap.test], 
      by=marker]$V1)
  # positive means flipped is better probs, neg means flipped is worse
  delta.prob <- subscore_weights[3] * (new.prob.major.hap-orig.prob.major.hap)
      
  # also save unweighted
  delta.haps <- orig_haps - new_haps
  delta.haps.wt <- subscore_weights[2] * (delta.haps.wt.pos - delta.haps.wt.neg)
  delta.total <- sum(
    subscore_weights * c(delta.haps.wt, delta.maj, delta.prob))
  
  # die if we get an NA
  stopifnot(!any(is.na(delta.total)))

  dog.block.test$hap.test <- NULL
  dog.block.test$hap.test.maj <- NULL
  dog.block.test$block.test <- NULL
  return(list(
    delta.maj=delta.maj,
    delta.haps.wt=delta.haps.wt,
    delta.prob=delta.prob,
    delta.total=delta.total))
}

# This is a debug plot for where to flip the chromosomal phasing and the
# associated scores. For info about the scores/colors, see delta.scores()
# function comment above.
plot.chr.phase.stats <- function(dog.block.chr, flip.mkr.set=c()) {
  
  flip.pos.set <- unique(
    dog.block.chr[
      mkr.end %in% flip.mkr.set, 
      pos.end])
  
  delta.cols <- names(dog.block.chr)[grep('delta',names(dog.block.chr))]
  
  haps.dt <- unique(dog.block.chr[, list(
      hap, hap.maj, homozyg, block, pos.start,
      pos.end, breed, mkr.end)])[,
        variable := 'Haplotype Calls']
  
  score.dt <- melt(
    unique(dog.block.chr[, c('pos.end', delta.cols), with=F]), 
    measure.vars=delta.cols)
  score.dt[, variable := factor(variable, 
    levels=levels(variable),
    labels=c('Major','Merge/Split','Posterior','Total'))]
  
  score.colors <- c('black','red','green','blue')
  
  phase.plot <- ggplot() + 
    geom_rect(data=haps.dt, aes(
      xmin=pos.start, xmax=pos.end-1, 
      ymin=hap-0.45, ymax=hap+0.45, 
      fill=breed,
      color=(hap != hap.maj & homozyg == F))) +
    geom_vline(xintercept=dog.block.chr[,
      pos.end[which.max(delta.total)]],
      linetype=3, color='gray50') +
    geom_vline(xintercept=flip.pos.set,linetype=2) +
    scale_color_manual(
      values=c('white','black'), 
      name='Haplotype\non Major Chr', guide=FALSE) +
    scale_fill_discrete(name='Breed') +
    facet_grid(variable~., scales='free_y', as.table=F) +
    scale_y_continuous(breaks=c(1,2), labels=c('1','2')) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill=NA),
      legend.box.just = "top",
      legend.justification = c(0, 1)
    )
  
  score.plot <- ggplot() + 
    geom_hline(yintercept=c(0)) +
    geom_vline(xintercept=dog.block.chr[,
      pos.end[which.max(delta.total)]],
      linetype=3, color='gray50') +
    geom_vline(xintercept=flip.pos.set,linetype=2) +
    geom_line(data=score.dt, aes(x=pos.end, y=value, color=variable)) +
    scale_color_manual(values=score.colors, name='Scores') +
    scale_x_continuous(, 
      limit=c(
        haps.dt[, min(pos.start)],
        score.dt[, max(score.dt$pos.end)])) +
    facet_grid(variable~., as.table=T) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill=NA)
    ) +
    xlab('Score') +
    ylab('Chr. Position')
  
  #http://stackoverflow.com/questions/36198451
  # /specify-widths-and-heights-of-plots-with-grid-arrange
  combined.plot <- plot_grid(
    phase.plot, score.plot, align = "v", nrow = 2, 
    rel_heights = c(1/3,2/3))

  return(combined.plot)
}

# helper function to flip a marker window in the tped based on:
# chr, pos.start, pos.end, and tped.col (for hap 1)
# (this is called at the very end to update the tped in place)
flip.tped.win <- function(mkr.win.start, mkr.win.end, tped.col) {

  new.hap.2 <- unlist(tped.phasefixed[
    mkr.win.start:mkr.win.end, paste0('V',tped.col), with=F])
  new.hap.1 <- unlist(tped.phasefixed[
    mkr.win.start:mkr.win.end, paste0('V',tped.col+1), with=F])
  
  tped.phasefixed[mkr.win.start:mkr.win.end,
    paste0('V',c(tped.col,tped.col+1)) := list(new.hap.1, new.hap.2)]
  return(NULL)
}