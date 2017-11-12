### R Functions to make per-dog plot pdfs.
library(ggplot2)

#set %ni% operator
`%ni%` <- Negate(`%in%`)

make_dog_plots <- function(
    this.dog, this.name, tped.probs.dt, 
    min.hap.len, min.seg, min.pct, min.prob) {
  
  # First we filter on haplotype length (min.hap.len)
  # Then we filter on either posterior win probability (min.prob) or the mean
  # posterior probability for the whole haplotype
  # Then we filter on breeds with at least min.seg segments and min.pct percentage
  # of the genome called.
  
  message(paste('Making plots for', this.dog,':', this.name))
  
  dog.block <- make.dog.block(tped.probs.dt)
  
  # rename breed names to remove underscores
  dog.block[, breed := factor(
    breed,
    levels=levels(breed),
    labels=gsub('_', ' ', as.character(levels(breed))))]

  # generate a per-haplotype version of dog.block called dog.hap
  dog.hap <- dog.block[, list(
      mkr.start=mkr.start[1], mkr.end=mkr.end[1], 
      prob=mean(prob), breed=breed[1], num.win=num.win[1],
      pos.start=pos.start[1], pos.end=pos.end[1], pos.len=pos.len[1]), 
    by=c('indiv', 'block','hap','chr')]
  
  seg.len.total <- dog.hap[indiv==this.dog, sum(as.numeric(pos.len))]
  
  # first, mark small segments with min.hap.len
  dog.hap <- dog.hap[, filter.out := pos.len < min.hap.len]
  dog.block <- dog.block[, filter.out := pos.len < min.hap.len]
  
  # do a per-haplotype cutoff of mean window posterior calls
  dog.hap <- dog.hap[prob < min.prob, filter.out := T]
  dog.block <- dog.block[prob < min.prob, filter.out := T]
  
  # use dog.hap to get breed.hapsize summary info
  this.dog.breed.hapsize.dt <- breed.breakdown.by.hapsize(
      this.dog, dog.hap, min.hap.len, min.prob)
  filter.seg.len.cat <- levels(this.dog.breed.hapsize.dt$seg.len.cat)[1]
  
  # use dog.blocks to get breed.posterior summary info
  this.dog.breed.posterior.dt <- breed.breakdown.by.posterior(
      this.dog, dog.block, min.prob, min.hap.len, seg.len.total)
  
  # perform breed-based segment and percentage cutoffs
  breeds.above.min.seg.and.pct <- this.dog.breed.hapsize.dt[
    seg.len.cat != filter.seg.len.cat, 
    list(num.seg=sum(num.seg), pct.breed.total=sum(pct.cat.total)), 
    by=breed][
      num.seg > min.seg & pct.breed.total > min.pct, breed]

  #save this dog's pcts
  kept.breed.pcts <- this.dog.breed.hapsize.dt[
    seg.len.cat != filter.seg.len.cat, 
    list(num.seg=sum(num.seg), pct.breed.total=sum(pct.cat.total)), 
    by=breed][
      num.seg > min.seg & pct.breed.total > min.pct]

  #save the mean posterior
  mean.post <- mean(dog.block[indiv==this.dog,prob])
  
  # do a per-window cutoff of breeds
  dog.block <- dog.block[breed %ni% breeds.above.min.seg.and.pct, filter.out := T]

  # do a per-haplotype cutoff of breeds
  dog.hap <- dog.hap[breed %ni% breeds.above.min.seg.and.pct, filter.out := T]
  
  #make multicolored bar graph
  this.bar <- multicolor.bar(this.dog.breed.posterior.dt[
    breed %in% breeds.above.min.seg.and.pct], this.name)
  pdf(file= paste(out.bar, this.name, '_bar','.pdf',sep=''),
      width=9.4, height=5)
  print(this.bar)
  dev.off()
  
  # make pie chart
  this.pie <- make.pie(this.dog.breed.posterior.dt[
    breed %in% breeds.above.min.seg.and.pct], this.name)
  pdf(file= paste(out.pie, this.name, '_pie','.pdf',sep=''),
      width=9, height=6)
  print(this.pie)
  dev.off()

  # paint chromosomes
  this.painting <- paint.chr(dog.hap, this.dog, this.name, breeds.above.min.seg.and.pct)
    pdf(file= paste(out.paint, this.name, '_chrpaint','.pdf',sep=''),
      width=10, height=8)
  print(this.painting)
  dev.off()

  #output statistics 
  write.table(kept.breed.pcts, paste(out.stats, this.name, '_breedPcts','.txt',sep=''), quote=F )
  write.table(mean.post, paste(out.stats, this.name, '_meanPosterior','.txt',sep=''), quote=F )
}

make.dog.block <- function(tped.probs.dt) {
  
  # get the last marker on every chr so we don't hop to the next to get
  # block/window ends
  chr.max.mkr <- markers[, marker[which.max(pos)], by=chr]$V1
  
  dog.block <- tped.probs.dt[order(indiv,chr,hap,marker), 
    # use rleid() to convert runs of breed windows into haplotype blocks
    block := rleid(as.vector(breed_guess)), by=c('indiv','hap','chr')][,
      `:=`(
        prob= mean(prob), 
        mkr.start= min(marker)-79,
        mkr.end= max(marker),
        num.win= .N),
      by=c('indiv','hap','chr','block')][
      # merge with markers to get block START positions
      markers[, list(
        mkr.start=marker,
        pos.start=pos, 
        mrg.start=morgans)], on='mkr.start', nomatch=0][
        # merge with markers to get block END positions
        markers[marker %ni% (chr.max.mkr+1), list(
          # minus 1 to go to the start of the next window, unless end of chr
          mkr.end=ifelse(marker %in% chr.max.mkr,marker,marker-1), 
          pos.end=pos, # start on last base of previous window 
          mrg.end=morgans)], on='mkr.end', nomatch=0][,
          # add position/morgan length to each block
          `:=`(pos.len = pos.end-pos.start, mrg.len = mrg.end-mrg.start)][]
  
  dog.block <- dog.block[, 
    `:=`(mkr.win.start=marker-79, mkr.win.end=marker)][
      # merge with markers to get WINDOW START positions
      markers[, list(
        mkr.win.start=marker,
        pos.win.start=pos, 
        mrg.win.start=morgans)], on='mkr.win.start', nomatch=0][
        # merge with markers to get WINDOW END positions
        markers[marker %ni% (chr.max.mkr+1), list(
          # minus 1 to go to the start of the next window, unless end of chr
          mkr.win.end=ifelse(marker %in% chr.max.mkr,marker,marker-1),
          pos.win.end=pos, # start on last base of previous window 
          mrg.win.end=morgans)], on='mkr.win.end', nomatch=0][,
          # add position/morgan length to each block
          `:=`(
            pos.win.len = pos.win.end-pos.win.start, 
            mrg.win.len = mrg.win.end-mrg.win.start)][]
  
  # rename breed_guess to breed
  names(dog.block)[which(names(dog.block)=='breed_guess')] <- 'breed'

  stopifnot(all(dog.block[, pos.len > 0]))
  stopifnot(all(dog.block[, mrg.len > 0]))
  return(dog.block)
}

# Return a breakdown of window counts per dog/breed/posterior cut category
breed.breakdown.by.posterior <- function(
    this.dog, dog.block, min.prob, min.hap.len, seg.len.total) {

  # generate probability breaks every 10% past the min prob
  prob.breaks <- c(0, seq(min.prob, 1, by=0.1))

  dog.win.dt <- dog.block[indiv==this.dog][, 
    prob.cat := cut(prob, breaks=prob.breaks)][,
      list(
        num.win= .N,
        # multiply the pos.len by the T/F value of 
        # if the pos.len is > than the min hap len
        pct.cat.total= (
          sum(as.numeric(pos.win.len) * (pos.len > min.hap.len))/seg.len.total)
      ),
      by=c('breed','prob.cat')][,
        pct.breed.total := sum(pct.cat.total),
        by=c('breed')]
    
  # flag window category below min prob
  dog.win.dt[, below.min.prob := prob.cat == levels(prob.cat)[1]]
  
  return(dog.win.dt)
}

breed.breakdown.by.hapsize <- function(this.dog, dog.hap, min.hap.len, min.prob) {
  seg.len.breaks <- c(0,min.hap.len,1.3e7,1.75e7,2.25e7,2.9e7,1.9e8)/1e6
  seg.len.total <- dog.hap[indiv==this.dog, sum(as.numeric(pos.len))]

  this.dog.breed.hapsize.dt <- dog.hap[indiv==this.dog & prob > min.prob][, 
    seg.len.cat := cut(pos.len/1e6, breaks=seg.len.breaks)][, 
      list(
        num.seg= .N, 
        pct.cat.total= sum(as.numeric(pos.len))/seg.len.total), 
      by=c('breed','seg.len.cat')][,
        list(
          num.seg=num.seg,
          pct.breed.total=sum(pct.cat.total),
          seg.len.cat,
          pct.cat.total=pct.cat.total), 
        by=c('breed')]
  return(this.dog.breed.hapsize.dt)
}

### MULTICOLOR BAR CHART ###

multicolor.bar <- function(this.dog.breed.posterior.dt, this.name) {
  
  ggplot(this.dog.breed.posterior.dt[below.min.prob==F]) +
    geom_bar(stat='identity', 
      aes(
        x=reorder(breed, -pct.cat.total, FUN= sum), 
        y=pct.cat.total), fill='white') +
    geom_bar(stat='identity', 
      aes(
        x=reorder(breed, -pct.cat.total, FUN= sum), 
        y=pct.cat.total,
        fill=breed,
        alpha=prob.cat)) +
    scale_y_continuous(labels = scales::percent) +
    scale_alpha_discrete(range=c(0.3, 1),
     labels=c(">50",">60",">70", ">80", ">90")) +
    scale_fill_discrete(guide=F) +
    labs(x="Breed", y="Percent of Genome", title=this.name) +
    guides(alpha=guide_legend(title="Estimated percent correct")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1))
}

### PIE CHART ###

make.pie <- function(this.dog.breed.posterior.dt, this.name) {
  
  pct.called <- this.dog.breed.posterior.dt[
      below.min.prob == F, sum(pct.cat.total)]

  not.called <- data.table(breed="no call", V1=1-pct.called)

  #set breed pcts
  this.dog.breed.pct <- this.dog.breed.posterior.dt[
      below.min.prob == F, list(sum(pct.cat.total)), by=breed]

  #set breed order by percentage
  this.dog.breed.pct[, breed := factor(breed,levels = rev(breed[order(V1)]),
      ordered = TRUE)]
  
  #set nocall
  this.dog.breed.pct <- rbind(this.dog.breed.pct, not.called)

  blank_theme <- theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      axis.text.x=element_blank())

  ggplot(this.dog.breed.pct, aes(x="", y=V1, fill=breed)) +
    geom_bar(width = 1, stat = "identity", position='stack', color='gray70', size=0.2) +
    scale_fill_manual(values=c(rainbow(nrow(this.dog.breed.pct)-1), "grey"), 
      name = "Breed") +
    coord_polar("y", start=0) +
    blank_theme +
    labs(title = this.name)
}

### CHROMOSOME PAINTING ###
paint.chr <- function(dog.hap, this.dog, this.name, my.breeds) {
  
  this.dog.block.nocall <- dog.hap[indiv==this.dog]
  this.dog.block.nocall <- this.dog.block.nocall[
    filter.out == T, breed := "No Call"]

  # breed order fix?
  #breed.order <- this.dog.data[, sum(pos.len), by=breed][, breed[order(V1)]]
  
  ggplot(this.dog.block.nocall,
    aes(x=pos.start, xend=pos.end, 
        y=hap, yend=hap,
        color=breed)) + 
    ylim(0.8,2.2) +
    xlab("Chromosome Position") +
    ylab("") +
    labs(title = this.name) +
    geom_segment(size=12, alpha=0.5) +
    scale_color_manual(
        #labels=c(gsub('_',' ',as.character(my.breeds)), 'No Call'),
        values=c(rainbow(length(my.breeds)), "dark grey"),
        name = "Breed") +
    guides(color=guide_legend(ncol=1, override.aes = list(size=4))) +
    facet_wrap(~chr) +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),
      legend.key = element_rect(fill='white'),
      strip.background = element_rect(fill='white'),
      strip.text = element_text(hjust=0.1),
      line = element_blank())
}
