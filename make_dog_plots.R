### 0. get arguments #######################################

args = commandArgs(trailingOnly=TRUE)

### 1. Load packages and functions  ########################

require(data.table)
require(ggplot2)
require(ggrepel)
require(base)
require(scales)

# get directory name

# this way of getting script dir was failing (sometimes?)
#script.dir <- dirname(sys.frame(1)$ofile)
# now passing it in 
script.dir <- args[13]

# load recombination functions from scripts dir
source(file.path(script.dir, 'recomb_fxn.R'))

# load dog plot functions from scripts dir
source(file.path(script.dir, 'dog_plot_fxn.R'))


### 2. Load Darwin's Dogs Data ########################

breed.leg <- data.table(read.table(args[1],header=T))

# get path for this output
out.data.path <- '/Users/boettger/postdoc/mutts/darwinsdogs_SMoutput/March2_12perbreed/'
        
# get tped, tfam, probs, and marker files
out.tped.file <- paste(args[2],'.tped',sep='')
out.tfam.file <- paste(args[2],'.tfam',sep='')
out.probs.file <- paste(args[2],'.Probs.tped',sep='')

# get markers
markers <- data.table(read.table(args[3], header=F)[,c(2,1,3:6)])

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

# get dog names
dog.names <- data.table(read.table(args[4], header=F))
colnames(dog.names) <- c("ID", "Name")

# Finally, load the tped.probs for this dataset.
tped.probs.dt <- load_supportmix_output(
  out.tped.file, out.tfam.file, out.probs.file, markers, breed.leg)

### 3. Call plot-making functions  ########################

# Cutoff for minimum number of segments
min.seg <- as.numeric(args[5])

# Cutoff for minimum breed percentage
min.pct <- as.numeric(args[6])

# Cutoff for min haplotype size (approx 1mb per window)
# a test to examine window size:
# dog.block.win.dt[pos.len < 2e6, table(num.win)]
# (dog.block.win.dt is inside count.dog.windows())
min.hap.len <- as.integer(args[7])

# Cutoff for posterior probability
min.prob <- as.numeric(args[8])

out.bar <- args[9]
out.pie <- args[10]
out.paint <- args[11]
out.stats <- args[12]

# get dog blocks (so we can look at it outside the fxns)
dog.block <- make.dog.block(tped.probs.dt)

# get all stats
dog.names[, make_dog_plots(ID, Name, tped.probs.dt, 
    min.hap.len, min.seg, min.pct, min.prob), by=ID]

# make all plots
breed.pct.table <- dog.names[, make_dog_plots(ID, Name, tped.probs.dt, 
    min.hap.len, min.seg, min.pct, min.prob), by=ID]

# rename breed names to add underscores
breed.pct.table[, breed := factor(
  breed,
  levels=levels(breed),
  labels=gsub(' ','_', as.character(levels(breed))))]

##output tally of percents
write.table(breed.pct.table, paste(out.stats,'breed_pct_table.txt', sep=''), row.names=F,
            quote=F)

#1: a550771-4306235-112317-918_B11.CEL    Boots



