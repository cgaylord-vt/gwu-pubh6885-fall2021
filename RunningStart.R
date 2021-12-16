# Using dada2 as preliminary pipeline to omeClust and MaAsLin2
# Clark Gaylord
# for PUBH 6885 -- Computational Biology
# and PUBH 6860 -- Principles of Bioinformatics
# Fall 2021

# now with timing!
TIMING = T
DEBUG  = F
SINK   = F

if (SINK) {
  sink_file = "RunningStart.out"
  if (file.exists(sink_file)) {
    file.rename(sink_file, paste0(sink_file,".bak"))
  }
  sink(file = sink_file, split = T)
}

# Uncomment any installs that are required
# These should be conditionals instead
# dada2 probably needs to be installed from source
# install.packages("dada2")

# library(devtools)
# devtools::install_github("benjjneb/dada2", ref="v1.18")
# change the ref argument to get other versions

# BiocManager::install("bioDist")
# BiocManager::install("Maaslin2")

# Ended up not using phyloseq
# BiocManager::install("phyloseq")

library(lubridate)
library(tools)
library(readxl)
library(dada2)
library(ggplot2)
library(bioDist)
library(Maaslin2)
# We don't use phyloseq but it is a good choice to explore with dada2
# library(phyloseq)

if (TIMING) now()

# We now assume ENAfiles, etc are subdirectories of cwd

# An early analysis used only the files that started with ERR146
# as a test run
# Infant_146_seq_files <-
#     list.files(path = "../ENAfiles/",
#                pattern = "ERR146",
#                full.names = T)
# plotQualityProfile(Infant_146_seq_files)

# Choose aa arbitrary file to view the quality scores for
# plotQualityProfile("../ENAfiles/ERR1464419.fastq")
plotQualityProfile("ENAfiles/ERR1464420.fastq.gz")

# This one is legitimately random (chosen out of band).
# The other two were just "first"
plotQualityProfile("ENAfiles/ERR3947849.fastq.gz")

Infant_ALL_seq_files <- list.files(
    "ENAfiles",
    pattern = "fastq.gz",
    full.names = T)

# Eyeballing the quality plots we choose to trim the ends
# One can also add filtering on quality score (say minQ=15 or 20)
# Filtered fastq files are placed into FilteredENAfiles directory
if (TIMING) now()

filtFastqReads <- filterAndTrim(
    "ENAfiles",
    "FilteredENAfiles",
    trimLeft = 25,
    trimRight = 10)

# Remove the .fastq.gz
# this is pretty fragile, relying on the "double extension"
# by trying to avoid an ugly string search/replace got something worse
accessionNames <- file_path_sans_ext(file_path_sans_ext(basename(Infant_ALL_seq_files)))
row.names(filtFastqReads) <- accessionNames
# get rid of this ugliness
# substr(row.names(filtFastqReads), start=1, stop=10)

# Some summarization of the filtered fastq reads
head(filtFastqReads)
summary(filtFastqReads)
print(paste(
    "Total reads ",
    sum(filtFastqReads[,1]),
    "and after filter",
    sum(filtFastqReads[,2])))
print("The proportion of reads retained is")
sum( filtFastqReads[,2]) / sum(filtFastqReads[,1] )
plot(filtFastqReads[,2]~filtFastqReads[,1],
     main = "Filtered Reads",
     ylab = "Total sequence filtered",
     xlab = "Total sequence raw")
# Add the reference line y=x
abline(a=0, b=1)

# The axes should be more dynamically formatted
plot(density(filtFastqReads[,2], adj=0.5),
     main = "Density estimate of number of reads per sample",
     xaxt = "n", yaxt = "n")
axis(side = 1,
     labels = c("0", "100,000", "200,000", "300,000", "400,000"),
     at = c(0,100000,200000,300000,400000))
rug(filtFastqReads[,2])

if (TIMING) now()

Infant_ALL_filt_files <- list.files(
    "FilteredENAfiles",
    pattern = "ERR",
    full.names = T)
errRates <- learnErrors(Infant_ALL_filt_files, multithread = T)
plotErrors(errRates, nominalQ = T)

##############################################################
### We found Q25 removed about half the data and error rates
### were not better
### Removed for final run
### filtFastqReads_q25 <- filterAndTrim("ENAfiles", 
###   "FilteredENAfiles_Q25", trimLeft = 25, trimRight = 10,
###   minQ = 25, multithread = T)
### head (filtFastqReads_q25)
### print(paste(
###     "Total reads ",
###     sum(filtFastqReads_q25[,1]),
###     "and after filter",
###     sum(filtFastqReads_q25[,2])))
### print("The proportion of reads retained is")
### sum(filtFastqReads_q25[,2])/sum(filtFastqReads_q25[,1])
### Infant_ALL_q25_files <- list.files(
###     "FilteredENAfiles_Q25/",
###     pattern = "ERR",
###     full.names = T)
### errRates_q25 <- learnErrors(Infant_ALL_q25_files, multithread = T)
### plotErrors(errRates_q25, nominalQ = T)

# q25 is very conservative and throws away about half the reads
# error plot with simple trimmed sequences looks ok (based on dada2
# tutorial)
# We do not use the q25 sample moving forward
# Maybe q15 would be better. Ad hoc investigations previously
# did not seem to make a compelling case for using minQ
##############################################################

#
# The final sequence identification
#
if (TIMING) now()

dada_Infant <- dada(
    Infant_ALL_filt_files,
    err = errRates,
    multithread = TRUE)

if (TIMING) now()

seqTable_Infant <- makeSequenceTable(dada_Infant)
row.names(seqTable_Infant) <- accessionNames
dim(seqTable_Infant)
table(nchar(getSequences(seqTable_Infant)))

seqTable_Infant.nochim <- removeBimeraDenovo(
    seqTable_Infant,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE)
dim(seqTable_Infant.nochim)
table(nchar(getSequences(seqTable_Infant.nochim)))
seqTable_Infant.nochim.df <- as.data.frame(seqTable_Infant.nochim)

# Now what would be really cool is to randomly select a couple
# sequences (or even just one), dynamically construct the plot
# and even do the BLAST query to get the best genus match.
# But nah.
#
# Create an example abundance plot
exampleSequence <- colnames(seqTable_Infant.nochim.df[11])
print(paste("Example sequence is ", exampleSequence))
subTitle <- paste("Nonzero abundance (",sum(seqTable_Infant.nochim.df[11]>0),")")
boxplot(sqrt(seqTable_Infant.nochim.df[11])[seqTable_Infant.nochim.df[11]>0],
        main="Example Bifidobacterium Abundance", sub = subTitle, ylab="sqrt(count)")


getN <- function(x) sum(getUniques(x))
track <- cbind(
    filtFastqReads,
    sapply(dada_Infant, getN), 
    # sapply(dadaRs, getN),
    # sapply(mergers, getN),
    rowSums(seqTable_Infant.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c(
    "input",
    "filtered",
    "denoisedF",
    # "denoisedR",
    # "merged",
    "nonchim")
rownames(track) <- accessionNames
head(track)
total_track <- colSums(track)
prop_track <- total_track/total_track[1]
total_track
prop_track

if (TIMING) now()

taxa <- assignTaxonomy(
    seqTable_Infant.nochim,
    "silva_nr_v132_train_set.fa.gz",
    multithread=TRUE)
# with 16S V4, there is probably no point to addSpecies but
# it is included in the dada2 tutorial
if (TIMING) now()

# Note at this time we don't actually do anything with the taxonomy
# You could save yourself a little time by removing this.
taxa <- addSpecies(
    taxa,
    "silva_species_assignment_v132.fa.gz")
taxa.df <- as.data.frame(taxa)
# should we chuck sequences with <50 bases? Unusual sequence lengths?

sequences <- (data.frame(cbind(
    rownames(taxa.df),
    as.numeric(nchar(rownames(taxa.df))))))
colnames(sequences) <- c('sequence','length')
sequences.df <- as.data.frame(sequences)

##############################################################
#
# selfConsistent ASV sounds like a great idea but when we brought this
# all the way through omeClust and MaAsLin2 the results were nearly
# identical. In some cases, the reproducibility of the dada process
# might not be this strong.
#
# Removed in interest of time from final run but left in for the reader.
#
# Compare with selfConsitent ASV. This takes a lot longer to run than without
# but all samples are used
# dada_Infant_selfConsist <- dada(
#	Infant_ALL_filt_files,
#	err = errRates,
#	multithread = TRUE,
#	selfConsist = TRUE)
# seqTable_Infant_selfConsist <- makeSequenceTable(dada_Infant_selfConsist)
# dim(seqTable_Infant_selfConsist)

# table(nchar(getSequences(seqTable_Infant_selfConsist)))

# seqTable_Infant_selfConsist.nochim <- removeBimeraDenovo(
#	seqTable_Infant_selfConsist, 
#	method = "consensus", 
#	multithread = TRUE, 
#	verbose = TRUE)
# table(nchar(getSequences(seqTable_Infant_selfConsist.nochim)))
# taxa.selfCons <- assignTaxonomy(
#	seqTable_Infant.nochim,
#	"silva_nr_v132_train_set.fa.gz",
#	multithread=TRUE)
# taxa.selfCons <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
# taxa.selfCons.df <- as.data.frame(taxa)
#
# sequences.selfCons <- (data.frame(cbind(
#	rownames(taxa.selfCons.df),
#	as.numeric(nchar(rownames(taxa.selfCons.df))))))
#
##############################################################

#
# Let's fold the metadata into the mix.
# There is much wrangling that has to happen here -- do what you need
# munge your metadata
# Ultimately omeClust and MaAsLin2 will need relatively simple arrays
#
#
# And creating "Possible" worksheet is just ugly beyond words
Metadata_Baby_Seeding_all_samples_final <-
    read_excel("Metadata_Baby_Seeding_all_samples_final.xlsx",
               sheet = "Possible")
metadata.sorted <-
    Metadata_Baby_Seeding_all_samples_final[
        order(Metadata_Baby_Seeding_all_samples_final$ena_accession),]
#metadata.reorder <- as.matrix((cbind(
#	metadata.sorted$ena_accession,
#	as.matrix(metadata.sorted)[,1:4])))
metadata.reorder <- as.matrix(cbind(
    metadata.sorted$ena_accession,
    metadata.sorted$baby_sex,
    metadata.sorted$birth_mode))
colnames(metadata.reorder) <- c(
    "ena_accession",
    "baby_sex",
    "birth_mode")
birth_mode <- as.data.frame(
    metadata.reorder[,3])
colnames(birth_mode) <- "birth_mode"
row.names(birth_mode) <- metadata.reorder[,1]

# I think this is fixed now from above basename
# XXX
#row.names( seqTable_Infant.nochim ) <- substr(
#    row.names(seqTable_Infant.nochim), start = 1, stop = 10)

#row.names( seqTable_Infant_selfConsist.nochim ) <- substr(
#	row.names(seqTable_Infant_selfConsist.nochim), start = 1, stop = 10)
# seqTable_Infant_selfConsist.nochim.df <- as.data.frame(seqTable_Infant_selfConsist.nochim)


# dist_seqTable_Infant.nochim <- dist(seqTable_Infant.nochim)
dist_seqTable_Infant.nochim <- as.matrix(dist(seqTable_Infant.nochim))
# dist_seqTable_Infant_selfConsist.nochim <- as.matrix(dist(seqTable_Infant_selfConsist.nochim))
# dist_seqTable_Infant.nochim.spearman <- as.matrix(spearman.dist(seqTable_Infant.nochim))

# Now we write out the data that omeClust will need
write.table(
    metadata.reorder,
    file="infant_metadata.txt",
    sep="\t",
    row.names = F,
    quote = F)
write.table(
    dist_seqTable_Infant.nochim,
    file="infant_distance.txt",
    row.names = T,
    col.names = T,
    sep="\t",
    quote=F)

# write.table(dist_seqTable_Infant_selfConsist.nochim, file="infant_distance_selfCons.txt", row.names = T, col.names = T, sep="\t", quote=F)
# write.table(dist_seqTable_Infant.nochim.spearman, file="../omeClust-analysis/infant_distance_spear.txt", row.names = T, col.names = T, sep="\t", quote=F)

# How to run omeClust (it isn't an R package)
#  omeClust
#      -i infant_distance.txt \
#      -o omeClust-results \
#      --metadata infant_metadata.txt \
#      --plot
#
# Wouldn't it be cool to go get the omeClust results and insert them here?
# That would sure make knitting a lot more complete.
if (TIMING) now()
system("omeClust -i infant_distance.txt -o omeClust-results --metadata infant_metadata.txt --plot")

if (TIMING) now()
Maaslin2(
    sqrt(seqTable_Infant.nochim.df),
    birth_mode,
    output="seqInfant-maaslin2-sqrt",
    reference = c("birth_mode", "CS"))

# Maaslin2(
#     sqrt(seqTable_Infant_selfConsist.nochim.df),
#     birth_mode, output="seqInfantselfCons-maaslin2-sqrt",
#     reference = c("birth_mode", "CS"))

if (TIMING) now()
if (SINK) sink()

