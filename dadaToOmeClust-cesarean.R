# Using dada2 as preliminary pipeline to omeClust and MaAsLin2
# Clark Gaylord
# for PUBH 6885 -- Computational Biology
# Fall 2021

# Simplified version for the tutorial
# This will be basis for R Markdown

# now with timing!
TIMING = T
DEBUG  = T
SINK   = T
sink_file = "pipeline.out"

# Where are the .fastq.gz files?
data_directory <- "./ENAfiles"
# Choose aa arbitrary file to view the quality scores for
example_sequence_file <- "ERR3947849.fastq.gz"
metadata_file <- "./fecal_onemonth-metadata.txt"

if (SINK) {
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

# Pick a file to eyeball the quality scores of. You can do more like this.
plotQualityProfile(paste0(data_directory, "/", example_sequence_file))

# raw data files are compressed
Infant_ALL_seq_files <- list.files(
    data_directory,
    pattern = "fastq.gz",
    full.names = T)

# Eyeballing the quality plots we choose to trim the ends
# One can also add filtering on quality score (say minQ=15 or 20)
# Filtered fastq files are placed into FilteredENAfiles directory
if (TIMING) now()

filtered_directory <- paste0(data_directory, ".filtered")
filtFastqReads <- filterAndTrim(
    data_directory,
    filtered_directory,
    trimLeft = 25,
    trimRight = 10)

# Remove the .fastq.gz
# This is fragile, relying on the "double extension"
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
    filtered_directory,
    pattern = "ERR",
    full.names = T)
errRates <- learnErrors(Infant_ALL_filt_files, multithread = T)
plotErrors(errRates, nominalQ = T)

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
	# we have single-ended reads; do these for double-ended
    # sapply(dadaRs, getN),
    # sapply(mergers, getN),
    rowSums(seqTable_Infant.nochim))
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

# We don't do taxonomy -- commenting out for reader's benefit
# taxa <- assignTaxonomy(
#     seqTable_Infant.nochim,
#     "silva_nr_v132_train_set.fa.gz",
#     multithread=TRUE)
# with 16S V4, there is probably no point to addSpecies but
# it is included in the dada2 tutorial
# if (TIMING) now()

# Note at this time we don't actually do anything with the taxonomy
# You could save yourself a little time by removing this.
# taxa <- addSpecies(
#     taxa,
#     "silva_species_assignment_v132.fa.gz")
# taxa.df <- as.data.frame(taxa)
# should we chuck sequences with <50 bases? Unusual sequence lengths?

# sequences <- (data.frame(cbind(
#     rownames(taxa.df),
#     as.numeric(nchar(rownames(taxa.df))))))
# colnames(sequences) <- c('sequence','length')
# sequences.df <- as.data.frame(sequences)

#
# Let's fold the metadata into the mix.
# There is much wrangling that has to happen here -- do what you need
# munge your metadata
# Ultimately omeClust and MaAsLin2 will need relatively simple arrays
#
#
# And creating "Possible" worksheet is just ugly beyond words
metadata <- read.delim("fecal_onemonth-metadata.txt")
# Metadata_Baby_Seeding_all_samples_final <-
#     read_excel("Metadata_Baby_Seeding_all_samples_final.xlsx",
#                sheet = "Possible")
metadata.sorted <-
    metadata[order(metadata$ena_accession),]
# The cool kids would use tidyverse
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

dist_seqTable_Infant.nochim <- as.matrix(dist(seqTable_Infant.nochim))

# Now we write out the data that omeClust will need
write.table(
    metadata.reorder,
    file="ome_infant_metadata.txt",
    sep="\t",
    row.names = F,
    quote = F)
write.table(
    dist_seqTable_Infant.nochim,
    file="ome_infant_distance.txt",
    row.names = T,
    col.names = T,
    sep="\t",
    quote=F)


# How to run omeClust (it isn't an R package)
#  omeClust
#      -i ome_infant_distance.txt \
#      -o ome_results \
#      --metadata ome_infant_metadata.txt \
#      --plot
#
# Wouldn't it be cool to go get the omeClust results and insert them here?
# That would sure make knitting a lot more complete.
if (TIMING) now()
system("omeClust -i ome_infant_distance.txt -o ome_results --metadata ome_infant_metadata.txt --plot")

if (TIMING) now()
Maaslin2(
    sqrt(seqTable_Infant.nochim.df),
    birth_mode,
    output="maaslin_infant-sqrt",
    reference = c("birth_mode", "CS"))

if (TIMING) now()
if (SINK) sink()

