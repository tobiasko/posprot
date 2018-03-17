# installation
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")

# I/O
p <- AAString("LGGNEQVTRYILAGVENSKGTFIIDPGGVIRGTFIIDPAAVIRGAGSSEPVTGLDAKTPVISGGPYEYRVEATFGVDESNAKTPVITGAPYEYRDGLDAASYYAPVRADVTPADFSEWSKLFLQFGAQGSPFLK")
Views(p)
p <- readAAStringSet(filepath = "irtfusion.fasta")

#read proteomes from fasta file
e.coli <- readAAStringSet(filepath = "uniprot-proteome_UP000000625.fasta")
summary(e.coli)
H.sapiens <- readAAStringSet(filepath = "fgcz_9606_cnl_20150422.fasta")
#peptides <- AAStringSet(c("DYWRALQNRIREGHVEDVYAYRRRQ", "QRWQSVLARDPNADGEFVFAV"))
summary(peptides)

# toy example
m <- vector("list", length(peptides)) #pre-allocate an empty list of right length
for (i in seq_along(peptides)){
  m[[i]] <- vcountPattern(pattern = peptides[[i]], subject = e.coli)
}

# test job
peptides <- unique(as.character(PepQuantDef$EG.StrippedSequence))
peptides <- AAStringSet(peptides[1:10])
m <- vector("list", length(peptides)) #pre-allocate an empty list of right length
system.time(for (i in seq_along(peptides)){
  m[[i]] <- vmatchPattern(pattern = peptides[[i]], subject = H.sapiens)
})

# can this code be vectorized?
foo <- function(x){vmatchPattern(pattern = x, subject = H.sapiens)}
system.time(
  l <- lapply(as.list(peptides), foo)
)
# yes - but with coercion to primitive data type
# Hmmm...let's use the lapply from BioC
foo <- function(x){vmatchPattern(pattern = x, subject = H.sapiens)}
system.time(
  l <- BiocGenerics::lapply(peptides, foo)
)
# :-) :-) :-)
#does vapply also work?
foo <- function(x){vmatchPattern(pattern = x, subject = H.sapiens)}
system.time(
  l <- BiocGenerics::vapply(peptides, foo)
)
#Nop!

# get ranges
extractAllMatches(subject = H.sapiens[[1]], mindex = l[[1]])
Views(subject = H.sapiens[[1]], l[[1]])

# more testing
# extract 10 substrings [1:10] from the first 10 e.coli proteins
pep <- BiocGenerics::lapply(e.coli[1:10], Views, start = 1, end = 10)
p <- BiocGenerics::sapply(pep, toString)
p <- AAStringSet(p)
# find those guys
foo <- function(x){vmatchPattern(pattern = x, subject = e.coli[1:10])}
system.time(
  index <- BiocGenerics::lapply(p, foo)
)

# plyr version
library("plyr")
system.time(index <-llply(as.list(p), foo))
#nice
peptides <- unique(as.character(PepQuantDef$EG.StrippedSequence))
peptides <- AAStringSet(peptides[1:10])
system.time(
  index <-llply(as.list(peptides), foo)
)

# AhoCorasickTrie
# O(n+m) baby!
library("AhoCorasickTrie")
# extract 10 substrings [1:10] from the first 10 e.coli proteins
pep <- BiocGenerics::lapply(e.coli[1:10], Views, start = 1, end = 10)
p <- BiocGenerics::sapply(pep, toString)
e  <- BiocGenerics::sapply(e.coli, toString)
#p <- AAStringSet(p)
system.time(
  l <- AhoCorasickTrie::AhoCorasickSearch(keywords = p, text = e, alphabet = "aminoacid")
)
system.time(
  l <- AhoCorasickTrie::AhoCorasickSearchList(keywords = p, textList = e, alphabet = "aminoacid")
)

#unlist and convert version


# sesssion
sessionInfo()
citation("Biostrings")