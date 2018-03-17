###########
# info ####
###########
# this script generates AssayLibs for DIA-based analysis in Spectronaut (Biognosys)
# from DDA data/Mascot search results (.dat) converted to a spectral library using the BiblioSpec implementation in skyline-dily


######################
# install package ####
######################

# specL
source("https://bioconductor.org/biocLite.R")
biocLite("specL")
library("specL")
vignette("specL")

# biocG
biocLite("BiocGenerics")
library("BiocGenerics")

###########
# MAIN ####
###########

# Which RT normalization peptides are used by default?
(iRTpeptides) #???

# Generate iRT-C18 scale
Sequence <- c("LGGNEQVTR", "GAGSSEPVTGLDAK", "VEATFGVDESNAK", "YILAGVENSK", "TPVISGGPYEYR", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "LFLQFGAQGSPFLK")
Name <- c("iRT-pep a", "iRT-pep b", "iRT-pep c", "iRT-pep d", "iRT-pep e", "iRT-pep f", "iRT-pep g", "iRT-pep h", "iRT-pep i", "iRT-pep k", "iRT-pep l")
iRT <- c(-24.92, 0.00, 12.39, 19.79, 28.71, 33.38, 42.26, 54.62, 70.52, 87.23, 100.00)
iRT.C18 <- data.frame(Sequence, iRT, row.names = Name)
# adjust column names
t <- data.frame(peptide = Sequence, rt = iRT, row.names = Name)
str(t)

# import blib files
# test case 1 - single .dat to single .blib
# SpecLib <- specL::read.bibliospec("F245177.blib")
# SpecLib.red <- specL::read.bibliospec("F245177.redundant.blib")
# plot(SpecLib)
# plot(SpecLib.red)

# test case 2 - many .dat to single .blib
# here I used skyline-daily to build the blib files (keep peptides with more than 1 best hit)
# A - two searches per .raw file (tryptic + semi ArgC)
# looks like BiblioSpec has problems creating these blib properly
SpecLib <- specL::read.bibliospec("joined_SpecLib.blib")
SpecLib.red <- specL::read.bibliospec("joined_SpecLib.redundant.blib")
plot(SpecLib)
# B - only semi ArgC search
SpecLib <- specL::read.bibliospec("dimethyl_GluC.blib")
SpecLib.red <- specL::read.bibliospec("dimethyl_GluC.redundant.blib")
plot(SpecLib)
# C - only Trp search

# generate SpecLib
# .blib files were filtered at 5% FDR, therefor raw mascot filtering is not needed
dia.windows <- seq(385, 1250, by=25) #classic 25 Da SWATH-MS
AssayLib <- genSwathIonLib(data=SpecLib, data.fit=SpecLib.red, iRT = t, ignoreMascotIonScore = TRUE) #use iRT-C18 std.
AssayLib <- genSwathIonLib(data=SpecLib, data.fit=SpecLib.red, iRT = t, ignoreMascotIonScore = TRUE, topN = 5, breaks = dia.windows) #use iRT-C18 std.
#AssayLib <- genSwathIonLib(data=SpecLib, data.fit=SpecLib.red)
summary(AssayLib)

# write AssayLib to disc
# A
write.spectronaut(AssayLib, file="specL2Spectronaut.top5.txt")
# B
write.spectronaut(AssayLib, file="specL2Spectronaut.semiArgC.txt")

#check assay lib from hard disc
specL2Spectronaut.top5 <- read.delim("~/Documents/RStudio/p2095_CasTC/specL2Spectronaut.top5.txt")
summary(specL2Spectronaut.top5)


# diagnostic plots
op<-par(mfrow=c(2,3))
plot(AssayLib)
par(op)
# output does not look as expected!

###################
# Session Info ####
###################
sessionInfo()

