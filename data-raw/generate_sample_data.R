## load Spectronaut report ------------------------------------------------------------------------
MSstats.report <- read.csv("GluC_endpoint_inclProtInf_ReportIV.csv", sep = ";")
MSstats.report$R.Condition <- relevel(MSstats.report$R.Condition, ref = "control")

## pre-process ------------------------------------------------------------------------------------
norm.input <- SpectronauttoMSstatsFormat(input = MSstats.report, intensity = "NormalizedPeakArea", filter_with_Qvalue = 0.01, useUniquePeptide = FALSE, fewMeasurements = FALSE)

## peptide level summary --------------------------------------------------------------------------
Protein <- norm.input$Protein
norm.input$ProteinName <- norm.input$PeptideSequence
s <- sample(norm.input$ProteinName, 10)
peptide.quant <- dataProcess(raw = norm.input[norm.input$ProteinName %in% s,], censoredInt = "0", normalization = FALSE, n_top_feature = 5, remove_proteins_with_interference = FALSE)

## two-group comparision ---------------------------------------------------------------------------
levels(quant$ProcessedData$GROUP_ORIGINAL)
c1 <- matrix(c(-1,1), nrow = 1)
c <- c1
row.names(c) <- c("control vs. GluC")
sample.data <- groupComparison(contrast.matrix = c, data = peptide.quant)
