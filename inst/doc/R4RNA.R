### R code from vignette source 'R4RNA.Rnw'

###################################################
### code chunk number 1: R4RNA.Rnw:41-55
###################################################
library(R4RNA)

message("TRANSAT prediction in helix format")
transat_file <- system.file("extdata", "helix.txt", package = "R4RNA")
transat <- readHelix(transat_file)

message("RFAM structure in dot bracket format")
known_file <- system.file("extdata", "vienna.txt", package = "R4RNA")
known <- readVienna(known_file)

message("Work with basepairs instead of helices for more flexibility")
message("Breaks all helices into helices of length 1")
transat <- expandHelix(transat)
known <- expandHelix(known)


###################################################
### code chunk number 2: R4RNA.Rnw:60-75
###################################################
#read fasta
fasta_file.mltp <- system.file("extdata", "fasta.mltp.txt", package = "R4RNA")
fasta.mltp <- readBStringSet(fasta_file.mltp)
#read helix
transat.file.mltp <- system.file("extdata", "transat.helix.mltp.txt" ,package = "R4RNA")
helix.transat.mltp <- readHelixMltp(transat.file.mltp)
known.file.mltp <- system.file("extdata", "known.helix.mltp.txt", package = "R4RNA")
helix.known.mltp <- readHelixMltp(known.file.mltp)

is.helix.mltp(helix.transat.mltp)
helix.transat.exp <- expandHelixMltp(helix.transat.mltp)
helix.known.exp <- expandHelixMltp(helix.known.mltp)

is.helix.mltp(helix.transat.exp)
is.helix.mltp(helix.known.exp)


###################################################
### code chunk number 3: R4RNA.Rnw:85-87
###################################################
plotHelix(known, line = TRUE, arrow = TRUE)
mtext("Known Structure", side = 3, line = -2, adj = 0)


###################################################
### code chunk number 4: R4RNA.Rnw:94-96
###################################################
plotHelixMltpSingleLine(copy(helix.known.exp),scale = FALSE)
mtext("Known Structure", side = 3, line = -2, adj = 0)


###################################################
### code chunk number 5: R4RNA.Rnw:100-102
###################################################
plotHelixMltpDoubleLine(copy(helix.known.exp),scale = FALSE)
mtext("Known Structure", side = 3, line = -10, adj = 0)


###################################################
### code chunk number 6: R4RNA.Rnw:111-114
###################################################
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)
mtext("TRANSAT\nPredicted\nStructure", side = 3, line = -5, adj = 0)
mtext("Known Structure", side = 1, line = -2, adj = 0)


###################################################
### code chunk number 7: R4RNA.Rnw:117-120
###################################################
plotDoubleHelixMltpSingleLine(helix.known.exp,helix.transat.exp,scale = FALSE)
mtext("Predicted\nStructure", side = 1, line = -5, adj = 0)
mtext("Known Structure", side = 3, line = -2, adj = 0)


###################################################
### code chunk number 8: R4RNA.Rnw:123-127
###################################################
plotDoubleHelixMltpDoubleLine(helix.known.exp,helix.transat.exp,
                              stable = TRUE,scale = FALSE)
mtext("Predicted\nStructure", side = 3, line = -5, adj = 1)
mtext("Known Structure", side = 3, line = -5, adj = 0)


###################################################
### code chunk number 9: R4RNA.Rnw:136-138
###################################################
message("Filter out helices above a certain p-value")
transat <- transat[which(transat$value <= 1e-3), ]


###################################################
### code chunk number 10: R4RNA.Rnw:144-152
###################################################
message("Assign colour to basepairs according to p-value")
transat$col <- col <- colourByValue(transat, log = TRUE)

message("Coloured encoded in 'col' column of transat structure")
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)

legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.05, bty = "n", border = NA, cex = 0.75, title = "TRANSAT P-values")


###################################################
### code chunk number 11: R4RNA.Rnw:157-165
###################################################
helix.transat.exp <- helix.transat.exp[,1:6]
helix.transat.exp$col <- col <- colourByValueMltp(helix.transat.exp, log = TRUE)

message("Coloured encoded in 'col' column of transat structure")
plotDoubleHelixMltpSingleLine(helix.known.exp,helix.transat.exp,scale = FALSE)

legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
	inset = 0.05, bty = "n", border = NA, cex = 0.75, title = "TRANSAT P-values")


###################################################
### code chunk number 12: R4RNA.Rnw:168-176
###################################################
helix.transat.exp <- helix.transat.exp[,1:6]
helix.transat.exp$col <- col <- colourByValueMltp(helix.transat.exp, log = TRUE)

message("Coloured encoded in 'col' column of transat structure")
plotDoubleHelixMltpDoubleLine(helix.known.exp,helix.transat.exp,scale = FALSE,stable = TRUE)

legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
	inset = 0.05, bty = "n", border = NA, cex = 0.75, title = "TRANSAT P-values")


###################################################
### code chunk number 13: R4RNA.Rnw:189-190
###################################################
plotOverlapHelix(transat, known, line = TRUE, arrow = TRUE, scale = FALSE)


###################################################
### code chunk number 14: R4RNA.Rnw:195-200
###################################################
plotComparisonHelixMltpSingleLine(helix1 = helix.known.exp[,1:6],
                                  helix2 = helix.transat.exp[,1:6],
                                  scale = FALSE)
mtext("Known Structure", side = 3, line = -5, adj = 0)
mtext("Predicted\nStructure", side = 1, line = -2, adj = 0)


###################################################
### code chunk number 15: R4RNA.Rnw:203-208
###################################################
plotComparisonHelixMltpDoubleLine(helix1 = helix.known.exp[,1:6],
                                  helix2 = helix.transat.exp[,1:6],
                                  scale = FALSE)
mtext("Known Structure", side = 3, line = -5, adj = 0)
mtext("Predicted\nStructure", side = 1, line = -2, adj = 0)


###################################################
### code chunk number 16: R4RNA.Rnw:225-232
###################################################
message("Multiple sequence alignment of interest")
library(Biostrings)
fasta_file <- system.file("extdata", "fasta.txt", package = "R4RNA")
fasta <- as.character(readBStringSet(fasta_file))

message("Plot covariance in alignment")
plotCovariance(fasta, known, cex = 0.5)


###################################################
### code chunk number 17: R4RNA.Rnw:237-245
###################################################
message("Multiple sequence alignment of interest")
fasta_file.mltp <- system.file("extdata", "fasta.mltp.txt", package = "R4RNA")
fasta.mltp <- readBStringSet(fasta_file.mltp)

message("Plot covariance in alignment")
plotCovarianceMltpSingleLine(helix = helix.transat.exp[,1:7],
                             msa = fasta.mltp,grid = FALSE,
                             scale = FALSE)


###################################################
### code chunk number 18: R4RNA.Rnw:249-257
###################################################
message("Multiple sequence alignment of interest")
fasta_file.mltp <- system.file("extdata", "fasta.mltp.txt", package = "R4RNA")
fasta.mltp <- readBStringSet(fasta_file.mltp)

message("Plot covariance in alignment")
plotCovarianceMltpDoubleLine(helix = helix.transat.exp[,1:7],
                             msa = fasta.mltp,grid = FALSE,
                             legend = FALSE,scale = FALSE)


###################################################
### code chunk number 19: R4RNA.Rnw:270-276
###################################################
plotCovarianceComparisonMltpSingleLine(msa = fasta.mltp,
                                       helix1 = helix.known.exp[,1:6],
                                       helix2 = helix.transat.exp[,1:6],
                                       scale = FALSE,legend = FALSE,
                                       species = 0)



###################################################
### code chunk number 20: R4RNA.Rnw:281-286
###################################################
plotDoubleCovarianceComparisonMltpSingleLine(msa = fasta.mltp,
                                             helix1 = helix.known.exp[,1:6],
                                             helix2 = helix.transat.exp[,1:6],
                                             scale = FALSE,legend = FALSE,
                                             species = 0,dist.y.between = 5 )


###################################################
### code chunk number 21: R4RNA.Rnw:292-297
###################################################
plotCovarianceComparisonMltpDoubleLine(msa = fasta.mltp,
                                       helix1 = helix.known.exp[,1:6],
                                       helix2 = helix.transat.exp[,1:6],
                                       scale = FALSE,legend = FALSE,
                                       species = 0)


###################################################
### code chunk number 22: R4RNA.Rnw:308-309
###################################################
plotCovariance(fasta, transat, cex = 0.5, conflict.col = "grey")


###################################################
### code chunk number 23: R4RNA.Rnw:314-321
###################################################
message("Multiple sequence alignment of interest")
fasta_file.mltp <- system.file("extdata", "fasta.mltp.txt", package = "R4RNA")
fasta.mltp <- readBStringSet(fasta_file.mltp)

message("Plot covariance in alignment")
plotCovarianceMltpSingleLine(helix = helix.transat.exp[,1:7],msa = fasta.mltp,
                             conflict.col = "grey",legend = TRUE,grid = FALSE,scale = FALSE)


###################################################
### code chunk number 24: R4RNA.Rnw:324-331
###################################################
message("Multiple sequence alignment of interest")
fasta_file.mltp <- system.file("extdata", "fasta.mltp.txt", package = "R4RNA")
fasta.mltp <- readBStringSet(fasta_file.mltp)

message("Plot covariance in alignment")
plotCovarianceMltpDoubleLine(helix = helix.transat.exp[,1:7],msa = fasta.mltp,
                             conflict.col = "grey",grid = FALSE,scale = FALSE)


###################################################
### code chunk number 25: R4RNA.Rnw:342-346
###################################################
col <- colourByCovariation(known, fasta, get = TRUE)
plotCovariance(fasta, col, grid = TRUE, legend = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.1, bty = "n", border = NA, cex = 0.37, title = "Covariation")


###################################################
### code chunk number 26: R4RNA.Rnw:349-355
###################################################

col <- colourByCovariationMltp(helix = helix.transat.exp,msa = fasta.mltp,get = TRUE)

plotCovarianceMltpSingleLine(helix = col,msa = fasta.mltp,grid = TRUE,scale = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
       inset = 0.1, bty = "n", border = NA, cex = 0.5, title = "Covariation")


###################################################
### code chunk number 27: R4RNA.Rnw:361-366
###################################################
custom_colours <- c("green", "blue", "cyan", "red", "black", "grey")
plotCovariance(fasta, col <- colourByConservation(known, fasta, get = TRUE),
    palette = custom_colours, cex = 0.5)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "Conservation")


###################################################
### code chunk number 28: R4RNA.Rnw:369-377
###################################################
custom_colours <- c("green", "blue", "cyan", "red", "black", "grey")
plotCovarianceMltpSingleLine(msa = fasta.mltp,
                             helix = col <- colourByConservationMltp(helix = helix.transat.exp,
                                                                     msa = fasta.mltp,get = TRUE,
                                                                     cols = custom_colours),
                             palette = custom_colours,scale = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
	inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "Conservation")


###################################################
### code chunk number 29: R4RNA.Rnw:383-387
###################################################
col <- colourByCanonical(known, fasta, custom_colours, get = TRUE)
plotCovariance(fasta, col, base.colour = TRUE, cex = 0.5)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
    inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "% Canonical")


###################################################
### code chunk number 30: R4RNA.Rnw:390-394
###################################################
col <- colourByCanonicalMltp(helix = helix.transat.exp,msa = fasta.mltp,get = TRUE)
plotCovarianceMltpSingleLine(helix = col,msa = fasta.mltp,scale = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
	inset = 0.15, bty = "n", border = NA, cex = 0.75, title = "% Canonical")


###################################################
### code chunk number 31: R4RNA.Rnw:400-404
###################################################
col <- colourByUnknottedGroups(known, c("red", "blue"), get = TRUE)
plotCovariance(fasta, col, base.colour = TRUE, legend = FALSE,
               species = 23, grid = TRUE, text = TRUE,
               text.cex = 0.2, cex = 0.5)


###################################################
### code chunk number 32: R4RNA.Rnw:407-409
###################################################
col <- colourByUnknottedGroupsMltp(helix = helix.transat.exp,cols = c("red","blue"),get = TRUE)
plotCovarianceMltpSingleLine(helix = col,msa = fasta.mltp,base.colour = TRUE,scale = FALSE)


###################################################
### code chunk number 33: R4RNA.Rnw:416-418
###################################################
hic_file <- system.file("extdata", "GSE63525.chr10.chr11.500kb.txt", package = "R4RNA")
helix_hic <- StrawToHelix(file = hic_file, chr1 = "chr10",chr2 = "chr11",scale = 500000)


###################################################
### code chunk number 34: R4RNA.Rnw:424-427
###################################################
helix_hic.trimmed <- helix_hic[value > 40]
attr(helix_hic.trimmed,"length") <- attr(helix_hic,"length")
helix_hic.trimmed <- colourByValueMltp(helix_hic.trimmed,get = TRUE)


###################################################
### code chunk number 35: R4RNA.Rnw:432-433
###################################################
plotHelixMltpSingleLine(helix = helix_hic.trimmed,shape = "triangle",scale = FALSE)


###################################################
### code chunk number 36: R4RNA.Rnw:438-439
###################################################
plotHelixMltpSingleLine(helix = helix_hic.trimmed,shape = "heatmap",scale = FALSE)


###################################################
### code chunk number 37: R4RNA.Rnw:446-448
###################################################
toLatex(sessionInfo())



