############################################################################################
# Get set up
############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~
# Load in libraries
#~~~~~~~~~~~~~~~~~~~~~~~~
# Required for segment_data()
library("aroma.affymetrix")	#bioconductor
library("aroma.cn")			#bioconductor
library("DNAcopy")			#bioconductor

#The first two packages may need to be installed using the following install script instead of (or in addition to?) the usual bioConductor method
#source("http://aroma-project.org/hbLite.R");
#hbInstall("aroma.affymetrix");
#hbInstall("aroma.cn");


# Required for get_gene_level_summaries()
library("CNTools")	#bioconductor
library("gdata")	# cran
data("geneInfo")	# load in gene info to use for gene-level summaries, based on build 36, HG18

# Optional settings
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
log <- verbose
options(digits=4)  # Dont display too many decimals.



############################################################################################
# segment_data(dataSet, chipType="GenomeWideSNP_6", dataSetControl="HapMap_Controls")
#
#	Function to perform CBS segmentation on raw CNV data
#
#	input:	dataSet = directory of samples to be analyzed
#			chipType = chip type being analyzed, defaults to "GenomeWideSNP_6"
#			dataSetControl = directory of normal control reference samples, defaults to "HapMap_Controls"
#
#	output:	pathname = path and file name where segmented file is written
#
#	NB: directory tree must conform to requirements of aroma.affymetrix (see http://www.aroma-project.org/setup)
#		The CDF GenomeWideSNP_6,Full.cdf must be available in the associated annotationData directory
############################################################################################

segment_data <- function(dataSet, chipType="GenomeWideSNP_6", dataSetControl="HapMap_Controls"){

#~~~~~~~~~~~~~~~~~~~~~~~~
# Verify annotation files
#~~~~~~~~~~~~~~~~~~~~~~~~
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

gi <- getGenomeInformation(cdf)
print(gi)

si <- getSnpInformation(cdf)
print(si)

acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
print(acs)


#~~~~~~~~~~~~~~~~~~~~~~~~
# Declare raw dataset
#~~~~~~~~~~~~~~~~~~~~~~~~
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
# Get cell line samples
dsR <- AffymetrixCelSet$byName(dataSet, chipType=chipType,checkChipType=FALSE, cdf=cdf)
# Get normal control samples
dsC <- AffymetrixCelSet$byName(dataSetControl, chipType=chipType,checkChipType=FALSE, cdf=cdf)



#~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up cel file names for control files
# Hacky long way to do this
#~~~~~~~~~~~~~~~~~~~~~~~~
setFullNamesTranslator(dsC, function(names, ...) {
# Remove SNP6 extra info info
	names <- gsub("\\_\\_", "_", names)  # replace double underscore with single underscore: __ ==> _
	names <- gsub("\\_\\SNP6\\.0", "", names)
	names <- gsub("\\_\\(GenomeWideSNP\\_6\\)", "", names)
	names <- gsub("\\_2", "", names)
	names <- gsub("2\\_", "", names)
	names <- gsub("[[:space:]]+", "", names)
	
# Turn into comma-separated tags
	names <- gsub("_", ",", names)
	
# Drop tags we dont need
	names <- gsub(",SNP6.0", "", names)
	names <- gsub("080122,", "", names)
	names <- gsub("080214,", "", names)
	names <- gsub("080219,", "", names)
	
# Remove well tag: well followed by comma
	names <- gsub("(A01|A02|A03|A04|A05|A06|A07|A08|A09|A10|A11|A12),", "", names)
	names <- gsub("(B01|B02|B03|B04|B05|B06|B07|B08|B09|B10|B11|B12),", "", names)
	names <- gsub("(C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12),", "", names)
	names <- gsub("(D01|D02|D03|D04|D05|D06|D07|D08|D09|D10|D11|D12),", "", names)
	names <- gsub("(E01|E02|E03|E04|E05|E06|E07|E08|E09|E10|E11|E12),", "", names)
	names <- gsub("(F01|F02|F03|F04|F05|F06|F07|F08|F09|F10|F11|F12),", "", names)
	names <- gsub("(G01|G02|G03|G04|G05|G06|G07|G08|G09|G10|G11|G12),", "", names)
	names <- gsub("(H01|H02|H03|H04|H05|H06|H07|H08|H09|H10|H11|H12),", "", names)
	
# Remove well tag: well preceeded by comma
	names <- gsub(",(A01|A02|A03|A04|A05|A06|A07|A08|A09|A10|A11|A12)", "", names)
	names <- gsub(",(B01|B02|B03|B04|B05|B06|B07|B08|B09|B10|B11|B12)", "", names)
	names <- gsub(",(C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12)", "", names)
	names <- gsub(",(D01|D02|D03|D04|D05|D06|D07|D08|D09|D10|D11|D12)", "", names)
	names <- gsub(",(E01|E02|E03|E04|E05|E06|E07|E08|E09|E10|E11|E12)", "", names)
	names <- gsub(",(F01|F02|F03|F04|F05|F06|F07|F08|F09|F10|F11|F12)", "", names)
	names <- gsub(",(G01|G02|G03|G04|G05|G06|G07|G08|G09|G10|G11|G12)", "", names)
	names <- gsub(",(H01|H02|H03|H04|H05|H06|H07|H08|H09|H10|H11|H12)", "", names)
	
#names <- gsub("(01|02|03|04|05|06|07|08|09|10|11|12),", "", names)
	
	names
})


print(getFullNames(dsR))  # CELL LINES
print(getFullNames(dsC))  # REFERENCE



# ~~~~~~~~~~~~~~~~~~~
# Process cell lines
# ~~~~~~~~~~~~~~~~~~~
# 1. Calibration for crosstalk btween allele probe pairs [CELL LINES]
acc <- AllelicCrosstalkCalibration(dsR, model="CRMAv2")
print(acc)

csC <- process(acc, verbose=verbose)
print(csC)


# 2. Normalization for nucleotide-position probe sequence effects [CELL LINES]
bpn <- BasePositionNormalization(csC, target="zero")
print(bpn)

csN <- process(bpn, verbose=verbose)
print(csN)


# 3. Probe summarization -- fit *total* CN signals [CELL LINES]
plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE)
print(plm)

if (length(findUnitsTodo(plm)) > 0) {
# Fit CN probes quickly (~5-10s/array + some overhead)
	units <- fitCnProbes(plm, verbose=verbose)
	str(units)
# int [1:945826] 935590 935591 935592 935593 935594 935595 ...
	
# Fit remaining units, i.e. SNPs (~5-10min/array)
	units <- fit(plm, verbose=verbose)
	str(units)
}

ces <- getChipEffectSet(plm)
print(ces)


# 4. Normalization for PCR fragment-length effects [CELL LINES]
fln <- FragmentLengthNormalization(ces, target="zero")
print(fln)

cesN <- process(fln, verbose=verbose)
print(cesN)



# ~~~~~~~~~~~~~~~~~~~
# Process Reference samples
# ~~~~~~~~~~~~~~~~~~~
# 1.1. Calibration for crosstalk between allele probe pairs [REFERENCE]
acc1 <- AllelicCrosstalkCalibration(dsC, model="CRMAv2")
csC1 <- process(acc1, verbose=verbose)

# 2.1. Normalization for nucleotide-position probe sequence effects [REFERENCE]
bpn1 <- BasePositionNormalization(csC1, target="zero")
print(bpn1)

csN1 <- process(bpn1, verbose=verbose)
print(csN1)

# 3.1 Probe summarization -- fit *total* CN signals [REFERENCE]
plm1 <- AvgCnPlm(csN1, mergeStrands=TRUE, combineAlleles=TRUE)
print(plm1)

if (length(findUnitsTodo(plm1)) > 0) {
# Fit CN probes quickly (~5-10s/array + some overhead)
	units1 <- fitCnProbes(plm1, verbose=verbose)
	str(units1)
# int [1:945826] 935590 935591 935592 935593 935594 935595 ...
	
# Fit remaining units, i.e. SNPs (~5-10min/array)
	units1 <- fit(plm1, verbose=verbose)
	str(units1)
}

ces1 <- getChipEffectSet(plm1)
print(ces1)


# 4.1 Normalization for PCR fragment-length effects [REFERENCE]
fln1 <- FragmentLengthNormalization(ces1, target="zero")
print(fln1)

cesN1 <- process(fln1, verbose=verbose)
print(cesN1)


#~~~~~~~~~~~~~~~~~~~~~~~~
# Segment the data
#~~~~~~~~~~~~~~~~~~~~~~~~
# Compute average of REF samples
cesN1.ref <- getAverageFile(cesN1, verbose=verbose)

# Set up CBS model with normals as reference
cbs <- CbsModel(cesN, cesN1.ref)
print(cbs);

# Run cbs segmentation
fit(cbs, verbose=verbose)

# Write cbs regions directly (without generating ChromosomeExplorer report)
pathname <- writeRegions(cbs, verbose=verbose) 

# Print path and name of segmented file
return(pathname)

}



############################################################################################
# get_gene_level_summaries(segmentedFile)
#
#	Function to generate gene-level measure of changes in copy number
#
#	input:	segmentedFile = file of segmented CN data (output from segment_data)
#	output:	geneSummaryFile = file name of gene-level CNVs (tab-delimited txt file)
############################################################################################

get_gene_level_summaries <- function(segmentedFile, geneMapInfo){

	# Read in segmented snp6 data
	snp6 = read.table(segmentedFile, sep="\t", header=TRUE)

	# Remove last column with url info
	snp6 = snp6[, c(1:4,6,5)]
	
	# Change column headers so that they are what CNTools expects
	colnames(snp6) = c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
	cnseg <- CNSeg(snp6)

	# Get gene-level calls
#	rdseg.gene <- getRS(cnseg, by="gene", imput=FALSE, XY=FALSE, geneMap=geneMapInfo)
	rdseg.gene <- getRS(cnseg, by="gene", imput=FALSE, XY=FALSE, geneMap=geneMapInfo, what="median")

	reducedseg.gene <- rs(rdseg.gene)

	# Create name of file to write gene level summaries
	geneSummaryFile = paste(segmentedFile, "_geneSum.txt", sep='')
	write.table(reducedseg.gene, geneSummaryFile, sep="\t", quote=FALSE, row.names=FALSE)
	
	return(geneSummaryFile)
}
