#!/usr/bin/env Rscript

# Copyright 2011 Anneleen Daemen, Obi Griffith, Laura Heiser

# DrugResponseToolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# DrugResponseToolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with DrugResponseToolbox. If not, see <http://www.gnu.org/licenses/>.



###Load raw data at individual patient-level (any combination of following data formats):
# Argument 1 = patient ID (string)
# Argument 2 = 'NA' or 1) Affymetrix GeneChip Human Genome U133A CEL-file
# Argument 3 = 'NA' or 2) Illumina Infinium Human Methylation27 BeadChip Kit tab-delimited file, with cg-probe ids in column 1 and values between 0 and 1 in column 2. The values represent the proportion of methylated DNA at each single CpG locus.
# Argument 4 = 'NA' or 3) Affymetrix Genome-Wide Human SNP Array 6.0 CEL-file

### Example run on command line:
### ./DrugResponseToolbox.R "FullCombination_test" "U133A_test_tumor.CEL" "Meth_test_tumor.txt" "SNP6_test_cellline.CEL"

#args=c('FullCombination_test','U133A_test_tumor.CEL','Meth_test_tumor.txt','SNP6_test_cellline.CEL')
args=(commandArgs(TRUE))
args.all=commandArgs(trailingOnly = FALSE)
PatientID=args[1] #'FullCombination_test'
AffyCELfile=args[2] #'U133A_test_tumor.CEL'
AffyCELfile_clean=gsub("_",x=AffyCELfile,replacement="\\\\_")#Escapes "_" in CELfile names
Methfile=args[3] #'Meth_test_tumor.txt'
Methfile_clean=gsub("_",x=Methfile,replacement="\\\\_")#Escapes "_" in CELfile names
SNPCELfile=args[4] #'SNP6_test_cellline.CEL'
SNPCELfile_clean=gsub("_",x=SNPCELfile,replacement="\\\\_")#Escapes "_" in CELfile names

#Cutoff/Threshold values
ModelAUCcutoff=0.6 #Require model training AUC of at least 0.6

if (Methfile=='NA'&SNPCELfile=='NA') {DataCombi=1} else #U133A
if (AffyCELfile=='NA'&SNPCELfile=='NA') {DataCombi=2} else #Methylation
if (AffyCELfile=='NA'&Methfile=='NA') {DataCombi=3} else #SNP6
if (SNPCELfile=='NA') {DataCombi=4} else #U133A+meth
if (Methfile=='NA') {DataCombi=5} else #U133A+SNP6
if (AffyCELfile=='NA') {DataCombi=6} else #meth+SNP6
{DataCombi=7} #U133A+meth+SNP6

###Determine path and name of script being run
R_script=sub("--file=", "", args.all[grep("--file=", args.all)])
R_script_name=strsplit(R_script,"/")[[1]][length(strsplit(R_script,"/")[[1]])] 
script_path=paste(getwd(),R_script_name,sep="/")
script_path=gsub("_",x=script_path,replacement="\\\\_") #Escapes "_" in path names

###Libraries
library(Biobase)
library(affy)
library(randomForest)
library(DNAcopy) #BioConductor
library(CNTools) #BioConductor
library("aroma.affymetrix") #hbInstall with source("http://aroma-project.org/hbLite.R")
library("aroma.cn") #hbInstall with source("http://aroma-project.org/hbLite.R")
library("gdata") # CRAN
library("R.cache")

###Set work directory
path='/gne/research/workspace/daemena/STM/Rtoolbox/'
setwd(path)
pathPatientFiles='TumorData/'

###Results files (note tex will be converted to pdf with pdflatex
txt_output=paste(PatientID,"_results.txt",sep="")
tex_output=paste(PatientID,"_results.tex",sep="")

###Load patient data and preprocess
# U133A - RMA, standard CDF (together with CEL files of the core cell line panel)
if (DataCombi==1|DataCombi==4|DataCombi==5|DataCombi==7) {
	affy.data1=ReadAffy(filenames=paste(pathPatientFiles,AffyCELfile,sep=""))
	setwd('CoreCELfiles/')
	affy.data2=ReadAffy()
	# Combine tumor CEL file and cell line CEL files
	affy.data=merge(affy.data1,affy.data2)
	rma.standard.data=rma(affy.data)
	exprSet=exprs(rma.standard.data)
	
	# Conversion of probe sets to genes (selection of most varying probe according to the cell line panel)
	GeneProbeMapping.rma=read.table(paste(path,'ExtraFiles/Neve_AffyRMA_gene_maxvarprobe_mapping_stringent.txt',sep=""),header=TRUE)
	data.U133A=exprSet[as.vector(GeneProbeMapping.rma[,2]),]
	data.U133A=as.data.frame(data.U133A[,1])
	rownames(data.U133A)=paste('U133A__',as.vector(GeneProbeMapping.rma[,1]),sep="")
	colnames(data.U133A)=PatientID
	setwd(path)
}
# Methylation - check values are within [0,1]
if (DataCombi==2|DataCombi==4|DataCombi==6|DataCombi==7) {
	MethData=read.table(paste(pathPatientFiles,Methfile,sep=""),sep="\t")
	rownames(MethData)=MethData[,1]
	
	# Name conversion of cg-probe ids to genes
	GeneProbeMapping.meth=read.csv(paste(path,'ExtraFiles/Methylation_annotation_stringent.csv',sep=""),sep=",",header=TRUE)
	data.Meth=as.data.frame(MethData[as.vector(GeneProbeMapping.meth[,1]),2])
	rownames(data.Meth)=GeneProbeMapping.meth[,1]
	rownames(data.Meth)=paste('Meth__',rownames(data.Meth),'__',as.vector(GeneProbeMapping.meth[,4]),sep="")
	colnames(data.Meth)=PatientID
}
# SNP6 processing data (CSB in DNAcopy), with summary at gene level (CNTools)
if (DataCombi==3|DataCombi==5|DataCombi==6|DataCombi==7) {
	#Load SNP6 processing functions
	source("SNP6_processing_functions.R")

	#Copy file from TumorData to SNP6 directory for processing
	SNPpath=paste(pathPatientFiles,SNPCELfile,sep="")
	file.copy(SNPpath, 'ExtraFiles/SNP6/rawData/patientdata/GenomeWideSNP_6/')
	setwd('ExtraFiles/SNP6/')

	#load in gene mapping info to use for SNP6 gene-level summaries, based on build 36, HG18
	#Originally obtained from data("geneInfo") in library('CNTools')
	geneInfo=read.table(file="geneInfo_hg18.txt", sep="\t", header=TRUE, row.names=1, colClasses = c("character","character","numeric","numeric","character","character"))

	#Perform CBS segmentation on raw SNP6 data
	segdata_file=segment_data(dataSet="patientdata")

	#Generate gene-level summaries on segmented data
	segdatagene_file=get_gene_level_summaries(segmentedFile=segdata_file, geneMapInfo=geneInfo)

	#Load gene-level SNP6 copy number values
	SNP6Data=read.table(file=segdatagene_file, header=TRUE, sep="\t")
	GeneProbeMapping.SNP=read.csv(paste(path,'ExtraFiles/SNP6/SNP6_genelevel_stringent_std0.7.csv',sep=""),sep=",",header=TRUE)
	data.SNP=as.data.frame(SNP6Data[which(SNP6Data[,5] %in% GeneProbeMapping.SNP[,"GeneSymbol"]),6])
	rownames(data.SNP)=paste('SNP__',SNP6Data[which(SNP6Data[,5] %in% GeneProbeMapping.SNP[,"GeneSymbol"]),5],sep="")
	colnames(data.SNP)=PatientID

	#Remove intermediary files from rawData/patientdata, probeData/, plmData/ and cbsData/ (but leave HapMap_Controls files)
	unlink("plmData/patientdata,ACC,ra,-XY,BPN,-XY,AVG,A+B", recursive=TRUE)
	unlink("plmData/patientdata,ACC,ra,-XY,BPN,-XY,AVG,A+B,FLN,-XY", recursive=TRUE)
	unlink("cbsData/patientdata,ACC,ra,-XY,BPN,-XY,AVG,A+B,FLN,-XY,paired", recursive=TRUE)
	unlink("probeData/patientdata,ACC,ra,-XY", recursive=TRUE)
	unlink("probeData/patientdata,ACC,ra,-XY,BPN,-XY", recursive=TRUE)
	unlink(paste("rawData/patientdata/GenomeWideSNP_6/",SNPCELfile,sep=""))
	unlink("cbsData/patientdata,ACC,ra,-XY,BPN,-XY,AVG,A+B,FLN,-XY,paired", recursive=TRUE)
	setwd(path)
}

setwd(paste(path,'Results/',sep=""))

### Testing of best predictor per drug compound on tumor sample (best predictor defined upon the inputed data combination) for set of 91 drug compounds
drugs_interest=read.table(paste(path,'ExtraFiles/DrugCompounds.txt',sep=""),header=TRUE,sep="\t")
predictor_perf = data.frame(cbind(drugs_interest, Model=NA, ModelAUC=NA, DataCombination=NA, Score=NA, Probability=NA, Sensitivity=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drugs_interest[[1]]

### Mapping table with best predictor per drug compound
### Drug x Data Combi matrix with best predictor given a data combi (might be based on only 1 data set, despite of available other data types)
predictor_mapping=read.table(paste(path,'ExtraFiles/PredictorsMapping.txt',sep=""),sep="\t",header=TRUE,row.names=1)
predictor_mapping_AUC=read.table(paste(path,'ExtraFiles/PredictorsMapping_AUC.txt',sep=""),sep="\t",header=TRUE,row.names=1)


for (drug_compound in drugs_interest[[1]]) {
	drug_compound=as.vector(drug_compound)
	model=as.character(predictor_mapping[drug_compound,DataCombi])
	modelAUC=predictor_mapping_AUC[drug_compound,DataCombi]
	if (regexpr('_U133A_Meth_SNP_',model)[1]>0) {
		CombiData_LSSVM=rbind(data.U133A,data.Meth,data.SNP)
		medianU133A=median(data.U133A[!is.na(data.U133A)])
		data.U133A[is.na(data.U133A)]=medianU133A
		medianMeth=median(data.Meth[!is.na(data.Meth)])
		data.Meth[is.na(data.Meth)]=medianMeth
		medianSNP=median(data.SNP[!is.na(data.SNP)])
		data.SNP[is.na(data.SNP)]=medianSNP
		CombiData_RF=rbind(data.U133A,data.Meth,data.SNP)
		DataCombiModel="U133A-Meth-SNP"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	} else if (regexpr('_U133A_Meth_',model)[1]>0) {
		CombiData_LSSVM=rbind(data.U133A,data.Meth)
		medianU133A=median(data.U133A[!is.na(data.U133A)])
		data.U133A[is.na(data.U133A)]=medianU133A
		medianMeth=median(data.Meth[!is.na(data.Meth)])
		data.Meth[is.na(data.Meth)]=medianMeth
		CombiData_RF=rbind(data.U133A,data.Meth)
		DataCombiModel="U133A-Meth"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(1,dim(data.Meth)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	} else if (regexpr('_U133A_SNP_',model)[1]>0) {
		CombiData_LSSVM=rbind(data.U133A,data.SNP)
		medianU133A=median(data.U133A[!is.na(data.U133A)])
		data.U133A[is.na(data.U133A)]=medianU133A
		medianSNP=median(data.SNP[!is.na(data.SNP)])
		data.SNP[is.na(data.SNP)]=medianSNP
		CombiData_RF=rbind(data.U133A,data.SNP)
		DataCombiModel="U133A-SNP"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(10,dim(data.SNP)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	} else if (regexpr('_Meth_SNP_',model)[1]>0) {
		CombiData_LSSVM=rbind(data.Meth,data.SNP)
		medianMeth=median(data.Meth[!is.na(data.Meth)])
		data.Meth[is.na(data.Meth)]=medianMeth
		medianSNP=median(data.SNP[!is.na(data.SNP)])
		data.SNP[is.na(data.SNP)]=medianSNP
		CombiData_RF=rbind(data.Meth,data.SNP)
		DataCombiModel="Meth-SNP"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)	} else if (regexpr('_U133A_',model)[1]>0) {
		CombiData_LSSVM=data.U133A
		medianU133A=median(data.U133A[!is.na(data.U133A)])
		data.U133A[is.na(data.U133A)]=medianU133A
		CombiData_RF=rbind(data.U133A)
		DataCombiModel="U133A"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	} else if (regexpr('_U133A_',model)[1]>0) {
                CombiData_LSSVM=data.U133A
                medianU133A=median(data.U133A[!is.na(data.U133A)])
                data.U133A[is.na(data.U133A)]=medianU133A
                CombiData_RF=rbind(data.U133A)
                DataCombiModel="U133A"
                rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
                rangeClKernel[1,]=c(rep(1,dim(data.U133A)[1]))
                colnames(rangeClKernel)=rownames(CombiData_LSSVM)
        } else if (regexpr('_Meth_',model)[1]>0) {
		CombiData_LSSVM=data.Meth
		medianMeth=median(data.Meth[!is.na(data.Meth)])
		data.Meth[is.na(data.Meth)]=medianMeth
		CombiData_RF=rbind(data.Meth)
		DataCombiModel="Meth"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(1,dim(data.Meth)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	} else if (regexpr('_SNP_',model)[1]>0) {
		CombiData_LSSVM=data.SNP
		medianSNP=median(data.SNP[!is.na(data.SNP)])
		data.SNP[is.na(data.SNP)]=medianSNP
		CombiData_RF=rbind(data.SNP)
		DataCombiModel="SNP"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(10,dim(data.SNP)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	}
	
	if (regexpr('LSSVM',model)[1]>0) {
		ModelDrug=read.table(paste(path,'ModelsPlatformSpecific/',model,'.txt',sep=""),as.is=1)
		TrainDrug=read.table(paste(path,'ModelsPlatformSpecific/Trainingdata_',model,'.txt',sep=""),header=TRUE,sep="\t",row.names=1)
		Asigm=ModelDrug[[1]][1]
		Bsigm=ModelDrug[[1]][2]
		bcoeff=ModelDrug[[1]][3]
		alpha=ModelDrug[[1]][-(1:3)]
		target=as.matrix(TrainDrug[1,])
		traindata=as.matrix(TrainDrug[-1,])
		FeatureSubset=rownames(TrainDrug[-1,])
		CombiSubData=as.data.frame(CombiData_LSSVM[FeatureSubset,])
		rownames(CombiSubData)=FeatureSubset
		rangeClKernelSub=rangeClKernel[1,FeatureSubset]
		
		### Reduction to model features that are present in the input data file(s)
		removal=which(is.na(CombiSubData[,1]))
		if (length(removal)>0) {	
			traindata=traindata[-removal,]
			FeatureSubset=FeatureSubset[-removal]
			CombiSubData=as.data.frame(CombiSubData[-removal,])
			rownames(CombiSubData)=FeatureSubset
			rangeClKernel=rangeClKernelSub[FeatureSubset]
		}
		
		Kfeature=0
		clinicalKernel <- function(i) {
			clinicalKernelIntern <- function(j) {
				Kfeature[j]=(rangeClKernel[i]-abs(traindata[i,j]-CombiSubData[i,]))
			}
			Kfeature=sapply(c(1:dim(traindata)[2]),clinicalKernelIntern)
			Kfeature=Kfeature/rangeClKernel[i]
		}
		KtestPerFeature=sapply(c(1:length(FeatureSubset)),clinicalKernel)
		Ktest=rowSums(KtestPerFeature)/length(FeatureSubset)
		#Ktest=t(traindata) %*% as.matrix(CombiSubData)
		latent=((alpha*target) %*% Ktest) + bcoeff
		predictor_perf[drug_compound,"Score"]=format(latent,digits=5)
		prob=1/(1+exp(Asigm*latent+Bsigm))
		predictor_perf[drug_compound,"Probability"]=format(prob,digits=5)
		predictor_perf[drug_compound,"Sensitivity"]=(prob>0.5)
		predictor_perf[drug_compound,"Model"]="LSSVM"
		predictor_perf[drug_compound,"DataCombination"]=DataCombiModel
		predictor_perf[drug_compound,"ModelAUC"]=modelAUC
	} else if (regexpr('RF',model)[1]>0) {
		RFdata=t(CombiData_RF)
		load(file=paste(path,'ModelsPlatformSpecific/',model,'.Rdata',sep=""))
		RF_predictions_response=predict(rf_model, RFdata, type="response")
		RF_predictions_prob=predict(rf_model, RFdata, type="prob")
		RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
		predictor_perf[drug_compound,"Score"]=RF_predictions_vote[1,"sensitive"]
		predictor_perf[drug_compound,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
		predictor_perf[drug_compound,"Sensitivity"]=RF_predictions_response=="sensitive"
		predictor_perf[drug_compound,"Model"]="RF"
		predictor_perf[drug_compound,"DataCombination"]=DataCombiModel
		predictor_perf[drug_compound,"ModelAUC"]=modelAUC
	}
}

### Convert Sensitivity=TRUE/FALSE to +/-
predictor_perf[which(predictor_perf[,"Sensitivity"]==TRUE),"Sensitivity"]="+"
predictor_perf[which(predictor_perf[,"Sensitivity"]==FALSE),"Sensitivity"]="-"

#Apply cutoffs to limit results to only those drugs with sufficient input data
predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Score"]=NA
predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Probability"]=NA
predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Sensitivity"]=NA

### Order results in decreasing predicted probability of response
prob_order=order(predictor_perf[,"Probability"],decreasing=TRUE)
predictor_perf=predictor_perf[prob_order,]

### Save results to txt file
write.table(predictor_perf, file=txt_output, sep="\t", row.names=FALSE)

### Determine "recommended treatment"
if (length(which(!is.na(predictor_perf[,"Probability"])))>0){
	recommended_treatment=rownames(predictor_perf)[which(predictor_perf[,"Probability"]==max(predictor_perf[,"Probability"], na.rm=TRUE))]
}else{
	recommended_treatment="NA"
}

# Format results for print to pdf
predictor_perf[,"Score"]=format(predictor_perf[,"Score"], digits=3)
predictor_perf[,"Probability"]=format(as.numeric(predictor_perf[,"Probability"]), digits=1)
predictor_perf[,"ModelAUC"]=format(as.numeric(predictor_perf[,"ModelAUC"]), digits=3)


### Create latex/pdf report
sink(tex_output)
cat ("\\documentclass{article}\n")
cat ("\\usepackage{amsmath}\n")
cat ("\\usepackage{amscd}\n")
cat ("\\usepackage[tableposition=top]{caption}\n")
cat ("\\usepackage{ifthen}\n")
cat ("\\usepackage[utf8]{inputenc}\n")
cat ("\\usepackage[hmargin=3cm,vmargin=3cm]{geometry}\n")
cat ("\\usepackage{hyperref}\n")
cat ("\\hypersetup{\n")
cat ("   colorlinks = true,\n")
cat ("   linkcolor = blue}\n")

cat ("\\begin{document}\n")

#Page 1
cat ("\\title{Treatment Response Prediction Report}\n")
cat ("\\author{Anneleen Daemen and Obi Griffith (PI: Spellman/Gray)}\n")
cat ("\\maketitle\n")

cat ("\\subsubsection*{Sample Information:}\n")
cat ("Affymetrix expression CEL File:",AffyCELfile_clean,"\n")
cat ("\n")
cat ("Methylation File:",Methfile_clean,"\n")
cat ("\n")
cat ("Affymetrix SNP6 CEL File:",SNPCELfile_clean,"\n")

cat ("\\subsubsection*{Recommended Treatment:}\n")
cat (recommended_treatment,"\n")

cat ("\\subsubsection*{Treatment Response Predictions:\\footnote[1]{Patient considered sensitive to treatment if probability is greater than 0.5}}\n")
cat ("\\begin{tabular}{| l || l | c | l | c | c | c |}\n")
cat ("\\hline\n")

#Create results table:
cat ("Treatment & Model & AUC & Preprocessing & Score & Probability & Sensitivity \\\\ \\hline \n")
#Loop through and print row for each drug
for (i in 1:25){
drug_clean=gsub("_",x=predictor_perf[i,1],replacement="\\\\_")#Escapes "_" in drug names
cat (drug_clean," & ",predictor_perf[i,2]," & ",predictor_perf[i,3]," & ",predictor_perf[i,4]," & ",predictor_perf[i,5]," & ",predictor_perf[i,6]," & ",predictor_perf[i,7]," \\\\ \\hline \n")
}
cat ("\\end{tabular}\n")
cat ("\\footnotetext[2]{Code location:",script_path,"}\n")

#Page 2
cat ("\\subsubsection*{Treatment Response Predictions (continued):\\footnote[1]{Patient considered sensitive to treatment if probability is greater than 0.5}}\n")
cat ("\\begin{tabular}{| l || l | c | l | c | c | c |}\n")
cat ("\\hline\n")

#Create results table:
cat ("Treatment & Model & AUC & Preprocessing & Score & Probability & Sensitivity \\\\ \\hline \n")
#Loop through and print row for each drug
for (i in 26:67){
drug_clean=gsub("_",x=predictor_perf[i,1],replacement="\\\\_")#Escapes "_" in drug names
cat (drug_clean," & ",predictor_perf[i,2]," & ",predictor_perf[i,3]," & ",predictor_perf[i,4]," & ",predictor_perf[i,5]," & ",predictor_perf[i,6]," & ",predictor_perf[i,7]," \\\\ \\hline \n")
}
cat ("\\end{tabular}\n")

#Page 3
cat ("\\subsubsection*{Treatment Response Predictions (continued):\\footnote[1]{Patient considered sensitive to treatment if probability is greater than 0.5}}\n")
cat ("\\begin{tabular}{| l || l | c | l | c | c | c |}\n")
cat ("\\hline\n")

#Create results table:
cat ("Treatment & Model & AUC & Preprocessing & Score & Probability & Sensitivity \\\\ \\hline \n")
#Loop through and print row for each drug
for (i in 68:length(rownames(predictor_perf))){
drug_clean=gsub("_",x=predictor_perf[i,1],replacement="\\\\_")#Escapes "_" in drug names
cat (drug_clean," & ",predictor_perf[i,2]," & ",predictor_perf[i,3]," & ",predictor_perf[i,4]," & ",predictor_perf[i,5]," & ",predictor_perf[i,6]," & ",predictor_perf[i,7]," \\\\ \\hline \n")
}
cat ("\\end{tabular}\n")
cat ("\\end{document}\n")

sink()
system(paste('pdflatex -interaction=batchmode',tex_output,'1>/dev/null')) 

#Delete unnecessary files. If pdf file not produced as you expect, comment these out
tex_aux=paste(PatientID,"_results.aux",sep="")
tex_out=paste(PatientID,"_results.out",sep="")
tex_log=paste(PatientID,"_results.log",sep="")
system(paste('rm',tex_aux))
system(paste('rm',tex_out))
system(paste('rm',tex_log))
system(paste('rm',tex_output))
