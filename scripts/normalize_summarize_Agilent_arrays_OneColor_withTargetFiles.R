# Process Agilent text files for one-color-array and do between array normalization.
# produce G-norm files
# Albert Cheng
#
#
#
# R --vanilla < normalize_summarize_Agilent_arrays_OneColor_withTargetFiles.R
#
# or in R prompt:
# 
# source('normalize_summarize_Agilent_arrays_OneColor_withTargetFiles.R')
#

origPath=getwd()
# Go to directory with *.txt files from scanner
setwd("../raw")


#spotTypeFileName="SpotTypes.txt"
targetFileName="targets.txt"

outputPrefix="../limmaOutput"
withControlPath="withControl"

# This package needed
library(limma)

targets=readTargets(targetFileName)
#spottypes <- readSpotTypes(spotTypeFileName)
#scanFiles = dir(pattern = ".*.txt$") #not using simple txt pattern match, use target file!~
num.arrays = length(targets$FileName)

# Only reads in the G values, fake the R values as G as well

rawObj <- read.maimages(targets,source="agilent",green.only=TRUE)
print("done reading")

print("background correct")
obj.corrected<-backgroundCorrect(rawObj,method="normexp",offset=1)
E <- normalizeBetweenArrays(obj.corrected,method="quantile")
colnames(E$E)=targets$Cy3

#does it have a batch vector?
if(length(targets$batch)>0){
	E.batchEffectRemoved<-removeBatchEffect(E,targets$batch)
	
}



# Go back to starting directory
unlink(outputPrefix,recursive=TRUE)


dir.create(outputPrefix)
# Go to output directory
setwd(outputPrefix)





#write tables
if(length(targets$batch)>0){
	write.table(cbind(E$genes, E.batchEffectRemoved), file="exp.norm.ber.txt", sep="\t", quote=FALSE, row.names=F)
	write.table(cbind(GeneName=E$genes$GeneName,E.batchEffectRemoved), file="exp.norm.ber.gn.txt", sep="\t", quote=FALSE, row.names=F)
}
write.table(cbind(E$genes, E$E), file="exp.norm.txt", sep="\t", quote=FALSE, row.names=F)
write.table(cbind(GeneName=E$genes$GeneName, E$E), file="exp.norm.gn.txt", sep="\t", quote=FALSE, row.names=F)

#go back to wherever we start off.
setwd(origPath)
