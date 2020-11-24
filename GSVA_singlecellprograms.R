#Script to perform single sample GSEA using GSVA in R
#Here I'm using ssGSEA to score human bulk RNAseq glioma samples for single cell derived signatures

library(tidyverse)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(stringr)
library(Biobase)
library(limma)


##Read in the single cell signatures that you're interested in and convert to human orthologs
#Each column is a signature
Tumor<-read.table(file="SingleCell_Signatures", header=TRUE, sep="\t", stringsAsFactors = FALSE)

#Convert genes to the human homologs, here using a table with 2 columns containing
#human and mouse gene symbols
anno <- "human.mouse.genenames.csv"
anno.table <- read.csv(anno, header=TRUE)
anno.table$mouse <- as.character(anno.table$mouse)
anno.table$human <- as.character(anno.table$human)
anno.table <- anno.table %>% dplyr::select(mouse,human)

#For loop to run through each signature and convert mouse gene symbols to human gene symbols
newtable<-list()
for (i in seq_len(ncol(Tumor))){
  sign<-names(Tumor)[i]
  genes<-as.character(Tumor[,i])
  human<-anno.table %>% filter(mouse %in% genes) %>% dplyr::select(human)
  vector<-human$human
  vector<-c(vector,rep(NA, 367-length(vector)))
  newtable[[sign]]<-vector
}

#Convert the individual lists (each signature) from for loop to dataframe:
Sign<-as.data.frame(newtable)

#Convert any columns that are factors into characters
Sign<-Sign %>% mutate_if(is.factor, as.character)


#Remove any signatures you don't want to use for ssGSEA step
#Here I removed the cyclying signatures and the mixed lineage signatures
Sign2<-Sign[, -c(3,4,6,9)]

##Load in the bulk RNAseq metadata and expression data

#Meta data for UCSC Treehouse data
Trk<-read.csv("TrkAlteredSamples.csv", sep=",", header=T, stringsAsFactors = FALSE)

#To match sample names in expression dataframe, replace "-" with "."
Trk$th_sampleid<-gsub("-", ".", Trk$th_sampleid)

#Remove samples with large scale CNAs that include Trk amplifications:
Trk<-Trk %>% filter(!(site_donor_id %in% c("TCGA-FG-A6J3", "TCGA-HT-7677", "TCGA-HT-A61B"))) 
#24 tumors remaining with Trk amplifications or fusions

#Read in expression (logCPM) data:
CPM<-read.csv("UCSC_TrkGlioma_CPM.csv", header=T, sep=",", stringsAsFactors = FALSE)
CPM[1:5,1:5]
Ensembl<-CPM$X

#Filter expression data to samples with Trk amp or fusion
CPM<-CPM[,Trk$th_sampleid]
CPM$Ens<-Ensembl

#Add the gene symbol to the Ensembl IDs
df<-read.csv("Gencode22_Ensembl_Symbol.csv", stringsAsFactors = FALSE, header=TRUE, row.names=1)
results1<-dplyr::left_join(CPM, df)

results1$Symbol<-make.unique(results1$Symbol)
row.names(results1)<-results1$Symbol
results1[1:5,1:5]
results1[1:5, 23:26]  #columns 25 adn 26 are the Ensembl IDs and gene symbol

#Get samples in same order as the metadata:
results<-results1[, match(Trk$th_sampleid, names(results1))]
#convert to matrix
expr.matrix<-as.matrix(results)


#Put data into an ExpressionSet:
#Phenotype data:
meta<-data.frame(Grade=Trk$Grade, stringsAsFactors = FALSE)
meta$Age<-Trk$Age
meta$Comp<-Trk$DiseaseType
meta$Comp<-ifelse(meta$Grade=="G3", "LGG", meta$Comp)
row.names(meta)<-Trk$th_sampleid

#Make sure the rownames of meta and column names of expression matrix match
all(row.names(meta)==colnames(expr.matrix))

phenos<-colnames(meta)
phenos


## descriptive data for phenotype columns, for now just put in labels of our meta data
pheno.desc<-phenos
metadata<-data.frame(labelDescription=
                         pheno.desc,
                       row.names=phenos)

phenoData<-new("AnnotatedDataFrame",
                 data=meta, varMetadata=metadata)

phenoData

#create the Expression Set using BioBase
ES<-ExpressionSet(assayData=expr.matrix,  
                         phenoData=phenoData)


##ssGSEA step:
#Input expression set (or input matrix of expression values with genes as rows
#and samples as columns), signatures of interest, method (options=gsva, ssgsea, zscore or plage),
#ssgsea.norm, here I'm skipping the last normlization step

ssGSEA <- gsva(ES, Sign2, method="ssgsea", ssgsea.norm=FALSE, verbose=FALSE) 

#dataframe of ssGSEA results
data<-data.frame(ssGSEA)
#Write results to csv file
write.csv(data, "ssGSEA_SingleCellSignatures_results.csv")


#Scores for each sample
Heat<-data[,1:5]

#Transpose df
Heat<-data.frame(t(Heat))

#Annotation for heatmap
Annotation<-data.frame(Grade=Trk$Grade)
row.names(Annotation)<-Trk$th_sampleid
Annotation<-Annotation %>% arrange(desc(Grade))
Order<-row.names(Annotation)

#Annotation colors
ann_color = list(Grade = c(G4 = "goldenrod", G3 = "purple", G2="darkgreen"))

#Rearrange columns, want columns from highest grade to lowest, use dput to first print current order as a vector:
dput(Order)
#neworder<-here I created a vector for the order of the samples in the heatmap

#Order dataframe 
Heat<-Heat[,(match(neworder, names(Heat)))]

pheatmap(Heat, scale="row",
         annotation_col = Annotation,
         annotation_colors=ann_color,
         cluster_cols = F,
         color=colorRampPalette(c("dodgerblue", "white", "red3"))(250),
         show_colnames = TRUE)

#Find significant differences between ssGSEA scores in HGG vs LGG:
adjPvalueCutoff <- 0.1
design <- model.matrix(~factor(ssGSEA$Comp))
colnames(design) <- c("Int", "LGGvHGG")
fit <- lmFit(ssGSEA, design)
fit <- eBayes(fit)

#Differential gene sets:
DEgeneSets <- topTable(fit, coef="LGGvHGG", number=Inf,
                        p.value=adjPvalueCutoff, adjust="BH")


#Visualize ssGSEA scores for the signatures with significant differences between LGG and HGG
#Violin plots:
ggplot(data, aes(x=Grade, y=SetA))+geom_violin()+geom_point()
ggplot(data, aes(x=Grade, y=SetB))+geom_violin()+geom_point()


#Ridge plots:
library(ggridges)
ggplot(data, aes(x=SetA, y=Grade, fill=Grade)) +
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")

ggplot(data, aes(x=SetB, y=Grade, fill=Grade)) +
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")

ggplot(data, aes(x=SetC, y=Grade, fill=Grade)) +
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")

