if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")

# Install ensembldb and EnsDb.Hsapiens.v79
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v79")

library(tidyverse)
library(dplyr) #DATA wrangling filtering, selecting and summarising
library(Seurat)
library(patchwork)#combines multiple ggplot2 plots into a single plot
library(GEOquery)
library(ggplot2)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ensembldb)
library(EnsDb.Hsapiens.v79)
library(SingleR)
library(celldex)
#install.packages('devtools') #try this after installing rtools
#devtools::install_github('immunogenomics/presto')

#bcells <- readRDS("C:/Users/Ananya/Desktop/SJSU/Stanford/SLE_data/seurat_B_PB_scRNA_SLE_HC_GSE174188.rds")
# bcells <- readRDS("/home/agupta1/SLE/seurat_B_PB_scRNA_SLE_HC_GSE174188.rds")

bcells <- readRDS("C:/Users/Ananya/Desktop/SJSU/Stanford/SLE/SLE_data/bcells_subset.rds")
print(bcells$RNA) # 21429 features for 20000 cells
VariableFeatures(bcells)
levels(bcells) # just 9 clusters
slotNames(bcells)


#just exploring it
as_tibble(bcells@meta.data)
#rownames(bcells@assays$RNA)genes



########### getting gene names ############


#to check for duplicate genes
which(duplicated(rownames(bcells@assays$RNA))) #none from here, duh
#going to make the genes unique


#barely any were labelled with Annotatedbi
#im tryinga new suggestion from https://support.bioconductor.org/p/87454/


# Load the EnsDb object
edb <- EnsDb.Hsapiens.v79
#str(edb, max.level = 3)

# Extract gene annotations
gene_annotations <- genes(edb, return.type = "data.frame")

#ENSEMBL ID and gene symbol
gene_mapping <- gene_annotations %>%
  dplyr::select(gene_id, gene_name)

head(gene_mapping) 

#getting current genes
features <- rownames(bcells@assays$RNA)

#getting genes with ensg
ensembl_ids <- features[grep("^ENSG", features)]

# Map ENSEMBL IDs to gene symbols using the EnsDb mapping
#find index of ensembl_id in the gene_id(ensembl) col of gene_mapping
# extract the gene_name at that index and assign to gene_symbol
gene_symbols <- gene_mapping$gene_name[match(ensembl_ids, gene_mapping$gene_id)]

# Create a named vector for mapping
feature_mapping <- features #keys and values both are original feature names
names(feature_mapping) <- features #name is key and val



#i want to update only the mapped Ensmbl ids
for (x in seq_along(ensembl_ids)) {#iterate along its length
  if (!is.na(gene_symbols[x])) { #at the index if the symbol is not na
    if (gene_symbols[x] %in% feature_mapping) {
      # Append the Ensembl ID to make it unique
      feature_mapping[ensembl_ids[x]] <- paste0(gene_symbols[x], "_", ensembl_ids[x])
    } else {
      # Use the gene symbol as is
      feature_mapping[ensembl_ids[x]] <- gene_symbols[x]
    }
  }
}

head(feature_mapping)
head(rownames(bcells@assays$RNA))
# Update the feature names in the Seurat object
rownames(bcells@assays$RNA) <- feature_mapping[rownames(bcells@assays$RNA)]

#checking it and it worked !!!
head(rownames(bcells@assays$RNA))
print(bcells$RNA)
VariableFeatures(bcells) #it got updated here too YAYYYY

# i was getting duplicate row name errors
# some ensembl ids mapped to the same gene symbols as the ones already existing
# in the data

#i dont want to remove them
#1. is it good to merge them? is that even possible?
#2. remove?
#3. Im opting to append it to the ensg val. will it cause problems in singler?


# Check raw counts
head(LayerData(bcells, assay = "RNA", layer = "counts")[, 1:5])


#########fixing barcodes col names ################


#weird ahh barcode names
head(colnames(bcells), 10) #they look so bad ill kms
clean_cols <- gsub("-.*$", "", colnames(bcells)) #pattern,replace,item
clean_cols <- make.unique(clean_cols)
#colnames(bcells) <- clean_cols sseurat 5 not letting me
bcells <- RenameCells(bcells, new.names = clean_cols)
head(colnames(bcells), 10)
any(duplicated(colnames(bcells))) #falseeee





#####QC#############

#so im doing this to see if i still need to perform QC
#looks like I do
bcells[["percent.mt"]]<-PercentageFeatureSet(bcells, pattern = "^MT-")
summary(bcells@meta.data$percent.mt)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1002  2.6603  3.4204  3.7700  4.3843 59.1613


#post the 20k subset
# summary(bcells@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1002  2.6636  3.4286  3.7780  4.4026 53.8231 

#ribosomal genes
bcells[["percent.rb"]]<- PercentageFeatureSet(bcells, pattern = "^RP[SL]")
#well idk what this plot is tbh so ill just leave it 
#will check later if its affecting clustering


# Plot percent.rb vs. PC1 and PC2
FeaturePlot(bcells, features = "percent.rb", reduction = "pca")
# Add number of genes per UMI for each cell to metadata
bcells$log10GenesPerUMI <- log10(bcells$nFeature_RNA) / log10(bcells$nCount_RNA)

#numi = nCount_RNA,ngene=nFeature_RNA

## visualising the nCount_RNA/cell

plot_nCountqc <- bcells@meta.data %>%
  ggplot(aes(color=SLE_status, x=nCount_RNA, fill= SLE_status)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) #most cells have 500 or greater
  
plot_nCountqc  


#visualising genes detected per cell

plot_nFeaturesqc <- bcells@meta.data %>% 
  ggplot(aes(color=SLE_status, x=nFeature_RNA, fill= SLE_status)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200)+
  geom_vline(xintercept = 2500)
  ggtitle("plot_nFeaturesqc")#looks pretty filtered

plot_nFeaturesqc

#everything other than mitochondrial genes seems filtered
bcells <- subset(x = bcells, 
                 subset= (percent.mt < 5))




########### Using SingleR ###########

# first check is the data normalised
# Check raw counts
head(LayerData(bcells, assay = "RNA", layer = "counts")[, 1:5])
head(LayerData(bcells, assay = "RNA", layer = "data")[, 18:5]) #sparse so i cant tell


# Raw counts should be integers
summary(as.vector(LayerData(bcells, assay = "RNA", layer = "counts")))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0     0.1     0.0  6023.0 

# Normalized data should be float/log-transformed
summary(as.vector(GetAssayData(bcells, assay = "RNA", layer = "data")))  #log normalized
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.06888 0.00000 8.57346 


#### singler applying #########

#ribosome remove
#cell cycle
#singleR datat load immune data reference
#then manual
#eg marker for PB





###########plots######################

#top 10 variable features
top10 <- head(VariableFeatures(bcells),10)
plot1var <- VariableFeaturePlot(bcells)
plot1var_labeled <-LabelPoints(plot = plot1var, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
print(plot1var_labeled)

# Examine and visualize PCA results a few different ways
print(bcells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(bcells, dims = 1:2, reduction = "pca")


## Plot PCA
pca<-PCAPlot(bcells,
        split.by = "SLE_status",pt.size = 1.5, alpha = 0.8)
        ggtitle("PCA Plot of B Cells by SLE Status")
pca
umap_og <- DimPlot(bcells, reduction = "umap", label =TRUE)+
        ggtitle("UMAP Plot Initial")#clear separated clusters, lets look at their specific types
umap_og



#by clust
##gene names to be changed in the beginning

bcells.markers <- FindAllMarkers(bcells, only.pos = TRUE,logfc.threshold = 0.25)
write.csv(bcells.markers, file = "bcells_markers_subset.csv", row.names = FALSE)
#write.csv(bcells.markers, file = "bcells_markers.csv", row.names = FALSE)
#bcells.markers <- read.csv("bcells_markers.csv")

#compare to https://www.genome.jp/kegg/pathway.html, 
#http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarkerSearch.jsp?quickSearchInfo=CPVL2&index_key=2#framekuang,
#ncbi



#manual markers
clust0.markers <- bcells.markers %>%
  filter(cluster == 0,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust0.markers


clust1.markers <- bcells.markers %>%
  filter(cluster == 1,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust1.markers


clust2.markers <- bcells.markers %>%
  filter(cluster == 2,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust2.markers


clust3.markers <- bcells.markers %>%
  filter(cluster == 3,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust3.markers


clust4.markers <- bcells.markers %>%
  filter(cluster == 4,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust4.markers

clust5.markers <- bcells.markers %>%
  filter(cluster == 5,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust5.markers

clust6.markers <- bcells.markers %>%
  filter(cluster == 6,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust6.markers


clust7.markers <- bcells.markers %>%
  filter(cluster == 7,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust7.markers

clust8.markers <- bcells.markers %>%
  filter(cluster == 8,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust8.markers

clust2.markers <- bcells.markers %>%
  filter(cluster == 9,
         pct.1 > 0.712) %>%  # pct.1 above threshold
  arrange(desc(avg_log2FC),desc(pct.1))
clust9.markers

#PB-Ki67
#PB
#DN1
#DN2
#CD27
#FolB
#IFNhigh
#Naiveii
#Naivei


### final labels          0               1           2            3      4       5
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                     "NK", "DC", "Platelet")
#                     6       7     8        9
#names(new.cluster.ids) <- levels(bcells)
#bcells <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(bcells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()