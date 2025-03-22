if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

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

#install.packages('devtools') #try this after installing rtools
#devtools::install_github('immunogenomics/presto')

bcells <- readRDS("C:/Users/Ananya/Desktop/SJSU/Stanford/SLE_data/seurat_B_PB_scRNA_SLE_HC_GSE174188.rds")
# bcells <- readRDS("/home/agupta1/SLE/seurat_B_PB_scRNA_SLE_HC_GSE174188.rds")


set.seed(42)
tot_cells <- ncol(bcells)
sample_b <- sample(colnames(bcells), size=20000, replace=FALSE)
bcells <- subset(bcells,cells=sample_b)

levels(bcells) # just 9 clusters
slotNames(bcells)

print(bcells$RNA)# 21429 features for 152962 (now 20k) cells

#just exploring it
as_tibble(bcells@meta.data)
rownames(bcells@assays$RNA)
VariableFeatures(bcells) #cleaned data i think


########### getting gene names ############



# # Get the current feature names (ENSEMBL IDs and gene symbols)
# features <- rownames(bcells@assays$RNA)
# 
# # Extract ENSEMBL IDs (those starting with "ENSG")
# ensembl_ids <- features[grep("^ENSG", features)]
# 
# # Map ENSEMBL IDs to gene symbols
# gene_symbols <- mapIds(
#   org.Hs.eg.db,
#   keys = ensembl_ids,
#   keytype = "ENSEMBL",
#   column = "SYMBOL"
# )
# 
# # View the mapping
# head(gene_symbols)
# sum(is.na(gene_symbols)) #1085
# 
# str(gene_symbols)# tot is 1156

#barely any were labelled?? 
#im tryinga new suggestion from https://support.bioconductor.org/p/87454/


# Load the EnsDb object
edb <- EnsDb.Hsapiens.v79
str(edb, max.level = 3)

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
    feature_mapping[ensembl_ids[x]] <- gene_symbols[x]
  }
  
}

head(feature_mapping)

# Update the feature names in the Seurat object
rownames(bcells@assays$RNA) <- feature_mapping[rownames(bcells@assays$RNA)]
#checking it and it worked
head(rownames(bcells@assays$RNA))

VariableFeatures(bcells)

# Check raw counts
head(LayerData(bcells, assay = "RNA", layer = "counts")[, 1:5])



bcells[["pca"]] #why does it look like that is that barcodes

# Examine and visualize PCA results a few different ways
print(bcells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(bcells, dims = 1:2, reduction = "pca")




#####QC#############

#so im doing this to see if i still need to perform QC
#looks like I do
bcells[["percent.mt"]]<-PercentageFeatureSet(bcells, pattern = "^MT-")
summary(bcells@meta.data$percent.mt)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1002  2.6603  3.4204  3.7700  4.3843 59.1613

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

###########plots######################
## Plot PCA
pca<-PCAPlot(bcells,
        split.by = "SLE_status",pt.size = 1.5, alpha = 0.8)
        ggtitle("PCA Plot of B Cells by SLE Status")
pca
umap_og <- DimPlot(bcells, reduction = "umap", label =TRUE)+
        ggtitle("UMAP Plot Initial")#clear separated clusters, lets look at their specific types
umap_og




install.packages("biomaRt")
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# Get the ENSG IDs from your Seurat object
ensg_ids <- rownames(bcells)

# Query biomaRt for gene symbols
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensg_ids,
  mart = ensembl
)

#check
head(gene_annotations)

#named vector for mapping
gene_symbols <- gene_annotations$hgnc_symbol
names(gene_symbols) <- gene_annotations$ensembl_gene_id

#updating VariableFeatures
variable_features_ensg <- VariableFeatures(bcells)
variable_features_symbols <- gene_symbols[variable_features_ensg]
VariableFeatures(bcells) <- variable_features_symbols#None of the features specified are present in this assay
head(VariableFeatures(bcells))

#by clust
##gene names to be changed in the beginning

bcells.markers <- FindAllMarkers(bcells, only.pos = TRUE,logfc.threshold = 0.25)
write.csv(bcells.markers, file = "bcells_markers.csv", row.names = FALSE)
#bcells.markers <- read.csv("bcells_markers.csv")

#compare to https://www.genome.jp/kegg/pathway.html, 
#http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarkerSearch.jsp?quickSearchInfo=CPVL2&index_key=2#framekuang,
#ncbi

#ribosome remove
#cell cycle
#singleR datat load immune data reference
#then manual
#eg marker for PB

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