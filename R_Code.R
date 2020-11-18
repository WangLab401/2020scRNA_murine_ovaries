
###Scripts of data progressing for single-cell transcriptome analyse during primordial follicle assembly in mice.

########  Single-Library Analysis with Cell Ranger in server
########  CellRanger version=3.1.0

nohup /home/.../10xgenomics/cellranger-3.1.0/cellranger count --id=SampleName \
--localcores 20 \
--transcriptome=/home/.../10xgenomics/refdata-cellranger-mm10-3.0.0 \
--fastqs=/home/.../Singlecell_Rawdata   \
--sample=SampleName \
--chemistry=threeprime \
--force-cells=6000 &


#############################################################################################
         #########///////=========Seurat version=3.1.5========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)

###Load Data
Seurat_Object <- Read10X(data.dir = "D:\\...\\SampleName\\filtered_feature_bc_matrix\\")

Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = "Seurat_Object", min.cells = 3, min.features = 200)

Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object, pattern = "^mt-")

VlnPlot(Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

###Value setting of filteration
Seurat_Object <- subset(Seurat_Object, subset = nFeature_RNA > value1 & nFeature_RNA < value2)

###======Perform integration
pancreas.anchors <- FindIntegrationAnchors(object.list = list(E16_5, PD0, PD3), dims = 1:30)

PF.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = PF.integrated) <- "integrated"

PF.integrated <- ScaleData(PF.integrated, verbose = FALSE)
PF.integrated <- RunPCA(PF.integrated, npcs = 30, verbose = FALSE)
PF.integrated <- RunUMAP(PF.integrated, reduction = "pca", dims = 1:20)

#### t-SNE and Clustering

PF.integration <- FindNeighbors(PF.integrated, reduction = "pca", dims = 1:20)

PF.integration <- FindClusters(PF.integration, resolution = 0.3)
table(PF.integration@meta.data$seurat_clusters)

p1 <- DimPlot(PF.integration, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(PF.integration, reduction = "umap", group.by = "seurat_clusters")
plot_grid(p1, p2)

###find markers
DefaultAssay(PF.integration) <- "RNA"
PF.integration.All.Markers <- FindAllMarkers(PF.integration, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25

###subset single cell type

cell_type <- subset(PF.integration, idents = c("cell_type"))
saveRDS(cell_type, file = "D:\\...\\cell_type.rds")


########====== single cell type with seurat

Seurat_Germs <- readRDS(file = "D:\\...\\Germ_cells.rds")
Seurat_Germs

Seurat_Germs[["percent.mt"]] <- PercentageFeatureSet(Seurat_Germs, pattern = "^MT-")
VlnPlot(Seurat_Granulosa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Seurat_Germs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_Germs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Germ_cells <- subset(Seurat_Germs, subset = nFeature_RNA > 200 & nFeature_RNA < ValueX )

######## Data processing

Germ_cells <- FindVariableFeatures(Germ_cells, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Germ_cells), 10)

all.genes <- rownames(Germ_cells)
Germ_cells <- ScaleData(Germ_cells, features = all.genes)

###Perform linear dimensional reduction
Germ_cells <- RunPCA(Germ_cells, features = VariableFeatures(object = Germ_cells))

print(Germ_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Germ_cells, dims = 1:2, reduction = "pca")

DimPlot(Germ_cells, reduction = "pca")
DimHeatmap(Germ_cells, dims = 1, cells = 500, balanced = TRUE)

###Determine the ‘dimensionality’ of the dataset
Germ_cells <- JackStraw(Germ_cells, num.replicate = 100)
Germ_cells <- ScoreJackStraw(Germ_cells, dims = 1:20)

JackStrawPlot(Germ_cells, dims = 1:15)
ElbowPlot(Germ_cells)

Germ_cells <- FindNeighbors(Germ_cells, dims = 1:10)
Germ_cells <- FindClusters(Germ_cells, resolution = 0.2)

#####Run non-linear dimensional reduction (UMAP/tSNE)

Germ_cells <- RunUMAP(Germ_cells, dims = 1:10)

##Marker
Germ_cells_cluster.markers <- FindAllMarkers(Germ_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



#############################################################################################
         #########///////=========Monocle version=2.10.1========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Monocle
### http://cole-trapnell-lab.github.io/monocle-release/docs/

library("Seurat")
library('monocle')
Subset_Germs <- readRDS(file = "D:\\...\\Subset_Germs_seurat.rds")
Subset_Germs

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Subset_Germs@assays$RNA@data), 'sparseMatrix')
dim(data)
head(data)[1:5,1:5]
pd <- new('AnnotatedDataFrame', data = Subset_Germs@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
subGerm_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


dim(subGerm_cds)
head(fData(subGerm_cds))
head(pData(subGerm_cds))

dim(exprs(subGerm_cds))

###Estimate size factors and dispersions

subGerm_cds <- estimateSizeFactors(subGerm_cds)
subGerm_cds <- estimateDispersions(subGerm_cds)


##Trajectory step 1: choose genes that define a cell's progress
## Select genes that differ between clusters/stages
diff_test_res <- differentialGeneTest(subGerm_cds,fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-5))

subGerm_cds <- setOrderingFilter(subGerm_cds, order_genes)
plot_ordering_genes(subGerm_cds)

###Trajectory step 2: reduce data dimensionality
subGerm_cds <- reduceDimension(subGerm_cds,max_components = 2,method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
subGerm_cds <- orderCells(subGerm_cds)
plot_cell_trajectory(subGerm_cds, color_by = "State")

##########BEAM Function
subGerm_cds <- orderCells(subGerm_cds,root_state = 2)
BEAM_res <- BEAM(subGerm_cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

library("colorRamps")
library("RColorBrewer")
hmcols<-colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(65)
plot_genes_branched_heatmap(subGerm_cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                         branch_point = 1,
                                         num_clusters = 4,
                                         cores = 1,
                                         hmcols = hmcols,
                                         use_gene_short_name = F,
                                         branch_colors = c("#bebebe","#009E73", "indianred2"),
                                         show_rownames = F,return_heatmap = T)



#############################################################################################
         #########///////=========SCENIC (version 0.9.1)========\\\\\\\#############
#############################################################################################
#####load the dataset 
library(Seurat)
library(monocle)

setwd("~/R")
subGerm_cds <- readRDS(file = "subGerm_cds.rds")
expr_Mat <- as.matrix(exprs(subGerm_cds))

cellInfo <- pData(subGerm_cds)
cellInfo$nGene <- colSums(expr_Mat>0)
cellInfo <- data.frame(cellInfo)
dim(cellInfo)


dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(State = c("1"="springgreen3",
                          "2"="goldenrod", 
                          "3"="violetred2"))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$State, legend=names(colVars$State))


##Initialize SCENIC settings
library(SCENIC)
org="mgi" # or hgnc, or dmel
dbDir="/home/.../crisTarget_databases"  # RcisTarget databases location
myDatasetTitle="SCENIC of mouse germ cell" # choose a name for your analysis

scenicOptions <- initializeScenic(org="mgi",dbDir="/.../crisTarget_databases",
                                  datasetTitle="SCENIC of mouse germ cell",
                                  nCores=30) 

mm10_dbs <- list('500bp' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather',
                 '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
db_mcVersion <- 'v9'
db_path <- '/home/.../crisTarget_db'
scenicOptions@settings$dbs <- mm10_dbs
scenicOptions@settings$dbDir <- db_path
scenicOptions@settings$db_mcVersion <- db_mcVersion

# Modify if needed
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#########Gene filter/selection

genesKept <- geneFiltering(expr_Mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(expr_Mat),
                           minSamples=ncol(expr_Mat)*.01)

######The following steps are executed according to SCENIC pipeline with default parameters. 
#####Please refer to the url below: 
### https://github.com/aertslab/SCENIC
### https://rawcdn.githack.com/aertslab/SCENIC/6aed5ef0b0386a87982ba4cc7aa13db0444263a6/inst/doc/SCENIC_Running.html



