<div style="text-align:center"><img style="width:60%; height: auto" src="https://github.com/yanailab/enhance/raw/master/images/splash.jpg"/></div>

## ENHANCE: Accurate denoising of single-cell RNA-Seq data

This repository contains an R implementation of the ENHANCE algorithm for denoising single-cell RNA-Seq data ([Wagner et al., 2019](https://www.biorxiv.org/content/10.1101/655365v1)).

The **python implementation** can be found in a [separate repository](https://github.com/yanailab/enhance).

### Running ENHANCE

Follow these instructions to run the R implementation

1. Install dependencies

Make sure you have `rsvd` and `Matrix` installed. 

2. Download the GitHub repository

[Download ENHANCE](https://github.com/yanailab/enhance-R/archive/master.zip), and extract the contents into a folder.

3. To run enhance on a matrix

enhance(data_raw)

data_raw should be your raw count expression matrix, with genes as rows and cells as columns. 
All other arguments are optional. 
The output is the denoised count expression matrix, with genes as rows and cells as columns. 

    Arguments
    data_raw: The raw count expression matrix (rows=genes, columns=cells).
    ratio_pcs (optional): Variance ratio between simulated and real matrix
    to use to determined the threshold at which to keep principal components
    Default: 2
    k_nn (optional): Number of neighbors to aggregate
    Default: Calculated so that aggregation yields ~ target_transcripts per cell
    target_transcripts (optional): Number of transcripts per cell to aim for in aggregation
    Default: 2*10^5
    percent_cells_max (optional): Maximum percentage of cells to aggregate
    Default: 2

4. To run enhance on a Seurat object

enhance_seurat(object)

object should be a Seurat object containing an 'RNA' assay. 
All other arguments are optional. 
The output is the Seurat object with a new assay called "enhance", in which the "counts" slot is the denoised counts matrix

    Arguments:
    object: A Seurat object
    setDefaultAssay: Whether to set enhance as the default assay
    Default: TRUE
    assay: Which assay to run enhance on
    Default: "RNA"
    ratio_pcs (optional): Variance ratio between simulated and real matrix
    to use to determined the threshold at which to keep principal components
    Default: 2
    k_nn (optional): Number of neighbors to aggregate                    
    Default: Calculated so that aggregation yields ~ target_transcripts per cell
    target_transcripts (optional): Number of transcripts per cell to aim for in aggregation
    percent_cells_max (optional): Maximum percentage of cells to aggregate
    Default: 2

### Example

The following example uses pbmc-4k-expression.tsv.gz. 
Once in the directory of the github repository, run the following lines in R which should generate the plots shown below. 

    file = 'data/pbmc-4k-expression.tsv.gz' # file containing the expression matrix
    markers = c('CCR7','CD3D') # vector containing marker genes to plot

    # Loading data

    library(Seurat)
    library(Matrix)
    source('enhance.R')
    counts = read.table(gzfile(file), header = TRUE, row.names = 1)

    # Making a Seurat object using the raw data

    srt = CreateSeuratObject(counts = counts)

    srt = NormalizeData(srt)
    srt = ScaleData(srt, features = rownames(srt))

    srt = FindVariableFeatures(srt, nfeatures = 1000)
    srt = RunPCA(srt, features = VariableFeatures(srt))
    srt = RunTSNE(srt, dims = 1:10)

    # Plotting the raw expression of marker genes

    DefaultAssay(srt) # check that this is "RNA"
    FeaturePlot(srt, markers, slot = 'counts', reduction = 'tsne', col = c('grey','darkblue'))

    # Running enhance

    srt = enhance_seurat(srt, 
    assay = 'RNA',
    ratio_pcs = 2,
    percent_cells_max = 2,
    setDefaultAssay = TRUE)

    srt = NormalizeData(srt)
    srt = ScaleData(srt, features = rownames(srt))

    # Plotting the enhance expression of marker genes

    DefaultAssay(srt) # check that this is "enhance"
    FeaturePlot(srt, markers, slot = 'counts', reduction = 'tsne', col = c('grey','darkred'))



