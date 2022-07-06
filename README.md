# awesome spatial omics

### review papers

* [Statistical and machine learning methods for spatially resolved transcriptomics data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02653-7). first author Zexian was my colleague when I was at DFCI.
* [Spatial omics and multiplexed imaging to explore cancer biology](https://www.nature.com/articles/s41592-021-01203-6)
* [Method of the Year: spatially resolved transcriptomics](https://www.nature.com/articles/s41592-020-01033-y)
* [Computational challenges and opportunities in spatially resolved transcriptomic data analysis](https://www.nature.com/articles/s41467-021-25557-9) by Jean Fan.
* [Spatial components of molecular tissue biology](https://www.nature.com/articles/s41587-021-01182-1)

### tutorial

* [Deconvolution vs Clustering Analysis for Multi-cellular Pixel-Resolution Spatially Resolved Transcriptomics Data](https://jef.works/blog/2022/05/03/deconvolution-vs-clustering/) A blog post by Jean Fan.

### benchmarking 

* [Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution](https://www.nature.com/articles/s41592-022-01480-9) We found that Tangram, gimVI, and SpaGE outperformed other integration methods for predicting the spatial distribution of RNA transcripts, whereas Cell2location, SpatialDWLS, and RCTD are the top-performing methods for the cell type deconvolution of spots.

### methods

[In situ polyadenylation enables spatial mapping of the total transcriptome](https://www.biorxiv.org/content/10.1101/2022.04.20.488964v1)

### tools

* [Giotto](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2) a toolbox for integrative analysis and visualization of spatial expression data
* [nnSVG](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1): scalable identification of spatially variable genes using nearest-neighbor Gaussian processes
* [DestVI](https://www.nature.com/articles/s41587-022-01272-8) identifies continuums of cell types in spatial transcriptomics data. DestVI is available as part of the open-source software package scvi-tools (https://scvi-tools.org).
* Here we present [spateo](https://spateo-release.readthedocs.io/en/latest/), a open source framework that welcomes community contributions for quantitative spatiotemporal modeling of spatial transcriptomics.
* [SpaGene]( https://biorxiv.org/content/10.1101/2022.04.20.488961v1.full.pdf
): Scalable and model-free detection of spatial patterns and colocalization
* [Palo](https://www.biorxiv.org/content/10.1101/2022.03.13.484080v1): Spatially-aware color palette optimization for single-cell and spatial data
* **squidpy** [paper](https://www.nature.com/articles/s41592-021-01358-2) - [code](https://squidpy.readthedocs.io/en/latest/): Squidpy: a scalable framework for spatial omics analysis
* **ncem** [paper](https://www.biorxiv.org/content/10.1101/2021.07.11.451750v1) - [code](https://ncem.readthedocs.io/en/latest/): Learning cell communication from spatial graphs of cells
* [Spatially informed cell-type deconvolution for spatial transcriptomics](https://www.nature.com/articles/s41587-022-01273-7) Here, we introduce a deconvolution method, conditional autoregressive-based deconvolution (CARD), that combines cell-type-specific expression information from single-cell RNA sequencing (scRNA-seq) with correlation in cell-type composition across tissue locations. https://github.com/YingMa0107/CARD
* Reconstruction of the cell pseudo-space from single-cell RNA sequencing data with [scSpace](https://www.biorxiv.org/content/10.1101/2022.05.07.491043v1)
* [SpatialCorr](https://www.biorxiv.org/content/10.1101/2022.02.04.479191v1.full): Identifying Gene Sets with Spatially Varying Correlation Structure
* [RCTD](https://www.nature.com/articles/s41587-021-00830-w): Robust decomposition of cell type mixtures in spatial transcriptomics
* [Supervised spatial inference of dissociated single-cell data with SageNet](https://www.biorxiv.org/content/10.1101/2022.04.14.488419v1): a graph neural network approach that spatially reconstructs dissociated single cell data using one or more spatial references. [code](https://github.com/MarioniLab/SageNet)
* [SpotClean adjusts for spot swapping in spatial transcriptomics data](https://www.biorxiv.org/content/10.1101/2021.06.11.448105v3.full): A quality issue in spatial transcriptomics data, and a statistical method to adjust for it. [R Package](https://github.com/zijianni/SpotClean).
* [Nonnegative spatial factorization](https://arxiv.org/abs/2110.06122)
* [SPICEMIX: Integrative single-cell spatial modeling of cell identity](https://www.biorxiv.org/content/10.1101/2020.11.29.383067v3)
* [De novo reconstruction of cell interaction landscapes from single-cell spatial transcriptome data with DeepLinc](https://pubmed.ncbi.nlm.nih.gov/35659722/)
* [Bayesian Modeling of Spatial Molecular Profiling Data via Gaussian Process](https://arxiv.org/abs/2012.03326)
* [Decoding functional cell-cell communication events by multi-view graph learning on spatial transcriptomics](https://www.biorxiv.org/content/10.1101/2022.06.22.496105v1)
* [BANKSY](https://www.biorxiv.org/content/10.1101/2022.04.14.488259v1) unifies cell-type clustering and domain segmentation by constructing a product space of cells' own and microenvironment transcriptomes. [R](https://github.com/prabhakarlab/Banksy) and [python](https://github.com/prabhakarlab/Banksy_py) code.  

### imputation

* [Accurate inference of genome-wide spatial expression with iSpatial](https://www.biorxiv.org/content/10.1101/2022.05.23.493144v2)

### Interactive tool

* [VITESSCE](https://github.com/vitessce/vitessce) Visual Integration Tool for Exploration of Spatial Single-Cell Experiments
