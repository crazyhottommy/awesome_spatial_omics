# awesome spatial omics

### Review papers

* [Spatial landscapes of cancers: insights and opportunities](https://doi.org/10.1038/s41571-024-00926-7)
* [The emerging landscape of spatial profiling technologies](https://www.nature.com/articles/s41576-022-00515-3)
* [The expanding vistas of spatial transcriptomics](https://www.nature.com/articles/s41587-022-01448-2)
* [Exploring tissue architecture using spatial transcriptomics](https://www.nature.com/articles/s41586-021-03634-9)
* [Statistical and machine learning methods for spatially resolved transcriptomics data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02653-7). first author Zexian was my colleague when I was at DFCI.
* [Spatial omics and multiplexed imaging to explore cancer biology](https://www.nature.com/articles/s41592-021-01203-6)
* [Method of the Year: spatially resolved transcriptomics](https://www.nature.com/articles/s41592-020-01033-y)
* [Computational challenges and opportunities in spatially resolved transcriptomic data analysis](https://www.nature.com/articles/s41467-021-25557-9) by Jean Fan.
* [Spatial components of molecular tissue biology](https://www.nature.com/articles/s41587-021-01182-1)
* [Methods and applications for single-cell and spatial multi-omics](https://www.nature.com/articles/s41576-023-00580-2)
* [The dawn of spatial omics](https://pubmed.ncbi.nlm.nih.gov/37535749/)
  
### Tutorials

* [Orchestrating Spatially-Resolved Transcriptomics Analysis with Bioconductor](https://lmweber.org/OSTA-book/)
* [Deconvolution vs Clustering Analysis for Multi-cellular Pixel-Resolution Spatially Resolved Transcriptomics Data](https://jef.works/blog/2022/05/03/deconvolution-vs-clustering/) A blog post by Jean Fan.
* [Analysis, visualization, and integration of spatial datasets with Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)

### Benchmarking 

* [Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution](https://www.nature.com/articles/s41592-022-01480-9) We found that Tangram, gimVI, and SpaGE outperformed other integration methods for predicting the spatial distribution of RNA transcripts, whereas Cell2location, SpatialDWLS, and RCTD are the top-performing methods for the cell type deconvolution of spots.
* [Robust alignment of single-cell and spatial transcriptomes with CytoSPACE](https://www.biorxiv.org/content/10.1101/2022.05.20.488356v1.full.pdf)
* [A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics](https://www.nature.com/articles/s41467-023-37168-7)
* [Comparative analysis of MERFISH spatial transcriptomics with bulk and single-cell RNA sequencing](https://doi.org/10.1101/2022.03.04.483068)

### Databases

* [SODB facilitates comprehensive exploration of spatial omics data](https://www.nature.com/articles/s41592-023-01773-7) [[website]](https://gene.ai.tencent.com/SpatialOmics/)
* [Museum of Spatial Transcriptomics](https://pachterlab.github.io/LP_2021/index.html)

### Methods
Generally grouped according to [Rao et al](https://www.nature.com/articles/s41586-021-03634-9), though several current in situ sequencing methods also rely on ex situ sequencing

#### Ex situ sequencing-based
* **Spatial Transcriptomics** - [Visualization and analysis of gene expression in tissue sections by spatial transcriptomics](https://www.science.org/doi/10.1126/science.aaf2403) (*Note* now commercialized as [Visium, 10x Genomics](https://www.10xgenomics.com/products/spatial-gene-expression))
* **HDST** - [High-definition spatial transcriptomics for in situ tissue profiling](https://www.nature.com/articles/s41592-019-0548-y)
* **STRS** - [In situ polyadenylation enables spatial mapping of the total transcriptome](https://www.nature.com/articles/s41587-022-01517-6) [[code]](https://github.com/mckellardw/STRS)
* **Expansion Spatial Transcritptomics** - [Expansion spatial transcriptomics](https://www.nature.com/articles/s41592-023-01911-1)
* **SlideSeq** - [Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution](https://www.science.org/doi/10.1126/science.aaw1219)
* **SlideSeq2** - [Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2](https://doi.org/10.1038/s41587-020-0739-1) (*Note* now commercialized as Seeker, [Curio Bioscience](https://curiobioscience.com/))
* **SlideTags** - [Slide-tags: scalable, single-nucleus barcoding for multi-modal spatial genomics](https://www.biorxiv.org/content/10.1101/2023.04.01.535228v1)
* **StereoSeq** - [Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays](https://www.sciencedirect.com/science/article/pii/S0092867422003993) (*Note* now commercialized as STOmics, [MGI/BGI/Complete Genomics](https://bgi-australia.com.au/stomics))
* **DNA Microscopy** - [DNA Microscopy: Optics-free Spatio-genetic Imaging by a Stand-Alone Chemical Reaction](https://www.sciencedirect.com/science/article/pii/S0092867419305471)
* **Volumetric DNA Microscopy** - [Volumetric imaging of an intact organism by a distributed molecular network](https://www.biorxiv.org/content/10.1101/2023.08.11.553025v1)

#### In situ sequencing-based
* **STARmap** - [Three-dimensional intact-tissue sequencing of single-cell transcriptional states](https://www.science.org/doi/10.1126/science.aat5691)
* **STARmap PLUS** - [Integrative in situ mapping of single-cell transcriptional states and tissue histopathology in a mouse model of Alzheimer’s disease](https://www.nature.com/articles/s41593-022-01251-x)
* **IGS** - [In situ genome sequencing resolves DNA sequence and structure in intact biological samples](https://www.science.org/doi/10.1126/science.aay3446)

#### Imaging/hybridization-based
* **MERFISH** - [Spatially resolved, highly multiplexed RNA profiling in single cells](https://www.science.org/doi/10.1126/science.aaa6090)
* **HiPR-FISH** - [Highly Multiplexed Spatial Mapping of Microbial Communities](10.1038/s41586-020-2983-4)
* **EEL-FISH** - [Scalable in situ single-cell profiling by electrophoretic capture of mRNA using EEL FISH](https://www.nature.com/articles/s41587-022-01455-3)

### normalization

* [Gene count normalization in single-cell imaging-based spatially resolved transcriptomic](https://www.biorxiv.org/content/10.1101/2023.08.30.555624v1)
* >Spatial omics data analysts sometimes use the “log1p” (y=log[1+x’]) transform incorrectly. Key fact: x’ represents normalized, not raw, umis/spot counts. And it really matters how you normalize! The figure shows three normalizations of the raw count
  
![image](https://github.com/crazyhottommy/awesome_spatial_omics/assets/4106146/c09126a3-9876-48cf-9984-1102e079c7c8)

https://twitter.com/shyam_lab/status/1698170321155850433?s=51&t=sLukUyq0ReWrcwOwgUR_XA

### Computational Tools

* [Sopa](https://github.com/gustaveroussy/sopa)  enables processing and analyses of image-based spatial-omics using a standard data structure and output. We currently support the following technologies: Xenium, MERSCOPE, CosMX, PhenoCycler, MACSIMA, Hyperion. Sopa was designed for generability and low-memory consumption on large images (scales to 1TB+ images).
* [Monkeybread](https://monkeybread.readthedocs.io/en/latest/notebooks/tutorial.html) A python package developed at Immunitas to do spatial analysis for Merfish data.
* [Giotto](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2) a toolbox for integrative analysis and visualization of spatial expression data
* [Voyager](https://pachterlab.github.io/voyager/index.html) is a package that facilitates exploratory spatial data analysis and visualization for spatial genomics data represented by SpatialFeatureExperiment objects.
* [nnSVG](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1): scalable identification of spatially variable genes using nearest-neighbor Gaussian processes
* [DestVI](https://www.nature.com/articles/s41587-022-01272-8) identifies continuums of cell types in spatial transcriptomics data. DestVI is available as part of the open-source software package scvi-tools (https://scvi-tools.org).
* Here we present [spateo](https://spateo-release.readthedocs.io/en/latest/), a open source framework that welcomes community contributions for quantitative spatiotemporal modeling of spatial transcriptomics.
* [SpaGene]( https://biorxiv.org/content/10.1101/2022.04.20.488961v1.full.pdf
): Scalable and model-free detection of spatial patterns and colocalization
* [Palo](https://www.biorxiv.org/content/10.1101/2022.03.13.484080v1): Spatially-aware color palette optimization for single-cell and spatial data
* **squidpy** - [paper](https://www.nature.com/articles/s41592-021-01358-2) - [code](https://squidpy.readthedocs.io/en/latest/): Squidpy: a scalable framework for spatial omics analysis
* **ncem** - [paper](https://www.biorxiv.org/content/10.1101/2021.07.11.451750v1) - [code](https://ncem.readthedocs.io/en/latest/): Learning cell communication from spatial graphs of cells
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
* **StereoCell** - [StereoCell enables high accuracy single cell segmentation for spatial transcriptomic dataset](https://www.biorxiv.org/content/10.1101/2023.04.01.535228v1)
* **cell2location**
* **STcEM** - [Spatial Transcriptomics-correlated Electron Microscopy maps transcriptional and ultrastructural responses to brain injury](https://www.nature.com/articles/s41467-023-39447-9)

#### Deconvolution
* spacedeconv is a unified interface to 31 deconvolution tools with a focus on spatial transcriptomics datasets. The package is able to directly estimate cell type proportions of immune cells and can deconvolute any cell type if an annotation single-cell reference dataset is available https://github.com/omnideconv/spacedeconv
  
#### Differential expression

* [A statistical method to uncover gene expression changes in spatial transcriptomics](https://www.nature.com/articles/s41592-022-01576-2) Cell type-specific inference of differential expression ([C-SIDE](https://github.com/dmcable/spacexr)) is a statistical model that identifies which genes (within a determined cell type) are differentially expressed on the basis of spatial position, pathological changes or cell–cell interactions.
* Niche differential gene expression analysis in spatial transcriptomics data identifies context-dependent cell-cell interactions https://www.biorxiv.org/content/10.1101/2023.01.03.522646v1

#### Spatial domain 

* BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis https://www.nature.com/articles/s41588-024-01664-3
* 
#### Integration 

* [Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST](https://www.nature.com/articles/s41467-023-35947-w)
* [High-resolution alignment of single-cell and spatial transcriptomes with CytoSPACE](https://www.nature.com/articles/s41587-023-01697-9)
* [Search and Match across Spatial Omics Samples at Single-cell Resolution](https://www.biorxiv.org/content/10.1101/2023.08.13.552987v1)
* [Alignment of spatial genomics data using deep Gaussian processes](https://www.nature.com/articles/s41592-023-01972-2)
* [Starfysh integrates spatial transcriptomic and histologic data to reveal heterogeneous tumor–immune hubs](https://www.nature.com/articles/s41587-024-02173-8)

#### 3D reconstruction

* [Alignment and integration of spatial transcriptomics data](https://www.nature.com/articles/s41592-022-01459-6)

#### Clustering

* [BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02734-7)

* [DeepST: A versatile graph contrastive learning framework for spatially informed clustering, integration, and deconvolution of spatial transcriptomics](https://www.biorxiv.org/content/10.1101/2022.08.02.502407v1)

#### Cell-cell interaction

* [De novo reconstruction of cell interaction landscapes from single-cell spatial transcriptome data with DeepLinc](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02692-0)

* [Modeling intercellular communication in tissues using spatial graphs of cells](https://www.nature.com/articles/s41587-022-01467-z)

#### Imputation

* [Accurate inference of genome-wide spatial expression with iSpatial](https://www.biorxiv.org/content/10.1101/2022.05.23.493144v2)

#### Interactive tool

* [VITESSCE](https://github.com/vitessce/vitessce) Visual Integration Tool for Exploration of Spatial Single-Cell Experiments
