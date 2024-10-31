# BIOF501A-Project
Repository for BIOF501 Project
- Goal is to take in Spatial Transcriptomic Data, and output enriched pathways by using DEGs from cell-enriched regions
- Cell-enriched regions could be defined by some threshold, or user-defined

## General steps
  1) Preprocessing
     - load in data with seurat
     - QC: remove capture spots that
         - i) are (heavily) contaminated with mitochondrial genes
         - ii) have "extreme" gene counts
  2) Cell Deconvolution with SpaceXR
     - user needs to provide reference
     - if no reference provided, do unsupervised clustering 
  3) Differential Gene Expression with Seurat Findmarkers (for 1 sample) or LIMMA-VOOM (for multiple samples)
     - Should be performed for top `n` enriched cells
     - find DEGs off of unsupervised clusters if no reference provided
  4) Gene Ontology with Enrichr
