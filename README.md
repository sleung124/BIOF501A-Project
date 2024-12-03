# BIOF501A-Project

## Current way to run pipeline:
- `nextflow run workflow.nf -profile docker --degs.QUICK_SAMPLE 100`

## Big Issues to fix:
- [x] Output of modules not correct, but we rolling with it
      - pipeline runs but only if `temp_output` folder with subdirectories exists
      - ideally, i'm able to store temporary `.rds` files and call them in other modules while saying images to publishDir  
- [x] Have to specify `-with-docker sleung124/spatial-pipeline:latest`, or else pipeline runs without docker image for whatever reason
      - i know it's a parameter i have to modify in `nextflow.config`, but the things i've tried have not worked so far 
      - SOLVED: we now just use `-with-profile docker`
- [ ] Need to see if user-inputted path to a Visium data directory breaks my pipeline or not
      - don't think it will but who knows
- [ ] Missing comprehensive comments for code and README

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

### Dummy Data
- Currently from [10x genomics](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard)
