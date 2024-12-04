# SOP Template (replace with name of pipeline here)

## Steps to highlight:
### 1) Background and Rationale
- Include the what's and why's, as well as aims
- Include package dependencies that are required (bullet points OK)

    - include the what’s and why’s; also your aims
    - include any package dependencies that are required (bullet points are ok for this)
    - You can include your DAG here
### 2) Usage
    Make sure you format everything so that step by step usage details are included. If we can’t run your pipeline then we can’t give you marks.

### 3) Installation (if necessary) 
- including any datasets that are to be used if they are not provided
(i.e. how to download them using wget or curl – exact paths need to be specified and the data
must be accessible)

- Exact step by step usage with descriptive comments on what action is being performed in each
step

### 4) Input
- Describe the format of the input data, explaining all fields.


### 5) Output
- Describe the format of the output including files and visualizations. Treat this section like the results of
a paper. You can look at readthedocs pages of popular bioinformatics tools to get inspired for this.

## Background

### Workflow Overview

```
flowchart TD
    input --> Preprocessing
    Preprocessing --> Cell_Deconvolution
    Cell_Deconvolution --> DEGs
    DEGs --> Pathways

    subgraph Preprocessing[**Preprocessing**]
        direction LR
        mito_filter(Seurat) --> mitoplot
    end

    subgraph Cell_Deconvolution[**Cell Deconvolution**]
        direction LR
        STdeconvolve --> spatial_celltype
    end

    subgraph DEGs[**Differential Gene Analysis**]
        direction LR
        FindMarkers(Seurat) --> degs_csv
    end

    subgraph Pathways[**Enriched Pathways**]
        direction LR
        Enrichr --> enrichr_plot
    end

    input[**Visium Spatial Transcriptomic Data**]
    mitoplot[Mitochondrial Spatial Plot]
    spatial_celltype[Cell Type Enrichment Spatial Plots]
    degs_csv[CSV of Differentially Expressed Genes]
    enrichr_plot[Pathway Enrichment Plots]

    classDef inputPath fill:salmon,stroke:#333,stroke-width:2px;
    class input inputPath;

    classDef output fill:lightgreen,stroke:#333,stroke-width:2px;
    class mitoplot,spatial_celltype,degs_csv,enrichr_plot output;

    classDef process fill:#ccf,stroke:#333,stroke-width:2px;
    class mito_filter,STdeconvolve,FindMarkers,Enrichr process;

```

### Core R Package Versions
TODO: add package versions for docker, git, nextflow
```r
# Format is [package] - [version]
enrichR - v3.2
STdeconvolve - v1.3.2
here - v1.0.1
forcats - v0.5.1
stringr - v1.5.1
dplyr - v1.1.4
purrr - v1.0.2
readr - v2.0.0
tidyr - v1.3.1
tibble - v3.2.1
ggplot2 - v3.5.1
tidyverse - v1.3.1
Seurat - v5.0.1
SeuratObject - v5.0.1
sp - v2.1-4
```
### Repository Structure
