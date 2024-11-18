# install devtools to install CRAN packages
install.packages("devtools")

# from https://github.com/satijalab/seurat-docker/blob/master/latest/Dockerfile
install.packages('BiocManager')
bioc_lst = c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi', 'variancePartition')
BiocManager::install(bioc_lst)

# Install CRAN suggests
# need statmod to fix error in running fitLM after duplicate correlation
cran_lst = c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'mixtools', 'ggrepel', 'harmony', 'patchwork', 'statmod', 'hdf5r', 'pals')
install.packages(cran_lst, repos = 'https://cloud.r-project.org')

# Install Seurat
remotes::install_version("Matrix", version="1.6-3")
remotes::install_version("sctransform", version="0.4.1")
remotes::install_version("SeuratObject", version="5.0.1")
remotes::install_version("Seurat", version="5.0.1")

# Install SeuratDisk
remotes::install_github('mojaveazure/seurat-disk')

## additoinal packages needed to run seurat visium vignette
devtools::install_github('satijalab/seurat-data')

options(timeout = 600000000)
# install spaceXR package
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# install enrichR package
devtools::install_github("wjawaid/enrichR")

# issue with using latest version of matrixStats, so we downgrading
remotes::install_version("matrixStats", version="1.1.0")
