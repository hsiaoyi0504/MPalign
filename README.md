# MPalign

## Working on GreatLake

### Installation

``` shell
conda deactivate
module load R/4.0.2 gcc/8.2.0 gdal/3.0.1 proj/6.1.1 geos/3.7.2
```

In R,

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
```
