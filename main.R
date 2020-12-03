library(arrow)
library(monocle3)
library(dplyr)
library(umap)

df <- read_parquet("./covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")
df_plate_map <- read_parquet("./plate_map_TS_202008.parquet")
df_image_scores <- read_parquet("./image_scores_CQ1_TS_202008.parquet")

# 429 is a column with NA values
cell_features <- cbind(df[,23:428], df[,430:824])

populate_cds <- function(
  cell_features,
  embedding_type = c("UMAP"),
  embedding = NULL,
  verbose = FALSE) {


  n_cells <- nrow(cell_features)
  n_features <- ncol(cell_features)
  cat("Loading cell dataset with dimensions [<feature>, <cell>] = [", n_features, ", ",  n_cells, "]\n", sep = "")

  expression_data <- cell_features %>%
    as.matrix() %>%
    t()

  gene_metadata <- colnames(cell_features)
  gene_metadata <- as.matrix(gene_metadata)
  colnames(gene_metadata) <- c("gene_short_name")
  rownames(gene_metadata) <- gene_metadata
  
  cell_metadata <- rownames(cell_features)
  cell_metadata <- as.matrix(cell_metadata)
  colnames(cell_metadata) <- c("id")
  rownames(cell_metadata) <- cell_metadata

  colnames(expression_data) <- cell_metadata
  
  cds <- new_cell_data_set(expression_data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)
  cds

  # if(verbose){
  #   cat("Creating a SingleCellExperiment object ...\n")
  # }
  # 
  # # unpack monocle3::new_cell_data_set(...)
  # # to not use dgCMatrix for the expression matrix they are dense feature matrices
  # sce <- SingleCellExperiment(
  #   list(counts = expression_data),
  #   rowData = gene_metadata,
  #   colData = cell_metadata)
  # 
  # if(verbose){
  #   cat("Creating a Cell Data Set object ...\n")
  # }
  # cds <- methods::new(
  #   Class = "cell_data_set",
  #   assays = SummarizedExperiment::Assays(list(counts = expression_data)),
  #   colData = colData(sce),
  #   int_elementMetadata = int_elementMetadata(sce),
  #   int_colData = int_colData(sce),
  #   int_metadata = int_metadata(sce),
  #   metadata = S4Vectors::metadata(sce),
  #   NAMES = NULL,
  #   elementMetadata = elementMetadata(sce)[,0],
  #   rowRanges = rowRanges(sce))
  # 
  # if(verbose){
  #   cat("Configuring the cell data set ...\n")
  # }
  # 
  # S4Vectors::metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  # clusters <- stats::setNames(S4Vectors::SimpleList(), character(0))
  # # cds <- monocle3::estimate_size_factors(cds)
  # 
  # row.names(SummarizedExperiment::colData(cds)) <- expression_data %>% ncol %>% seq_len
  # if (!is.null(embedding)) {
  #   SingleCellExperiment::reducedDims(cds)[[embedding_type]] <- embedding
  # }
  # cds
}

cds <- populate_cds(cell_features, verbose=TRUE)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

jpeg('plot_pc_variance.jpg', width=2000, height=2000, res=300)
plot_pc_variance_explained(cds)
dev.off()

jpeg('plot_cell.jpg', width=2000, height=2000, res=300)
plot_cells(cds)
dev.off()

# different time series -> different batch
# one option: shape feature only?
# use shape feature to visualize?
