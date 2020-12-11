library(arrow)
library(monocle3)
library(dplyr)


# on local
# df <- read_parquet("./covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")
# df_plate_map <- read_parquet("./plate_map_TS_202008.parquet")
# df_image_scores <- read_parquet("./image_scores_CQ1_TS_202008.parquet")
df <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")
df_plate_map <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/plate_map_TS_202008.parquet")
df_image_scores <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/image_scores_CQ1_TS_202008.parquet")

# df_uninfected = df %>% dplyr::select( based on time_point)
cell_features_uninfected <- df %>% dplyr::filter(time_point == "Uninfected" )
cell_features_8hrs <- df %>% dplyr::filter(time_point == "8 hours")
cell_features_24hrs <- df %>% dplyr::filter(time_point == "24 hours")
cell_features_30hrs <- df %>% dplyr::filter(time_point == "30 hours")
cell_features_36hrs <- df %>% dplyr::filter(time_point == "36 hours")
cell_features_48hrs <- df %>% dplyr::filter(time_point == "48 hours")

prepare_cds <- function(
  cell_features,
  verbose = FALSE) {

  # 429 is a column with NA values
  metadata <- cell_features[,1:22]
  cell_features <- cbind(cell_features[,23:428], cell_features[,430:824])

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
  cell_metadata <- cbind(cell_metadata, metadata)

  cds <- new_cell_data_set(expression_data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)
  cds
}

set.seed(1)


# first half
sampled_cell_features <- cell_features_uninfected %>% slice(seq(0.5 * n()))
cds <- prepare_cds(sampled_cell_features , verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds_reduced <- reduce_dimension(cds)
jpeg('plot_cell_first_half.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds_reduced, color_cells_by="Image_Metadata_WellID") + theme(legend.position = "right"))
dev.off()

# second half
sampled_cell_features <- cell_features_uninfected %>% slice(-seq(0.5 * n()))
cds <- prepare_cds(sampled_cell_features , verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds_reduced <- reduce_dimension(cds)
jpeg('plot_cell_second_half.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds_reduced, color_cells_by="Image_Metadata_WellID") + theme(legend.position = "right"))
dev.off()

# for (val in 1:10) {
  # for debug purpose
  # cell_features <- cell_features[1:1000,]
  # sampled_cell_features <- sample_n(cell_features_uninfected, size=1000)

  # first half
  # sampled_cell_features <- cell_features_uninfected %>% slice(seq(0.5 * n()))
  # second half
  # sampled_cell_features <- cell_features_uninfected %>% slice(-seq(0.5 * n()))
  # cds <- prepare_cds(sampled_cell_features , verbose=TRUE)
  # cds <- preprocess_cds(cds, num_dim = 100)
  # cds_reduced <- reduce_dimension(cds)
  # jpeg(paste0('plot_cell_', val,'.jpg'), width=2000, height=2000, res=300)
  # print(plot_cells(cds_reduced, color_cells_by="Image_Metadata_WellID"))
  # dev.off()
# }


cds <- cluster_cells(cds_reduced, resolution=1e-5)



# jpeg('plot_pc_variance.jpg', width=2000, height=2000, res=300)
# plot_pc_variance_explained(cds)
# dev.off()

# jpeg('plot_cell.jpg', width=2000, height=2000, res=300)
# plot_cells(cds)
# dev.off()

# different time series -> different batch
# one option: shape feature only?
# use shape feature to visualize?
