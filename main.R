library(arrow)
library(monocle3)
library(dplyr)
library(philentropy)

# on local
# df <- read_parquet("./covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")
# df_plate_map <- read_parquet("./plate_map_TS_202008.parquet")
# df_image_scores <- read_parquet("./image_scores_CQ1_TS_202008.parquet")
df <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet")
df_plate_map <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/plate_map_TS_202008.parquet")
df_image_scores <- read_parquet("/nfs/turbo/umms-jzsexton/MPalign/image_scores_CQ1_TS_202008.parquet")

cell_metadata_columns <- tibble::tibble(
    feature = df_image_scores %>% names())

objects <- c("Cells", "Nuclei", "Cytoplasm", "InfectedCells")
dyes <- c("ConA", "Hoe", "NP", "Spike")
coordinates <- c("X", "Y", "Z")

cell_feature_columns <- tibble::tibble(
    feature = df %>% names(),
    transform = "identity") %>%
    dplyr::anti_join(cell_metadata_columns, by = "feature") %>%
    dplyr::anti_join(
        expand.grid(child = objects, parent = objects) %>%
        dplyr::mutate(feature = paste0(child, "_Parent_", parent)),
        by = "feature") %>%
    dplyr::filter(!(feature %in% c("Nuclei_Distance_Centroid_InfectedCells"))) %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            feature = c(
                "Location",
                "AreaShape"),
            coordinate = coordinates) %>%
        dplyr::mutate(
            feature = paste(object, feature, "Center", coordinate, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            statistic = c(
                "MaxIntensity",
                "CenterMassIntensity"),
            coordinate = coordinates,
            dye = dyes) %>%
        dplyr::mutate(
            feature = paste(object, "Location", statistic, coordinate, dye, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            object = objects,
            class = c("Positive", "Negative")) %>%
        dplyr::mutate(
            feature = paste(object, "Classify", class, sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        expand.grid(
            parent = objects,
            child = objects) %>%
        dplyr::mutate(
            feature = paste(parent, "Children", child, "Count", sep = "_")),
        by = "feature") %>%
    dplyr::anti_join(
        data.frame(objects) %>%
        dplyr::mutate(feature = paste0(objects, "_Number_Object_Number")),
        by = "feature")


# df_uninfected = df %>% dplyr::select( based on time_point)
cell_features_uninfected <- df %>% dplyr::filter(time_point == "Uninfected" )
cell_features_8hrs <- df %>% dplyr::filter(time_point == "8 hours")
cell_features_24hrs <- df %>% dplyr::filter(time_point == "24 hours")
cell_features_30hrs <- df %>% dplyr::filter(time_point == "30 hours")
cell_features_36hrs <- df %>% dplyr::filter(time_point == "36 hours")
cell_features_48hrs <- df %>% dplyr::filter(time_point == "48 hours")

prepare_cds <- function(
  cell_features,
  cell_metadata_columns,
  cell_feature_columns,
  verbose=FALSE) {

  metadata <- cell_features[, cell_metadata_columns]
  cell_features <- cell_features[, cell_feature_columns]
  # if (artificial_batch){
  #   # 429 is a column with NA values
  #  metadata <- cell_features[,1:22], cell_features[,ncol(cell_features)])
  #  cell_features <- cbind(cell_features[,23:428], cell_features[,430:(ncol(cell_features)-1)])
  #} else {
  #  # 429 is a column with NA values
  #  metadata <- cell_features[,1:22]
  #  cell_features <- cbind(cell_features[,23:428], cell_features[,430:ncol(cell_features)])
  #}

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

# exp 0
# first half
sampled_cell_features <- cell_features_uninfected %>% slice(seq(0.5 * n()))
cds <- prepare_cds(sampled_cell_features, cell_metadata_columns$feature, cell_feature_columns$feature, verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds_reduced <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
jpeg('plot_cell_first_half.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds_reduced, reduction_method="PCA", color_cells_by="Image_Metadata_FieldID", label_cell_groups=FALSE))
dev.off()

# second half
sampled_cell_features <- cell_features_uninfected %>% slice(-seq(0.5 * n()))
cds <- prepare_cds(sampled_cell_features, cell_metadata_columns$feature, cell_feature_columns$feature, verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds_reduced <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
jpeg('plot_cell_second_half.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds_reduced, reduction_method="PCA", color_cells_by="Image_Metadata_FieldID", label_cell_groups=FALSE))
dev.off()

# As a proof of concept, try to cluster field 1 and field 2 independently of uninfecte_cell_features
field_1 <- cell_features_uninfected %>% dplyr::filter(Image_Metadata_FieldID== "0001")
cds <- prepare_cds(field_1, cell_metadata_columns$feature, cell_feature_columns$feature, verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
# seems not work
cds <- cluster_cells(cds, reduction_method="PCA", resolution=1e-5)
jpeg('plot_cell_first_field.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds, reduction_method="PCA"))
dev.off()

field_2 <- cell_features_uninfected %>% dplyr::filter(Image_Metadata_FieldID== "0002")
cds <- prepare_cds(field_2, cell_metadata_columns$feature, cell_feature_columns$feature, verbose=TRUE)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
cds <- cluster_cells(cds, reduction_method="PCA", resolution=1e-5)
jpeg('plot_cell_second_field.jpg', width=2000, height=2000, res=300)
print(plot_cells(cds, reduction_method="PCA"))
dev.off()

# exp 1.
# each time point -> divide it into two batches -> align these batches within one time point (using the removing the batches effect function)
exp1 <- function(cell_features, filename_prefix){
  cell_features$artificial_batch <- sample(c("batch 1", "batch 2"), length(cell_features$Image_Metadata_PlateID), replace = TRUE)
  cds <- prepare_cds(cell_features, c(cell_metadata_columns$feature, "artificial_batch"), cell_feature_columns$feature, verbose=TRUE)
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
  jpeg(paste0(filename_prefix, '_before_remove_batch_effect.jpg'), width=2000, height=2000, res=300)
  print(plot_cells(cds, reduction_method="PCA", color_cells_by="artificial_batch", label_cell_groups=FALSE))
  dev.off()
  cds <- align_cds(cds, num_dim = 100, alignment_group = "artificial_batch")
  cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
  jpeg(paste0(filename_prefix, '_after_remove_batch_effect.jpg'), width=2000, height=2000, res=300)
  print(plot_cells(cds, reduction_method="PCA", color_cells_by="artificial_batch", label_cell_groups=FALSE))
  dev.off()
  cds <- cluster_cells(cds, reduction_method="PCA", resolution=1e-5)
  temp <- data.frame(colData(cds)$artificial_batch, clusters(cds, reduction_method="PCA"))
  colnames(temp) <- c("artificial.batch", "cluster")
  count_batch_1 <- temp %>% filter(artificial.batch == "batch 1") %>% count(cluster, sort=TRUE)
  count_batch_2 <- temp %>% filter(artificial.batch == "batch 2") %>% count(cluster, sort=TRUE)
  p <- count_batch_1["n"] / sum(count_batch_1["n"])
  q <- count_batch_2["n"] / sum(count_batch_2["n"])
  cat("KL 1: ", KL(t(cbind(p,q))))
  cat("KL 2: ", KL(t(cbind(q,p))))
  # use JSD instead because KL is not symemtric
  cat("JSD: ", JSD(t(cbind(p,q))))
}

cat("\nexp1:")
cat("Uninfected:")
exp1(cell_features_uninfected, "uninfected")
cat("8hrs:")
exp1(cell_features_8hrs, "8hrs")
cat("24hrs:")
exp1(cell_features_24hrs, "24hrs")
cat("30hrs")
exp1(cell_features_30hrs, "30hrs")
cat("36hrs")
exp1(cell_features_36hrs, "36hrs")
cat("48hrs")
exp1(cell_features_48hrs, "48hrs")

# exp 2
# pairs of time points (uninfected, 8 hrs) -> cluster together -> align alignment group
exp2 <- function(cell_features_time_1, cell_features_time_2, time_1, time_2){
  cell_features <- rbind(cell_features_time_1, cell_features_time_2)
  cds <- prepare_cds(cell_features, cell_metadata_columns$feature, cell_feature_columns$feature, verbose=TRUE)
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
  jpeg(paste0(time_1, "_", time_2, '_before_remove_batch_effect.jpg'), width=2000, height=2000, res=300)
  print(plot_cells(cds, reduction_method="PCA", color_cells_by="time_point", label_cell_groups=FALSE))
  dev.off()
  cds <- align_cds(cds, num_dim = 100, alignment_group = "time_point")
  cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="PCA")
  jpeg(paste0(time_1, "_", time_2, '_after_remove_batch_effect.jpg'), width=2000, height=2000, res=300)
  print(plot_cells(cds, reduction_method="PCA", color_cells_by="time_point", label_cell_groups=FALSE))
  dev.off()
  cds <- cluster_cells(cds, reduction_method="PCA", resolution=1e-5)
  temp <- data.frame(colData(cds)$time_point, clusters(cds, reduction_method="PCA"))
  colnames(temp) <- c("time_point", "cluster")
  count_time_1 <- temp %>% filter(time_point == time_1) %>% count(cluster, sort=TRUE)
  count_time_2 <- temp %>% filter(time_point == time_2) %>% count(cluster, sort=TRUE)
  p <- count_time_1["n"] / sum(count_time_1["n"])
  q <- count_time_2["n"] / sum(count_time_2["n"])
  cat("KL 1: ", KL(t(cbind(p,q))))
  cat("KL 2: ", KL(t(cbind(q,p))))
  # use JSD instead because KL is not symemtric
  cat("JSD: ", JSD(t(cbind(p,q))))
}

cat("\nexp2:")
cat("uninfected -> 8hrs")
exp2(cell_features_uninfected, cell_features_8hrs, "uninfected", "8hrs")
cat("uninfected -> 24hrs")
exp2(cell_features_uninfected, cell_features_24hrs, "uninfected", "24hrs")
cat("uninfected -> 30hrs")
exp2(cell_features_uninfected, cell_features_30hrs, "uninfected", "30hrs")
cat("uninfected -> 36hrs")
exp2(cell_features_uninfected, cell_features_36hrs, "uninfected", "36hrs")
cat("uninfected -> 48hrs")
exp2(cell_features_uninfected, cell_features_48hrs, "uninfected", "48hrs")
cat("8hrs -> 24hrs")
exp2(cell_features_8hrs, cell_features_24hrs, "8hrs", "24hrs")
cat("8hrs -> 30hrs")
exp2(cell_features_8hrs, cell_features_30hrs, "8hrs", "30hrs")
cat("8hrs -> 36hrs")
exp2(cell_features_8hrs, cell_features_36hrs, "8hrs", "36hrs")
cat("8hrs -> 48hrs")
exp2(cell_features_8hrs, cell_features_48hrs, "8hrs", "48hrs")
cat("24hrs -> 30hrs")
exp2(cell_features_24hrs, cell_features_30hrs, "24hrs", "30hrs")
cat("24hrs -> 36hrs")
exp2(cell_features_24hrs, cell_features_36hrs, "24hrs", "36hrs")
cat("24hrs -> 48hrs")
exp2(cell_features_24hrs, cell_features_48hrs, "24hrs", "48hrs")
cat("30rs -> 36hrs")
exp2(cell_features_30hrs, cell_features_36hrs, "30hrs", "36hrs")
cat("30rs -> 48hrs")
exp2(cell_features_30hrs, cell_features_48hrs, "30hrs", "48hrs")
cat("36rs -> 36hrs")
exp2(cell_features_36hrs, cell_features_48hrs, "36hrs", "48hrs")

# exp 3: do what we have in exp 1 and exp 2 100 times

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


# cds <- cluster_cells(cds_reduced, resolution=1e-5)

# jpeg('plot_pc_variance.jpg', width=2000, height=2000, res=300)
# plot_pc_variance_explained(cds)
# dev.off()

# jpeg('plot_cell.jpg', width=2000, height=2000, res=300)
# plot_cells(cds)
# dev.off()

# different time series -> different batch
# one option: shape feature only?
# use shape feature to visualize?

