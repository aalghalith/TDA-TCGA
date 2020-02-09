library(TDAmapper)
library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(umap)
library(bioDist)

#create_mut_object <- function(exp_table, mut_table, num_intervals=c(20,20), percent_overlap=35, var_threshold, num_cores=1)


# INPUT AND CLEAN ---------------------------------------------------------

#exp_table <- readline("Input gene expression table: ")
#exp_table <- (read.csv("C:/Users/Adam/Desktop/Camara Lab/TDA_TCGA_Project/Test Data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #Laptop
exp_table <- (read.csv("/home/rstudio/documents/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #Lab
#mut_table <- (read.csv("blah", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))

exp_table <- exp_table[!duplicated(rownames(exp_table)),!duplicated(colnames(exp_table))] #No duplicated genes/samples
#exp_table <- exp_table[!duplicated(as.list(exp_table))]
mut_table <- mut_table[!duplicated(rownames(mut_table)),!duplicated(colnames(mut_table))]

if(all(colnames(exp_table) %in% colnames(mut_table))) { #Gene names are same CHANGE
  if(all(colnames(mut_table) %in% colnames(exp_table))) {
    print("Matching genes between expression and mutation tables")
  }
} else {
  print("Unmatched genes between expression and mutation tables")
}

exp_table_top <- exp_table[,order(-apply(exp_table,2,var))][,1:5000] #5000 genes with highest var.. CHANGE

# MAPPER ------------------------------------------------------------------

plot.umap = function(emb) {
  ggplot(emb) +
    geom_point(aes(x=V1, y=V2, color=0), size=0.9)
}

build.mapper = function(dist, umap) {
  mapper2D(dist, umap, c(20, 20), percent_overlap=50, num_bins_when_clustering=10)
}

dist_matrix <- as.matrix(cor.dist(as.matrix(exp_table_top))) #from bioDist library in Bioconductor

umap_emb = umap(dist_matrix)$layout %>% as.data.frame
plot.umap(umap_emb)

mapperObj <- build.mapper(dist_matrix,umap_emb)
adj <- graph.adjacency(mapperObj$adjacency, mode="undirected")
plot(adj, layout = layout.auto(adj))


#### TRY OUT FROM STEVEN

build.cover = function(filter_values=umap_emb, num_intervals=c(20, 20), percent_overlap=30) {
  filter_min_1 <- min(filter_values[[1]])
  filter_max_1 <- max(filter_values[[1]])
  filter_min_2 <- min(filter_values[[2]])
  filter_max_2 <- max(filter_values[[2]])

  interval_length_1 <- (filter_max_1 - filter_min_1) / (num_intervals[1] - (num_intervals[1] - 1) * percent_overlap/100 )
  interval_length_2 <- (filter_max_2 - filter_min_2) / (num_intervals[2] - (num_intervals[2] - 1) * percent_overlap/100 )

  step_size_1 <- interval_length_1 * (1 - percent_overlap/100)
  step_size_2 <- interval_length_2 * (1 - percent_overlap/100)

  num_levels <- num_intervals[1] * num_intervals[2]

  cover = list()

  level_indices_1 <- rep(1:num_intervals[1], num_intervals[2])
  level_indices_2 <- rep(1:num_intervals[2], each=num_intervals[1])
  # begin mapper main loop
  for (level in 1:num_levels) {
    level_1 <- level_indices_1[level]
    level_2 <- level_indices_2[level]

    min_value_in_level_1 <- filter_min_1 + (level_1 - 1) * step_size_1
    min_value_in_level_2 <- filter_min_2 + (level_2 - 1) * step_size_2
    max_value_in_level_1 <- min_value_in_level_1 + interval_length_1
    max_value_in_level_2 <- min_value_in_level_2 + interval_length_2

    cover[[level]] = c(min_value_in_level_1, max_value_in_level_1, min_value_in_level_2, max_value_in_level_2)
  } # end mapper main loop

  cover
}

