library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(purrr)
library(bioDist)
library(dimRed)
library(TDAmapper)
library(umap)
library(dbscan)
library(RayleighSelection)

create_mut_object <- function(exp_table, mut_table, filter_method, max_interval, max_percent_overlap, var_threshold, num_cores=1) {


  # FUNCTIONS ---------------------------------------------------------------

  plot.umap = function(emb) {
    ggplot(emb) +
      geom_point(aes(x=V1, y=V2, color=0), size=0.9)
  }

  build.mapper = function(dist, umap_emb, num_intervals, percent_overlap) {
    mapper2D(dist, umap_emb, num_intervals, percent_overlap)
  }

  plot.mapper = function(mapper_obj) {

    adj_graph <- graph.adjacency(mapper_obj$adjacency, mode="undirected")

    V(adj_graph)$size <- log2(as.numeric(lapply(mapper_obj$points_in_vertex, length))+1)
    E(adj_graph)$size <- rep(0, gsize(adj_graph))
    V(adj_graph)$label <- ""
    V(adj_graph)$frame.color <- "black"

    plot(adj_graph, layout.auto(adj_graph))

  }

  average.expression.for.node = function(points_in_vertex, expression) {
    M <- Matrix(0, nrow=length(points_in_vertex), ncol=ncol(expression), dimnames=list(NULL, colnames(expression)))

    for (i in 1:length(points_in_vertex)) {
      points = points_in_vertex[[i]]
      if (length(points) == 1) {
        M[i, ] <- expression[points, ]
      }
      else {
        M[i, ] <- apply(expression[points, ], 2, mean)
      }
    }

    M
  } # Avg data across all points in node and assigns node avg value

  collapse = function(adj,mapper_obj) {
    collapsed_adj <- adj
    for (i in 1:nrow(adj)) {
      for (j in 1:ncol(adj)) {
        if (adj[i, j] == 1 && length(mapper_obj$points_in_vertex[[j]]) == 1) {
          collapsed_adj[i, j] <- 0 # moving to here to always remove
          for (k in 1:ncol(adj)) {
            if (adj[j, k] == 1) {
              # ok. so remove i -> j and j -> k and add i -> k?
              collapsed_adj[j, k] <- 0
              if (i != k) {
                collapsed_adj[i, k] <- 1
              }
            }
          }
        }
      }
        # if (adj[i, j] == 1 && length(mapperobj$points_in_vertex[[i]]) == 1) { # need to do the other way too right
        #   collapsed_adj[i, j] <- 0 # moving to here to always remove
        #
        #   # this is new and i'm not sure it's correct
        #   for (k in 1:ncol(adj)) {
        #     if (adj[k, i] == 1) {
        #       collapsed_adj[k, i] <- 0
        #       if (j != k) {
        #         collapsed_adj[k, j] <- 1
        #
        #       }
        #     }
        #   }
        # }

    }

    g_collapsed <- graph_from_adjacency_matrix(collapsed_adj)

    g_collapsed$level_of_vertex <- mapper_obj$level_of_vertex
    g_collapsed$points_in_vertex <- mapper_obj$points_in_vertex
    g_collapsed$points_in_level <- mapper_obj$points_in_level
    g_collapsed$level_of_vertex <- mapper_obj$level_of_vertex
    g_collapsed$vertices_in_level <- mapper_obj$vertices_in_level

    g_collapsed
  } # Remove redundancies in connectivity due to singletons

  build.cover = function(umap_emb, num_intervals, percent_overlap) {
    filter_min_1 <- min(umap_emb[[1]])
    filter_max_1 <- max(umap_emb[[1]])
    filter_min_2 <- min(umap_emb[[2]])
    filter_max_2 <- max(umap_emb[[2]])

    interval_length_1 <- (filter_max_1 - filter_min_1) / (num_intervals[1] - (num_intervals[1] - 1) * percent_overlap/100 )
    interval_length_2 <- (filter_max_2 - filter_min_2) / (num_intervals[2] - (num_intervals[2] - 1) * percent_overlap/100 )

    step_size_1 <- interval_length_1 * (1 - percent_overlap/100)
    step_size_2 <- interval_length_2 * (1 - percent_overlap/100)

    num_levels <- num_intervals[1] * num_intervals[2]

    cover = list()

    level_indices_1 <- rep(1:num_intervals[1], num_intervals[2])
    level_indices_2 <- rep(1:num_intervals[2], each=num_intervals[1])

    for (level in 1:num_levels) {
      level_1 <- level_indices_1[level]
      level_2 <- level_indices_2[level]

      min_value_in_level_1 <- filter_min_1 + (level_1 - 1) * step_size_1
      min_value_in_level_2 <- filter_min_2 + (level_2 - 1) * step_size_2
      max_value_in_level_1 <- min_value_in_level_1 + interval_length_1
      max_value_in_level_2 <- min_value_in_level_2 + interval_length_2

      cover[[level]] = c(min_value_in_level_1, max_value_in_level_1, min_value_in_level_2, max_value_in_level_2)
    }

    cover.df = as.data.frame(do.call(rbind, cover))
    plot.umap(umap_emb) + geom_rect(data=cover.df, mapping=aes(xmin=V1, xmax=V2, ymin=V3, ymax=V4),
                                    color="black", alpha=0.05)

  } # Visualize grid of mapper parameters on embedding

  # INPUT AND CLEAN ---------------------------------------------------------

    # Cleaning tables
  #exp_table <- (read.csv("C:/Users/Adam/Desktop/Camara Lab/TDA_TCGA_Project/Test Data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #Laptop
  exp_table <- (read.csv("/home/rstudio/documents/Test Data/LGG_Full_TPM_matrix.csv", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?"))) #Lab
  #mut_table <- (read.csv("blah", row.names=1, header=T, stringsAsFactors=F, na.strings=c("NA","NaN", " ", "?")))

  exp_table <- exp_table[!duplicated(rownames(exp_table)),!duplicated(colnames(exp_table))] # No duplicated genes/samples

    # Only if mutation data given
  # mut_table <- mut_table[!duplicated(rownames(mut_table)),!duplicated(colnames(mut_table))]
  #
  # if(all(colnames(exp_table) %in% colnames(mut_table))) { #Gene names are same CHANGE
  #   if(all(colnames(mut_table) %in% colnames(exp_table))) {
  #     print("Matching genes between expression and mutation tables")
  #   }
  # } else {
  #   print("Unmatched genes between expression and mutation tables")
  # }

  exp_table_top <- exp_table[,order(-apply(exp_table,2,var))][,1:var_threshold]
  dist_matrix <- as.matrix(cor.dist(as.matrix(exp_table_top)))

    # MAPPER PLOTS ---------------------------------------------------------

    # Initialize parameters for several Mapper complexes

  num_intervals <- max_interval/10
  #num_y_intervals <- max_y_intervals/10
  start_interval <- max_interval - 10*(num_intervals - 1)
  #start_y_interval <- max_y_interval - 10*(num_y_intervals - 1)
  interval_range <- seq(start_interval, max_interval, length.out = num_intervals)
  #y_interval_range <- seq(start_y_interval, max_y_interval, length.out = num_y_intervals)

  num_percents <- max_percent_overlap/10
  start_percent <- max_percent_overlap - 10*(num_percents - 1)
  percent_range <- seq(start_percent, max_percent_overlap, length.out = num_percents)

  count = 0

  mapper_complexes <- vector('list',length(num_percents*num_intervals))


    # Filter Functions

  if (filter_method == "UMAP") {

    emb = umap(dist_matrix)$layout %>% as.data.frame
    plot.umap(emb)

  }

    else if (filter_method == "PCA") {
      emb <- autoplot(prcomp(dist_matrix))$data[,1:2]
      autoplot(prcomp(dist_matrix))
    }

      else {
        print("Not a valid filter method")
      }

  for (i in interval_range){
    for (p in percent_range){
      count = count + 1
      print(paste0("Creating Mapper Complex ", count, " of ", num_intervals*num_percents))

      mapperObj <- capture.output ({ build.mapper(dist_matrix, emb, c(i,i), p) })
      mapper_complexes[[count]] <- mapperObj

    }
  }

  class(mapper_complexes) <- "DriverMut_TDA"

  #return(capture.output({ mapper_complexes }))

  #   # Laplacian eigenmap
  # leim <- LaplacianEigenmaps()
  # lap_emb <- leim@fun(as((exp_table_top), "dimRedData"), leim@stdpars)
  # m_lap <- build.mapper(dist_matrix,list(lap_emb@data@data[,1], lap_emb@data@data[,2]), c(20,20), 50)
  # g_lap <- nerve_complex(m_lap$points_in_vertex)
  # plot_skeleton(g_lap)

}


suppressWarnings({
  suppressMessages ({
    library(TDAmapper, quietly=T, warn.conflicts=FALSE)
  })
})
