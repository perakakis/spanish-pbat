# Calculate link communities and produce network plots

# Load necessary libraries
library(igraph)
library(linkcomm)
library(RColorBrewer)
library(reshape2)

# Set seed for reproducibility
set.seed(1)

# Load Boruta results
df <- read.csv("./results/Boruta_results_es.csv")
plot_title <- "Spain"
filename <- "./figures/LinkCommunities_es.pdf"

long_df <- melt(df, id.vars = "Predictor")
colnames(long_df)[c(2, 3)] <- c("Outcome", "Rank")

# Set importance threshold
rank_threshold <- 3
predictors_as_pies <- FALSE
size_as_degree <- TRUE
lc_cutat <- FALSE
# Keep important predictors
all_results <- long_df[long_df$Rank <= rank_threshold & !is.na(long_df$Rank), ]

# Invert the rank to make it a weight
all_results$Rank <- rank_threshold - all_results$Rank + 1

# Create a graph from the edge list
g <- graph_from_data_frame(all_results, directed = FALSE)
V(g)$type <- V(g)$name %in% unique(all_results$Outcome)

# Convert graph to edge list
edgelist <- as.data.frame(as_edgelist(g))

# Compute network centrality statistics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, directed = FALSE)
V(g)$closeness <- closeness(g, mode = "all", normalized = TRUE)
V(g)$eigenvector <- eigen_centrality(g)$vector

# Print centrality measures
centrality_df <- data.frame(
  Node = V(g)$name,
  Degree = V(g)$degree,
  Betweenness = V(g)$betweenness,
  Closeness = V(g)$closeness,
  Eigenvector = V(g)$eigenvector
)
# print(centrality_df)

# Run Link Communities (detect overlapping communities)
lc <- getLinkCommunities(edgelist,
  hcmethod = "average",
  bipartite = TRUE,
  use.all.edges = FALSE,
  plot = FALSE,
  verbose = FALSE
)

# plotLinkCommGraph(lc)

# Create edge labels from the edgelist
edge_labels <- paste(edgelist[, 1], edgelist[, 2], sep = " – ")

# Assign the edge labels to lc$hclust$labels
lc$hclust$labels <- edge_labels

# par(mar = c(1, 5, 1, 1)) # Adjust margins (bottom, left, top, right)
# plot(lc$hclust, hang = -1, main = "Dendrogram of Link Communities")
# rect.hclust(lc$hclust, h = lc$pdmax, border = "red")

# Cut the dendrogram if lc_cutat is set to a value between 0 and 1
if (lc_cutat) {
  lc <- newLinkCommsAt(lc, cutat = lc_cutat)
}

# Get community centrality
# print(getCommunityCentrality(lc))

# --- Community Index Remapping ---
edge2communities <- lapply(1:ecount(g), function(e) {
  which(sapply(lc$clusters, function(comm) e %in% comm))
})

# Get unique community indices and create sequential mapping
all_comm_indices <- unique(unlist(edge2communities))
if (length(all_comm_indices) > 0) {
  comm_mapping <- setNames(seq_along(all_comm_indices), all_comm_indices)
  edge2communities <- lapply(edge2communities, function(comms) {
    as.numeric(comm_mapping[as.character(comms)])
  })
  num_communities <- length(all_comm_indices)
} else {
  stop("No communities detected. Check network connectivity.")
}

# Initialize matrices with correct dimensions
outcome_comm_dist <- matrix(0, nrow = vcount(g), ncol = num_communities)
predictor_comm_dist <- matrix(0, nrow = vcount(g), ncol = num_communities)
all_nodes_contributions <- matrix(0, nrow = vcount(g), ncol = num_communities)

# --- Calculate Outcome Communities ---
for (node in V(g)$name[V(g)$type]) {
  node_id <- which(V(g)$name == node)
  edges <- as.numeric(incident(g, node_id))

  valid_edges <- edges[sapply(edge2communities[edges], length) > 0]
  if (length(valid_edges) == 0) next

  for (e in valid_edges) {
    comms <- edge2communities[[e]]
    if (any(comms < 1 | comms > num_communities)) {
      stop(paste("Invalid community index in edge", e))
    }
    w_contribution <- 1 / length(valid_edges)
    outcome_comm_dist[node_id, comms] <- outcome_comm_dist[node_id, comms] +
      w_contribution
  }
}

# --- Calculate Predictor Communities ---
for (predictor in V(g)$name[!V(g)$type]) {
  pred_id <- which(V(g)$name == predictor)
  edges <- as.numeric(incident(g, pred_id, mode = "out"))

  valid_edges <- edges[sapply(edge2communities[edges], length) > 0]
  predictor_comm_dist[pred_id, ] <- 0

  for (e in valid_edges) {
    outcome <- ends(g, e)[2]
    outcome_id <- which(V(g)$name == outcome)
    comms <- edge2communities[[e]]

    # Validate dimensions
    if (ncol(outcome_comm_dist) != ncol(predictor_comm_dist)) {
      stop("Community dimension mismatch
      between outcome and predictor matrices")
    }
    if (predictors_as_pies) {
      comm_contrib <- outcome_comm_dist[outcome_id, ]
      predictor_comm_dist[pred_id, ] <- predictor_comm_dist[pred_id, ] +
        comm_contrib
    } else {
      w_contribution <- 1 / length(valid_edges)
      predictor_comm_dist[pred_id, comms] <- predictor_comm_dist[pred_id, comms] +
        w_contribution
    }
  }
}

all_nodes_contributions[V(g)$type, ] <- outcome_comm_dist[V(g)$type, ]
all_nodes_contributions[!V(g)$type, ] <- predictor_comm_dist[!V(g)$type, ]

# --- Visualization ---
community_colors <- brewer.pal(
  max(num_communities, 3), "Set2"
)[1:num_communities]

# Get community indices for each node's non-zero contributions
vertex_communities <- lapply(1:vcount(g), function(i) {
  which(all_nodes_contributions[i, ] > 0)
})

# Get pie chart colors for each node
vertex_pie_colors <- lapply(vertex_communities, function(comms) {
  community_colors[comms] # Preserves community order
})

# Create pie values
vertex_pie <- lapply(1:vcount(g), function(j) {
  comms <- vertex_communities[[j]] # Get communities for node j
  all_nodes_contributions[j, comms] # Get contributions for node j
})

# Determine dominant color for non-pie nodes
vertex_colors <- sapply(vertex_communities, function(comms) {
  if (length(comms) == 0) {
    "lightgray"
  } else if (length(comms) == 1) {
    community_colors[comms] # Direct color match
  } else {
    NA # Pies will use pie.colors
  }
})

# Determine pie vs. circle nodes
is_pie <- sapply(vertex_communities, length) > 1

# Determine node size based on degree
if (size_as_degree) {
  vertex_sizes <- degree(g) * 4
} else {
  vertex_sizes <- rep(10, vcount(g))
}

# Determine frame color for outcomes and predictors
# vertex_frame_colors <- ifelse(V(g)$type, "black", "white")

pdf(filename, width = 15, height = 10)
par(mar = c(4, 4, 2, 2), cex.main = 3)
lay <- layout_with_kk(g)

# Plot network
plot(g,
  layout = lay,
  vertex.shape = ifelse(is_pie, "pie", "circle"),
  vertex.size = vertex_sizes,
  vertex.label = V(g)$name,
  vertex.color = vertex_colors, # NA lets pie colors show
  # vertex.frame.color = vertex_frame_colors,
  # vertex.frame.color = "white",
  # vertex.frame.width = 3,
  vertex.pie = vertex_pie,
  vertex.label.cex = 1.5,
  vertex.pie.color = vertex_pie_colors,
  edge.width = all_results$Rank * 1.5,
  edge.color = sapply(1:ecount(g), function(e) {
    comms <- edge2communities[[e]]
    if (length(comms) > 0) community_colors[comms[1]] else "gray"
  }),
  main = plot_title
)

# Add legend
legend("topright",
  legend = paste("Community", 1:num_communities),
  pch = 21,
  pt.bg = community_colors,
  pt.cex = 2,
  bty = "n"
)

dev.off()