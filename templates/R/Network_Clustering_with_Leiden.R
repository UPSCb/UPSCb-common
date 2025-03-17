# Load the libraries we need to calculate the clusters and plot
library(igraph)
library(ggplot2)

# Read in the edge list you created:
edgelist <- read.table("your_edgelist.tsv", header = TRUE, sep = "\t")

edgelist_standardized <- edgelist %>%
  rowwise() %>%
  mutate(
    Source = min(Source, Target),
    Target = max(Source, Target)
  ) %>%
  ungroup()

# Remove duplicate rows
edgelist_unique <- edgelist_standardized %>%
  distinct()

edgelist <- edgelist_unique
# The algorithm below produces a result based on a "seed" in the background. To maintain reproducibility we
# can set a seed so we get the exact same result for each dataset every time we run it. If we do not do this
# the seed will be randomly generate and the results could slightly vary each run.
set.seed(123) 

# Create an igraph graph from the sampled edge list
g <- graph_from_data_frame(edgelist, directed = FALSE)

# Apply the Leiden algorithm
leiden_clusters <- cluster_leiden(g, objective_function = "modularity")

# Access the clusters
clusters <- membership(leiden_clusters)

# Count the number of nodes in each cluster
cluster_counts <- table(clusters)

# Create a barplot for checking the overall cluster distribution
barplot(cluster_counts, 
        main = "Cluster Composition", 
        xlab = "Cluster", 
        ylab = "Number of Genes", 
        col = rainbow(length(cluster_counts)))

# Get the node names by accessing the network attribute with V
gene_names <- V(g)$name

# Convert clusters to a vector
cluster_vector <- as.vector(clusters)

# Create a data frame with node names and cluster memberships
gene_cluster_df <- data.frame(Gene = gene_names, Cluster = cluster_vector)
gene_cluster_df$Cluster <- as.character(gene_cluster_df$Cluster)

# Print the first 6 lines of the data frame to check the output
print(head(gene_cluster_df))
