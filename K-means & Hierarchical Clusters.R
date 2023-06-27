# Clustering analysis

# Example 1
# Simulated data:
n <- 200
sig <- 3
# x: sample 100 from mean of 0, 100 from mea of 10, both sd = sig
x <- rnorm(n, mean = rep(c(0, 10), each = n/2), sd = sig)
# y: sample 200 from mean 0, sd = sig
y <- rnorm(n, mean = 0, sd = sig)
# Creates two gaussian clusters centered at x = 0 and x = 10, about y = 0

dat1 <- data.frame(x = x, y = y)
plot(dat1, asp = 1)


# /// K-means clustering ///
km1 <- kmeans(dat1, centers = 2, instart = 50)
km1
# centers = chosen number of clusters
# instart = how many iterations to do, keeping the best one

# cluster assignment for each data point
km1$cluster

# plot, coloring points by cluster assignment
plot(dat1, asp = 1, col = km1$cluster)


# /// Hierarchical clustering ///
hc1 <- hclust(d = dist(dat1), method = "complete")

# dendrogram (just plot)
plot(hc1)

# cluster assignment
cutree(tree = hc1, k = 2) #takes tree input, then assigns points to k clusters via tree

# visualisation
plot(dat1, asp = 1, col = cutree(tree = hc1, k = 2))


# K selection
factoextra::fviz_nbcluster(...)
