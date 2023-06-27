# -------------------------- Packages/setting up data -------------------------
library("tictoc")
library("symbols")
library("aplpack")
library("scatterplot3d")
library("corrplot")
library("plotly")
library("factoextra")
library("GGally")
library("ggcorrplot")

setwd("~/R/Multivariate Analysis/Country Performance and Region Correlation (PCA&Clusters)")

# read in data
raw_data <- read.csv("data.csv")
colnames(raw_data) <- c("Country", "Region","Surface area", "Population", 
                        "Population density", "Sex ratio", "GDP growth rate", 
                        "GDP per capita","Agriculture %", "Industry %", 
                        "Services %", "Unemployment", "Labour force sex ratio",
                        "Urban population %", "Fertility rate", 
                        "Av. life expectancy", "Age distribution", 
                        "Infant mortality rate", "CO2 emissions")
head(raw_data)


# remove missing values
which_na <- sapply(1:nrow(raw_data), function(i) any(raw_data[i,] == -99))
data <- raw_data[!which_na,]

which_na <- sapply(1:nrow(data), function(i) any(data[i,] == "..."))
data <- data[!which_na,]

which_na <- sapply(1:nrow(data), function(i) any(data[i,] == ".../..."))
data <- data[!which_na,]

head(data)


# convert string fractions (in data) to numbers
eval_data1 <- data.frame(sapply(data$`Labour force sex ratio`, 
                                function(i) eval(parse(text=i)))) # text->expression->eval

data$`Labour force sex ratio` <- eval_data1[,1]

eval_data2 <- data.frame(sapply(data$`Age distribution`, 
                                function(i) eval(parse(text=i))))
data$`Age distribution` <- eval_data2[,1]


head(data)

# store country and region of data as labels
country <- data$Country
country
region <- data$Region
region

# separate numeric data (used for analysis)
numeric_data <- data.frame(data[, -c(1,2)])

# remove unwanted variables
numeric_data$Population <- NULL
numeric_data$CO2.emissions <- NULL
numeric_data$GDP.growth.rate <- NULL

# name rows by country
row.names(numeric_data) <- country
head(numeric_data)


# check all is numeric
for (i in 1:ncol(numeric_data)) {
  numeric <- is.numeric(numeric_data[, i])
} 

if (numeric == TRUE) {
  print("All data is numeric")
}


head(numeric_data)

# GENERALISING REGIONS
possible_regions <- unique(region)
possible_regions


main_region <- region

# REGION ASSIGNMENT BY NAME (For Habillage)
# Asia
main_region <- replace(main_region, main_region=="WesternAsia", "Asia")
main_region <- replace(main_region, main_region=="EasternAsia", "Asia")
main_region <- replace(main_region, main_region=="NorthernAsia", "Asia")
main_region <- replace(main_region, main_region=="SouthernAsia", "Asia")
main_region <- replace(main_region, main_region=="South-easternAsia", "Asia")
main_region <- replace(main_region, main_region=="CentralAsia", "Asia")

# Europe
main_region <- replace(main_region, main_region=="NorthernEurope", "Europe")
main_region <- replace(main_region, main_region=="EasternEurope", "Europe")
main_region <- replace(main_region, main_region=="SouthernEurope", "Europe")
main_region <- replace(main_region, main_region=="WesternEurope", "Europe")

# Africa
main_region <- replace(main_region, main_region=="NorthernAfrica", "Africa")
main_region <- replace(main_region, main_region=="EasternAfrica", "Africa")
main_region <- replace(main_region, main_region=="SouthernAfrica", "Africa")
main_region <- replace(main_region, main_region=="WesternAfrica", "Africa")
main_region <- replace(main_region, main_region=="MiddleAfrica", "Africa")

# Oceania
main_region <- replace(main_region, main_region=="Micronesia", "Oceania")
main_region <- replace(main_region, main_region=="Melanesia", "Oceania")
main_region <- replace(main_region, main_region=="Polynesia", "Oceania")

# North America
main_region <- replace(main_region, main_region=="NorthernAmerica", "North America")
main_region <- replace(main_region, main_region=="CentralAmerica", "North America")
main_region <- replace(main_region, main_region=="Caribbean", "North America")

# South America
main_region <- replace(main_region, main_region=="SouthAmerica", "South America")







# ////////// ANALYSIS //////////

# PAIRS PLOT AND CORRELOGRAM
pairs(numeric_data, lower.panel = NULL)

correlation <- cor(numeric_data)
ggcorrplot(correlation)


# APPLYING PCA/Classical MDS and COMPARING TIMINGS (Note: results are identical: Z & mds)
tic("PCA")
X <- scale(numeric_data, center = TRUE, scale = TRUE)
pca <- prcomp(X)
Z <- X %*% pca$rotation
toc()

tic("CMDS")
D <- dist(scale(numeric_data, center = TRUE, scale = TRUE))
mds <- cmdscale(D, k = 14)
toc()


# CONTINUING WITH PCA (quicker in general)
summary(pca)
pca
# scree plot
plot(pca)

# 2D visualisation (general)
fviz_pca_biplot(pca, repel = TRUE, label = "ind")
fviz_pca_biplot(pca, repel = TRUE, label = "var")

# 2D visualisation (main regions)
fviz_pca_biplot(pca, repel = TRUE, label = "ind", habillage = main_region)

# 2D visualisation (specific regions)
fviz_pca_biplot(pca, repel = TRUE, label = "none", habillage = region)
fviz_pca_ind(pca, repel = TRUE, label = "none", habillage = km_s_clustering)


# principal component contributions
# first two
fviz_contrib(pca, choice = "var", axes = c(1, 2))
# first three
fviz_contrib(pca, choice = "var", axes = c(1, 2, 3))


# VISUALISING ONLY SPECIFIC REGIONS
plot(Z[main_region=="Asia",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "Asia"]))
plot(Z[main_region=="Europe",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "Europe"]))
plot(Z[main_region=="Africa",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "Africa"]))
plot(Z[main_region=="North America",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "North America"]))
plot(Z[main_region=="South America",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "South America"]))
plot(Z[main_region=="Oceania",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(region[main_region %in% "Oceania"]))


# 3D VISUALISATION
# not used because it does not give enough extra information, hard to represent
# in a report as well
Z <- X %*% pca$rotation[, c(1:3)] # n x v %*% v x 3 --> PCs are unit vectors so these are projections onto PC axes


plotly::plot_ly(data=as.data.frame(Z), x = ~PC1, y= ~PC2, z = ~PC3, 
                color = factor(main_region))
plotly::plot_ly(data=as.data.frame(Z), x = ~PC1, y= ~PC2, z = ~PC3, 
                color = factor(region))





# // CLUSTERING ANALYSIS//
# data for clustering analysis
Z <- X %*% pca$rotation[, c(1:3)]


# Removing outliers?
outliers <- c("Qatar", "United Arab Emirates")
X_o <- X[!(row.names(X) %in% outliers), ]
main_region_o <- main_region[!(row.names(X) %in% outliers)]
pca_o <- prcomp(X_o)

fviz_pca_biplot(pca_o, repel = TRUE, habillage = main_region_o)
                  
Z_o <- X %*% pca_o$rotation[, c(1, 2)]


# APPLYING K-Means
# K = number of general regions (6)
K <- length(unique(main_region))
km_g <- kmeans(Z, centers = K, nstart = 1000)
km_g_clustering <- km_g$cluster

# visualising
fviz_pca_biplot(pca, repel = TRUE, label = "none", habillage = km_g_clustering)

# K = number of specific regions (22)
K <- length(unique(region))
km_s <- kmeans(Z, centers = K, nstart = 1000)
km_s_clustering <- km_s$cluster

# visualising together
fviz_pca_biplot(pca, repel = TRUE, label = "none", habillage = km_s_clustering)

# VISUALISING SPECIFIC REGIONS SEPARATELY
plot(Z[main_region=="Asia",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "Asia"]))
plot(Z[main_region=="Europe",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "Europe"]))
plot(Z[main_region=="Africa",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "Africa"]))
plot(Z[main_region=="North America",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "North America"]))
plot(Z[main_region=="South America",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "South America"]))
plot(Z[main_region=="Oceania",], xlim=c(-6, 6), ylim=c(-3, 8), col = factor(km_s_clustering[main_region %in% "Oceania"]))

