# --------------------------- Get data/Initialise ------------------------------

# load data
load("~/R/win-library/4.1/Projects/Assignment 2/assignment-2-redpanda/data/FishLengths.RData")
data_real <- x

# Assigning initial groupings to data: (based on data collected for 100 fish)
# initial grouping of fish from 0-30 (group 1), 30-55 (group 2), 55+ (group 3)
for (i in 1:1000) {
  if (data_real$Length[i] < 30) {
    data_real$Age[i] <- 1
  } else if ((data_real$Length[i] >=30) & (data_real$Length[i] < 55)) {
    data_real$Age[i] <- 2
  } else if ((data_real$Length[i] >=55) & (data_real$Length[i] < 100)) {
    data_real$Age[i] <- 3
  }
}
