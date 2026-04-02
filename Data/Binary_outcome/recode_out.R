# Read continuous data
pheno <- read.table("trait1_analysis_out.txt", header = TRUE)

# Determine the threshold for the 90th percentile
threshold_90 <- quantile(pheno$trait1, probs = 0.7)

# Recode data into categorical variable
pheno$trait1 <- ifelse(pheno$trait1 >= threshold_90, 2, 1)

# Check the distribution of the new categorical variable
table(pheno$trait1)

# If you need to ensure exactly 10% are 1 and 90% are 0:
num_1 <- sum(pheno$trait1 == 2)
num_0 <- sum(pheno$trait1 == 1)

# Calculate the difference between actual and desired counts
difference <- num_1 - 0.3 * (num_1 + num_0)

# If there are more 1s than desired, randomly select some to change to 0
if (difference > 0) {
  indices_1 <- which(pheno$trait1 == 2)
  indices_to_change <- sample(indices_1, size = difference)
  pheno$phent[indices_to_change] <- 1
}

# If there are fewer 1s than desired, randomly select some to change to 1
if (difference < 0) {
  indices_0 <- which(pheno$trait1 == 1)
  indices_to_change <- sample(indices_0, size = abs(difference))
  pheno$phent[indices_to_change] <- 2
}

# Check the distribution again to ensure it's approximately 10% 1s and 90% 0s
table(pheno$trait1)
write.table(pheno, "trait1_analysis_out_b.txt", col.names = TRUE, row.names = F, quote = F)
#########################################################################
########################################################################
############## recoding the target dataset #############################
# Read continuous data
pheno <- read.table("trait1_PRStrained_out.txt", header = TRUE)

# Determine the threshold for the 90th percentile
threshold_90 <- quantile(pheno$trait1, probs = 0.7)

# Recode data into categorical variable
pheno$trait1 <- ifelse(pheno$trait1 >= threshold_90, 2, 1)

# Check the distribution of the new categorical variable
table(pheno$trait1)

# If you need to ensure exactly 10% are 1 and 90% are 0:
num_1 <- sum(pheno$trait1 == 2)
num_0 <- sum(pheno$trait1 == 1)

# Calculate the difference between actual and desired counts
difference <- num_1 - 0.3 * (num_1 + num_0)

# If there are more 1s than desired, randomly select some to change to 0
if (difference > 0) {
  indices_1 <- which(pheno$trait1 == 2)
  indices_to_change <- sample(indices_1, size = difference)
  pheno$phent[indices_to_change] <- 1
}

# If there are fewer 1s than desired, randomly select some to change to 1
if (difference < 0) {
  indices_0 <- which(pheno$trait1 == 1)
  indices_to_change <- sample(indices_0, size = abs(difference))
  pheno$phent[indices_to_change] <- 2
}

# Check the distribution again to ensure it's approximately 10% 1s and 90% 0s
table(pheno$trait1)
write.table(pheno, "trait1_PRStrained_out_b.txt", col.names = TRUE, row.names = F, quote = F)