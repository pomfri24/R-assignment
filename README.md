# R-assignment

library(tidyverse)

#Dataset-----------------------------------------------------------------
dvst <- read_csv("https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-08-r/Dataset_S1.txt")

#Data inspection---------------------------------------------------------
# Inspect the structure of the dataset
str(dvst)

# View first few rows
head(dvst)

#View last few rows
tail(dvst)

# Summary statistics
summary(dvst)

# Check dimensions
dim(dvst)

# Data processing ---------------------------------------------------------
# Rename relevant columns if needed
colnames(dvst) <- c("SNP_ID", "Chromosome", "Position", "Genotype")

# Convert chromosome column to a factor for better handling
dvst$Chromosome <- as.factor(dvst$Chromosome)

# Adjust based on actual column count
colnames(dvst) <- c("SNP_ID", "Chromosome", "Position", "Genotype")

#If you’re unsure how many columns there are, use
names(dvst) <- make.names(names(dvst), unique = TRUE)

#Verify Structure Again
str(dvst)
head(dvst)

#Convert to Tibble (Fix Hidden Data Issues)
dvst <- as_tibble(dvst)

# Sort by Chromosome and Position
dvst_sorted <- dvst %>% arrange(Chromosome, Position)

# Check for missing values
sum(is.na(dvst))

# Replace missing genotypes with '?'
dvst_cleaned <- dvst %>%
  mutate(Genotype = ifelse(is.na(Genotype), "?", Genotype))

# Data visualization ------------------------------------------------------

#SNPs per Chromosome:
  
  # Plot total SNPs across categories (assuming this represents SNP counts for different categories)
  ggplot(dvst, aes(x = factor(1), y = `total SNPs`, fill = factor(1))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Total SNPs Across Categories", x = "Category", y = "Total SNPs") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Homozygosity vs. Heterozygosity:
  
  # Assuming you can create a column for Heterozygosity
  dvst <- dvst %>%
  mutate(Heterozygosity = ifelse(`Heterozygosity` > 0.5, "Heterozygous", "Homozygous"))

# Plot the distribution of homozygous vs heterozygous sites
ggplot(dvst, aes(x = Heterozygosity, fill = Heterozygosity)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Homozygous vs Heterozygous Sites", x = "Genotype", y = "Proportion")

#Missing Data:
  
  # Create a new column to indicate missing data
  dvst <- dvst %>%
  mutate(Missing = ifelse(is.na(`SNPs`), "Missing", "Not Missing"))

# Plotting missing data for each category
ggplot(dvst, aes(x = factor(1), fill = Missing)) +
  geom_bar(stat = "count", position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Missing Data", x = "Category", y = "Proportion") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Distribution of Heterozygosity:
  
  # Plotting the distribution of heterozygosity
  ggplot(dvst, aes(x = `Heterozygosity`, fill = factor(1))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of Heterozygosity", x = "Heterozygosity", y = "Density")

#Your Own Visualization:
  
  # Plotting Theta across samples
  ggplot(dvst, aes(x = `Theta`, fill = factor(1))) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Theta (Nucleotide Diversity)", x = "Theta", y = "Count")
