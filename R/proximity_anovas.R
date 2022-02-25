# IluAg 2022: Proximity ANOVAs

library(dplyr)
library(rstatix)

df <- read.csv("D:/IluAg/IluAg/means/proximiy.txt", row.names=NULL, sep="")

# Rename the columns
names <- c("experiment", "condition", "channel", "window", "mean", "subject")
colnames(df) <- names

# Clean the content of the last column
df$subject = substr(df$subject,1,2)

# Rename time windows to component names
df$window[df$window == "+100..+120"] <- "N1"
df$window[df$window == "+180..+220"] <- "P2"
df$window[df$window == "+310..+370"] <- "P3"

# Reshape the data
df[, "proximity"] <-
  NA

for (iRow in 1:length(df$condition)) {
  if (df[iRow, "condition"] == "TWonec") {
    df[iRow, "proximity"] <- 0
  } else if (df[iRow, "condition"] == "TWtwoc") {
    df[iRow, "proximity"] <- 1
  } else if (df[iRow, "condition"] == "TWthreec") {
    df[iRow, "proximity"] <- 1
  } else if (df[iRow, "condition"] == "TWfourc") {
    df[iRow, "proximity"] <- 0
  } 
}    

# Drop unnecessary columns
keepers <- c("condition", "channel", "window", "mean", "subject", "proximity")
df <- df[keepers]

# Convert to appropriate data types
df$condition <- as.factor(df$condition)
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.factor(df$subject)
df$proximity <- as.factor(df$proximity)
df$mean <- as.numeric(df$mean)

# Run four ANOVAs for the different components of interest
N1 <- subset(df, channel == "Fz" & window == "N1") 
P2<- subset(df, channel == "Fz" & window == "P2")
P3a <- subset(df, channel == "Fz" & window == "P3")
P3b <- subset(df, channel == "Pz" & window == "P3")

# N1 
data <- N1

# Summary statistics
summary <- data %>%
  group_by(subject, proximity) %>%
  get_summary_stats(mean, type = "mean")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV
aov <- anova_test(data = data, dv = mean, wid = subject, within = proximity)

# Get results
N1_aov <- get_anova_table(aov, correction = "none")
