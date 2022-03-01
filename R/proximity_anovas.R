# IluAg 2022: Proximity ANOVAs

library(dplyr)
library(rstatix)
library(ez)

# df <- read.csv("D:/IluAg/IluAg/means/proximity-congruency.txt", row.names=NULL, sep="")
df <- read.csv("~/Uni/IluAg/means/proximity-congruency.txt", row.names=NULL, sep="")

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
df[, "congruency"] <-
  NA

for (iRow in 1:length(df$condition)) {
  if (df[iRow, "condition"] == "CLOSECONc") {
    df[iRow, "proximity"] <- 1
    df[iRow, "congruency"] <- 1
  } else if (df[iRow, "condition"] == "CLOSEINCc") {
    df[iRow, "proximity"] <- 1
    df[iRow, "congruency"] <- 0
  } else if (df[iRow, "condition"] == "FARCONc") {
    df[iRow, "proximity"] <- 0
    df[iRow, "congruency"] <- 1
  } else if (df[iRow, "condition"] == "FARINCc") {
    df[iRow, "proximity"] <- 0
    df[iRow, "congruency"] <- 0
  } 
}    

# Drop unnecessary columns
keepers <- c("condition", "channel", "window", "mean", "subject", "proximity", "congruency")
df <- df[keepers]

# Convert to appropriate data types
df$condition <- as.factor(df$condition)
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.factor(df$subject)
df$proximity <- as.factor(df$proximity)
df$congruency <- as.factor(df$congruency)
df$mean <- as.numeric(df$mean)

# Run four ANOVAs for the different components of interest
N1_Fz <- subset(df, channel == "Fz" & window == "N1") 
N1_Cz <- subset(df, channel == "Cz" & window == "N1") 

P2_Fz<- subset(df, channel == "Fz" & window == "P2")
P2_Cz<- subset(df, channel == "Cz" & window == "P2")


P3a_Fz <- subset(df, channel == "Fz" & window == "P3")
P3a_Cz <- subset(df, channel == "Cz" & window == "P3")

P3b <- subset(df, channel == "Pz" & window == "P3")

# N1 
data <- N1_Cz

# Summary statistics
summary <- data %>%
  group_by(subject, proximity, congruency) %>%
  get_summary_stats(mean, type = "mean")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV

anova = ezANOVA(
  data = data
  , dv = .(mean)
  , wid = .(subject)
  , within = .(proximity, congruency)
)

# P2 
data <- P2_Cz

# Summary statistics
summary <- data %>%
  group_by(subject, proximity, congruency) %>%
  get_summary_stats(mean, type = "mean")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV

anova = ezANOVA(
  data = data
  , dv = .(mean)
  , wid = .(subject)
  , within = .(proximity, congruency)
)

# P3a
data <- P3a_Cz

# Summary statistics
summary <- data %>%
  group_by(subject, proximity, congruency) %>%
  get_summary_stats(mean, type = "mean")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV

anova = ezANOVA(
  data = data
  , dv = .(mean)
  , wid = .(subject)
  , within = .(proximity, congruency)
)

# P3b
data <- P3b

# Summary statistics
summary <- data %>%
  group_by(subject, proximity, congruency) %>%
  get_summary_stats(mean, type = "mean")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV

anova = ezANOVA(
  data = data
  , dv = .(mean)
  , wid = .(subject)
  , within = .(proximity, congruency)
)

library(ggplot2)

ggplot(P3a, aes(x=subject, y=mean)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
