# IluAg 2022: Order ANOVAs

library(dplyr)
library(rstatix)
library(ez)

df <- read.csv("D:/IluAg/IluAg/means/means10.txt", row.names=NULL, sep="")

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
df[, "order"] <-
  NA
df[, "congruency"] <-
  NA

for (iRow in 1:length(df$condition)) {
  if (df[iRow, "condition"] == "POSCONc") {
    df[iRow, "order"] <- 1
    df[iRow, "congruency"] <- 1
  } else if (df[iRow, "condition"] == "POSINCc") {
    df[iRow, "order"] <- 1
    df[iRow, "congruency"] <- 0
  } else if (df[iRow, "condition"] == "NEGCONc") {
    df[iRow, "order"] <- 0
    df[iRow, "congruency"] <- 1
  } else if (df[iRow, "condition"] == "NEGINCc") {
    df[iRow, "order"] <- 0
    df[iRow, "congruency"] <- 0
  } 
}    

# Drop unnecessary columns
keepers <- c("condition", "channel", "window", "mean", "subject", "order", "congruency")
df <- df[keepers]

# Convert to appropriate data types
df$condition <- as.factor(df$condition)
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.factor(df$subject)
df$order <- as.factor(df$order)
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

# Write here which test you want to do 
data <- P3b

# Summary statistics
summary <- data %>%
  group_by(order, congruency) %>%
  get_summary_stats(mean, type = "mean_se")

# Run ANOVA
data <- ungroup(data) # Ungroup so that aov can work# Run AOV

anova = ezANOVA(
  data = data
  , dv = .(mean)
  , wid = .(subject)
  , within = .(order, congruency)
)
anova

# Plot interaction of order and congruency in P3b

ggplot(summary, aes(y=mean, x=order, fill = congruency)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) + 
  scale_x_discrete(name = "Order", labels=c("sound before BP", "BP before sound")) + 
  scale_y_continuous((name = "Mean amplitude of P3b")) +
  scale_fill_discrete(name = "Congruency", labels=c("incongruent", "congruent")) + 
  theme_bw()
