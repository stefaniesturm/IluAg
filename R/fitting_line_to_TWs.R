# Fit a line to TW means

library(ggplot2)

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

# Drop unnecessary columns
keepers <- c("condition", "channel", "window", "mean", "subject")
df <- df[keepers]

# Convert to appropriate data types
df$condition <- as.factor(df$condition)
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.factor(df$subject)
df$mean <- as.numeric(df$mean)

P3a <- subset(df, channel == "Fz" & window == "P3")

ggplot(P3a, aes(x=condition, y=mean, group=1)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)

# Normalise the data per subject
P3a$subject <- as.numeric(P3a$subject)

df_norm <- data.frame()
for (iSub in 1:26) {
  tmp <- subset(P3a, subject == iSub)
  tmp$mean_norm <- (tmp$mean - mean(tmp$mean)) / sd(tmp$mean)
  df_norm <- rbind(df_norm, tmp)
}

ggplot(df_norm, aes(x=condition, y=mean_norm, group=1)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
