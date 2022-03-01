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
df$condition <- factor(df$condition, levels=c("TWonec", "TWtwoc", "TWthreec", "TWfourc"))
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.factor(df$subject)
df$mean <- as.numeric(df$mean)

P3b <- subset(df, channel == "Pz" & window == "P3")

ggplot(P3b, aes(x=condition, y=mean, group=1)) +
  geom_point() +
  geom_smooth(se=FALSE, method="loess")


ggplot(P3b, aes(x=condition, y=mean, group=1)) +
  geom_point() +
  geom_smooth()

# Normalise the data per subject
P3a$subject <- as.numeric(P3a$subject)

df_norm <- data.frame()
for (iSub in 1:26) {
  tmp <- subset(P3b, subject == iSub)
  tmp$mean_norm <- (tmp$mean - mean(tmp$mean)) / sd(tmp$mean)
  df_norm <- rbind(df_norm, tmp)
}

ggplot(df_norm, aes(x=condition, y=mean_norm, group=1)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(name = "Mean amplitude normalized by z-score") +
  scale_x_discrete(name = "Distance sound-BP in ms", labels = c("-120 to -60", "-60 to 0", "0 to 60", "60 to 120"))

# Different time window (based on proximity difference in ERP)

df <- read.csv("D:/IluAg/IluAg/means/tws-p2.txt", row.names=NULL, sep="")

# Rename the columns
names <- c("experiment", "condition", "channel", "window", "mean", "subject")
colnames(df) <- names

# Clean the content of the last column
df$subject = substr(df$subject,1,2)

# Rename time windows to component names
df$window[df$window == "+210..+290"] <- "P2"

# Drop unnecessary columns
keepers <- c("condition", "channel", "window", "mean", "subject")
df <- df[keepers]

# Convert to appropriate data types
df$condition <- factor(df$condition, levels=c("TWonec", "TWtwoc", "TWthreec", "TWfourc"))
df$channel <- as.factor(df$channel)
df$window <- as.factor(df$window)
df$subject <- as.numeric(df$subject)
df$subject <- as.factor(df$subject)
df$mean <- as.numeric(df$mean)

P2 <- subset(df, channel == "Cz" & window == "P2")

ggplot(P2, aes(x=condition, y=mean, group=1)) +
  geom_point() +
  geom_smooth(method='loess')

P2_norm <- data.frame()
for (iSub in 1:26) {
  tmp <- subset(P2, subject == iSub)
  tmp$mean_norm <- (tmp$mean - mean(tmp$mean)) / sd(tmp$mean)
  P2_norm <- rbind(P2_norm, tmp)
}

ggplot(P2_norm, aes(x=condition, y=mean_norm, group=1)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(name = "Mean amplitude normalized by z-score") +
  scale_x_discrete(name = "Distance sound-BP in ms", labels = c("-120 to -60", "-60 to 0", "0 to 60", "60 to 120"))
