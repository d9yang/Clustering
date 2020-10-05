library(readxl)

#number of samples
reads <- apply(OTUfile[,-1],1,sum)
length(reads)
summary(reads)
##calculate the proportion of 0
pr0 <- lapply(OTUfile[,-1], function(x) length(which(x==0))/length(x))

notu <- length(which(pr0<0.8))

#### merge data
PD1 <- OTUfile[, c(1,c(which(pr0<0.8))+1)]
PD2 <- merge(PD1, Covfile, by.x = "ID", by.y = "ID")

PD <- PD2[, c(2:length(PD2))]
colnames(PD)[1:length(PD1)] <- c(paste0("OTU",which(pr0<0.8)),"class")
PD$class <- as.numeric(as.factor(PD$class))

