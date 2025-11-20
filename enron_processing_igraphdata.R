library(dplyr)
library(igraph)
library(igraphdata)
data(enron)

elist <- as.data.frame(as_edgelist(enron, names = TRUE))
colnames(elist) <- c("source", "target")

elist$time <- E(enron)$Time
elist$year <- as.integer(substring(elist$time, 1, 4))
elist$reciptype <- E(enron)$Reciptype

elist <- elist[(elist$year>=1998),]
elist <- elist[elist$source!=elist$target,]
elist <- elist[elist$reciptype=="to",]

elist <- elist %>% 
  group_by(source, target) %>%
  summarise(weight = n(), .groups = "drop")
write.csv(elist, "./dat/enron_elist.csv", row.names = F)

email <- V(enron)$Email
role <- V(enron)$Note
nlist <- data.frame(id=c(1:length(V(enron))), email=email, role=role)
write.csv(nlist, "./dat/enron_nlist.csv", row.names = F)

