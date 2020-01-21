require(quanteda)
require(data.table)
require(stringr)
require(igraph)
require(plotly)
require(RColorBrewer)
require(ggplot2)
require(dplyr)
require(ggrepel)
require(ggraph)
require(leiden)

############################################
# DATA
############################################

df <- read.csv("patents_first_clean.csv", stringsAsFactors = F, strip.white = T) #194322 rows
meta <- fread("patents_mergedclaims.csv")
colnames(meta) <- c("row", as.character(meta[1,2:14]))
meta <- meta[2:nrow(meta),]
meta$patent_id <- meta$id
df <- merge(df, meta, by = "patent_id", all.x = T)
df <- df[!is.na(df$`filing/creation date`),] # 194284 rows
df$year <- str_extract(df$`filing/creation date`, "\\d\\d$")
table(df$year)
df$year <- as.numeric(df$year)
df$year <- ifelse(df$year < 50, df$year+2000, df$year+1900)


############################################
# PROCESSING
############################################

claim.dtm <- dfm(corpus(df$claim, docvars = df)) # 194,284 documents, 16,751 features
bin.claim.dtm <- dfm_weight(claim.dtm, scheme = "boolean") # binarize counts
bin.claim.dtm <- dfm_remove(bin.claim.dtm, stopwords(), valuetype = "fixed") #194,284 documents, 16,636 features
bin.claim.dtm <- dfm_keep(bin.claim.dtm, pattern = ".{3,}", valuetype = "regex") #194,284 documents, 16,520 features

bin.claim.dtm <- dfm_trim(bin.claim.dtm, min_docfreq = 2) # 194,284 documents, 13,392
bin.claim.dtm <- dfm_subset(bin.claim.dtm, rowSums(bin.claim.dtm) > 1) # 190,767 documents, 13,392 features
bin.claim.dtm <- dfm_trim(bin.claim.dtm, min_docfreq = 2) # 190,767 documents, 13,387 features
bin.claim.dtm <- dfm_subset(bin.claim.dtm, rowSums(bin.claim.dtm) > 1)  # 190,767 documents, 13,387 features

table(docvars(bin.claim.dtm,"year")) #there are year gaps before 1968
docvars(bin.claim.dtm,"year") <- ifelse(docvars(bin.claim.dtm,"year") <= 1968,
                                        1968, docvars(bin.claim.dtm,"year"))
year.claim.dtm <- dfm_group(bin.claim.dtm, "year")
year.claim.term.freqs <- textstat_frequency(year.claim.dtm )
table(year.claim.term.freqs$docfreq) # there are 1996 words that appear in a single year
bin.claim.dtm <- dfm_remove(bin.claim.dtm,
                            year.claim.term.freqs[year.claim.term.freqs$docfreq == 1,]$feature,
                            valuetype = "fixed") # 190,767 documents, 11,391 features

bin.claim.dtm <- dfm_trim(bin.claim.dtm, min_docfreq = 2) # 190,767 documents, 11,391 features
bin.claim.dtm <- dfm_subset(bin.claim.dtm, rowSums(bin.claim.dtm) > 1) # 190,735 documents, 11,391 features
bin.claim.dtm <- dfm_trim(bin.claim.dtm, min_docfreq = 2) # 190,735 documents, 11,391 features
bin.claim.dtm <- dfm_subset(bin.claim.dtm, rowSums(bin.claim.dtm) > 1)  # 190,735 documents, 11,391 features

year.claim.dtm <- dfm_group(bin.claim.dtm, "year")
year.claim.term.freqs <- textstat_frequency(year.claim.dtm )
table(year.claim.term.freqs$docfreq) # there are no single year words

# All words appear in more than one year, and in more than one document.
# All documents have more than one word

claim.term.freqs <- textstat_frequency(bin.claim.dtm )

save(bin.claim.dtm, file = "~/Dropbox/patents/data/dtm_nov.rda")

############################################
# NETWORKS
############################################
# TODO:
# MOVE FROM CO-OCURR TO CLARK_DE

year.network <- function(dtm, yr){
  temp.dtm <- dfm_subset(dtm, year == yr)
  temp.dtm <- dfm_trim(temp.dtm, min_docfreq = 2)
  temp.dtm <- dfm_subset(temp.dtm, rowSums(temp.dtm) > 1)
  temp.fcm <- fcm(temp.dtm, tri = F)
  ####################################################
  # PMI implementation, see supplementary materials in "Lexical Shifts"
  # look for 'weighted pmi'
  tcmrs = Matrix::rowSums(temp.fcm)
  tcmcs = Matrix::colSums(temp.fcm)
  N = sum(tcmrs)
  colp = tcmcs/N
  rowp = tcmrs/N
  pp = temp.fcm@p+1
  ip = temp.fcm@i+1
  tmpx = rep(0,length(temp.fcm@x)) # new values go here, just a numeric vector
  # iterate through sparse matrix:
  for(i in 1:(length(temp.fcm@p)-1) ){ 
    ind = pp[i]:(pp[i+1]-1)
    not0 = ip[ind]
    icol = temp.fcm@x[ind]
    tmp = log( (icol/N) / (rowp[not0] * colp[i] )) # PMI
    tmpx[ind] = tmp    
  }
  temp.fcm@x = tmpx
  james <- function(x) (abs(x)+x)/2
  nafn <- function(x) ifelse(is.na(x),0,x)
  temp.fcm <- james(temp.fcm)
  temp.fcm <- Matrix::drop0(temp.fcm)
  ############################################################################
  g <- graph_from_adjacency_matrix(adjmatrix = temp.fcm,
                                   weighted = T,
                                   diag = F, mode = "undirected")
  g <- delete.vertices(simplify(g), degree(g)==0)
  cl <- clusters(g)
  g <- induced_subgraph(g, which(cl$membership == which.max(cl$csize)))
  
  partition <- leiden(as_adjacency_matrix(g),
                        seed = 123,
                        n_iterations = -2)
  V(g)$leiden.comm <- partition
  c <- cluster_louvain(g)
  V(g)$louvain.comm <- c$membership
  
  return(g)
}




############################################
# NODE CONTRIBUTIONS
############################################
# how much does each node contribute to its corresponding community?

#TODO:
# IMPLEMENT CLARK_DE?
# PARALLELIZE contrib.to.comm

contrib.to.comm <- function(g, comms = "leiden"){
  # Rule et al. p.2, "Network Mapping"
  # how important is a node to a community
  if(comms == "leiden"){
    comm.kind = "leiden.comm"
  }else{
    comm.kind = "louvain.comm"
  }
  comm.df <- data.frame(name= V(g)$name,
                        contrib.to.comm = NA,
                        comm = eval(expr = parse(text = (paste("V(g)$",comm.kind, sep = "")))),
                        stringsAsFactors = F)
  for(node in 1:nrow(comm.df)){
    node.comm <- comm.df[comm.df$comm == comm.df[node, "comm"],]$name
    temp.subgraph <- induced_subgraph(g, node.comm)
    node.temp.index <- which(V(temp.subgraph)$name == comm.df[node,"name"])  #temp.subgraph has new indices
    node.strength <- as.numeric(strength(temp.subgraph, node.temp.index))
    total <- as.numeric(sum(strength(temp.subgraph)))
    comm.df[node,"contrib.to.comm"] <- (node.strength/total)
  }
  return(comm.df)
}

contrib.to.node <- function(g, comms = "leiden"){
  # Rule et al. p.2, "Network Mapping"
  # how important is a community to a node
  if(comms == "leiden"){
    comm.kind = "leiden.comm"
  }else{
    comm.kind = "louvain.comm"
  }
  node.contrib <- data.frame(name= V(g)$name,
                             contrib.to.node = NA,
                             comm = eval(expr = parse(text = (paste("V(g)$",comm.kind, sep = "")))),
                             stringsAsFactors = F)
  for(node in 1:nrow(node.contrib)){
    node.comm <- node.contrib[node.contrib$comm == node.contrib[node, "comm"],]$name
    node.neigh <- neighbors(g, node)
    
    numerator.nodes <- intersect(node.comm, node.neigh$name)
    numerator.vp <- as.vector(t(as.matrix(expand.grid(c(node.contrib[node, "name"]), numerator.nodes))))
    numerator.ei <- get.edge.ids(g,numerator.vp)
    
    denominator.vp <- as.vector(t(as.matrix(expand.grid(c(node.contrib[node, "name"]), node.neigh$name))))
    denominator.ei <- get.edge.ids(g,denominator.vp)
    
    contrib.to.comm <- sum(E(g)[numerator.ei]$weight)/sum(E(g)[denominator.ei]$weight)
    node.contrib[node,"contrib.to.node"] <- contrib.to.comm 
  }
  return(node.contrib)
}



############################################
#  COMMUNITY ALIGNMENTS
############################################

r.dist <- function(v.t2,v.t1){
  # Rule et al. p.2, "River Networks"
  return(1- (  ((sqrt(2))^-1) * sqrt( sum( ( (sqrt(v.t2) - sqrt(v.t1))^2 )) )))
}

align <- function(path.t1, path.t2){
  comm.df.1 <- read.csv(path.t1)
  comm.df.2 <- read.csv(path.t2)
  colnames(comm.df.1) <- paste(colnames(comm.df.1), "1", sep=".")
  colnames(comm.df.2) <- paste(colnames(comm.df.2), "2", sep=".")
  
  clust.df <- merge(comm.df.1 , comm.df.2, by.x = "name.1", by.y= "name.2", all=T)
  clust.df[is.na(clust.df$contrib.to.comm.1),]$contrib.to.comm.1 <- 0
  clust.df[is.na(clust.df$contrib.to.comm.2),]$contrib.to.comm.2 <- 0
  
  out.df <- expand.grid(unique(clust.df$comm.1),
                        unique(clust.df$comm.2))
  out.df <- out.df[complete.cases(out.df),]
  out.df$alignment <- NA
  for(comm.1 in unique(out.df$Var1)){
    for(comm.2 in unique(out.df$Var2)){
      temp <- clust.df
      temp[is.na(temp$comm.1),]$comm.1 <- 0 
      temp[is.na(temp$comm.2),]$comm.2 <- 0 
      temp$contrib.to.comm.1 <- ifelse(temp$comm.1 != comm.1, 0, temp$contrib.to.comm.1)
      temp$contrib.to.comm.2 <- ifelse(temp$comm.2 != comm.2, 0, temp$contrib.to.comm.2)
      condition <- ((out.df$Var1 == comm.1) & (out.df$Var2 == comm.2))
      out.df[condition,]$alignment <- r.dist(temp$contrib.to.comm.2, 
                                             temp$contrib.to.comm.1)
    }
  }
  colnames(out.df) <- c("comm.t1", "comm.t2", "alignment")
  out.df$t1 <- str_extract(path.t1, "\\d{4}")
  out.df$t2 <- str_extract(path.t2, "\\d{4}")
  return(out.df)
}


############################################
#  INIT
############################################
#TODO: PARALLELIZE 

#load("~/Dropbox/patents/data/dtm_nov.rda")

years <- sort(unique(docvars(bin.claim.dtm,"year")))

for(year in years){
  print(year)
  g <- year.network(bin.claim.dtm, year)
  save(g, file = paste("~/Dropbox/patents/data/results/graphs/",
                       "graph_",
                       as.character(year),
                       ".rda",
                       sep = ""))
}


for(year in years){
  load(paste("~/Dropbox/patents/data/results/graphs/",
             "graph_",
             as.character(year),
             ".rda",
             sep = ""))
  comm.df <- contrib.to.comm(g)
  write.csv(comm.df, paste("~/Dropbox/patents/data/results/leiden_commdfs/",
                           "leiden_commdf_",
                           as.character(year),
                           ".csv", sep = ""))
  }

for(year in years[1:(length(years)-1)]){
  print(year)
  path1 <- paste("~/Dropbox/patents/data/results/leiden_commdfs/",
                 "leiden_commdf_",
                 as.character(year),
                 ".csv", sep = "")
  
  path2 <- paste("~/Dropbox/patents/data/results/leiden_commdfs/",
                 "leiden_commdf_", 
                 as.character(year+1),
                 ".csv", sep = "")
  temp <- align(path1, path2)
  write.csv(temp, paste("~/Dropbox/patents/data/results/leiden_aligned_comms/","leiden_aligndf",
                        as.character(year),"_",
                        as.character(year+1),
                        ".csv", sep = ""))
}


############################################
#  VISUALIZE
############################################

#TODO:
# ADD YEAR INTERVAL FUNCTIONALITY
# ADD LABELLING FUNCTIONALITY

#TODO:
# Study threshold sensitivity, look at supplemental materials in "Lexical Shifts"

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=T)
  do.call(rbind, lapply( filenames, function(nam){ 
    cbind(name=nam, read.csv(paste0(nam), header=TRUE))}))}

plot.rivernet <- function(path, title = "", thresh = 0, labs = NULL){
  #TODO: add community labels functionality
  edge.list <- multmerge(path)
  labels <- sort(unique(c(paste(edge.list$t2,edge.list$comm.t2, sep = "|"),
                          paste(edge.list$t1,edge.list$comm.t1, sep = "|"))))
  timepoints <- table(str_extract(labels, "^\\d{4}"))
  
  colors <- mapply(rep, brewer.pal(length(timepoints), "Spectral"), as.vector(timepoints))
  colors <- Reduce(c, colors)
  edge.list$source <- paste(edge.list$t1,edge.list$comm.t1, sep = "|")
  edge.list$target <- paste(edge.list$t2,edge.list$comm.t2, sep = "|")
  
  temp.edge.list <- edge.list[edge.list$alignment >= thresh,]
  
  plot_ly(
    type = "sankey",
    orientation = "h",
    hoverinfosrc = "source",
    node = list(
      label = labels,
      color = colors,
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    
    link = list(
      source = match(temp.edge.list$source,labels) - 1,
      target = match(temp.edge.list$target,labels) - 1,
      value =  temp.edge.list$alignment
    )
  ) %>% 
    layout(
      title = title,
      font = list(
        size = 10
      )
    )
}

plot.rivernet("~/Dropbox/patents/data/results/aligned_comms/toy", thresh = .1)


############################################
#  ANALYSIS
############################################

# Revise paths

path <- "~/Dropbox/patents/data/results/leiden_aligned_comms/"
edge.list <- multmerge(path)

edge.list$source <- paste(edge.list$t1, edge.list$comm.t1, sep = "|")
edge.list$target <- paste(edge.list$t2, edge.list$comm.t2, sep = "|")
edge.list$weight <- edge.list$alignment
g <- graph_from_data_frame(edge.list[,c("source", "target", "weight")])

psi.plus <- c()
for(node in V(g)$name){
  out.edges <- incident(g, node, "out")
  out.strength <- as.numeric(strength(g, node, "out"))
  if(length(out.edges) < 1){
    # psi is not defined for new communities
    psi.plus <- c(psi.plus,NA)
  }else{
    l.two <- sum(((out.edges$weight)/(out.strength))^2)^.5 
    psi.plus <- c(psi.plus,l.two)
  }
}

psi.minus <- c()
for(node in V(g)$name){
  in.edges <- incident(g, node, "in")
  in.strength <- as.numeric(strength(g, node, "in"))
  if(length(in.edges) < 1){
    # psi is not defined for new communities
    psi.minus <- c(psi.minus,NA)
  }else{
    l.two <- sum(((in.edges$weight)/(in.strength))^2)^.5 
    psi.minus <- c(psi.minus,l.two)
  }
}

result.df <- data.frame(comm = V(g)$name,
                        psi.plus = psi.plus,
                        psi.minus = psi.minus)
result.df$year <- as.numeric(str_extract(result.df$comm, "\\d{4}"))



p0 <- ggplot(result.df, aes(psi.minus, psi.plus)) + 
  geom_hline(yintercept = .5, lwd = 1) +
  geom_vline(xintercept = .5, lwd = 1) +
  scale_y_continuous(breaks = c(0,1)) +
  scale_x_continuous(breaks = c(0,1)) +
  labs(x = bquote(psi[c]^"-"),
       y = bquote(psi[c]^"+"),
       title = "Relation between Outgoing and Received Coherence") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))
  

p1 <- ggplot(result.df, aes(psi.minus, psi.plus)) +
  geom_point(aes(colour = year)) +
  labs(x = bquote(psi[c]^"-"),
       y = bquote(psi[c]^"+"),
       title = "Most Communities are Reorganizers") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))


p2 <- ggplot(result.df, aes(psi.minus, psi.plus)) +
  geom_point(color = "grey", alpha = .7) +
  labs(x = bquote(psi[c]^"-"),
       y = bquote(psi[c]^"+"),
       title = "Average Outgoing Coherence Falls Below the Diagonal") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0))) +
  geom_abline(intercept = 0, slope = 1, lwd = 1) +
  geom_smooth(method = "lm", lwd = 1, color =  "red", se = F)


centroids <- aggregate(cbind(psi.minus, psi.plus)~year,
                       result.df,mean)
p2.1 <- ggplot(centroids, aes(psi.minus, psi.plus)) +
  geom_point(aes(colour = year)) +
  labs(x = bquote(psi[c]^"-"),
       y = bquote(psi[c]^"+"),
       title = "Year Centroids Cluster in Reorganizer Quadrant") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0))) +
  xlim(c(.4,1)) + ylim(c(.4,1))


year.df <- result.df[complete.cases(result.df),] %>%
  group_by(year) %>%
  summarise(psi.dot = 1-((1/n())*(psi.plus%*%psi.minus)))

p3 <- ggplot(year.df, aes(year, psi.dot)) +
  labs(x = "t",
       y = bquote(psi[t]),
       title = "Discontinuities Decrease Over Time") +
  geom_line(lwd = 1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))

ggsave(filename = "~/Dropbox/patents/figs/Fig0_TextXD.png",
       plot = p0, dpi = 400)
ggsave(filename = "~/Dropbox/patents/figs/Fig1_TextXD.png",
       plot = p1, dpi = 400)
ggsave(filename = "~/Dropbox/patents/figs/Fig2_TextXD.png",
       plot = p2, dpi = 400)
ggsave(filename = "~/Dropbox/patents/figs/Fig2_1_TextXD.png",
       plot = p2.1, dpi = 400)
ggsave(filename = "~/Dropbox/patents/figs/Fig3_TextXD.png",
       plot = p3, dpi = 400)



plot(g, 
     vertex.size = 2, 
     vertex.label = NA,
     edge.arrow.width = .5,
     edge.arrow.size = .1,
     edge.arrow.mode = "-")


V(g)$Year <- as.numeric(str_extract(V(g)$name, "^\\d{4}"))

p4 <- ggraph(g) + 
  geom_edge_link(color = "gray", alpha = .5) + 
  geom_node_point(aes(colour = Year)) + 
  labs(title = "Community Alignment Graph") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))

ggsave(filename = "~/Dropbox/patents/figs/Fig4_TextXD.png",
       plot = p4, dpi = 400)


freq.df <- data.frame(table(docvars(bin.claim.dtm, "year")))
colnames()

p5 <- docvars(bin.claim.dtm, c("patent_id", "year")) %>%
  group_by(year) %>%
  summarize(patent.count = n()) %>%
  ggplot(aes(x = year, y = patent.count)) +
  geom_line()  + 
  labs(title = "Number of Patents over Time",
       x = "Year",
       y= "Num. Patents") +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))

ggsave(filename = "~/Dropbox/patents/figs/Fig5_TextXD.png",
       plot = p5, dpi = 400)


load("~/Dropbox/patents/data/results/graphs/graph_2010.rda")
comm.df <- read.csv("~/Dropbox/patents/data/results/leiden_commdfs/leiden_commdf_2010.csv")

table(V(g)$leiden.comm)
subg <- induced_subgraph(g, which(V(g)$leiden.comm == 4))
comm.4 <- comm.df[comm.df$comm == 4,]
comm.4 <- comm.4[order(comm.4$contrib.to.comm, decreasing = T),]

set.seed(101)
p6 <- ggraph(subg) + 
  geom_edge_link(color = "gray", alpha = .8) +
  geom_node_point(color = "darkgray") +
  labs(title = "Math Toys in 2010") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))
set.seed(101)
p6.1 <- ggraph(subg) + 
  geom_edge_link(color = "gray", alpha = .8) +
  geom_node_point(color = "darkgray") + 
  geom_node_label(aes(label =  name), size = 3, repel = T) +
  labs(title = "Math Toys in 2010") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, margin=margin(0,0,15,0)))

ggsave(filename = "~/Dropbox/patents/figs/Fig6_TextXD.png",
       plot = p6, dpi = 400)
ggsave(filename = "~/Dropbox/patents/figs/Fig6_1_TextXD.png",
       plot = p6.1, dpi = 400)

############################################
#  CONTRIB TO NODE
############################################
#TODO


############################################
#  CUMULATIVE NETWORK
############################################
#TODO

load("~/Dropbox/patents/data/dtm_nov.rda")
years <- sort(unique(docvars(bin.claim.dtm,"year")))

for(year in years){
  print(year)
  load(paste("~/Dropbox/patents/data/results/graphs/",
             "graph_",
             as.character(year),
             ".rda",
             sep = ""))
  comm.df <- contrib.to.node(g)
  write.csv(comm.df, paste("~/Dropbox/patents/data/results/leiden_nodedfs/",
                           "leiden_nodedf_",
                           as.character(year),
                           ".csv", sep = ""))
  }

