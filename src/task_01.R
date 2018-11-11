# Load and install necessary packages
requiredPackages <- c("igraph", "data.table")

for (pac in requiredPackages) {
    if(!require(pac,  character.only=TRUE)){
        install.packages(pac, repos="http://cran.rstudio.com")
        library(pac,  character.only=TRUE)
    } 
}
rm(pac)
rm(requiredPackages)



## Working directory
wd = getwd()
if(grepl("nora", wd)) {
    setwd("~/Documents/18-19/ASM/HW/04")
} else {
    ## Put working path for Carolina
}
rm(wd)

####################

karate <- graph.famous("Zachary")
wc <- walktrap.community(karate)
(modularity(wc))
(membership(wc))
plot(wc, karate)

computeDiffCommunities <- function(graph) {
    return (list(edge.betweenness.community(graph),
                 fastgreedy.community(graph),
                 label.propagation.community(graph),
                 leading.eigenvector.community(graph),
                 multilevel.community(graph),
                 spinglass.community(graph),
                 walktrap.community(graph),
                 infomap.community(graph)))
}

communitiesNames = c("edge.betweenness",
                    "fastgreedy",
                     "label",
                     "leading",
                     "multilevel",
                     "spinglass",
                     "walktrap",
                     "infomap")


computeTrianglePartitionRatio <- function(graph){
    return(0)
}

computeExpansion <- function(graph){
    return(0)
}

computeConductance <- function(graph){
    return(0)
}


computeModularity <- function(graph){
    return(0)
}


communities = computeDiffCommunities(karate)
resTable <- data.table("Method" = character(),
                       "TBT" = numeric(),
                       "Conductance" = numeric(),
                       "Modularity" = numeric(),
    stringsAsFactors = FALSE)

for(i in seq(length(communities))){
    
    c <- communities[i]
    
    name <- communitiesNames[i]
    TPT <- computeTrianglePartitionRatio(c)
    expansion <- computeExpansion(c)
    mod <- computeModularity(c)
    
    resTable <- rbind(resTable, list(name, TPT, expansion, mod))
}
