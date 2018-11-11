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
    setwd("~/Documents/18-19/CSN/LABS/04/src")
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


computeTrianglePartitionRatio <- function(graph, communityMethod){
    
    numSubComm <- max(communityMethod$membership)
    verticesMembership <- communityMethod$membership
    verticesPos <- seq(length(V(graph)))
    res = c()
    for(subGIdx in seq(numSubComm)){
        # Take the vertices on that subcommunity
        vOfSubComm <- verticesPos[verticesMembership == subGIdx]
        # Create subgraph of subcommunity subGIdx
        subG = induced_subgraph(graph, vids = vOfSubComm)
        tri = sum(count_triangles(subG)) # Sum of all the vertices triangles
        metric = tri / numSubComm
        res = append(res, metric) 
    }
    
    # TODO: Check how to adapt metric for whole network
    res <- sum(res)
    return(res)
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


computeTableForGraph <- function(graph){

    communities = computeDiffCommunities(graph)
    resTable <- data.table("Method" = character(),
                           "TBT" = numeric(),
                           "Expansion" = numeric(),
                           "Conductance" = numeric(),
                           "Modularity" = numeric(),
        stringsAsFactors = FALSE)
    
    verticesPos <- seq(length(V(graph)))
    for(i in seq(length(communities))){
        
        c <- communities[i]; c <- c[[1]] # Dealing with the list messing up w/ the structure
        
        name <- communitiesNames[i]
        TPT <- computeTrianglePartitionRatio(graph, c)
        expansion <- computeExpansion(c)
        conduct <- computeConductance(c)
        mod <- computeModularity(c)
        
        resTable <- rbind(resTable, list(name, TPT, expansion, conduct, mod))
    }
}

graph = karate
computeTableForGraph(karate)