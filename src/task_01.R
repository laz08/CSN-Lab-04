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


computeMetrics <- function(graph, communityMethod){
    m <- length(E(graph))
    n <- length(V(graph))
  
    numSubComm <- max(communityMethod$membership)
    verticesMembership <- communityMethod$membership
    verticesPos <- seq(n)
    degrees <- degree(graph)

    
    triangles = c()
    conductances <- c()
    modularities <- c()
    expansions <- c()
    for(subGIdx in seq(numSubComm)){
        # Take the vertices on that subcommunity
        vOfSubComm <- verticesPos[verticesMembership == subGIdx]
        # Create subgraph of subcommunity subGIdx
        subG = induced_subgraph(graph, vids = vOfSubComm)
        mc <- length(E(subG))
        nc <- length(V(subG))
        
        #TRIANGLES
        tri = sum(count_triangles(subG)) # Sum of all the vertices triangles
        tpt = tri / numSubComm
        weightedTPT = tpt * (nc/n)
        triangles = append(triangles, weightedTPT) 
        
        #EXPANSION
        subdegree <- degree(subG)
        subsetOrgDegrees <- degrees[verticesMembership == subGIdx]
        difdegree <- subsetOrgDegrees - subdegree
        fc <- sum(difdegree)
        expansion <- fc/numSubComm
        weightedExpansion <- expansion * (nc/n)
        expansions <- append(expansions, weightedExpansion)
        
        
        #MODULARITY
        p <- m / ((1/2) * n * (n - 1))
        expected <- p*((1/2) * nc * (nc - 1))
        modularity <- (1/(4*m))*(mc - expected)
        weightedModularity <- modularity * (nc/n)
        modularities <- append(modularities, weightedModularity)
        
        #CONDUCTANCE
        conductance <- fc/(fc + 2*(mc))
        weightedConductance <- conductance*(nc/n)
        conductances <- append(conductances, weightedConductance)
    }
    
    # TODO: Check how to adapt metric for whole network
    res <- c(sum(triangles), sum(modularities), sum(conductances), sum(expansions))
    return(res)
}


computeTableForGraph <- function(graph){

    communities = computeDiffCommunities(graph)
    resTable <- data.table("Method" = character(),
                           "TPT" = numeric(),
                           "Expansion" = numeric(),
                           "Conductance" = numeric(),
                           "Modularity" = numeric(),
        stringsAsFactors = FALSE)
    
    verticesPos <- seq(length(V(graph)))
    for(i in seq(length(communities))){
        
        c <- communities[i]; c <- c[[1]] # Dealing with the list messing up w/ the structure
        
        name <- communitiesNames[i]
        metrics <- computeMetrics(graph, c)
        TPT <- metrics[1]
        expansion <- metrics[4]
        conduct <- metrics[3]
        mod <- metrics[2]
        
        resTable <- rbind(resTable, list(name, TPT, expansion, conduct, mod))
    }
    return (resTable)
}

graph = karate
(t = computeTableForGraph(karate))