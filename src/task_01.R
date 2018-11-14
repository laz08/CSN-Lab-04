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
    setwd("~/Documents/18-19/CSN/LABS/04/data")
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

computeDiffCommunities <- function(graph, communitiesNames) {
    
    times = list()
    startTime = as.numeric(Sys.time())
    
    eb = edge.betweenness.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    fastgreedy = fastgreedy.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    labelProp = label.propagation.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    leadingEigen = leading.eigenvector.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    multilevel = multilevel.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    spinglass = spinglass.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    walktrap = walktrap.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    startTime = as.numeric(Sys.time())
    
    infomap = infomap.community(graph)
    times = append(times, as.numeric(Sys.time()) - startTime )
    
    df = data.table("Method" = communitiesNames,
                    "Elapsed time" = times,
                    stringsAsFactors = FALSE)
    
    communitiesList = list(eb, fastgreedy, labelProp, leadingEigen, multilevel, spinglass, walktrap, infomap)
    
    return(list(df, communitiesList))
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
    
    res <- c(sum(triangles), sum(modularities), sum(conductances), sum(expansions))
    return(res)
}


computeTableForGraph <- function(graph){

    results = computeDiffCommunities(graph, communitiesNames)
    resTable = results[1][[1]]
    #resTable <- data.table("Method" = character(),
    #                       "TPT" = numeric(),
    #                       "Expansion" = numeric(),
    #                       "Conductance" = numeric(),
    #                       "Modularity" = numeric(),
    #    stringsAsFactors = FALSE)
    
    
    metricsTPT = c()
    metricsExpansion = c()
    metricsConduct = c()
    metricsModularity = c()
    
    for(i in seq(length(communitiesNames))){
        
        #c <- communities[i]; c <- c[[1]] # Dealing with the list messing up w/ the structure
        c <- results[2][[1]][[i]] # Dealing with the list messing up w/ the structure
        # name <- communitiesNames[i]
        metrics <- computeMetrics(graph, c)
        metricsTPT = append(metricsTPT, metrics[1])
        metricsExpansion = append(metricsExpansion, metrics[4])
        metricsConduct = append(metricsConduct, metrics[3])
        metricsModularity = append(metricsModularity, metrics[2])
    }
    
    resTable = cbind(resTable, "TPT" = metricsTPT)
    resTable = cbind(resTable, "Expansion" = metricsExpansion)
    resTable = cbind(resTable, "Conduct" = metricsConduct)
    resTable = cbind(resTable, "Modularity" = metricsModularity)
    return (resTable)
}

graph = karate

results = computeDiffCommunities(graph, communitiesNames)

(t = computeTableForGraph(karate))

wikiG <- read.graph("wikipedia.gml", format="gml")

start = Sys.time()
t = walktrap.community(wikiG)
end = Sys.time()

delta = end - start
cat("Elapsed time: ", delta, "\n")
