# Load and install necessary packages
requiredPackages <- c("igraph", "data.table", "knitr")

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
    setwd("~/Documents/FEUP/5A/1S/CSN/Lab/Lab04/CSN-Lab-04/data")
}
rm(wd)


####################
###### FLAGS #######
####################
PRINT_PLOTS = FALSE
TEST_ZACHARY = FALSE
COMPUTE_WIKI_COMMS = FALSE
###################
###################

################
### FUNCTIONS
### ############
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
    times = round(as.numeric(times), 5)
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
    expansions <- c()
    
    for(subGIdx in seq(numSubComm)){
        # Take the vertices on that subcommunity
        vOfSubComm <- verticesPos[verticesMembership == subGIdx]
        # Create subgraph of subcommunity subGIdx
        subG = induced_subgraph(graph, vids = vOfSubComm)
        mc <- length(E(subG))
        nc <- length(V(subG))
        
        #TRIANGLES
        alltriad = count_triangles(subG)
        tri = length(alltriad[alltriad > 0]) # count vertices that belong to triangles
        tpt = tri / nc
        weightedTPT = tpt * (nc/n)
        triangles = append(triangles, weightedTPT) 
        
        #EXPANSION
        # Number of edges outside the graph
        subdegree <- degree(subG)
        subsetOrgDegrees <- degrees[verticesMembership == subGIdx]
        difdegree <- subsetOrgDegrees - subdegree
        fc <- sum(difdegree)
        expansion <- fc/nc
        weightedExpansion <- expansion * (nc/n)
        expansions <- append(expansions, weightedExpansion)
      
        
        #CONDUCTANCE
        conductance <- fc/(fc + 2*(mc))
        weightedConductance <- conductance*(nc/n)
        conductances <- append(conductances, weightedConductance)
    }
    
    res <- c(sum(triangles), modularity(communityMethod), sum(conductances), sum(expansions), numSubComm)
    return(res)
}


computeTableForGraph <- function(graph){

    results = computeDiffCommunities(graph, communitiesNames)
    resTable = results[1][[1]] # Dealing with the list messing up w/ the structure
    
    metricsTPT = c()
    metricsExpansion = c()
    metricsConduct = c()
    metricsModularity = c()
    nrSubCommsFound = c()
    
    for(i in seq(length(communitiesNames))){
        
        c <- results[2][[1]][[i]] # Dealing with the list messing up w/ the structure
        
        metrics <- computeMetrics(graph, c)
        metricsTPT = append(metricsTPT, metrics[1])
        metricsExpansion = append(metricsExpansion, metrics[4])
        metricsConduct = append(metricsConduct, metrics[3])
        metricsModularity = append(metricsModularity, metrics[2])
        nrSubCommsFound = append(nrSubCommsFound, metrics[5])
    }
    
    resTable = cbind(resTable, "# C" = nrSubCommsFound)
    resTable = cbind(resTable, "TPT" = round(metricsTPT, 4))
    resTable = cbind(resTable, "Expansion" = round(metricsExpansion, 4))
    resTable = cbind(resTable, "Conduct" = round(metricsConduct, 4))
    resTable = cbind(resTable, "Modularity" = round(metricsModularity, 4))
    return (resTable)
}


computeSummaryTable <- function(graphs){
    
    graphsNames = c("Zachary", "Tutte", "Coxeter", "Modified Thomassen", "3 full conn. graphs")
    table <- data.table("Graph" = character(),
                        "N" = numeric(),
                        "E" = numeric(),
                        "k" = numeric(),
                        "delta" = numeric(),
                        stringsAsFactors = FALSE)
    
    for (x in 1:length(graphsNames)){
        
        g = graphs[[x]]
        gName = graphsNames[x]
        
        E = length(E(g))
        N = length(V(g))
        k = 2*E/N
        delta = 2*E/(N * (N-1))
        
        table <- rbind(table, list(gName, N, E, round(k, 2), round(delta, 2)))
    }
    return(table)
}

computeSummaryTableForWiki <- function(wikiG){
    

    table <- data.table("Graph" = character(),
                        "N" = numeric(),
                        "E" = numeric(),
                        "k" = numeric(),
                        "delta" = numeric(),
                        stringsAsFactors = FALSE)
    
    g = wikiG
    gName = "Wikipedia"
    
    E = length(E(g))
    N = length(V(g))
    k = 2*E/N
    delta = 2*E/(N * (N-1))
    
    table <- rbind(table, list(gName, N, E, round(k, 2), round(delta, 6)))
        
    return(table)
}

plotGraphBeautifully <- function(graph){
    plot.igraph(graph,layout=layout.auto,vertex.size=23,vertex.label.color="yellow",vertex.shape = "sphere", vertex.label.font=2,vertex.color="darkblue",edge.color="black")
}
  


##############################
##############################

# Graphs to study

karate <- graph.famous("Zachary")
tutte <- graph.famous("Tutte")
coxeter <- graph.famous("Coxeter")

Thomassen <- graph.famous("Thomassen")
Thomassen <- Thomassen + make_full_graph(5) + make_full_graph(10)
Thomassen <- add_edges(Thomassen, c(34, 35, 37, 45))

full3Graphs <- make_full_graph(10) + make_full_graph(10) +  make_full_graph(10) 
full3Graphs <- add_edges(full3Graphs, c(10, 11, 20, 21))


listOfG <- list(karate, tutte, coxeter, Thomassen, full3Graphs)


# Making the metrics table
summaryTable <- computeSummaryTable(listOfG)


####################

karateMetricsTable = computeTableForGraph(karate)
tutteMetricsTable = computeTableForGraph(tutte)
coxeterMetricsTable = computeTableForGraph(coxeter)
modThomassenMetricsTable = computeTableForGraph(Thomassen)
full3MetricsTable = computeTableForGraph(full3Graphs)

#################################
#################################

if(TEST_ZACHARY){
    wc <- walktrap.community(karate)
    
    (modularity(wc))
    (membership(wc))
    plot(wc, karate)
}


#### Task 02
computeWalktrap <- function(graph){
    start = Sys.time()
    walktrapComms = walktrap.community(graph)
    end = Sys.time()
    
    delta = end - start
    cat("Elapsed time: ", delta, "\n")
    return(walktrapComms)
}


printNumberOfCommunitiesFound <- function(communities){
    numSubComm <- max(communities$membership)
    cat("Total number of communities found: ", numSubComm, "\n")
}

plotGraphFirstXCommunities <- function(communities, graph, x) {
  
    #par(mfrow=c(1,2))
    
  verticesMembership <- communities$membership
  verticesPos <- seq(length(communities$membership))
  
  for(subGIdx in seq(x)){
      # Take the vertices on that subcommunity
      vOfSubComm <- verticesPos[verticesMembership == subGIdx]
      # Create subgraph of subcommunity subGIdx
      subG = induced_subgraph(graph, vids = vOfSubComm)
      plot(subG, main=paste("Plot of community", subGIdx))
      box(which="plot")
  }

}



wikiG <- read.graph("wikipedia.gml", format="gml")

wikiSummaryTable <- computeSummaryTableForWiki(wikiG)
if(COMPUTE_WIKI_COMMS){
    walktrapCommsWiki = computeWalktrap(wikiG)    
}


if(PRINT_PLOTS){
    plotGraphFirstXCommunities(walktrapCommsWiki, wikiG, 4)
}
