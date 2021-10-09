####### filter the cells in ArchR project########

Filter_proj = function(x,cellnames){
	k = which(x$cellNames %in% cellnames == T)
	cellsPass <- x$cellNames[k]
	x_cl = x[cellsPass]
	return(x_cl)
}


###### add peak matrix #######

addPeakSet_Matrix <- function(x,Peak_set_GR){
	x = addPeakSet(ArchRProj = x, peakSet = Peak_set_GR, force = T)
	x = addPeakMatrix(
  		ArchRProj = x,
  		ceiling = 4,
  		binarize = TRUE,
  		verbose = TRUE,
  		threads = 10,
  		parallelParam = NULL,
  		force = TRUE
	)
	return(x)
}

##### demension reduction analysis with PeakMatrix ####

Process_project <- function(Proj){
Proj <- addIterativeLSI(
    ArchRProj = Proj,
    useMatrix = "PeakMatrix", 
    name = "IterativeLSI", 
    iterations = 3, 
    clusterParams = list( #See Seurat::FindClusters)
        resolution = c(0.3), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 2:15,
    force=T
)
return(Proj)
}
