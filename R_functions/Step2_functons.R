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
    dimsToUse = 2:30,
    force=T
)
return(Proj)
}


#### extrxt p2g from ArchR project ####
Get_p2g_fun <- function(x){
	corCutOff = 0.20
	FDRCutOff = 1e-6
	varCutOffATAC = 0.7
	varCutOffRNA = 0.3
	p2g <- metadata(x@peakSet)$Peak2GeneLinks
	p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
	if(!is.null(varCutOffATAC)){
    	p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
	}
	if(!is.null(varCutOffRNA)){
    	p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
	}
	mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
	mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
	p2g$peak <- paste0(rowRanges(mATAC))
	p2g$gene <- rowData(mRNA)$name
	return(p2g)
}

