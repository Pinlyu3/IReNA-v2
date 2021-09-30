### Read fragments ####
Read_fragment_to_GR <- function(x){
    tab <- read_tsv(x,col_names = F)
    GR <- GRanges(seqnames=tab$X1,ranges=IRanges(start=tab$X2,end=tab$X3),index=tab$X4,dup=tab$X5)
    return(GR)
}


### GR list to GR ####
unlist_fun <- function(GR_list){
	GR_Ori = GR_list[[1]]
	for(i in 2:length(GR_list)){
		GR_Ori = c(GR_Ori,GR_list[[i]])
	}
	return(GR_Ori)
}

