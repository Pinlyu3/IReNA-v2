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


### 
fragments_to_bam_GR <- function(fragments_cl){
    chr = seqnames(fragments_cl)
    s = start(fragments_cl)
    e = end(fragments_cl)
    GR_plus = GRanges(seqnames=chr,IRanges(s,s+1),strand='+',Reads=fragments_cl$Reads)
    GR_minus = GRanges(seqnames=chr,IRanges(e-1,e),strand='-',Reads=fragments_cl$Reads)
    GR <- list(GR_plus,GR_minus)
    return(GR)
}


###
Convert_to_bedpe <- function(gr_list){
	library(tibble)
	######
	gr_left = gr_list[[1]]
	gr_right = gr_list[[2]]
	######
	chr1 = as.character(seqnames(gr_left))
	s1 = start(gr_left)
	e1 = end(gr_left)
	######
	chr2 = as.character(seqnames(gr_right))
	s2 = start(gr_right)
	e2 = end(gr_right)
	######
	name = gr_left$Reads
	scores=c(rep(".", length(gr_left)))
	######
	strands1=as.character(strand(gr_left))
	strands2=as.character(strand(gr_right))
	######
    df <- tibble(chr1,s1,e1,chr2,s2,e2,name,scores,strands1,strands2)
    return(df)
}

