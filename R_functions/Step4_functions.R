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


###
Output_to_bed_files = function(Fragment_list_cl,folder){
	tags = names(Fragment_list_cl)
	### to bam files ###
	fragments_cl_bamGR_list = list()
	for(i in 1:length(Fragment_list_cl)){
    	print(i)
    	temp_Fragment = Fragment_list_cl[[i]]
    	temp_Fragment$Reads = paste(tags[i],1:length(temp_Fragment),sep='_')
    	temp_Fragment_bamGR = fragments_to_bam_GR(temp_Fragment)
    	fragments_cl_bamGR_list = c(fragments_cl_bamGR_list,list(temp_Fragment_bamGR))
	}
	###
	for(i in 1:length(fragments_cl_bamGR_list)){
    	print(i)
    	library('readr')
    	library(tibble)
    	setwd(folder)
    	FN = paste(tags[i],"fragments_cl_bamGR_pe.bed",sep='_')
    	print(FN)
    	fragments_cl_bamGR_tab = Convert_to_bedpe(fragments_cl_bamGR_list[[i]])
    	print(head(fragments_cl_bamGR_tab))
    	write_tsv(fragments_cl_bamGR_tab, path=FN, col_names= FALSE)
	}
}


#####

Check_normalized_Signal <- function(file,savefile){
	library(rtracklayer)
	temp_bw = import.bw(file)
	print(summary(width(temp_bw)))
	print(summary(temp_bw$score))
	setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/All_normalized')
	saveRDS(temp_bw,file=savefile)
}


#####
Str_to_GR <- function(x){
	sp = strsplit(x,split='-')
	chr = sapply(sp,function(x) x[[1]])
	s = as.numeric(sapply(sp,function(x) x[[2]]))
	e = as.numeric(sapply(sp,function(x) x[[3]]))
	GR_out = GRanges(chr,IRanges(s,e))
	return(GR_out)
}

Seqnames <- function(x){
	out= as.character(seqnames(x))
	return(out)
}

Start <- function(x){
	out= as.numeric(start(x))
	return(out)
}

End <- function(x){
	out= as.numeric(end(x))
	return(out)
}

Must_to_GR <- function(x){
	library('pbapply')
	chr_all = pblapply(x,Seqnames)
	start_all = pblapply(x,Start)
	end_all = pblapply(x,End)
	len_all = pblapply(x,function(x) length(x))
	#####
	chr_all = as.character(unlist(chr_all))
	start_all = as.numeric(unlist(start_all))
	end_all = as.numeric(unlist(end_all))
	##### ##### 
	names_all = rep(names(x),len_all)
	GR_out = GRanges(chr_all,IRanges(start_all,end_all),motifs=names_all)
	return(GR_out)
}




Calculate_signal_bw <- function(GR,Signal){
	####
	res = findOverlaps(GR,Signal)
	res_score = Signal$score[subjectHits(res)]
	####
	res_score_merge = tapply(res_score,queryHits(res),sum)
	####
	GR$score = 0
	####
	GR$score[as.numeric(names(res_score_merge))] = as.numeric(res_score_merge)
	####
	return(GR$score)
}


Calculate_footprint <- function(footprint_GR,Signal){
	width = width(footprint_GR)
	#####
	left_s = start(footprint_GR)-width*3
	left_e = start(footprint_GR)-1
	#####
	left_GR = footprint_GR
	start(left_GR) = left_s
	end(left_GR) = left_e
	######
	right_s = end(footprint_GR)+1
	right_e = end(footprint_GR)+width*3
	######
	right_GR = footprint_GR
	start(right_GR) = right_s
	end(right_GR) = right_e
	######
	center_GR = footprint_GR
	######
	######
	print('Calculate_footprint_center')
	footprint_GR$center = Calculate_signal_bw(center_GR,Signal)
	print('Calculate_footprint_left')
	footprint_GR$left = Calculate_signal_bw(left_GR,Signal)
	print('Calculate_footprint_right')
	footprint_GR$right = Calculate_signal_bw(right_GR,Signal)
	######
	return(footprint_GR)
}


Calculate_footprint_celltypes <- function(footprint_GR,Signal,TFs,out_all_ext){
	motifs_need = out_all_ext$Motif[which(out_all_ext$TFs %in% TFs == T)]
	####
	print(length(motifs_need))
	####
	footprint_GR_cl = footprint_GR[which(footprint_GR$motifs %in% motifs_need == T)]
	####
	footprint_GR_cl_out = Calculate_footprint(footprint_GR_cl,Signal)
	####
	return(footprint_GR_cl_out)
}




Filter_footprints <- function(footprint_GR_cl_out,delta=0.1){
	####
	left_delta = (footprint_GR_cl_out$left)/3 - footprint_GR_cl_out$center
	right_delta = (footprint_GR_cl_out$right)/3 - footprint_GR_cl_out$center
	####
	k = which(footprint_GR_cl_out$left > delta*3 & footprint_GR_cl_out$right > delta*3 & footprint_GR_cl_out$center < -delta)
	####
	footprint_GR_cl_out_filtered = footprint_GR_cl_out[k]
	####
	return(footprint_GR_cl_out_filtered)
}



