#### combine footprint and cis-regulatory elements #####

Reg_one_cells_RPC_MG <- function(footprint_res,peak_gene_list,out_all_ext,diff_list){
	###### process the list ##########
	###### generate peak-gene tables  ######
	peak_gene_list$Gene1_TSS$Correlation = 'TSS'
	peaks = c()
	genes = c()
	PtoG = c()
	for(i in 1:length(peak_gene_list)){
		peaks = c(peaks,peak_gene_list[[i]]$peaks)
		print(length(peak_gene_list[[i]]$peaks))
		genes = c(genes,peak_gene_list[[i]]$gene_name)
		print(length(peak_gene_list[[i]]$gene_name))
		PtoG = c(PtoG,peak_gene_list[[i]]$Correlation)
	}
	#### ###### ###### ###### ###### #####
	### k = which(genes %in% rownames(diff_list) == T)
	#print(length(k)/length(genes))
	#print(length(genes[!duplicated(genes)]))
	####
	Out_Dat = data.frame(Peaks =peaks ,Target=genes,PtoG=PtoG)
	dim(Out_Dat)
	####
	#### Next filter the motifs according to gene expression ############
	print(length(diff_list))
	out_all_ext_cl = out_all_ext[which(out_all_ext$TFs %in% diff_list == T),]
	Motif_need = out_all_ext_cl$Motif
	#### filter the footprint ###################
	k = which(footprint_res$motifs %in% Motif_need == T)
	footprint_res_cl = footprint_res[k]
	#### merge Out_Dat with footprint_res_cl ######
	Out_Dat_GR = GRanges(Out_Dat$Peaks,Target=Out_Dat$Target,PtoG=PtoG)
	#### remove their not overlapped regions ####################
	k1 = which(countOverlaps(Out_Dat_GR,footprint_res_cl) >0)
	k2 = which(countOverlaps(footprint_res_cl,Out_Dat_GR) >0)
	#### keep overlapped ######
	Out_Dat_GR_Overlap = Out_Dat_GR[k1]
	footprint_res_cl_Overlap = footprint_res_cl[k2]
	#### merge Out_Dat_GR_Overlap and footprint_res_cl_Overlap ####
	Res_find = data.frame(findOverlaps(footprint_res_cl_Overlap,Out_Dat_GR_Overlap))
	####
	Res_find$Motifs = footprint_res_cl_Overlap$motifs[Res_find$queryHits]
	Res_find$footprint = as.character(footprint_res_cl_Overlap)[Res_find$queryHits]
	####
	Res_find$Target = Out_Dat_GR_Overlap$Target[Res_find$subjectHits]
	Res_find$peaks = as.character(Out_Dat_GR_Overlap)[Res_find$subjectHits]
	Res_find$PtoG = Out_Dat_GR_Overlap$PtoG[Res_find$subjectHits]
	##### add TF names to it #######
	##### one motif may corresponding to multiple TF genes #####
	colnames(out_all_ext) = c('Motifs','TFs')
	out_all_ext_cl = out_all_ext[which(out_all_ext$TFs %in% diff_list == T),]
	#####
	Res_find_merge = merge(Res_find,out_all_ext_cl)
	####
	Res_find_merge = Res_find_merge[,c(8,1,4,6,7,5)]
	####
	return(Res_find_merge)
}


#### 
#### filter the GRNs with TF-target gene expression corr ####
####

Add_Cor_to_GRN_network_and_Filter <- function(Cor_Res,Reg,diff_list){
	###### filter Cor_Res gene pairs #####
	k = which((Cor_Res$Var1 %in% diff_list ==T) & (Cor_Res$Var2 %in% diff_list ==T))
	Cor_Res = Cor_Res[k,]
	####### prepare the gene pair index in correlation of gene pairs ####
	Cor_Res_cl_1 = Cor_Res
	Cor_Res_cl_2 = Cor_Res
	Cor_Res_cl_1$index = paste(Cor_Res$Var1,Cor_Res$Var2)
	Cor_Res_cl_2$index = paste(Cor_Res$Var2,Cor_Res$Var1)
	Cor_Res_cl_combine = rbind(Cor_Res_cl_1,Cor_Res_cl_2)
	####### prepare the gene pair index in GRN networks ######
	Reg$index = paste(Reg$TFs,Reg$Target)
	####### added correlation and regulatory types #####
	m = match(Reg$index,Cor_Res_cl_combine$index)
	Reg$Cor = Cor_Res_cl_combine$value[m]
	Reg$Tag = Cor_Res_cl_combine$tag[m]
	###### filtered out the regulatory gene pairs if their correlation are not signifiant #########
	Reg_cl = Reg[which(Reg$Tag != 'No'),]
	#####
	return(Reg_cl)
}


