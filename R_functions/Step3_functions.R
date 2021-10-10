##### get potential cis-regulatory elements for each candidate gene #####
#####

Selection_peaks_for_one <- function(All_peaks_list,All_genes,p2g,distance_F=100000,mm10_TSS_GR_all){
		All_genes = All_genes
		#### first to find the TSS Peaks ##########
		Genes1 = All_genes
		TSS_peak_list = All_peaks_list$TSS
		TSS_peak_list$peaks = as.character(TSS_peak_list)
		### filter ###
		Genes1_TSS1 = TSS_peak_list[which(TSS_peak_list$gene_name %in% Genes1 == T)]
		m = match(Genes1_TSS1$gene_name,Genes1)
		Genes1_TSS1 = Genes1_TSS1[order(m)]
		print(Genes1_TSS1$gene_name[!duplicated(Genes1_TSS1$gene_name)])
		############# Second ######
		##### Gene Body list ##############
		Body_peak_list = All_peaks_list$GeneBody
		Body_peak_list$peaks = as.character(Body_peak_list)####
		### need p2g correlation ###############
		p2g$G_P = paste(p2g$gene,p2g$peak,sep='@')
		Body_peak_list$G_P = paste(Body_peak_list$gene_name,Body_peak_list$peaks,sep='@')
		#### add P2G corr #####
		k = which(Body_peak_list$G_P %in% p2g$G_P == T)
		Body_peak_list_cl = Body_peak_list[k]
		m = match(Body_peak_list_cl$G_P,p2g$G_P)
		Body_peak_list_cl$Correlation = p2g$Correlation[m]
		#####
		##### filter #####
		Genes1_Body1 = Body_peak_list_cl[which(Body_peak_list_cl$gene_name %in% Genes1 == T)]
		m = match(Genes1_Body1$gene_name,Genes1)
		Genes1_Body1 = Genes1_Body1[order(m)]
		print(Genes1_Body1$gene_name[!duplicated(Genes1_Body1$gene_name)])
		########### Postive ########## Good !!! ###
		###### integenetic peaks ################
		Inter_peak_list = All_peaks_list$Intergenic
		Inter_peak_list$peaks = as.character(Inter_peak_list)
		#### first filter p2g data.frame ######
		p2g_cl = p2g[which(p2g$gene %in% c(Genes1) == T),]
		p2g_clcl = p2g_cl[which(p2g_cl$peak %in% Inter_peak_list$peaks == T),]
		### p2g to GRanges #######
		p2g_clcl_GR = GRanges(p2g_clcl$peak)
		p2g_clcl_GR$peaks = p2g_clcl$peak
		p2g_clcl_GR$gene_name = p2g_clcl$gene
		p2g_clcl_GR$Correlation = p2g_clcl$Correlation
		#### filter with the distance ###########
		mm10_TSS_GR_all$peaks = as.character(mm10_TSS_GR_all)
		#######
		m = match(p2g_clcl_GR$gene_name,mm10_TSS_GR_all$gene_name)
		p2g_clcl_GR$TSS_region = mm10_TSS_GR_all$peaks[m]
		p2g_clcl_GR$distance = distance(GRanges(p2g_clcl_GR$peaks),GRanges(p2g_clcl_GR$TSS_region))
		p2g_clcl_GR = p2g_clcl_GR[which(p2g_clcl_GR$distance < distance_F)]
		########  #####
		Gene1_inter1 = p2g_clcl_GR[which(p2g_clcl_GR$gene_name %in% Genes1 == T)]
		m = match(Gene1_inter1$gene_name,Genes1)
		Gene1_inter1 = Gene1_inter1[order(m)] #######
		print(Gene1_inter1$gene_name[!duplicated(Gene1_inter1$gene_name)])
		#### Output the list !!!! ########
		Out_list = list(Genes1_TSS1,Genes1_Body1,Gene1_inter1)
		names(Out_list) <- c('Gene1_TSS','Gene1_Body','Gene1_Inter')
		return(Out_list)
}