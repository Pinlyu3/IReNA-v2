

Process_DEGs_to_Celltypes <- function(x,index=c('RPC_S3','L_N','AC/HC','Rod')){
	####
	All_Genes = x$gene[!duplicated(x$gene)]
	print(length(All_Genes))
	####
	####
	tab = data.frame(matrix(0,nrow=length(All_Genes),ncol=length(index)))
	####
	tab$genes = All_Genes
	####
	colnames(tab)[1:length(index)] = index
	####
	for(i in 1:length(index)){
		sub_x = x[which(x$cluster == index[i]),]
		m = match(tab$genes,sub_x$gene)
		value = sub_x$avg_logFC[m]
		value[is.na(value)] <- 0
		tab[,i] = value
	}
	####
	tmp_mat = tab[,c(1:length(index))]
	tmp_num = apply(tmp_mat,1,function(x) which(x == max(x)))
	tmp_max = apply(tmp_mat,1,function(x) max(x))
	####
	tab$enrich = index[tmp_num]
	tab$max = tmp_max
	####
	return(tab)
}