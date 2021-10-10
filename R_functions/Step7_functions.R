#### found feedback pairs ####

FoundFeedBackPairsOne <- function(tmp_motif1,tmp_motif2){
	##### Pos ###############
	k1 = which(tmp_motif1$Tag == 'pos')
	k2 = which(tmp_motif1$Tag == 'neg')
	tmp_motif1_pos = tmp_motif1[k1,]
	tmp_motif1_neg = tmp_motif1[k2,]
	##### Neg ###############
	k3 = which(tmp_motif2$Tag == 'pos')
	k4 = which(tmp_motif2$Tag == 'neg')
	tmp_motif2_pos = tmp_motif2[k3,]
	tmp_motif2_neg = tmp_motif2[k4,]
	##### Pos pairs ##########
	Pos_index1_rev = paste(tmp_motif1_pos$Target,tmp_motif1_pos$TFs)
	##### Overlapped only #########
	Pos_index1_rev_index = which(Pos_index1_rev %in% tmp_motif2_pos$index == T)
	print(length(Pos_index1_rev_index))
	if(length(Pos_index1_rev_index) > 0){
		Pos_index1_rev_res = tmp_motif1_pos[Pos_index1_rev_index,]
		Pos_index1_rev_res = Pos_index1_rev_res[!duplicated(Pos_index1_rev_res$index),]
		Pos_index1_rev_res = Pos_index1_rev_res[,c(1,6,9,8)]
	}else{
		Pos_index1_rev_res = data.frame(TFs = 'ND',Target='ND',Tag='ND',Cor=0)
	}
	##### Neg Pairs ##########
	Neg_index1_rev = paste(tmp_motif1_neg$Target,tmp_motif1_neg$TFs)
	Neg_index1_rev_index = which(Neg_index1_rev %in% tmp_motif2_neg$index == T)
	print(length(Neg_index1_rev_index))
	if(length(Neg_index1_rev_index) > 0){
		Neg_index1_rev_res = tmp_motif1_neg[Neg_index1_rev_index,]
		Neg_index1_rev_res = Neg_index1_rev_res[!duplicated(Neg_index1_rev_res$index),]
		Neg_index1_rev_res = Neg_index1_rev_res[,c(1,6,9,8)]
	}else{
		Neg_index1_rev_res = data.frame(TFs = 'ND',Target='ND',Tag='ND',Cor=0)
	}
	######
	Res = rbind(Pos_index1_rev_res,Neg_index1_rev_res)
	Res_nd = which(Res$TFs == 'ND')
	if(length(Res_nd) > 0){
		Res = Res[-Res_nd,]
	}
	return(Res)
}



#### find feedback pairs from the GRNs list ####
FoundFeedBackPairs_new <- function(Motif_list){
	######
	Out_list = list()
	######
	for(i in 1:length(Motif_list)){
		for(j in 1:length(Motif_list)){
			print(c(i,j))
			Names_1 = names(Motif_list)[i]
			Names_2 = names(Motif_list)[j]
			print(c(Names_1,Names_2))
			tmp_motif1 = Motif_list[[i]]
			tmp_motif2 = Motif_list[[j]]
			#####
			tmp_out = FoundFeedBackPairsOne(tmp_motif1,tmp_motif2)
			tmp_out$celltypes = paste(Names_1,Names_2,sep=':')
			tmp_out$TFs = paste(Names_1,tmp_out$TFs,sep=':')
			tmp_out$Target = paste(Names_2,tmp_out$Target,sep=':')
			Out_list = c(Out_list,list(tmp_out))
		}
	}
	#######
	Out_list_out = do.call("rbind",Out_list)
}


### add annotations ####
Process_the_Feedback_res = function(RPCMG_Feedback_res){
	#######
	sp_TFs = strsplit(RPCMG_Feedback_res$TFs,split=':',fixed=T)
	#######
	sp_Target = strsplit(RPCMG_Feedback_res$Target,split=':',fixed=T)
	#######
	RPCMG_Feedback_res$TF_index = sapply(sp_TFs,function(x) x[[1]])
	RPCMG_Feedback_res$TF_gene = sapply(sp_TFs,function(x) x[[2]])
	RPCMG_Feedback_res$Target_index = sapply(sp_Target,function(x) x[[1]])
	RPCMG_Feedback_res$Target_gene = sapply(sp_Target,function(x) x[[2]])
	#######
	RPCMG_Feedback_res$TF_target = paste(RPCMG_Feedback_res$TF_gene,RPCMG_Feedback_res$Target_gene,sep='::')
	######
	return(RPCMG_Feedback_res)
}

