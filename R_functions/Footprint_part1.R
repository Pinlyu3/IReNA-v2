#####
#####
##### fragments to bam files #####
#####
#####

conda activate Signac2

R

library('readr')
library('GenomicRanges')


Read_fragment_to_GR <- function(x){
    tab <- read_tsv(x,col_names = F)
    GR <- GRanges(seqnames=tab$X1,ranges=IRanges(start=tab$X2,end=tab$X3),index=tab$X4,dup=tab$X5)
    return(GR)
}


Fragment_list = list()
tags = c('E11','E12','E14','E16','E18','P0','P2','P5','P8','P11','P14')
for(i in tags){
    tag = i
    file = 'fragments_cl.tsv'
    print(i)
    folder = paste('/zp1/data/plyu3/scATACseq_new_clean/',tag,sep='')
    print(folder)
    setwd(folder)
    tab_GR <- Read_fragment_to_GR(file)
    Fragment_list = c(Fragment_list,tab_GR)
}
names(Fragment_list) = tags

#########################
#########################

#### 输出所有细胞系的 fragments ####
#### ##### #### ##### ##### #####
#### load the celltype annotation #####

ATAC_meta <- read.delim('/zp1/data/plyu3/Arrow_Project/New_Figure2/ATAC_Meta_tab.txt',sep='\t',header=T)
table(ATAC_meta$celltypes)
k1 = which(ATAC_meta$celltypes == 'Early_Cone')
k2 = which(ATAC_meta$celltypes == 'Early_Rod')
ATAC_meta$celltypes[k1] = 'Cone'
ATAC_meta$celltypes[k2] = 'Rod'
ATAC_meta$celltypes = as.character(ATAC_meta$celltypes)
table(ATAC_meta$celltypes)


####
####
Cell_names_list = split(ATAC_meta$cell_id,ATAC_meta$celltypes)


#####
#####
#####
unlist_fun <- function(GR_list){
	GR_Ori = GR_list[[1]]
	for(i in 2:length(GR_list)){
		GR_Ori = c(GR_Ori,GR_list[[i]])
	}
	return(GR_Ori)
}

Fragments_filter = Fragment_list
Total_insertion_list = list()

for(i in 1:length(Cell_names_list)){
	print(paste('i=',i))
	###
	tmp_name = names(Cell_names_list)[i]
	###
	tmp_name_list = Cell_names_list[[i]]
	###
	print(tmp_name)
	### how many times #####
	tmp_times = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[1]])
	tmp_indexs = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[2]])
	### how many times #####
	print(levels(as.factor(tmp_times)))
	### split cell names by time ######
	print(tmp_name)
	print(table(tmp_times))
	tmp_indexs_list = split(tmp_indexs,tmp_times)
	###
	tmp_GR_list = list()
	###
	for(j in 1:length(tmp_indexs_list)){
		print(paste('j=',j))
		####
		tmp_times_tmp = names(tmp_indexs_list)[j]
		####
		print(names(tmp_indexs_list)[j])
		####
		tmp_indexs_tmp = tmp_indexs_list[[j]]
		####
		k_tmptmp = which(names(Fragments_filter) == tmp_times_tmp)
		####
		tmp_fragments = Fragments_filter[[k_tmptmp]]
		####
		tmp_GR = tmp_fragments[which(tmp_fragments$index %in% tmp_indexs_tmp == T)]
		####
		print(length(tmp_GR))
		####
		tmp_GR_list = c(tmp_GR_list,list(tmp_GR))
		####
		####
	}
	tmp_GR_list = unlist_fun(tmp_GR_list)
	###
	Total_insertion_list = c(Total_insertion_list,list(tmp_GR_list))
}

####
####
names(Total_insertion_list) = names(Cell_names_list)

####
####
####
setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/All')

save(Total_insertion_list,file='Total_insertion_list')


#### early ####
library('Seurat')

load('/zp1/data/plyu3/Arrow_Project/New_Figure5_202009/E14_E16_ATAC_seurat')

head(E14_E16_ATAC_seurat@meta.data)

Cell_names_list = split(E14_E16_ATAC_seurat$cell_id,E14_E16_ATAC_seurat$New_celltypes)

Total_insertion_list = list()

for(i in 1:length(Cell_names_list)){
	print(paste('i=',i))
	###
	tmp_name = names(Cell_names_list)[i]
	###
	tmp_name_list = Cell_names_list[[i]]
	###
	print(tmp_name)
	### how many times #####
	tmp_times = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[1]])
	tmp_indexs = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[2]])
	### how many times #####
	print(levels(as.factor(tmp_times)))
	### split cell names by time ######
	print(tmp_name)
	print(table(tmp_times))
	tmp_indexs_list = split(tmp_indexs,tmp_times)
	###
	tmp_GR_list = list()
	###
	for(j in 1:length(tmp_indexs_list)){
		print(paste('j=',j))
		####
		tmp_times_tmp = names(tmp_indexs_list)[j]
		####
		print(names(tmp_indexs_list)[j])
		####
		tmp_indexs_tmp = tmp_indexs_list[[j]]
		####
		k_tmptmp = which(names(Fragments_filter) == tmp_times_tmp)
		####
		tmp_fragments = Fragments_filter[[k_tmptmp]]
		####
		tmp_GR = tmp_fragments[which(tmp_fragments$index %in% tmp_indexs_tmp == T)]
		####
		print(length(tmp_GR))
		####
		tmp_GR_list = c(tmp_GR_list,list(tmp_GR))
		####
		####
	}
	tmp_GR_list = unlist_fun(tmp_GR_list)
	###
	Total_insertion_list = c(Total_insertion_list,list(tmp_GR_list))
}

names(Total_insertion_list) = names(Cell_names_list)

Early_Total_insertion_list = Total_insertion_list

setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/Early')

save(Early_Total_insertion_list,file='Early_Total_insertion_list')



#### late1 ####

library('Seurat')

load('/zp1/data/plyu3/Arrow_Project/New_Figure5_202009/E18_P2_ATAC_seurat')

head(E18_P2_ATAC_Seurat@meta.data)

Cell_names_list = split(E18_P2_ATAC_Seurat$cell_id,E18_P2_ATAC_Seurat$New_celltypes)

Fragments_filter = Fragment_list

Total_insertion_list = list()

for(i in 1:length(Cell_names_list)){
	print(paste('i=',i))
	###
	tmp_name = names(Cell_names_list)[i]
	###
	tmp_name_list = Cell_names_list[[i]]
	###
	print(tmp_name)
	### how many times #####
	tmp_times = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[1]])
	tmp_indexs = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[2]])
	### how many times #####
	print(levels(as.factor(tmp_times)))
	### split cell names by time ######
	print(tmp_name)
	print(table(tmp_times))
	tmp_indexs_list = split(tmp_indexs,tmp_times)
	###
	tmp_GR_list = list()
	###
	for(j in 1:length(tmp_indexs_list)){
		print(paste('j=',j))
		####
		tmp_times_tmp = names(tmp_indexs_list)[j]
		####
		print(names(tmp_indexs_list)[j])
		####
		tmp_indexs_tmp = tmp_indexs_list[[j]]
		####
		k_tmptmp = which(names(Fragments_filter) == tmp_times_tmp)
		####
		tmp_fragments = Fragments_filter[[k_tmptmp]]
		####
		tmp_GR = tmp_fragments[which(tmp_fragments$index %in% tmp_indexs_tmp == T)]
		####
		print(length(tmp_GR))
		####
		tmp_GR_list = c(tmp_GR_list,list(tmp_GR))
		####
		####
	}
	tmp_GR_list = unlist_fun(tmp_GR_list)
	###
	Total_insertion_list = c(Total_insertion_list,list(tmp_GR_list))
}

names(Total_insertion_list) = names(Cell_names_list)

Late1_Total_insertion_list = Total_insertion_list

setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/Late1')

save(Late1_Total_insertion_list,file='Late1_Total_insertion_list')



#### late2 ####

library('Seurat')

load('/zp1/data/plyu3/Arrow_Project/New_Figure5_202009/P5_P8_ATAC_Seurat')

head(P5_P8_ATAC_Seurat@meta.data)

k1 = which(P5_P8_ATAC_Seurat$New_celltypes %in% c('RPC_S3','MG'))

P5_P8_ATAC_Seurat$New_celltypes[k1] = 'RPCS3_MG'

Cell_names_list = split(P5_P8_ATAC_Seurat$cell_id,P5_P8_ATAC_Seurat$New_celltypes)

Fragments_filter = Fragment_list

Total_insertion_list = list()

for(i in 1:length(Cell_names_list)){
	print(paste('i=',i))
	###
	tmp_name = names(Cell_names_list)[i]
	###
	tmp_name_list = Cell_names_list[[i]]
	###
	print(tmp_name)
	### how many times #####
	tmp_times = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[1]])
	tmp_indexs = sapply(strsplit(as.character(tmp_name_list),split=':'),function(x) x[[2]])
	### how many times #####
	print(levels(as.factor(tmp_times)))
	### split cell names by time ######
	print(tmp_name)
	print(table(tmp_times))
	tmp_indexs_list = split(tmp_indexs,tmp_times)
	###
	tmp_GR_list = list()
	###
	for(j in 1:length(tmp_indexs_list)){
		print(paste('j=',j))
		####
		tmp_times_tmp = names(tmp_indexs_list)[j]
		####
		print(names(tmp_indexs_list)[j])
		####
		tmp_indexs_tmp = tmp_indexs_list[[j]]
		####
		k_tmptmp = which(names(Fragments_filter) == tmp_times_tmp)
		####
		tmp_fragments = Fragments_filter[[k_tmptmp]]
		####
		tmp_GR = tmp_fragments[which(tmp_fragments$index %in% tmp_indexs_tmp == T)]
		####
		print(length(tmp_GR))
		####
		tmp_GR_list = c(tmp_GR_list,list(tmp_GR))
		####
		####
	}
	tmp_GR_list = unlist_fun(tmp_GR_list)
	###
	Total_insertion_list = c(Total_insertion_list,list(tmp_GR_list))
}

names(Total_insertion_list) = names(Cell_names_list)

Late2_Total_insertion_list = Total_insertion_list

setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/Late2')

save(Late2_Total_insertion_list,file='Late2_Total_insertion_list')




/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/
/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/
/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/
/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/



setwd('/zp1/data/plyu3/Human_scATACseq_new/scATAC/Hint/All')












