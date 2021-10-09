#### select random cells for each cell type ####

random_cells_by_celltypes = function(x,celltypes){
	x$cell_id = colnames(x)
	####
	x_cl = subset(x,subset= New_celltypes %in% celltypes == T)
	####
	min_numbers = min(table(x_cl$New_celltypes))
	#### sample for each cell types #####
	cell_list = c()
	####
	for(i in 1:length(celltypes)){
		tmp_celltypes = celltypes[i]
		print(tmp_celltypes)
		tmp_cells = colnames(x_cl)[which(x_cl$New_celltypes == celltypes[i])]
		print(length(tmp_cells))
		index_sample = sample(1:length(tmp_cells),min_numbers,replace=F)
		###
		tmp_cells_choose = tmp_cells[index_sample]
		##
		cell_list = c(cell_list,tmp_cells_choose)
		###
		print(length(cell_list))
	}
	print(head(cell_list))
	#k = which(x_cl$cell_id %in% cell_list == T)
	x_cl_choose = subset(x_cl,subset = cell_id %in% cell_list)
	print(table(x_cl_choose$New_celltypes))
	####
	return(x_cl_choose)
	####
}


#### calculate the corelation between genes ####

RNA_Corr_Add_cutoff <- function(smooth_matrix){
	smooth_matrix_ori = smooth_matrix
	####
	k_index = which(smooth_matrix@x > 0)
	smooth_matrix_k = smooth_matrix
	smooth_matrix_k@x[k_index] = 1
	####
	k1 = which(Matrix::rowSums(smooth_matrix_k) < 50)
	k2 = which(is.na(Matrix::rowSums(smooth_matrix_k)) ==T)
	####
	genes_k1 = rownames(smooth_matrix_k)[k1]
	genes_k2 = rownames(smooth_matrix_k)[k2]
	####
	k_Cells = which(rownames(smooth_matrix) %in% c(genes_k1,genes_k2) == T)
	smooth_matrix = smooth_matrix[-k_Cells,]
	####
	print(dim(smooth_matrix))
	####
	Input = t(as.matrix(smooth_matrix))
	####
	Corr_res = sparse.cor3(Input)
	print(dim(Corr_res))
	####
	Corr_res[lower.tri(Corr_res,diag=T)]= 2
	####
	library(reshape2)
	Corr_res_out = melt(Corr_res)
	####
	k3 = which(Corr_res_out$value == 2)
	Corr_res_out_cl = Corr_res_out[-k3,]
	####
	####
	cutoff = quantile(Matrix::rowSums(smooth_matrix),0.3)
	k3 = which(Matrix::rowSums(smooth_matrix) > cutoff)
	genes_k3 = rownames(smooth_matrix)[k3]
	####
	k4 = which((Corr_res_out_cl$Var1 %in% genes_k3 == T) & (Corr_res_out_cl$Var2 %in% genes_k3 == T))
	Corr_res_out_clcl = Corr_res_out_cl[k4,]
	#########
	left = quantile(Corr_res_out_clcl$value,0.025)
	right = quantile(Corr_res_out_clcl$value,0.975)
	####
	Corr_res_out_cl$tag = 'No'
	k_pos = which(Corr_res_out_cl$value > right)
	k_neg = which(Corr_res_out_cl$value < left)
	Corr_res_out_cl$tag[k_pos] = 'pos'
	Corr_res_out_cl$tag[k_neg] = 'neg'
	####
	return(Corr_res_out_cl)
} 


#### calculate the correlations #####

sparse.cor3 <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    cSums <- colSums(x)
    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp
    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}
