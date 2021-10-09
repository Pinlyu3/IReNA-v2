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