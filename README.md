# IReNA-v2
### Constructing gene regulatory networks by integrating scRNA-seq and scATAC-seq data


 <div align="center">
 <img src="Summary.png" width="450" height = "450"/>
 </div>


## STEP 1:Selecting candidate genes
The DEGs were used as candidate genes for GRNs construction. For each developmental process which we aim to investigate in mouse and human, we identified the enriched genes for each cell type using the function ‘FindMarkers’ in Seurat. In constructing the GRNs of progenitors transition, the following parameters of ‘FindMarkers’ were used: min.pct = 0.05, logfc.threshold = 0.20, only.pos = TRUE, p-adjust < 0.01. In constructing GRNs regulating neurogenesis, the following parameters of ‘FindMarkers’ were used: min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE and p-adjust < 0.01.


## STEP 2:Identifying significant peak-to-gene links
We used the ArchR package to identify the significant peak-to-gene links. First, we integrated the age-matched scRNA-seq and scATAC-seq datasets for each time point using unconstrained Integration method with the function ‘addGeneIntegrationMatrix’. Then, using the function ‘addPeak2GeneLinks’, we calculated the correlation between accessibility peak intensity and gene expression. Finally, we identified the significant peak-to-gene links with the following cutoff: abs(correlation) > 0.2 and fdr < 1e-6.

## STEP 3:Identifying the potential cis-regulatory elements for each candidate gene
We identified potential cis-regulatory elements for each candidate gene based on their location and the peak-to-gene links from Step2. We first classified all peaks into three categories according to their genomic location related to their potential target genes: 1) Promoter. 2) Gene body. 3) Intergenic. For the peaks in the promoter region,we treated all of them as correlated accessible chromatin regions (CARs) of their overlapping target genes. For the peaks in the gene body region, we defined them as CARs of their overlapping genes if they met the following criteria: 1) the distance between the peak and the TSS of its overlapping gene is < 100kb. 2) the links between the peak and its overlapping gene is significant.  For the peaks in the intergenic region, we first find their target genes and construct the peak-gene pairs if the target genes’ TSS are located within the upstream 100kb or downstream 100 kb of the intergenic peaks. Then we keep the peak-gene pairs if their peak-to-gene links are significant in step2. These peaks were identified as CARs of their gene pairs.


##


##


##


##


##
