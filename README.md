# IReNA-v2
### Constructing gene regulatory networks by integrating scRNA-seq and scATAC-seq data


 <div align="center">
 <img src="Summary.png" width="450" height = "450"/>
 </div>

## Before following the IReNA-v2 analysis pipline:
We use E14-E16 scRNAseq/scATACseq datasets as example datasets. Seurat objects and ArchR objects can be downloaded by google drive: [Example datasets](https://drive.google.com/drive/folders/1BMwEuVM72ThIJj5MwqUAmGuhcvN-WChF?usp=sharing)

## STEP 1:Selecting candidate genes
The DEGs were used as candidate genes for GRNs construction. For each developmental process which we aim to investigate in mouse and human, we identified the enriched genes for each cell type using the function ‘FindMarkers’ in Seurat. In constructing the GRNs of progenitors transition, the following parameters of ‘FindMarkers’ were used: min.pct = 0.05, logfc.threshold = 0.20, only.pos = TRUE, p-adjust < 0.01. In constructing GRNs regulating neurogenesis, the following parameters of ‘FindMarkers’ were used: min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE and p-adjust < 0.01.


## STEP 2:Identifying significant peak-to-gene links
We used the ArchR package to identify the significant peak-to-gene links. First, we integrated the age-matched scRNA-seq and scATAC-seq datasets for each time point using unconstrained Integration method with the function ‘addGeneIntegrationMatrix’. Then, using the function ‘addPeak2GeneLinks’, we calculated the correlation between accessibility peak intensity and gene expression. Finally, we identified the significant peak-to-gene links with the following cutoff: abs(correlation) > 0.2 and fdr < 1e-6.

## STEP 3:Identifying the potential cis-regulatory elements for each candidate gene
We identified potential cis-regulatory elements for each candidate gene based on their location and the peak-to-gene links from Step2. We first classified all peaks into three categories according to their genomic location related to their potential target genes: 1) Promoter. 2) Gene body. 3) Intergenic. For the peaks in the promoter region,we treated all of them as correlated accessible chromatin regions (CARs) of their overlapping target genes. For the peaks in the gene body region, we defined them as CARs of their overlapping genes if they met the following criteria: 1) the distance between the peak and the TSS of its overlapping gene is < 100kb. 2) the links between the peak and its overlapping gene is significant.  For the peaks in the intergenic region, we first find their target genes and construct the peak-gene pairs if the target genes’ TSS are located within the upstream 100kb or downstream 100 kb of the intergenic peaks. Then we keep the peak-gene pairs if their peak-to-gene links are significant in step2. These peaks were identified as CARs of their gene pairs.


## STEP 4:Predicting cell-type specific TFs binding in cis-regulatory elements
With the cis-regulatory elements identified in Step 3, we next predicted the TF binding in these elements for each cell type with the PWMs extracted from TRANSFAC database. Firstly, we searching the motifs in all the cis-regulatory elements with the function ‘matchMotifs (p.cutoff = 5e-05)’ from the motifmatchr package. Then we filtered these motif regions according to their footprint score and their corresponding TF’s expression for each cell type.

To calculate the footprint score for each motif region in each cell type, we re-grouped the insertion fragments based on their origin of cell type and converted these cell-type-specific fragments into bam files using a custom script. Then we fed the bam files to TOBIAS software and obtained the bias-corrected Tn5 signal (log2(obs/exp)) with the default parameters except: ATACorrect --read_shift 0 0. Next, we calculated footprint scores including  NC, NL and NR for each motif's binding region. NC indicated the average bias-corrected Tn5 signal in the center of the motif. NL and NR indicated the average bias-corrected Tn5 signal in the left and right flanking regions of the motif, respectively. The flanking region is triple the size of the center region. We kept the motifs with the following criteria: NC < -0.1 and NL > 0.1 and NR > 0.1.

We further removed the motifs binding region for each cell type if the expression level of their corresponding TFs are not enriched in that cell type (from Step1).



## STEP 5:Calculating gene-gene correlation
We calculated the expression correlations between all the expressed genes at the single-cell level. First, we extracted the cell-by-matrix from Seurat objects and filtered out the non-expressed genes in the matrix (rowSums < 10). Then we applied the MAGIC software to impute missing values and recover the gene interactions with the cell-by-gene matrix. The output matrix from MAGIC was used to calculate gene-gene correlation using the function ‘cor’ in R.  To identify the significant gene-gene correlations, we ranked all the gene-gene correlations (~1X10e8). The top 2.5% correlations were treated as significant positive correlations (p < 0.025) and the bottom 2.5% correlations were treated as significant negative correlations  (p < 0.025).



## STEP 6:Constructing gene regulatory networks
By integrating data from Step1-Step5, We constructed cell-type specific GRNs with the following procedure:
We first obtained the peak-target links from Step 3, and cell-type specific TF-peak links from Step 4.  We then merged these 2 types of links to the cell-type specific TF-peak-target relationships. Next, we classified these TF-peak-target relationships into activation or repression relationships based on the sign of the expression correlation between TF and target from Step 5. The significant positive/negative correlated TF-targets were selected as the active/repressive regulations respectively.
Finally, we removed all the duplicated TF-target regulatory relationships for each cell type and merged them to the final GRNs which were used for the downstream analysis.



## STEP 7:Identifying and visualizing feedback TF pairs
With the GRNs constructed in the previous steps, we searched for TF pairs connected by either positive or negative feedback regulatory relationships. The TF pairs that activated each other were identified as ‘double positive’ pairs and the TF pairs repressed each other were identified as ‘double negative’ pairs. We visualized these feedback TFs pairs using Cytoscape software.
