# 3. singleR / celltypist

![WechatIMG365.jpeg](3%20singleR%20celltypist%201c035e5385be4d9692c1ee10dcd2b422/WechatIMG365.jpeg)

## 1. singleR

```r
# Clusters can be annotated by comparing markergenes from the dataset 
# and marker genes from reference dataset viaenrichment tests
# Jaccard index , 又称为Jaccard相似系数（Jaccard similarity coefficient）用于比较有限样本集之间的相似性与差异性。
# Jaccard系数值越大，样本相似度越高。
# 给定两个集合A,B，Jaccard 系数定义为A与B交集的大小与A与B并集的大小的比值

## build-in reference

# BiocManager::install("SingleR")
# conda install -c bioconda bioconductor-singler

library(celldex)
library(SingleR)

hpca.se = HumanPrimaryCellAtlasData()
Blue.se = BlueprintEncodeData() 
Immune.se = DatabaseImmuneCellExpressionData()
Nover.se = NovershternHematopoieticData()
MonacoIm.se = MonacoImmuneData()
save(hpca.se, Blue.se, Immune.se, Nover.se, MonacoIm.se, file='/glusterfs/home/local_pan_yuwen/resource/annotation/scRNA/singleR.anno.human.RData')

Mouse.se = MouseRNAseqData() #(鼠)
ImmGen.se = ImmGenData() #(鼠)
save(ImmGen.se, Mouse.se, file='/glusterfs/home/local_pan_yuwen/resource/annotation/scRNA/singleR.anno.mouse.RData')
```

```r
## human

library(celldex)
library(SingleR)

load('/glusterfs/home/local_pan_yuwen/resource/annotation/scRNA/singleR.anno.human.RData')

sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters = sce@meta.data$seurat_clusters

pred.hpca <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.main,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.hpca1 <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.fine,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

pred.Blue <- SingleR(test = sce_for_SingleR, ref = Blue.se, labels = Blue.se$label.main,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.Blue1 <- SingleR(test = sce_for_SingleR, ref = Blue.se, labels = Blue.se$label.fine,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

pred.Immune <- SingleR(test = sce_for_SingleR, ref = Immune.se, labels = Immune.se$label.main,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.Immune1 <- SingleR(test = sce_for_SingleR, ref = Immune.se, labels = Immune.se$label.fine,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

pred.Nover <- SingleR(test = sce_for_SingleR, ref = Nover.se, labels = Nover.se$label.main,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.Nover1 <- SingleR(test = sce_for_SingleR, ref = Nover.se, labels = Nover.se$label.fine,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

pred.MonacoIm <- SingleR(test = sce_for_SingleR, ref = MonacoIm.se, labels = MonacoIm.se$label.main,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.MonacoIm1 <- SingleR(test = sce_for_SingleR, ref = MonacoIm.se, labels = MonacoIm.se$label.fine,
                          clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType = data.frame(
	ClusterID = levels(sce@meta.data$seurat_clusters),
	hpca.main = pred.hpca$labels, 
	hpca.fine = pred.hpca1$labels, 
	Blue.main = pred.Blue$labels, 
	Blue.fine = pred.Blue1$labels, 
	Immune.main = pred.Immune$labels, 
	Immune.fine = pred.Immune1$labels, 
	Nover.main = pred.Nover$labels, 
	Nover.fine = pred.Nover1$labels, 
	MonacoIm.main = pred.MonacoIm$labels, 
	MonacoIm.fine = pred.MonacoIm1$labels
)
write.csv(cellType, "celltype_anno.csv", row.names=F)

sce@meta.data$singleR.main = factor(cellType[match(clusters, cellType$ClusterID),'hpca.main'])
sce@meta.data$singleR.fine = factor(cellType[match(clusters, cellType$ClusterID),'hpca.fine'])

sce@meta.data$singleR.immu_main = factor(cellType[match(clusters, cellType$ClusterID),'Immune.main'])
sce@meta.data$singleR.immu_fine = factor(cellType[match(clusters, cellType$ClusterID),'Immune.fine'])

save(sce, file='')

# Annotation diagnostics
# scores for all cells across all reference labels, which allows users to inspect the confidence of the predicted labels across the dataset. 
# Ideally, each cell (i.e., column of the heatmap) should have one score that is obviously larger than the rest, indicating that it is unambiguously assigned to a single label.
pdf('hpca.main.pdf')
plotScoreHeatmap(pred.hpca)
dev.off()
```

```r
## mouse

refdata = get(load("/glusterfs/home/local_deng_dongjie/Project/9HNSC/Download/SingleR_ref/ref_Mouse_all.RData"))
immudata = get(load("/glusterfs/home/local_deng_dongjie/Project/9HNSC/Download/SingleR_ref/ref_Mouse_imm.RData"))

sce = ifnb.combined
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters = sce@meta.data$seurat_clusters
save.image(file = "cluster1/sce_ifnb.combined.RData")

pred.mouseImmu <- SingleR(test = sce_for_SingleR, ref = immudata, labels = immudata$label.main,
                          clusters = clusters,
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.mouseImmu1 <- SingleR(test = sce_for_SingleR, ref = immudata, labels = immudata$label.fine,
                          clusters = clusters,
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.mouseRNA <- SingleR(test = sce_for_SingleR, ref = refdata, labels = refdata$label.main,
                         clusters = clusters,
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.mouseRNA1 <- SingleR(test = sce_for_SingleR, ref = refdata, labels = refdata$label.fine,
                         clusters = clusters,
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
 
cellType = data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    mouseImmu.main=pred.mouseImmu$labels,
                    mouseImmu.fine=pred.mouseImmu1$labels,
                    mouseRNA.main=pred.mouseRNA$labels,
                    mouseRNA.fine=pred.mouseRNA1$labels )
write.csv(cellType,"celltype_anno.csv",row.names=F)

new.cluster.ids <- c("Macrophages","Stem cells", "Monocytes","Macrophages","Stem cells",
		     "DCs","Fibroblasts","T cells","B cells","Macrophages","Neutrophils",
		     "Monocytes","Macrophages","NK","Monocytes","Endothelial cells","B cells",
		     "Fibroblasts")

cellType = data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels,
		    newcluster=new.cluster.ids)

sce@meta.data$singleR = factor(cellType[match(clusters, cellType$ClusterID), 'newcluster'])
```

![Untitled](3%20singleR%20celltypist%201c035e5385be4d9692c1ee10dcd2b422/Untitled.png)

![Untitled](3%20singleR%20celltypist%201c035e5385be4d9692c1ee10dcd2b422/Untitled%201.png)

## 2. celltypist

| model | description |
| --- | --- |
| Immune_All_Low.pkl | immune sub-populations combined from 20 tissues of 18 studies |
| Immune_All_High.pkl | immune populations combined from 20 tissues of 18 studies |
| Adult_Mouse_Gut.pkl | cell types in the adult mouse gut combined from eight datasets |
| Autopsy_COVID19_Lung.pkl | cell types from the lungs of 16 SARS-CoV-2 infected COVID-19 autopsy adult donors |
| COVID19_HumanChallenge_Blood.pkl | detailed blood cell states from 16 individuals after being challenged with SARS-CoV-2 |
| COVID19_Immune_Landscape.pkl | immune subtypes from lung and blood of COVID-19 patients and healthy controls |
| Cells_Fetal_Lung.pkl | cell types from human embryonic and fetal lungs |
| Cells_Intestinal_Tract.pkl | intestinal cells from fetal, pediatric (healthy and Crohn's disease) and adult human gut |
| Cells_Lung_Airway.pkl | cell populations from scRNA-seq of five locations of the human lungs and airways |
| Developing_Human_Brain.pkl | cell types from the first-trimester developing human brain |
| Developing_Human_Thymus.pkl | cell populations in embryonic, fetal, pediatric, and adult stages of the human thymus |
| Developing_Mouse_Brain.pkl | cell types from the embryonic mouse brain between gastrulation and birth |
| Healthy_COVID19_PBMC.pkl | peripheral blood mononuclear cell types from healthy and COVID-19 individuals |
| Human_IPF_Lung.pkl | cell types from idiopathic pulmonary fibrosis, chronic obstructive pulmonary disease, and healthy lungs of adult humans |
| Human_Lung_Atlas.pkl | integrated Human Lung Cell Atlas (HLCA) combining 46 datasets of the human respiratory system |
| Human_PF_Lung.pkl | cell types from different forms of pulmonary fibrosis lungs of adult humans |
| Lethal_COVID19_Lung.pkl | cell types from the lungs of individuals who died of COVID-19 and control individuals |
| Nuclei_Lung_Airway.pkl | cell populations from snRNA-seq of five locations of the human lungs and airways |
| Pan_Fetal_Human.pkl | stromal and immune populations from the human fetus |

```python
#### 文献信息
#https://mp.weixin.qq.com/s/LsDYAjKQ2QAeo8hUlm9YAQ
#Cross-tissue immune cell analysis reveals tissue-specific features in humans
#https://github.com/Teichlab/celltypist

conda create -n celltypist -y
conda activate celltypist
# conda install leidenalg
# conda install -c conda-forge python-annoy -y
pip install celltypist
#conda install -c bioconda -c conda-forge celltypist
```

```python
import celltypist
from celltypist import models

print(sys.path)

#Show all available models that can be downloaded and used.
models.models_description()
#Download available models
models.download_models() 

#Select the model from the above list
model = models.Model.load(model = 'Healthy_COVID19_PBMC.pkl')
#The model summary information.

############## Celltyping based on the input of count table
#demo test data. This is a UMI count csv file with cells as rows and gene symbols as columns.
#raw count matrix (reads or UMIs) 
input_file = celltypist.samples.get_sample_csv()

#Predict the identity of each input cell.
#predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl')
#Alternatively, the model argument can be a previously loaded `Model` as in 1.4.
#predictions = celltypist.annotate(input_file, model = model) #这个model也可以自己train
#In case your input file is a gene-by-cell mtx file.
#predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', transpose_input = True, gene_file = '/path/to/gene/file.txt', cell_file = '/path/to/cell/file.txt')

#In case your input file is a gene-by-cell table.
predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', transpose_input = True)

#Summary information for the prediction result.
predictions
#Examine the predicted cell type labels.
predictions.predicted_labels
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix

#Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
#Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', mode = 'prob match', p_thres = 0.5)

#Get an `AnnData` with predicted labels and confidence scores embedded into the observation metadata columns.
adata = predictions.to_adata(insert_labels = True, insert_conf = True)
#Get an `AnnData` with predicted labels, confidence scores, and decision matrix.
adata = predictions.to_adata(insert_labels = True, insert_conf = True, insert_decision = True)
#Get an `AnnData` with predicted labels, confidence scores, and probability matrix (recommended).
adata = predictions.to_adata(insert_labels = True, insert_conf = True, insert_prob = True)
adata.obs

#CellTypist provides a quick function to_plots to visualise your AnnotationResult and store the figures 
#without the need of explicitly transforming it into an AnnData. 可以直接做图而不用转换成AnnData的显式矩阵形式
#Visualise the predicted cell types overlaid onto the UMAP.
predictions.to_plots(folder = './', prefix = 'sampleTest',plot_probability=True) #图非常难看和乱，需要去R绘制

####################### 前面分析默认每个细胞独立，为了考虑相同细胞类型的细胞会聚类在一起，使用pass in majority_voting = True to the annotate function.
#Use a majority voting classifier combined with celltyping
#Turn on the majority voting classifier as well.
predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', majority_voting = True) #这个功能最好一开始传递一个聚类好的AnnData做input，否则会自己聚类

#Add your own over-clustering result.
#over_clustering argument. This argument can be specified in several ways:
#1. an input plain file with the over-clustering result of one cell per line.
#2. a string key specifying an existing cell metadata column in the AnnData (pre-created by the user).
#3. a list-like object (such as a numpy 1D array) indicating the over-clustering result of all cells.
#4. if none of the above is provided, will use a heuristic over-clustering approach, noted above.
predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', majority_voting = True, over_clustering = '/path/to/over_clustering/file')
predictions.to_plots(folder = './', prefix = 'sampleTest2') #图非常难看和乱，需要去R绘制

#There is also a min_prop parameter (defaults to 0) which controls the minimum proportion of cells from the dominant cell type required to name a given subcluster by this cell type. 
#Subcluster that fails to pass this proportion threshold will be assigned Heterogeneous.
predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', majority_voting = True, over_clustering = '/path/to/over_clustering/file',min_prop=10)

#Examine the predicted cell type labels.
predictions.predicted_labels
#Examine specifically the majority-voting results.
predictions.predicted_labels.majority_voting
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix
```

```r
############################# R下面运行

.libPaths()

library(Seurat)
library(dplyr)
library(ggplot2)
library(celldex)
library(SingleR)
library(future)

library(reticulate)
#use_python('/glusterfs/home/local_pan_yuwen/miniconda3/bin/python')

scanpy <- import("scanpy")
celltypist <- import("celltypist")
pandas <- import("pandas")
numpy <- import("numpy")

#library(cowplot)
#library(tidyverse)
#library(Matrix)
#library(ggpubr)

options(future.globals.maxSiz = 50 * 1024^3) #将全局变量上限调至50G，锚点整合很占内存
#导入并引用python包

#如果第一次运行需要下载训练好的模型数据
#celltypist$models$download_models(force_update = F) 

# 准备了一个seurat对象，这里是sce
sce <- readRDS("/titan3/local_yin_liufan/Project/BDGENE_10xSingleCell/analysis_merge/02_merge_seuratAnchor_GSE168732/filtered2/raw_CellType_obj.rds")

#使用scanpy构建一个AnnData
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(sce@assays$RNA@counts))),
                   obs = pandas$DataFrame(sce@meta.data),
                   var = pandas$DataFrame(data.frame(gene = rownames(sce@assays$RNA@counts),
                                                     row.names = rownames(sce@assays$RNA@counts)))
                   )
model = celltypist$models$Model$load(model = 'Healthy_COVID19_PBMC.pkl')
model$cell_types

adata = adata$copy() # to solve a bug, https://github.com/scverse/scanpy/issues/1183
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)

#细胞类型注释，使用seurat的聚类结果，按照majority_voting进行注释
predictions = celltypist$annotate(adata, model = 'Healthy_COVID19_PBMC.pkl', majority_voting = T, over_clustering = 'seurat_clusters')
sce = AddMetaData(sce, predictions$predicted_labels$majority_voting, "majority_voting")
saveRDS(sce,"./raw_CellType_obj_annoWithcelltypist.rds") #保存聚类后的seurat_obj格式文件，供下次读取

#注释完成，统计和做图
tempdf <- sce@meta.data[,c("seurat_clusters","celltype","majority_voting")]
tempdf <- tempdf[order(tempdf$seurat_clusters),]
tempdf <- tempdf[!duplicated(tempdf$seurat_clusters),]
write.csv(tempdf,"./celltype_annoWithCelltypist.csv",row.names=F)

#细胞注释作图
p3 = DimPlot(sce, reduction = "umap", pt.size=0.6,label=T, label.size=7,repel=T,)+
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18),
    axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
    legend.title=element_text(size=18),legend.text=element_text(size=18))
cellTypeList <- c("gdT","CD4.Naive","CD4.CM","CD8.TE","CD8.Naive","NK_56hi","NK_16hi","DC3","pDC","B_naive","B_switched_memory","Plasma_cell_IgA","CD14_mono","CD16_mono","Platelets","RBC")
sce@meta.data$majority_voting=factor(sce@meta.data$majority_voting,levels=as.character(cellTypeList))
colsList <- c("yellow","#EA867D","pink","#6DB9CD","#74C580","#CB8339","#8E8EC8","grey","#10AEE5","#7DCDAC","light blue","#DD89C6","#AF913C","6AA7CF","#CA65B6","#869C23")

plot = DimPlot(sce,reduction = "umap",group.by = 'majority_voting',pt.size=0.6,label=T, label.size=7,repel=T,
  cols=colsList) + 
theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18),
    axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
    legend.title=element_text(size=18),legend.text=element_text(size=18))
plot=p3+plott
ggsave(filename = './pbmc_umap_celltypistIdentify.pdf', height = 8, width = 19, plot = plot)
```