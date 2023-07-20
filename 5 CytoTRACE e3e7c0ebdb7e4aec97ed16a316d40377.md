# 5. CytoTRACE

**More diversity at the top**

> CytoTRACE, a computational framework based on the simple observation that transcriptional diversity—the number of genes expressed in a cell—decreases during differentiation.
> 

代码解析：

[https://www.jianshu.com/p/bff720322815](https://www.jianshu.com/p/bff720322815)

## 0. install

```bash
# shell
conda create -n cytotrace -y  # python3.7
conda activate cytotrace

conda install numpy
conda install -c conda-forge python-annoy -y
pip install scanoramaCT

## fix the bug related to the changes in python3.10
## no need for python3.7
# vim /glusterfs/home/local_pan_yuwen/miniconda3/lib/python3.10/site-packages/intervaltree/intervaltree.py
# collections.MutableSet -> collections.abc.MutableSet

conda install R -y  # R 4.2
conda install -c conda-forge r-devtools -y

conda install -c conda-forge r-hiclimr -y
conda install -c r r-ccapp -y
conda install -c conda-forge r-reticulate -y
conda install -c conda-forge r-ggpubr -y
conda install -c conda-forge r-plyr -y
conda install -c conda-forge r-egg -y
conda install -c bioconda bioconductor-sva -y
conda install -c bioconda r-seurat -y

# R
#install.packages("devtools")
#install.packages('BiocManager')

# https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
devtools::install_local("/glusterfs/home/local_pan_yuwen/software/CytoTRACE_0.3.3.tar.gz")

# library(reticulate)
# Sys.setenv(RETICULATE_PYTHON="/PATH/TO/PYTHON/BIN")
```

## 1. cytotrace

```bash
conda activate cytotrace
```

```r
library(Seurat)
library(CytoTRACE)

# CytoTRACE: 对自定义 scRNA-seq 数据集进行CytoTRACE 分析的函数
# iCytoTRACE: 跨多个 scRNA-seq 批次/数据集运行 CytoTRACE 的功能
# plotCytoTRACE: 用于生成 CytoTRACE、表型和基因表达的 2D 可视化的函数

table(Idents(sce))
p1 = DimPlot(sce, label = T)

####提取表型 & 坐标文件
table(pbmc3k.final$seurat_annotations)
phe <- pbmc3k.final$seurat_annotations
phe = as.character(phe)
names(phe) <- rownames(pbmc3k.final@meta.data)

emb <- data.frame(Embeddings(pbmc3k.final, reduction = "umap"))

####提取表达矩阵
mat_3k <- as.matrix(pbmc3k.final@assays$RNA@counts)
mat_3k[1:4,1:4]

# 使用“ncores”（默认值 = 1）进行多线程
# 使用“subsamplingsize”（默认值 = 1,000 个单元格）指示子采样大小
results <- CytoTRACE(mat = mat_3k, ncores = 8, subsamplesize = 1000)

## 在多个 scRNA-seq 批次/数据集上运行
# datasets <- list(marrow_10x_expr, marrow_plate_expr)
# results <- iCytoTRACE(datasets)

###############################################################
## (1) 第一项表示：每个细胞的轨迹分数，取值范围在0~1。越接近0，分化潜能越大；反之越小。
head(results[[1]])
# 10X_P7_2_AAACCTGCAGTAACGG 10X_P7_2_AAACGGGAGGACGAAA 10X_P7_2_AAACGGGAGGTACTCT
#                0.62346760                0.35901926                0.73555166
# 10X_P7_2_AAACGGGAGGTGCTTT 10X_P7_2_AAACGGGAGTCGAGTG 10X_P7_2_AAAGATGAGCTTCGCG
#                0.09457093                0.23467601                0.42644483

## (2) 第二项表示：每个细胞的分化潜能的排名
head(results[[2]])
# 10X_P7_2_AAACCTGCAGTAACGG 10X_P7_2_AAACGGGAGGACGAAA 10X_P7_2_AAACGGGAGGTACTCT
#                       713                       411                       841
# 10X_P7_2_AAACGGGAGGTGCTTT 10X_P7_2_AAACGGGAGTCGAGTG 10X_P7_2_AAAGATGAGCTTCGCG
#                       109                       269                       488

## (3) 第三项表示：每个基因表达与细胞群分化轨迹相关性
head(results[[3]])
#      Rpl4    Eef1a1      Rps5     Rps3a    Rpl13a     Rps4x
# 0.9234229 0.9225254 0.9222174 0.9216470 0.9196403 0.9161128
###############################################################

# phenotype 细胞类型注释
# gene 是否映射特定基因表达
#	emb 是否提供细胞降维坐标
#	outputDir 图片储存路径
plotCytoTRACE(results, emb = emb, phenotype = phe, gene = "CCR7", outputDir = "pbmc_3k")

plotCytoGenes(results, numOfGenes = 10,outputDir = "pbmc_3k")
```