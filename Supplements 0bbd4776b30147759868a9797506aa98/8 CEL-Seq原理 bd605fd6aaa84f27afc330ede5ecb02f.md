# 8. CEL-Seq原理

CEL-Seq（Cell expression by linear amplification and sequencing）是和SMART-seq同年发表的单细胞测序技术。在2012年发表在Cell Reports上。

![https://upload-images.jianshu.io/upload_images/15771939-373dc6e26cf4e020.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-373dc6e26cf4e020.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

2016年作者对CEL-Seq进行了改进，CEL-Seq2发表在Genome Biology上面。

![https://upload-images.jianshu.io/upload_images/15771939-93d1cb0e4046f838.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-93d1cb0e4046f838.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

与[SMARTseq](https://www.jianshu.com/p/6c5d663433a4)必须是以单细胞为体系进行反应不同，CEL-seq在逆转录中的oligo(dT)VN Primer增加了条形码标签。

![https://upload-images.jianshu.io/upload_images/15771939-d231fe7e47604f6e.png?imageMogr2/auto-orient/strip|imageView2/2/w/793/format/webp](https://upload-images.jianshu.io/upload_images/15771939-d231fe7e47604f6e.png?imageMogr2/auto-orient/strip|imageView2/2/w/793/format/webp)

步骤：

- 细胞裂解：和SMARTseq一样，CEL-seq需要首先分离获得单个细胞，并将单个细胞转移至低渗裂解缓冲溶液（裂解液中含有RNase抑制剂）中进行裂解，将转录组mRNA释放到裂解液中。
- 逆转录合成第一链：加入逆转录酶、dNTP、含有Barcode的oligo(dT)VN Primer等，对mRNAs进行反转录，获得cDNA第一条链。（由于第一链上带有barcode标签，这一步之后就可以将样品混合处理）
- 线性扩增。利用T7 RNA聚合酶对cDNA进行转录反应;
- 片段化及添加接头。对RNA进行片段化，连接Illumina 3测序接头；
- 反转录RNA形成DNA
- 筛选文库。利用PCR反应和Illimina 3‘接头及Illumina 5’接头对文库进行筛选;
- 获得的文库进行双端测序。其中R1序列包含barcode序列，而R2包含mRNA序列。

CEL-Seq2在CEL-Seq的基础上增加了UMI序列标记转录本，提高了准确性，并显著提高了RT效率，从而提高了检测灵敏度。