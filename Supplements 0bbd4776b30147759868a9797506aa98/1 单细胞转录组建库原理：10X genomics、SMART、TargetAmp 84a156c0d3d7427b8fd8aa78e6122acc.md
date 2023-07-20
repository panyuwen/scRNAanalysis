# 1. 单细胞转录组建库原理：10X genomics、SMART、TargetAmp

> 要实现单细胞mRNA测序，需要解决2个难题：
> 
> - 1. 一个人类细胞中，RNA总量大约只有10pg左右（1 pg = 10^-12 g），其中mRNA的量大约只有0.2pg。要把这么少的mRNA转变成越零点几个ug（1 ug = 10^-6 g）以上的核酸文库，意味着核酸的扩增量要达到几百万倍以上。 如何在核酸扩增过程中不引入太多的PCR偏差，一直是个大问题。所谓PCR偏差，就是在PCR扩增过程中，某些片段被大量扩增，而大部分片段被扩增的量很少，甚至没有被扩增。这就导致高通量测序只能测到所有样本中很少的一部分片段序列。PCR偏差会随着PCR循环次数的增多而指数放大。那么在这种情况下，一方面要把核酸扩增几百万倍甚至更多的倍数，另一方面又想得到均一覆盖的文库，这就是单细胞mRNA建库当中，所要解决的第一个大难题。
> - 2. 第二个难题是如何尽可能高效地得到mRNA文库，而不是含了大量rRNA序列的文库。因为rRNA在总RNA当中占了95%甚至更高的比例，而 mRNA在总RNA中只占了2%-3%的比例。如果不加区分的进行逆转录再扩增建库，很可能测序得到的绝大部分序列都是rRNA的序列。但是rRNA序列不能给我们带来有效的生物信息，只有mRNA序列，才是我们想要的信息。因此，如何能够选择性地把mRNA转化成测序文库，并且避免把rRNA带到测序文库中来，这就是单细胞mRNA建库当中，要解决的第二个大难题。

# 1**. 10X genomics**

## 1**.1 简介**

10X genomics是把微珠加DNA标签、微滴发生、酶反应和高通量测序后的数据分析这一系列的技术整合在一起的一个基于油包水乳浊液酶反应原理的分子生物学分析系统。该方法基于`微流控技术`（Microfluidics-based approaches）[4]，与SMART有相似的分子生物学原理，运用了模板转换技术，但与SMART的细胞捕获和通量不同。Droplet-based方法是将单个细胞包裹在一个小油滴中（含有barcode和RT primer）反转录成cDNA，然后油滴破裂释放cDNA，统一进行文库构建，增大了实验通量，但需要专门的实验设备。

![https://upload-images.jianshu.io/upload_images/15771939-514e6c6d744a01b5.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-514e6c6d744a01b5.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

## 1**.2 10X工作原理**

### **Gel beads及其核酸序列构成**

Gel beads，即凝胶微珠。每个凝胶微珠上有40-80万特定的核酸引物序列，该序列由以下几部分构成：1）Read 1 测序引物；2）10X Barcode序列：16碱基，一个Gel bead对应一种10X Barcode，共有~350W种10X Barcode，用于区分细胞。任意两个barcode之间至少差2个或2个以上的碱基（避免误读）；3）UMI（unique molecular identifier）：12碱基，随机序列，作用是在经过PCR扩增再深度测序得到的Reads，可以看出哪些reads是来自于一个原始的cDNA分子，用于区分同一细胞的不同转录本。可以排除各种cDNA因为PCR扩增效率的不同而导致的reads数的偏差（PCR bias）；4）poly dT反转录引物：30nt，作用是与mRNA的Poly(A)尾巴结合，作为逆转录的引物，逆转录出cDNA来。

![https://upload-images.jianshu.io/upload_images/15771939-a672912d0d8471b0.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-a672912d0d8471b0.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

一个凝胶微珠上的40-80万核酸引物序列的barcode是一样的，但UMI不同

### **芯片上的液流管路**

![https://upload-images.jianshu.io/upload_images/15771939-89ebe57e89f29bca.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-89ebe57e89f29bca.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

![https://upload-images.jianshu.io/upload_images/15771939-d3094dfe2a5df0c2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-d3094dfe2a5df0c2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

经过这个系统，制备出油包水小液滴的乳浊液。这些小液滴里面是水相，外面包裹的是油相。

![https://upload-images.jianshu.io/upload_images/15771939-db1629bd0ccaf5f4.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-db1629bd0ccaf5f4.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

细胞混悬液中约65%的细胞会被包到有微珠的小液滴当中。这些液滴中包含细胞的数目是符合泊松分布的，大部分细胞会被单独包裹在一个小液滴中。

### **测序文库制备**

1. 在得到乳浊液之后，将细胞膜破掉，让细胞当中的mRNA游离出来。游离出来的mRNA与小液滴中的水相混合，也就是和逆转录酶、结合在凝胶微珠上的核酸引物、以及dNTP底物相接触。

![https://upload-images.jianshu.io/upload_images/15771939-0e740522a5126a49.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-0e740522a5126a49.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

2. 接着发生逆转录反应，mRNA与凝胶微珠上带标签的DNA分子相结合，在逆转录酶的作用下，逆转录出cDNA `第一链`（下图紫色序列）。

![https://upload-images.jianshu.io/upload_images/15771939-3604855f85cce9d0.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-3604855f85cce9d0.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这样得到的cDNA分子第一链是带有这个微珠所特有的Barcode标签和各自特定的UMI标签的。有了这两个标签，cDNA分子就可以互相区分开来。

3. 以SMART方式使用TSO引物完成第二链合成

![https://upload-images.jianshu.io/upload_images/15771939-bcde41bb3cbb5d97.png?imageMogr2/auto-orient/strip|imageView2/2/w/1044/format/webp](https://upload-images.jianshu.io/upload_images/15771939-bcde41bb3cbb5d97.png?imageMogr2/auto-orient/strip|imageView2/2/w/1044/format/webp)

4. 油滴破碎，磁珠纯化cDNA一链，然后PCR扩增cDNA。

![https://upload-images.jianshu.io/upload_images/15771939-8538decd246547d1.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-8538decd246547d1.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

5. cDNA扩增完成后酶切片段化并磁珠筛选最适片段，通过末端修复、加A、接头连接Read2测序引物，再以PCR方式构建含有P5和P7接头的cDNA文库即可。

![https://upload-images.jianshu.io/upload_images/15771939-4ee7ce34d56bdbbf.png?imageMogr2/auto-orient/strip|imageView2/2/w/1148/format/webp](https://upload-images.jianshu.io/upload_images/15771939-4ee7ce34d56bdbbf.png?imageMogr2/auto-orient/strip|imageView2/2/w/1148/format/webp)

随后就可以进行测序和数据分析了。

> 10X genomics优点：
> 
> - 简单便捷：集单细胞分选、扩增、建库于一体。
> - 细胞通量高：每个样本细胞数可以达到5000-10000个。
> - 建库周期短：1天可以完成单细胞悬液制备、单细胞捕获、扩增及建库。
> - 捕获效率高：单个液滴捕获效率高达65%
> - 真正意义的单细胞：单个液滴捕获到多细胞概率极低（0.9% / 1000cells）
> 
> **10X genomics缺点：**
> 
> - 非全长信息：只能获得3‘端转录本信息
> - 样本要求高：单个样本细胞起适量达5x10^4 - 5x10^5 个，活细胞数目需要超过80%，最好在90%以上。

From chatGPT:

In 10X scRNA sequencing, only the 3' end of the RNA molecule is sequenced because the technique uses a method called "5' end tagging" or "oligo-dT priming" to selectively capture and amplify the mRNA molecules.

During library preparation, a unique barcode and a poly-dT oligonucleotide with a 10X-specific sequence are added to the 3' end of the mRNA molecules. This poly-dT tail binds to the poly-A tail found at the 3' end of most mRNA molecules, allowing for selective capture and amplification of the mRNA molecules that are polyadenylated.

When the cDNA is synthesized from the mRNA using reverse transcription, the oligo-dT primers anneal to the poly-A tail and extend the cDNA in the 3' direction. This results in the cDNA fragments that only contain the 3' end of the mRNA molecule. These cDNA fragments are then sequenced using next-generation sequencing (NGS) technologies, such as Illumina.

By sequencing only the 3' end of the mRNA molecules, 10X scRNA sequencing allows for high-throughput sequencing of single cells, while still providing valuable information on gene expression levels and cell type identification. Additionally, sequencing only the 3' end reduces the sequencing cost and increases the throughput of the sequencing run.

The **`--r1-length`** option in the Cell Ranger Count analysis specifies the length of Read 1 in the sequencing data. By default, Cell Ranger assumes that Read 1 contains the 16-base 10x barcode and the 10-base Unique Molecular Identifier (UMI), so the default value for **`--r1-length`** is 26.

The 10x barcode is used to identify which cell the mRNA transcript originated from, while the UMI is used to remove PCR duplicates that may arise during library preparation and sequencing. The remaining bases in Read 1 correspond to the cDNA fragment that was synthesized from the 3' end of the mRNA molecule.

The **`--r1-length`** option allows users to specify a different length for Read 1 if their sequencing data has a different configuration or if they have trimmed their reads for quality control purposes. However, if the length of Read 1 is shorter than the expected length of the barcode and UMI, the barcodes and UMIs may not be correctly identified, leading to a reduction in the accuracy of cell and molecular counting. On the other hand, if the length of Read 1 is longer than the expected length, it may increase the sequencing cost without improving the quality of the results.

Therefore, the default value of 26 for **`--r1-length`** in Cell Ranger Count is based on the typical 10x sequencing configuration and is optimized for efficient and accurate cell and molecular counting.

the length of the UMI in 10x Genomics sequencing is typically 12 bases, not 10. Therefore, the expected length of Read 1 in 10x Genomics sequencing data is 28 bases (16 for the barcode + 12 for the UMI).

In the Cell Ranger Count analysis, the **`--r1-length`** option specifies the length of Read 1 in the sequencing data. Therefore, in order to account for the expected length of the barcode and UMI, **`--r1-length`** should be set to 28 for 10x Genomics sequencing data with a 12-base UMI.

It's worth noting that the default value of 26 for **`--r1-length`** in Cell Ranger Count is actually used for older versions of 10x Genomics libraries, which used a 10-base UMI. However, for newer 10x Genomics libraries with a 12-base UMI, **`--r1-length`** should be set to 28.

# 2**. SMART (Clontech)**

## 2**.1 简介**

SMART方法的全称是 Switching Mechanism at 5' End of RNA Template，该方法发表于于2012年[1]，2013年发表了其改进技术的应用Smart-Seq2 [2]，2014年Smart-Seq2 protocol发表[3]。Smart-Seq2对原始的Smart-Seq实验流程进行了多项改进优化，它不再需要纯化步骤，可大大提高产量，最重要的改进是下面两项：（1）TSO 3'端最后一个鸟苷酸替换为锁核酸LNA(locked nucleic acid)。LNA单体的热稳定性增强，其退火温度增强非模板cDNA的3'延伸能力。（2）甜菜碱（一种具有两个重要作用的甲基供体：它会增加蛋白质的热稳定性，并通过破坏DNA螺旋来降低甚至消除了DNA热融变对碱基对组成的依赖性）与较高的MgCl2浓度结合使用。解决某些RNA形成二级结构（例如发夹或环）由于空间位阻，可能导致酶终止链延长的问题。

Smart技术是基于高保真的反转录酶、模板转换和前置放大来增加cDNA得率，实验流程2天，得到的是**全长转录本**。该方法有较好的覆盖范围，可检测到稀有转录本，因此应用范围较广。

## 2**.2 建库原理**

**建库流程图总览**

![https://upload-images.jianshu.io/upload_images/15771939-92cb209d31a5d9a9.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-92cb209d31a5d9a9.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

**Smart-Seq2建库原理**

1. 单细胞分选：使用流式细胞仪或显微操作进行细胞分选，体积不超过0.5 ul。
2. 细胞裂解：将分离细胞直接转移到细胞裂解液中进行细胞裂解。
3. 反转录（ 一链合成 ）：使用Oligo(dT) primer 对带有polyA尾的RNA( 主要mRNA )进行反转录。

![https://upload-images.jianshu.io/upload_images/15771939-088ad8bec93e9e9b.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-088ad8bec93e9e9b.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这个逆转录的起始引物先是一段通用序列，会用作PCR扩增的引物识别序列。中间是一长串的T，这些T专门识别mRNA的3‘末端的Poly(A)尾巴序列，与Poly(A)尾互补结合。引物最末端有一个定位结构，在3‘末端的倒数第二个碱基是一个非T的简并碱基（V表示A/C/G）。最后一个碱基则是简并碱基N（A/C/T/G都有可能）。引物的这个末端结构，就是让它正好结合在mRNA的3‘端连到Poly(A)尾巴的连接处，而不会结合到mRNA别的地方。这样就保证了逆转录的起始位置正好是mRNA的3'端的序列终止位置。

![https://upload-images.jianshu.io/upload_images/15771939-17f64b8a3e272550.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-17f64b8a3e272550.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

由于使用了特殊活性反转录酶（Moloney Murine Leukemia Virus , `MMLV`,莫洛尼鼠白血病病毒反转录酶）进行反转录，所以在它转录到mRNA的5'末端的时候，**会在cDNA链3'端加上几个不依赖于模版的C碱基**。

1. 模板置换（ 二链合成 ）：该步使用`TSO`（template-switching oligo, 特异性模板转换引物）引物合成了cDNA的二链，从而置换了与一链cDNA互补的RNA。

![https://upload-images.jianshu.io/upload_images/15771939-7b48a921e09897a3.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-7b48a921e09897a3.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

要注意的是TSO引物的 3'端有三个非脱氧的G碱基（RNA的G碱基），能与一链3'端MMLV多合成的几个C碱基互补，而最末端的+G是一个修饰过的G，能增加TSO的热稳定性，以及其与一链cDNA游离的3’端的互补的能力。互补杂交之后，可以引导MMLV酶再次发挥聚合作用，以刚才那条新合成的cDNA为模版来复制得到双链cDNA。

![https://upload-images.jianshu.io/upload_images/15771939-c80814ae868915bd.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-c80814ae868915bd.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这个双链cDNA两端都已经接好了我们人工设计的PCR引物序列（红圈），然后就加入常规PCR引物，进行常规PCR扩增。

1. PCR扩增：该步进行轻度的cDNA富集，将cDNA扩增至ng级即可。
2. 标记：利用改造后的高活性Tn5转座酶对DNA进行打断的同时将接头添加到cDNA的两端。标记完成后的DNA片段通常在200-600bp。
3. PCR富集及上机测序：在进行最后一次PCR扩增后，即可上机测序。

> 3个巧妙点：1. 先用一个定位引物，保证cDNA的合成是从mRNA的3'最末端开始的。同时让合成的cDNA在下游连上了一个通用PCR序列。2. 利用MMLV逆转录酶在新合成cDNA的3'端多加几个C碱基的特点，再用有3个G碱基的上游引物进行第二链的合成。这也就保证了只有完整的cDNA也就是那些带多个 C的cDNA（第一链）才能合成出cDNA第二链。这就保证了双链cDNA是全长的cDNA。3. 保证了PCR扩增效率的一致性。PCR扩增效率的最主要的影响因素是引物的序列，现在因为cDNA的5‘端和3’端都分别引入了统一的引物序列，就去除了因为引物序列的不同而引起PCR效率不同这个最主要的偏差因素。也就在较大程度上保证了PCR扩增效率的一致性，减少了PCR偏差。
> 

经过实验发现，用SMART方法，对一个细胞也就是10pg总RNA进行建库测序，RPKM为10的这些基因，有60%是被测序测到的。对RPKM为100的这些基因，有90%是可以被测序测到的，而且被测到的几率波动很小。

![https://upload-images.jianshu.io/upload_images/15771939-e83e5985b73c76db.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-e83e5985b73c76db.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这说明SMART方法是一个有效的单细胞mRNA测序的方法。

> Smart-Seq2优点：
> 
> - 相比于截短的cDNAs，MMLV逆转录酶更倾向于选择全长cDNAs作为其末端转移酶活性的底物。因此每个转录本的所有外显子都能被检测到。这使它可以用于检测`可变剪切`，还能在转录层面进行全面的`SNP`和`突变分析`，扩大了其应用范围。
> - 不同I5、I7 Index组合使其能够进行多样本混合测序。
> - 方案和原理公开，让研究者可以进一步对其进行改良。目前在这个方案的基础上涌现了很多单细胞测序的新成果。
> - 和10x Genomics相比，单细胞检测到的转录本更多。
> 
> **Smart-Seq2缺点：**
> 
> - 由于对聚腺苷酸化的RNA具有选择性，所以不能分析非poly(A)的RNA。
> - 测序reads不带有mRNA链特异性。

# 3**. TargetAmp ()**

## 3**.1 简介**

TargetAmp方法由Illumina公司旗下的EpiCentre公司开发，可以把少量RNA（最少到1个细胞，约10pg）扩增到ng级，以达到可以进行高通量测序所需的核酸量。

## 3**.2 建库原理**

**建库流程图总览**

![https://upload-images.jianshu.io/upload_images/15771939-0ea706be77c368ab.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-0ea706be77c368ab.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

**建库原理**

1. 第一链cDNA合成：用T7-Oligo(dT)的引物进行cDNA合成

![https://upload-images.jianshu.io/upload_images/15771939-d034c544ef81a87e.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-d034c544ef81a87e.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这个引物在5'端设计了一个T7启动子序列，3'端是多个T碱基，可以与mRNA的poly(A)尾巴相结合，作为逆转录的起始引物。

![https://upload-images.jianshu.io/upload_images/15771939-615143d329f07d3e.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-615143d329f07d3e.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

逆转录得到的第一链cDNA被引入了一个T7启动子。

1. cDNA第二链的合成：使用RNase H酶特异性降解RNA和DNA杂交链中的RNA链，剩下带有T7启动子的cDNA单链，再合成出第二条cDNA链来。

![https://upload-images.jianshu.io/upload_images/15771939-ee5f20c588e39b81.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-ee5f20c588e39b81.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

1. 得到的双链cDNA可以作为转录的模板，利用链上的T7启动子，启动体外转录生成大量反义aRNA（antisense-RNA）

![https://upload-images.jianshu.io/upload_images/15771939-87f7e710a174da5b.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-87f7e710a174da5b.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

1. 纯化aRNA
2. 第二轮的第一链cDNA合成：使用随机引物进行逆转录，得到第二轮cDNA

![https://upload-images.jianshu.io/upload_images/15771939-9df08a21fc1d3e7d.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-9df08a21fc1d3e7d.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

1. 第二轮中的第二链cDNA合成：用RNase H把DNARNA杂交产物中的RNA消化掉。用T7-Oligo(dT)引物粘到第二轮cDNA的poly(A)尾巴上，合成出cDNA双链。

![https://upload-images.jianshu.io/upload_images/15771939-37d890695d730059.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-37d890695d730059.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

1. 这个双链cDNA再经过第二轮的转录，又得到第二轮的反义RNA。

![https://upload-images.jianshu.io/upload_images/15771939-06ce3a73abc1a6d7.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-06ce3a73abc1a6d7.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

这些第二轮的反义RNA的量，足以达到微克级。再经过一轮逆转录，就可以得到几个微克的cDNA，足以进行建库测序之用。

> 本方法的巧妙之处：1. 它不是用PCR来扩增核酸，而是用转录的方法来增加核酸的量。因为扩增那么多倍的核酸，如果用PCR，需要几十个循环。那么PCR不同的扩增子的扩增效率，即使一开始是很小的差异，也会在几十个循环中被指数放大，变成一个很大的差异。TargetAmp用转录的方法，统一都用T7这个启动子，它转录的起始效率大体上就保持了一致。它的每一轮转录，都把核酸的量扩大了几千倍，经过两轮的扩增，就把核酸的量扩大了百万倍。这样一方面得到了足以用来建库的高达几微克的核酸，另一方面又避免了PCR过程，也就避免了PCR扩增偏差。2. 第一轮与第二轮都是线性扩增，大大减少了PCR反应的指数效应所引起的Bias。3. 高效扩增，一轮扩增可以扩增几千倍，把10pg级的Total RNA中的mRNA扩增到几个ng，达到二代测序的样本量要求。如果经过两轮扩增，就可以达到生物芯片所需的ug级的核酸量。
> 

# **4. 多种方法的比较**

10X和SMART的比较，出自张泽民团队[5]

![https://upload-images.jianshu.io/upload_images/15771939-43ac6a1ec7e4e089.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-43ac6a1ec7e4e089.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

2017年Molecular Cell文章，对6种单细胞转录组技术的比较[6]

![https://upload-images.jianshu.io/upload_images/15771939-e95af7ccc7eb00ae.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-e95af7ccc7eb00ae.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

2019年Nature Communications文章，对7种单细胞RNA测序方法进行比较[7]

![https://upload-images.jianshu.io/upload_images/15771939-3f694df16245dcd2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-3f694df16245dcd2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

![https://upload-images.jianshu.io/upload_images/15771939-91e02ea6cea5ff6f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-91e02ea6cea5ff6f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

7种单细胞RNA测序方法比较

**参考文献：**

1. Ramskold, D., et al., Full-length mRNA-Seq from single-cell levels of RNA and individual circulating tumor cells. Nat Biotechnol, 2012. 30(8): p. 777-82.
2. Picelli, S., et al., Smart-seq2 for sensitive full-length transcriptome profiling in single cells. Nat Methods, 2013. 10(11): p. 1096-8.
3. Picelli, S., et al., Full-length RNA-seq from single cells using Smart-seq2. Nat Protoc, 2014. 9(1): p. 171-81.
4. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets. Cell. 2015;161(5):1202-1214.
5. Direct Comparative Analyses of 10X Genomics Chromium and Smart-seq2. Genomics Proteomics Bioinformatics. 2021 Mar 1:S1672-0229(21)00048-6.
6. Comparative Analysis of Single-Cell RNA Sequencing Methods. Mol Cell. 2017;65(4):631-643.e4.
7. A systematic evaluation of single cell RNA-seq analysis pipelines. Nat Commun. 2019;10(1):4667.