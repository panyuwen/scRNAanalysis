# 3. 单细胞免疫组库：TCR基因重排原理和TCR测序建库方法

免疫组库（immune repertoire, IR）是指某个体在特定时间点其循环系统中所有功能多样性B细胞和T细胞的总和，在自身免疫疾病、传染病、癌症等多种疾病中扮演者至关重要的角色，也正因此免疫组库一直是各领域的研究重点。

# **TCR的结构**

![https://upload-images.jianshu.io/upload_images/15771939-496acb25dd14a6e2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-496acb25dd14a6e2.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

免疫是人类和脊椎动物最重要的防御机制，在暴露于外来抗原后产生特异性抗体。T细胞介导的抗原识别取决于T细胞受体（TCR）与抗原主要组织相容性复合体（MHC）分子的相互作用。TCR是高度多样化的异源二聚体，由大多数T细胞表达的α和β链（αβ TCR）的组成，或由外周血中的T细胞（1-5%）和在粘膜部位发现的T细胞表达的γδ链（γδ TCR）组成。每个T细胞表面约有3000～30000个TCR分子。

所有T细胞的TCR总和成为TCR profiling（TCR谱）。随着疾病不断恶化，TCR profiling会发生很大的变化。于是研究人员越来越关注在不同的疾病条件下免疫谱的状态，如癌症、自身免疫、炎症、传染病等。

![https://upload-images.jianshu.io/upload_images/15771939-a2dfca4726afa48f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp](https://upload-images.jianshu.io/upload_images/15771939-a2dfca4726afa48f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp)

TCR链由抗原识别的可变区域(V区)和恒定区域(C区)组成。编码人的TCR α和δ链的基因分别定位于第14号和第7号染色体。α链由多个可变的V、连接J基因和C基因(编码恒定C区)编码，而TCR β链由V、J、C和多样性D基因编码。V、D、J、C等基因座位又各自有不同的等位基因。

![https://upload-images.jianshu.io/upload_images/15771939-9e945ac045ae4a97.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp](https://upload-images.jianshu.io/upload_images/15771939-9e945ac045ae4a97.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp)

# **TCR基因重排**

TCR基因重排又称DNA重排或体细胞重组，指在T细胞在胸腺中分化成熟过程中，胚系状态的V区基因由分隔的、无转录活性的基因片段在特异性重组酶的作用下连接成一个完整的、有转录功能的活性基因的过程。重排时，Vβ基因先进行D、J连接(某―D片段与某―J片段相连)，再进行V、DJ连接(相连的DJ片段与某―V片段相连)形成`V|DJ`片段。由此产生一个有转录活性的Vβ基因(V|DJ基因)，后者再与C区基因相连，形成一个完整的β链功能基因。Vβ基因的成功重排，可诱导Vα基因的重排。Vα基因无D片段，直接进行V、J连接(某―V片段和某一J片段相连)形成`VJ`片段。重排后有转录活性的Vα基因(VJ基因)再与C区基因相连，形成一个完整的α链功能基因(TCR的特异性分别由α链和β链的V-J及V-D-J片段决定)。

![https://upload-images.jianshu.io/upload_images/15771939-7c31a5f3f7d6f529.png?imageMogr2/auto-orient/strip|imageView2/2/w/1128/format/webp](https://upload-images.jianshu.io/upload_images/15771939-7c31a5f3f7d6f529.png?imageMogr2/auto-orient/strip|imageView2/2/w/1128/format/webp)

TCR的V区基因重排只发生在T细胞分化的早期，因为特异性重组酶只存在早期，其活性具有严格的时限性和兰织细胞特异性。故T细胞在分化成熟过程中只能进行一次有效的基因重排，保证了一个T细胞克隆只能表达一种特异性TCR，只显示一种抗原识别特异性。

两条链基因重排后可形成千万种不同特异性的TCR分子，可识别环境中多种多样的抗原。正是这种随机重排导致的TCR多样化使免疫系统在识别和杀伤特异外源性抗原时更为快速、准确和有效。

过于多样化也是目前免疫组学测序最大的问题。理论上VDJ重排后可产生10^15~10^20种不同的克隆类型，但实际上只有约10^13种不同的克隆类型。这就意味着看似随机的VDJ重排其实是有规律的，并且受到各种条件的限制。

# **CDR区**

用恒定基因区段重组可变区产生功能性TCR链转录物。这个过程导致强组合（取决于哪个基因区域将重组）和连接多样性（将添加/删除多少核苷酸），从而产生大量且高度可变的TCR profiling，最终鉴定到大量的TCR profiling抗原。通过配对α和β或γ和δ链以形成功能性TCR来实现多样性。

每个TCR链包含三个高可变环区（hypervariable loops）,称之为互补决定区CDR1-3（complementarity determining region 1-3）。CDR1和CDR2由V基因编码，并且对TCR与MHC复合物互作至关重要。然而CDR3由V和J或D和J之间连接区编码，因此CDR3变化程度较大。由于**CDR3是与抗原直接接触的TCR区域**，因此CDR3在TCR与肽-MHC复合物的相互作用中起到了十分重要的作用。

![https://upload-images.jianshu.io/upload_images/15771939-d630d4f4811197b6.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-d630d4f4811197b6.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

Vβ基因重排

所以CDR3是T细胞常见的克隆区域，除非T细胞来源于相同的形式扩增出来的，一般来说**T细胞几乎不太可能表达相同的CDR3序列**。

# **文库构建方法**

[Chromium Next GEM Single Cell 5' Kit v2](https://links.jianshu.com/go?to=https%3A%2F%2Fassets.ctfassets.net%2Fan68im79xiti%2F4oB71TeT0kDoIHhfq9dPxd%2F05ce9121d027715321d2a9765b1e9b70%2FCG000331_ChromiumNextGEMSingleCell5_v2_UserGuide_RevA.pdf)

10X V(D)J测序与其转录组测序建库流程基本相同，都要经历微流控构建单细胞油包水反应体系，mRNA逆转录为cDNA并扩增的过程。与常规10X单细胞转录组测序把barcode和UMI序列放在3'端不同，V(D)J测序要把barcode和UMI序列放在5'端。样本逆转录的cDNA经过扩增后可以一分为二，一部分做5'端scRNA测序，另一部分用富集引物扩增V(D)J序列测序。

### **Step 1: GEM Generation & Barcoding**

![https://upload-images.jianshu.io/upload_images/15771939-a5426ef838fa4b78.png?imageMogr2/auto-orient/strip|imageView2/2/w/836/format/webp](https://upload-images.jianshu.io/upload_images/15771939-a5426ef838fa4b78.png?imageMogr2/auto-orient/strip|imageView2/2/w/836/format/webp)

![https://upload-images.jianshu.io/upload_images/15771939-b128ac36471c3646.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-b128ac36471c3646.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

### **Step 2: Post GEM-RT Cleanup & cDNA Amplification**

在GEM-RT reaction之后，使用磁珠富集纯化带有10X barcode的第一链cDNA，然后PCR扩增cDNA。这一步扩增的产物可以同时用来构建T细胞、B细胞建库（Step 3- 4）和5'基因表达文库（Step 5）

![https://upload-images.jianshu.io/upload_images/15771939-3ac294d3ba500d08.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-3ac294d3ba500d08.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

### **Step 3 & 4: V(D)J Amplification from cDNA and V(D)J Library Construction**

![https://upload-images.jianshu.io/upload_images/15771939-385325be11721efd.png?imageMogr2/auto-orient/strip|imageView2/2/w/892/format/webp](https://upload-images.jianshu.io/upload_images/15771939-385325be11721efd.png?imageMogr2/auto-orient/strip|imageView2/2/w/892/format/webp)

扩增得到的全长V(D)J序列会被酶切为长短不一的片段。因为每条V(D)J序列都有很多个PCR扩增的拷贝，并且这些拷贝有相同的UMI标签，所以每条V(D)J序列都会得到一组测序reads（通过UMI确定同一来源）。

### **Step 5: 5ʹ Gene Expression (GEX) Library Construction**

![https://upload-images.jianshu.io/upload_images/15771939-f5150288a9b15c9e.png?imageMogr2/auto-orient/strip|imageView2/2/w/868/format/webp](https://upload-images.jianshu.io/upload_images/15771939-f5150288a9b15c9e.png?imageMogr2/auto-orient/strip|imageView2/2/w/868/format/webp)

### **Step 6: Sequencing**

![https://upload-images.jianshu.io/upload_images/15771939-7a88eddaefa4d00c.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-7a88eddaefa4d00c.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

# **10X V(D)J和其他方法相比的优势**

- 多重PCR法：扩增范围局限在CDR3区
- 5'RACE法：只能获得单链全长信息
- 10X V(D)J1）可实现成对的重链和轻链（B细胞）或α和β链（T细胞）的全长测序2）精细到单细胞水平，可同时获得大量单细胞的免疫组库数据3）同批样本可同时测得mRNA表达谱数据4）灵活的获取量，可在7分钟内封装100-80000个单细胞，进而建库测序。