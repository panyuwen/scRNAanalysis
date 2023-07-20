# 2. DNA文库构建和Illumina测序化学原理

# **1. DNA文库构建**

所谓DNA文库，实际上是许多个DNA片段，在两端接上了特定的DNA接头，形成的DNA混合物。文库有2个特点：1. 当中这一段插入的DNA，它的序列是各种各样的。2. 它的两头的接头序列，是人工特异加上去的，是已知的。要构建文库，首先需要把基因组DNA用超声波打断，之后把两端用酶补平。再用Klenow酶在3’端加上一个A碱基，然后再用连接酶把接头给连上去。连好了接头的DNA混合物，我们就称为一个文库。

![https://upload-images.jianshu.io/upload_images/15771939-4a8793719d13b468.png?imageMogr2/auto-orient/strip|imageView2/2/w/1176/format/webp](https://upload-images.jianshu.io/upload_images/15771939-4a8793719d13b468.png?imageMogr2/auto-orient/strip|imageView2/2/w/1176/format/webp)

# **2. Illumina测序化学原理**

Illumina仪器对比。从最早的Miseq一天测三千万条read，到Hiseq一天测30亿条reads，再到Novaseq一天可以测130亿条read，通量还是有一个非常大的提升。

![https://upload-images.jianshu.io/upload_images/15771939-9abfabb134103fc5.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-9abfabb134103fc5.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

![https://upload-images.jianshu.io/upload_images/15771939-445db40532d52882.png?imageMogr2/auto-orient/strip|imageView2/2/w/880/format/webp](https://upload-images.jianshu.io/upload_images/15771939-445db40532d52882.png?imageMogr2/auto-orient/strip|imageView2/2/w/880/format/webp)

不同仪器芯片比较

### **2.1 桥式PCR**

文库构建好之后，后续就是做桥式PCR。桥式PCR是把文库种到芯片上去然后进行扩增的这样一个过程。1）首先要把文库加到芯片上。芯片的内表面种满两种不同类型的oligo(寡核苷酸序列)。因为文库两头的DNA序列和芯片上的引物是互补的，就可以发生互补杂交。2）随后加入dNTP和聚合酶，聚合酶会从引物开始，延着模版合成出一条全新的DNA链来。新的这条链和原来的链是完全互补的。3）接下来加入NaOH碱溶液，DNA在NaOH碱溶液存在的情况下，就解链了。液流一冲，原来的模版链（没有和芯片共价连接的链）就会被冲走，和芯片共价连接的 链就会被保留。4） 再往液流池中加入中性液体（中和前面加入的碱液），这时DNA链上的另外一端就会和玻璃板上的第二种引物发生互补杂交。5）加入酶和dNTP，聚合酶就沿着第二个引物合成出一条新的链来。6）然后再加碱，把两条链解链开，再加入新的中和液，这时候DNA链就会和新的引物杂交。再加酶，再加dNTP，又从新的引物上合成出新的链来。连续重复这一过程，DNA链的数量就会以指数方式增长。

![https://upload-images.jianshu.io/upload_images/15771939-9ea21e26572f6e84.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-9ea21e26572f6e84.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

### **2.2 双链变单链**

桥式PCR完成之后，接下来需要把合成的双链变成可以测序的单链。办法是通过一个化学反应，把一个引物上的一个特定基团给切断掉，然后再用碱溶液来洗芯片。碱让DNA双链解链，那根被切断了根的DNA链就被水冲掉了，留下那根共价键连在芯片上的链。接下来加入中性溶液，再在这个中性溶液里加入测序引物，随后就可以开始正式的测序工作了。

![https://upload-images.jianshu.io/upload_images/15771939-44a625a59903f04f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-44a625a59903f04f.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

### **2.3 测序**

在测序的时候加进去的主要是两个东西，一是带荧光标记的dNTP（3‘末端是被一个叠氮基堵住的），二是聚合酶。聚合酶就会选择哪个dNTP是和原来位置上的那个碱基是互补的，根据互补性原理，把这个dNTP合成到新的链上去。

![https://upload-images.jianshu.io/upload_images/15771939-0a9fec5bbc01fad8.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp](https://upload-images.jianshu.io/upload_images/15771939-0a9fec5bbc01fad8.png?imageMogr2/auto-orient/strip|imageView2/2/w/1042/format/webp)

因为dNTP的3‘端是被一个叠氮基团堵住的，所以它一个循环只能延长一个碱基。合成之后，用水把多余的dNTP和酶给冲掉，放到显微镜下去进行激光扫描，根据发出来的荧光判断它是哪个碱基。因为4种dNTP上面标记的荧光素都不一样，根据荧光就可以判断新合成的碱基是什么碱基。因为新合成的碱基和原来位置的碱基是互补的，就可以知道模版链的碱基是什么。

![https://upload-images.jianshu.io/upload_images/15771939-6cd3d1e7dd1050b8.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp](https://upload-images.jianshu.io/upload_images/15771939-6cd3d1e7dd1050b8.png?imageMogr2/auto-orient/strip|imageView2/2/w/1200/format/webp)

一个循环完成之后，就加入一些化学试剂，把叠氮基团和旁边标记的荧光基团切掉。切掉之后，3‘端的羟基就暴露出来，接下来加入新的dNTP和新的酶，就又延长一个碱基。之后把多余的酶和dNTP冲掉，再进行一轮显微的激光扫描，判断碱基是什么。重复这个过程，就可以把上百个甚至更多个碱基的序列读出来。

![https://upload-images.jianshu.io/upload_images/15771939-3d60873219cd74af.png?imageMogr2/auto-orient/strip|imageView2/2/w/744/format/webp](https://upload-images.jianshu.io/upload_images/15771939-3d60873219cd74af.png?imageMogr2/auto-orient/strip|imageView2/2/w/744/format/webp)

### **2.4 读取Index (Barcode)**

因为illumina的测序量很大，但一个样本往往用不了几亿条DNA。所以科学家就想了一个办法，在文库的接头上做了一些标记，每一个样本有一个特定的接头，每个接头里面有一段特定的序列，这段特定的序列，我们就称为Index/Barcode（特定序列标记了特定样本的来源）。

要读Index序列，先用碱把上面这跟测完‘Read 1’的序列上面的DNA链解链掉，加入中性液，再加入‘Read 2’的测序引物。Read 2的结合位点就在Index序列的旁边，接下来进行第二轮测序。一般是读6-8个碱基。读完以后就可以知道这某一个具体的一段DNA，它来自原始的哪个样本。

![https://upload-images.jianshu.io/upload_images/15771939-6360d8b2d713d4db.png?imageMogr2/auto-orient/strip|imageView2/2/w/556/format/webp](https://upload-images.jianshu.io/upload_images/15771939-6360d8b2d713d4db.png?imageMogr2/auto-orient/strip|imageView2/2/w/556/format/webp)

### **2.5 双端测序**

这是Illumina的最核心的另外一个技术 。双端测序就是一根DNA链，除了从正向读一遍，还可以从DNA的负向再读一遍。这样子就把Illumina测序的有效长度加了一倍。这个倒链的过程，是先让DNA合成，得到互补链。之后用化学试剂切断模版链根部，加入碱溶液洗掉，接下来就进行第2端的测序。原理和第1端是一样的。

最重要的是，我们可以理解，一个点经过几百个循环得到一条链几百个碱基的信息。但实际上这个芯片可以有上亿个点，也就是上亿个cluster（簇）。上亿个链同时在合成，因此每一个循环都可以读出上亿个序列，这就得到了很大的一个测序数据量。

### **2.6 Limits**

![https://upload-images.jianshu.io/upload_images/15771939-d7c423881329633a.png?imageMogr2/auto-orient/strip|imageView2/2/w/640/format/webp](https://upload-images.jianshu.io/upload_images/15771939-d7c423881329633a.png?imageMogr2/auto-orient/strip|imageView2/2/w/640/format/webp)

reads越长，出错的越多，真实信号会越来越弱。因此illumiina的测序读长被限制在300bp以内。