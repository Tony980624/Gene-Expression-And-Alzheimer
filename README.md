HBTRC公布了一组数据，包含230个去世的老人，差不多一半是带着阿兹海默症（老年痴呆）去世的。数据记录了每个对象的年龄，性别，是否患病，以及大脑基因表达序列。接下来我会尝试用PCA分析，单变量回归，Boosting模型，Lasso Logostic 回归，以及Graphical Lasso来分析 究竟在患者和正常人之间，存在哪些基因的表达差异？
数据获取地址: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE44772  (HBTRC)

# PCA

PCA的第一步是标准化数据，一般来说对每个变量就是减去均值，除以标准差:

$$
Z = \frac{X-\bar{X}}{sd(X)}
$$

第二步就是计算协方差矩阵Covariance matrix

$$
Cov = \frac{1}{n-1}Z^TZ
$$

第三步就是把协方差矩阵特征分解

$$
Cov*v = \lambda v
$$

算出来的 v 代表PC（主成分）的方向,代表各个变量的线性组合。 $\lambda$ 指示了主成分的重要性（方差解释能力），帮助我们选择主成分。我们选择  $\lambda$ 大的PC,因为它解释了最多的Variance,也就是数据的差异性。也就意味着含有更多信息。

这是对基因数据的PCA,画出了前30个PC解释的方差比例：

![pca](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/000010.png)

然后画出PC的热力图，颜色深的行（变量），代表对于该PC权重更大的变量。

![pca2](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/000011.png)

# 单变量回归

线性回归：拟合一条残差最小的直线

$$
min\sum_{i=1}(y-\beta X_i-\epsilon)^2  
$$

$$
\frac{\partial f}{\partial x}||y-X\beta||^2=0
$$

化解后求导得:

$$
\beta = (X^TX)^{-1}X^Ty
$$

数据包含了三千多个基因序列，这就意味着要进行三千多次回归和t-test。假设我们进行了 100 次 t 检验，每个检验的显著性水平设为 0.05。即使所有假设都是真阴性（即没有效应），我们仍然期望有大约 100×0.05=5100×0.05=5 次假阳性结果。随着检验次数的增加，误判的数量也会增加。这时候我们就需要用到FDR了，Benjamini-Hochberg (BH) 方法是控制错误发现率（False Discovery Rate, FDR）的经典方法之一。FDR 是多重假设检验中“错误发现”的比率，表示在所有被判定为显著的结果中，有多少比例是错误的（即假阳性）。BH 方法通过对 p 值进行排序和阈值调整，来控制 FDR 的上限。

1：将所有p值从小到大排列

2：为排序后得p值计算一个阈值 (i/m)*$\alpha$,其中i是p的排序位置，m是检验总数，$\alpha$  是想要控制的FDR水平，比如0.05

3：找到最大的i使得$p_i <= \frac{i}{m} \alpha$, 这个p值就是我们判断显著性的阈值，保留所有小于阈值的p，拒绝他们的零假设。

