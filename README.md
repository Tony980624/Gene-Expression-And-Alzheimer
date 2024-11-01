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

结果如下：

![pca](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/000014.png)

最后在剩下的显著变量里记录下系数绝对值大的基因

## Lasso Logistic回归

Lasoo其实是一种L1正则法，L1 正则化是指在模型的损失函数中添加一个基于系数绝对值的惩罚项。之所以在这里用到Lasso，是因为变量有几千个，而样本只有230个，这种情况下矩阵不是满序的，不可逆。所以我们要加入惩罚项，把不那么重要的系数逼为0.

损失方程：

$$
-\sum_{i=1}^n (y_ilog(p_i)+(1-y_i)log(1-p_i))
$$

优化目标函数：

$$
-\sum_{i=1}^n (y_ilog(p_i)+(1-y_i)log(1-p_i))+\lambda\sum^p_{j=1}|\beta_j|
$$

![pca2](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/000012.png)


最终模型的交叉验证错误率为10%左右，训练集错误率为4%，不算严重过拟合。

## Boosting 模型

Boosting 是一种提升模型性能的集成学习方法，通过将一系列弱学习器（通常是准确率略高于随机猜测的模型）组合成一个强学习器来提高预测效果。Boosting 的核心思想是在每一轮训练中关注错误分类的数据样本

之所以在这个数据采用Boosting模型，是因为Boosting模型是基于决策树的，它比较适用于高维数据，不太需要降维。

决策树训练依据:熵Entropy

$$
H(s) = -\sum^c_{i=1}p_ilog_2 (p_i)
$$

其中 $p_i$ 代表类别i的比例

基因重要度排列，前十名：
![ad](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/output.png)

最终模型交叉验证错误率约为10%,训练集错误率为0，存在过拟合问题。

## Graphical Lasso

$$
cor_{x_ix_j|x_{-ij}}=-\frac{\omega_{ij}}{\sqrt{\omega_{ii}\omega_{jj}}}
$$

$\Omega = \Sigma^{-1}$, $\Sigma$ 是协方差矩阵, $\Omega$ 被称为精确度矩阵。

$\omega_{ii}$ 以及 $\omega_{jj}$ 代表 精度矩阵的对角元素，代表其条件方差的倒数。

$\omega_{ij}$ 代表精度矩阵的第i行第j列。


```{r}
library(huge)
library(qgraph)
glasso.out = huge(no_group,method = 'glasso',lambda = seq(0.001,0.5,length.out=100))
glasso.sel = huge.select(glasso.out,criterion = 'ebic')
bestind = glasso.sel$opt.index
precm = glasso.sel$icov[[bestind]]
filtered_precm = precm
filtered_precm[abs(filtered_precm) < 0.1] = 0
qgraph(round(filtered_precm,2),graph = 'cor',edge.labels = T,filetype = "png", filename = "partial_correlation_graph2.png")
```

![ad](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/output2.png)

![ad](https://github.com/Tony980624/Gene-Expression-And-Alzheimer/blob/main/file01/output3.png)

对比发现老年痴呆患者有两组基因之间失去了关联。
