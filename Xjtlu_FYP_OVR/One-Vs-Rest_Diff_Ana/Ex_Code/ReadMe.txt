代码说明：
1. Data_Pre 对于学姐给的原始的数据补充基因值。然后做初步的筛选，选择tcga barcode为1和3的作为原位癌数据，删除掉超过25%表达值为0的基因。
2.Data_clean 对于上一步生成的数据集，选择基因的交集或者并集两种方法生成数据的矩阵
3.clinical 对临床的数据进行观察，一些原位癌数据分布之类的
4. Diff_groups 生成差异分析需要用到的one-vs-rest的表达矩阵和标签矩阵
5.Diff_select 差异分析结束后，我们对每一个癌症都可以获得一个基因列表，包含了一些评分例如log2Foldchange值，我们使用这些标准对基因再一次筛选，然后生成用于模型构建的表达矩阵。
6 model 使用一些机器学习方法建模 例如xgboost， catboost lr，绘制roc 计算fnr
7. move_example 根据训练的模型构建转移癌的表达矩阵，然后预测模型在转移癌上的预测能力。
8. MLP MLP训练的模型
9. TSP_Origin.R   TSP的训练模型
10. Diff.R 差异分析的R分析模型，