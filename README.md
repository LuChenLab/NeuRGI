# Neutrophil Regulatory Gene Identifier (NeuRGI)

hezhifeng 2025/01/07



## What is NeuRGI?

NeuRGI is a machine-learning pipeline based on random forest, which also uses PU learning, GMM and neural networks to find genes that affect neutrophil differeatiation and function. This model is developed based on neutrophils, but it can also be expanded to different fields according to different data, such as the hematopoietic differentiation process (myeloid, lymphoid, erythroid) and so on.

You only need to provide a dataset containing  positive set and an unlabeled set, and the NeuRGI process can be used to predict all genes in the unlabeled set.

![img](https://github.com/LuChenLab/MRGI/blob/main/MRGI.jpg?raw=true)

## Requirements

The pipeline has been tested in `R 4.0.0` and the following packages are needed to be installed in `R`

```R
dplyr(>= 1.1.3)
e1071(>= 1.7-13)
pROC(>= 1.18.5)
randomForest(>= 4.6-14)
caret(>= 6.0-94)
PRROC(>= 1.3.1)
mclust(>= 6.0.0)
```



## Tutorial

### 1.Prepare

Loading pipeline and example data set

```R
source("../data/NeuRGI.R")
## NeuRGI has successfully been loaded.
InputData <- readRDS("../data/InputData_example.rds") #example data
```

`InputData` is a dataframe, each row is a gene and the rowname is gene symbol (all caps), each column is a feature, and the last column is `Label`, which is `factor` class.

In our example InputData, positive set was `Label` = 1,unlabel set was `Label` = 2.

```R
dim(InputData)
## [1] 19874    28
table(InputData$Label)
##  1			2
## 411	19463 
```

```R
head(InputData)
```

```R
## Feature1 Feature2 Feature3 Feature4   Feature5 Feature6    Feature7     Feature8  Feature9   Feature10   Feature11    Feature12
## A1BG           0        0    0.000    0.000 0.0000e+00        0 -0.06127166 7.974781e-01 0.7663269  20.4633333 -0.23486545 0.7481653905
## A1CF           0        0    0.604    0.692 5.9218e-10       71 -0.36010805 1.188509e-01 1.0000000   0.2400000 -0.10472665 0.1005016969
## A2M            0        0    0.529    0.769 4.5229e-11      147  0.39735134 8.276791e-02 0.8577036 520.5550000  5.19989430 0.0020573143
## A2ML1          0        0    0.590    0.692 1.6109e-34       94 -0.38161285 9.686237e-02 0.9271207   1.3950000 -0.63259736 0.0047968588
## A3GALT2        0        0    0.000    0.000 0.0000e+00        0  0.38129316 9.716572e-02 0.9550781   0.1066667  0.01114488 0.7274622865
## A4GALT         0        0    0.751    0.346 7.6668e-06       29  0.87535118 4.351285e-07 0.8974895   2.3900000  0.96993644 0.0002381155
##          Feature13   Feature14 Feature15 Feature16   Feature17   Feature18  Feature19 Feature20 Feature21 Feature22 Feature23 Feature24
## A1BG     0.0000000 0.000000000 0.0000000     0.000  0.00000000 0.000000000 0.07423696         0         1         0         0         0
## A1CF    -0.6849616 0.003413323 0.8333333     0.030 -0.03099908 0.009778623 0.07563611         0         1         0         0         0
## A2M      0.0000000 0.000000000 0.0000000     0.000  0.00000000 0.000000000 0.00000000         0         1         0         0         0
## A2ML1    0.0000000 0.000000000 0.0000000     0.000  0.00000000 0.000000000 0.08414691         0         1         0         0         0
## A3GALT2 -0.1010351 0.709660153 0.5306122     0.030 -0.01984006 0.203556996 0.09147422         1         0         0         0         0
## A4GALT  -0.3985131 0.126288754 0.7380952     0.265 -0.18155648 0.026188603 0.00000000         0         1         0         0         0
##         Feature25 Feature26   Feature27 Label
## A1BG    0.9747373 0.4963035 0.006955117     2
## A1CF    0.9747373 0.4963035 0.006955117     2
## A2M     0.9747373 0.4963035 0.006955117     2
## A2ML1   0.0000000 0.0000000 0.000000000     2
## A3GALT2 0.0000000 0.0000000 0.000000000     2
## A4GALT  0.9747373 0.4963035 0.006955117     2
```

### 2.PU-learning for reliable negative genes

In this step, we used PU-learning spy algorithm to get reliable negative genes.

```R
RN <- PUlearningForNeg(InputData)
## Reliable Negatives has been successfully found.
```

### 3.Create a balance training data

In this step, users can create a balance training data with gene type annotation using down-sampling method, and we offered gene type file `gene_type.rds`.

```R
gene_type <- readRDS("../data/01.gene_type.rds")
TrainData <- TrainingSetDownSample(InputData = InputData,RN = RN, gene_type = gene_type)
## Training set has been created.
table(TrainData$Label)
##  0   1 
## 411 411 
```

### 4.Training

NeuRGI use random forest algorithm to train model with `k` fold cross validation.

```R
RFmodels <- RFModelTrain(TrainData = TrainData, k = 10)
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Setting direction: controls < cases
## Training has been done.
```

user can check the performance of model using mean AUC, AUC-PR, ACC, MCC and F1 Score.

```R
RFmodels$Performances
## AUC mean   AUC-PR mean      ACC mean      MCC mean F1 Score mean 
## 0.9749717     0.9795671     0.9234352     0.8493263     0.9253931 
```

### 5. Prediction

Predict all genes not in training set. 

```R
Pred <- Prediction(InputData = InputData, TrainData = TrainData, RFModels = RFmodels,gene_type = gene_type)
## Prediction has been completed.
head(Pred)
##      symbol Functional Score  type
## 586    ACTN4                1 Other
## 600    ACVR1                1    MP
## 603   ACVRL1                1    MP
## 614   ADAM10                1    MP
## 722  ADIPOR1                1 Other
## 1009   ALCAM                1    MP
```

### 6. GMM classification

In this step, we use Gaussion Mixture Model to claasify **positive , negative and uncertain** genes depending on the prediction.

You can first check the distribution of the predicted function scores to determine the number of Gaussian mixture model fits `G`.

```R
Pred <- GMM(Pred = Pred,G = 3)
```

```R
table(Pred$classification)
##  Negative Uncertain  Positive 
##     1469     14047      3536
```

### (Optional) Feature Importance

```R
Importance <- VarImportance(RFmodels = RFmodels)
```

```R
head(Importance$df)
##            features MeanDecreaseGini normalized
## Feature4   Feature4         70.53377  1.0000000
## Feature6   Feature6         61.31508  0.8693010
## Feature21 Feature21         43.50472  0.6167928
## Feature5   Feature5         30.73842  0.4357972
## Feature3   Feature3         19.23685  0.2727325
## Feature16 Feature16         15.39690  0.2182912
```

```R
Importance$p
```

![img](https://github.com/LuChenLab/NeuRGI)

### (Optional) *In silico* knockout

You can further use [OntoVAE](https://github.com/hdsu-bioquant/onto-vae) to perform *in silico* knockout of predicted functional genes based on the own GMM classification results. We have placed the script at `OntoVAE.ipynb`



## SessionInfo

```R
sessionInfo()
```

```R
##R version 4.0.2 (2020-06-22)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.5 LTS

## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
##  [1] mclust_6.1.1        PRROC_1.3.1         caret_6.0-94        lattice_0.22-5      ggplot2_3.4.4       randomForest_4.6-14
##  [7] pROC_1.18.5         e1071_1.7-13        ROSE_0.0-4          dplyr_1.1.3        

## loaded via a namespace (and not attached):
##  [1] jsonlite_1.8.7       QuickJSR_1.1.3       splines_4.0.2        foreach_1.5.2        prodlim_2023.08.28   RcppParallel_5.1.7  
##  [7] StanHeaders_2.32.6   stats4_4.0.2         yaml_2.3.8           globals_0.16.2       ipred_0.9-14         pillar_1.9.0        
## [13] glue_1.6.2           digest_0.6.33        hardhat_1.3.0        colorspace_2.1-0      recipes_1.0.8        Matrix_1.6-3        
## [19] plyr_1.8.9           timeDate_4022.108    pkgconfig_2.0.3      rstan_2.32.6         listenv_0.9.0        purrr_1.0.2         
## [25] scales_1.2.1         processx_3.8.2       openxlsx_4.2.5.2     gower_1.0.1          lava_1.7.3           timechange_0.2.0    
## [31] tibble_3.2.1         proxy_0.4-27         farver_2.1.1         generics_0.1.3       withr_2.5.2          nnet_7.3-14         
## [37] cli_3.6.1            crayon_1.5.2         survival_3.5-7       magrittr_2.0.3       ps_1.7.5             future_1.33.0       
## [43] fansi_1.0.5          parallelly_1.36.0    nlme_3.1-164         MASS_7.3-60          class_7.3-17         pkgbuild_1.4.2      
## [49] loo_2.7.0            prettyunits_1.2.0    tools_4.0.2          data.table_1.16.0    matrixStats_1.1.0    lifecycle_1.0.4     
## [55] stringr_1.5.0        munsell_0.5.0        zip_2.3.0            callr_3.7.3          compiler_4.0.2       tinytex_0.50        
## [61] rlang_1.1.2          grid_4.0.2           iterators_1.0.14     rstudioapi_0.15.0    labeling_0.4.3       gtable_0.3.3        
## [67] ModelMetrics_1.2.2.2 codetools_0.2-16     inline_0.3.19        reshape2_1.4.4       R6_2.5.1             gridExtra_2.3       
## [73] lubridate_1.9.3      knitr_1.45           future.apply_1.11.0  utf8_1.2.4           stringi_1.8.2        parallel_4.0.2      
## [79] Rcpp_1.0.11          vctrs_0.6.3          rpart_4.1-15         tidyselect_1.2.0     xfun_0.43
```



## Contact

Please contact Lu Chen ([luchen@scu.edu.cn](mailto:luchen@scu.edu.cn)) or Zhifeng He ([hzf4304@163.com](mailto:hzf4304@163.com)).



## Citation

If you use NeuRGI pipeline in your publication, please cite by

