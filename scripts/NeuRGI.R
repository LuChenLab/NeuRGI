#---
#title: "NeuRGI"
#author: "Hezhifeng"
#date: "2025/1/7"
#---

library(dplyr)
library(ROSE)
library(e1071)
library(pROC)
library(randomForest)
library(caret)
library(PRROC)
library(ggplot2)

mcc <- function(conf_matrix){
  TP <- conf_matrix[1, 1]
  TN <- conf_matrix[2, 2]
  FP <- conf_matrix[1, 2]
  FN <- conf_matrix[2, 1]
  
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  return(numerator / denominator)
}

F1_Score <- function(y_true, y_pred) {
  cm <- confusionMatrix(as.factor(y_pred), as.factor(y_true))
  precision <- cm$byClass['Precision']
  recall <- cm$byClass['Recall']
  return(2 * ((precision * recall) / (precision + recall)))
}

PUlearningForNeg <- function(InputData){
  groudTruth <- InputData[InputData$Label == 1, ]
  Unlabeled <- InputData[InputData$Label != 1, ]
  indx <- sample(2, nrow(groudTruth), replace = T, prob = c(0.9, 0.1))
  s <- groudTruth[indx == 2,]
  Ps <- groudTruth[indx == 1,]
  Us <- rbind(Unlabeled, s )
  Ps$NewLabel <- 1
  Us$NewLabel <- 0
  rn_train <- rbind(Ps, Us)
  balance.over <- rn_train[-(ncol(rn_train)-1)]
  classifier <- naiveBayes(NewLabel ~., balance.over)
  rn_pred_u <- predict(classifier, Unlabeled[, -c(ncol(rn_train)-1)], type = "raw" )
  rn_pred_s <- predict(classifier, s[,-c(ncol(rn_train)-1)], type = "raw" )
  tr <- quantile(rn_pred_s[,2], probs = seq(0,1,0.1))[2][[1]]
  index_RN <- which(rn_pred_u[,2] < tr)
  RN <- Unlabeled[index_RN,]
  RN$Label <- 0
  cat("Reliable Negatives has been successfully found.")
  return(RN)
}

TrainingSetDownSample <- function(InputData, RN, gene_type){
  RN_genes <- data.frame(symbol = rownames(RN))
  RN_genes$type[RN_genes$symbol%in%gene_type$Enz]<-"Enz"
  RN_genes$type[RN_genes$symbol%in%gene_type$Mp] <- "MP"
  RN_genes$type[RN_genes$symbol%in%gene_type$TF]<-"TF"
  RN_genes$type[RN_genes$symbol%in%gene_type$RBP]<-"RBP"
  RN_genes$type[!RN_genes$type%in%c("Enz","MP","TF","RBP")]<-"Other"
  
  groudTruth <- InputData[InputData$Label == 1, ]
  groudTruth$type[rownames(groudTruth)%in%gene_type$Enz]<-"Enz"
  groudTruth$type[rownames(groudTruth)%in%gene_type$Mp] <- "MP"
  groudTruth$type[rownames(groudTruth)%in%gene_type$TF]<-"TF"
  groudTruth$type[rownames(groudTruth)%in%gene_type$RBP]<-"RBP"
  groudTruth$type[!groudTruth$type%in%c("Enz","MP","TF","RBP")]<-"Other"
  
  sample_sizes = groudTruth$type %>% table()
  
  RN_genes <- RN_genes %>%
    group_by(type) %>%
    filter(type %in% names(sample_sizes)) %>%  
    group_modify(~ .x %>% slice_sample(n = sample_sizes[.y$type], replace = FALSE)) %>%  # 使用 group_modify 进行按类型采样
    ungroup()
  
  balance.over_2 <- rbind(InputData[InputData$Label == 1, ],InputData[rownames(InputData)%in%RN_genes$symbol,])
  balance.over_2$Label <- as.numeric(balance.over_2$Label)
  balance.over_2$Label[balance.over_2$Label != 1] <- 0
  balance.over_2$Label <- as.factor(balance.over_2$Label)
  cat("Training set has been created.")
  return(balance.over_2)
}

RFModelTrain <- function(TrainData,k){
  n <- length(names(TrainData))
  rates <- rep(NA, n-1) 
  for(i in 1:(n-1)){
    rf_train<-randomForest(TrainData$Label~.,data=TrainData,mtry=i,na.action = na.roughfix,ntree=1000)
    rates[i]<-mean(rf_train$err.rate)  
    # print(rates[i])
  }
  mtry_num = which(rates == min(rates)) 

  folds <- createFolds(y=TrainData$Label,k=k)
  RF_classifier <- list()
  roc_obj <- list()
  auc_value <- numeric(k)
  pr_curve_data <- list()
  pr_auc_values <- numeric(k)
  f1_score <- numeric(k)
  mcc_value <- numeric(k)
  acc_value <- numeric(k)
  for(i in 1:10){
    fold_test <- TrainData[folds[[i]],] 
    fold_train <- TrainData[-folds[[i]],]
    
    RF_classifier_v1 <- randomForest(Label ~ ., data = fold_train, mtry = mtry_num, ntree = 1000, importance = T)
    RF_classifier[[i]] <- RF_classifier_v1
    
    #####  AUROC
    pred_rf_v1 <- predict(RF_classifier_v1, newdata = fold_test[,-ncol(fold_test)], type = "prob")
    obs_p_rf_v1 <- data.frame(pred = pred_rf_v1[,2],
                              obs = fold_test$Label )
    rf_roc_v1 <- roc(obs ~ pred, obs_p_rf_v1, levels = c("0", "1"))
    auc_value_v1 <- as.numeric(auc(rf_roc_v1))
    
    
    # PRROC
    pr <- pr.curve(scores.class0 = obs_p_rf_v1$pred[obs_p_rf_v1$obs == 1],
                   scores.class1 = obs_p_rf_v1$pred[obs_p_rf_v1$obs == 0],
                   curve = T)
    pr_area_v1 <- pr$auc.integral  
    
    # F1 Score
    threshold <- 0.5  # You can adjust the threshold if necessary
    preds_class <- ifelse(obs_p_rf_v1$pred >= threshold, 1, 0)
    f1_score_v1 <- F1_Score(obs_p_rf_v1$obs, preds_class)
    
    # MCC
    cm <- confusionMatrix(as.factor(preds_class), as.factor(obs_p_rf_v1$obs))
    mcc_value_v1 <- mcc(cm$table)
    
    # Accuracy
    acc_value_v1 <- cm$overall['Accuracy']
    
    # Store results
    auc_value[i] <- auc_value_v1
    f1_score[i] <- f1_score_v1
    mcc_value[i] <- mcc_value_v1
    acc_value[i] <- acc_value_v1
    pr_auc_values[i] <- pr$auc.integral
  }
  
  result <- c(mean(auc_value),
              mean(pr_auc_values),
              mean(acc_value),
              mean(mcc_value),
              mean(f1_score)) 
  names(result) <- c("AUC mean","AUC-PR mean","ACC mean", "MCC mean", "F1 Score mean")
  cat("Training has been done.")
  return(list(Models = RF_classifier, Performances = result))
}

VarImportance <- function(RFmodels){
  Gini <- list()  
  for (i in 1:length(RFmodels)){
    tmp1 <- RFmodels$Models[[i]]$importance %>% as.data.frame()
    tmp2 <- tmp1[,4] %>% as.data.frame()
    rownames(tmp2) <- rownames(tmp1)
    Gini[[i]] <- tmp2
  }
  Gini <- do.call(cbind,Gini)
  Gini <- apply(Gini,1,mean)
  
  dat <- data.frame(features = names(Gini),MeanDecreaseGini = Gini)
  dat <- dat[order(dat$MeanDecreaseGini,decreasing = TRUE),]
  dat$features <- factor(dat$features, levels = dat$features)
  dat$normalized <- dat$MeanDecreaseGini/max(dat$MeanDecreaseGini)
  
  
  p_var <- ggplot(dat,aes(x = features,y = normalized))+
    geom_point(shape = 21, size = 4)+
    theme_classic()+
    theme(axis.text.x = element_text(size = 12,color = "black",angle = 45,hjust = 1),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(.8,.4))+
    labs(y = "Relative Importance score", x = "")+
    scale_x_discrete(labels=function(x) stringr::str_wrap(x, width=45))
  return(list(df = dat, p = p_var))
}

Prediction <- function(InputData,TrainData,RFModels,gene_type){
  testdata <- InputData[!rownames(InputData)%in% rownames(TrainData),]
  RF_pred_list <- list()
  for (i in 1:length(RFModels$Models)){
    RF_pred_v1 <- predict( Models$Models[[i]],testdata[,-ncol(testdata)], type = "prob") %>% as.data.frame()
    colnames(RF_pred_v1) <- c( "NonFunc", "Functional Score" )
    RF_pred_list[[i]] <- RF_pred_v1[,-1]
  }
  Func_Score <- as.data.frame(do.call(cbind,RF_pred_list)) %>% apply(.,1,mean)
  RF_mean_all <- as.data.frame(cbind(symbol = rownames(testdata),`Functional Score` = Func_Score))
  RF_mean_all$`Functional Score` <- as.numeric(RF_mean_all$`Functional Score`)
  RF_mean_all <- RF_mean_all[order(RF_mean_all$`Functional Score`,decreasing = TRUE),]
  RF_mean_all$type[RF_mean_all$symbol%in%gene_type$Enz]<-"Enz"
  RF_mean_all$type[RF_mean_all$symbol%in%gene_type$Mp] <- "MP"
  RF_mean_all$type[RF_mean_all$symbol%in%gene_type$TF]<-"TF"
  RF_mean_all$type[RF_mean_all$symbol%in%gene_type$RBP]<-"RBP"
  RF_mean_all$type[!RF_mean_all$type%in%c("Enz","MP","TF","RBP")]<-"Other"
  
  cat("Prediction has been completed.")
  return(RF_mean_all)
}

GMM <- function(Pred,G){
  gmm <- Mclust(Pred$`Functional Score`, G = G, modelNames = c("E","V"))
  if (G == 3){
    Pred$classification <- gmm$classification
    Pred$classification <- as.factor(Pred$classification)
    Pred$classification <- plyr::mapvalues(from = c("1", "2", "3"), 
                                           to = c("Negative","Uncertain","Positive"), 
                                           x = Pred$classification)
  }else{
    Pred$classification <- gmm$classification
  }
  
  return(Pred)
}

cat("NeuRGI has successfully been loaded.")