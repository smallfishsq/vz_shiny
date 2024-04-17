library(survey)
library(dplyr)
library(mice)
library(xgboost)
library(splitstackshape)
library(survival)
library(survminer)
library(ggplot2)
library(mice)
library(mlexperiments)
library(mlsurvlrnrs)
library(ParBayesianOptimization)

xgboostcox <- function(data, features, test='Train'){
  set.seed(123)
  params <- list(eta = 0.1, 
                 alpha = 0.75, 
                 lambda = 1, 
                 max_depth = 5,
                 max_child_weight = 30,
                 subsample = 0.8,
                 colsample_bytree = 0.8,
                 colsample_bylevel = 0.7,
                 colsample_bynode = 0.7,
                 verbose = FALSE,
                 objective = "survival:cox",
                 eval_metric = 'cox-nloglik')
  num_boost_round = 450
  early_stopping_round = 43
  
  features <- c("c_Age",
                "c_Gender",
                "c_SMQ",
                "c_BMI",
                "c_FHHtak",
                "c_FHAthma",
                "c_FHDibetes",
                "c_HighBP",
                "c_Insurance",
                "c_HeavyDrinker",
                "c_Anemia",
                "c_Asthma",
                "c_ConfusionMemory",
                "c_CHF",
                "c_CHD",
                "c_Htak",
                "c_Kidney",
                "c_WorkType",
                "c_TimeSinceHealcare", 
                "c_NumHopital",
                "c_Arthritis",
                "c_CancerType",
                "c_Diabetes")
  features_mort <- c(features, 'mortstat', 'ftime', 'wt')
  
  train_valid_split <- function(data){
    train_ind <- runif(nrow(data))
    data[, 'Split'] <- ifelse(train_ind<0.8,0,ifelse(train_ind<0.9,1,2))
    return(data)
  }
  
  create_dmatrix <- function(data, subset,test=test){
    if (test == 'Train'){
      datamat <- as.matrix(data[subset, features])
      mortstat <- data[subset, ]$mortstat
      ftime <- data[subset, ]$ftime+0.1
      wt <-data[subset, ]$wt
      return(xgb.DMatrix(data = datamat,
                         label = ifelse(mortstat==1,
                                        ftime,
                                        -ftime),
                         # label = ifelse(datamat[, 'mortstat']==1,
                         #                datamat[, 'ftime']+0.1,
                         #                -datamat[, 'ftime']),
                         #feature_names = features,
                         weight = wt))
    }else{
      wt <-data[subset, ]$wt
      return(xgb.DMatrix(data = data[subset, features],
                         #feature_names = features,
                         weight = wt))
    }
  }
  
  train_model <- function(data){
    data_train <- create_dmatrix(data, subset = which(data[, 'Split']==0), test='Train')
    data_validation <- create_dmatrix(data, which(data[, 'Split']==1), test='Train')
    watchlist <- list(train = data_train, 
                      eval = data_validation)
    trainmodel <- xgb.train(params = params,
                            data = data_train,
                            nrounds = num_boost_round,
                            early_stopping_rounds =early_stopping_round,
                            watchlist = watchlist)
    importance_matrix = xgb.importance(colnames(data_train), model = trainmodel)
    save(importance_matrix, file = 'featureImportance.RDATA')
    return(trainmodel)
  }
  
  predictfun <- function(model, data, test='Train'){
    if(test=='Train'){
      data_all <- create_dmatrix(data, subset = seq(1, nrow(data)), test='Train')
    }else{
      data_all <- create_dmatrix(data, subset = seq(1, nrow(data)), test='Test')
    }
    
    data[, 'Prediction'] <- predict(model, data_all)#, iterationrange = c(0, model$best_iteration)
    return(data)
  }
  
  cohort <- function(data){
    age_cohort <- cut(data$c_Age, breaks =c(17.9, 30, 40, 50, 60, 86) ,right = T)
    data$cohort <- paste0(age_cohort,'*',data$c_Gender,'*',data$c_SMQ)
    cohortmean <- data[, c('cohort', 'Prediction')]%>%group_by(cohort)%>%
      summarise(mean=mean(Prediction), sd=sd(Prediction))
    data_c <- merge(data, data.frame(cohortmean), by='cohort', all.x = T)
    data$Prediction_center <- (data_c$Prediction-data_c$mean)/data_c$sd
    saveRDS(data.frame(cohortmean), file='cohortmean.RDS')
    return(data)
  }
  
  data <- train_valid_split(data)
  
  filename <- 'model_coxph.model'
  
  if(test == 'Train'){
    trainmodel <- train_model(data)
    xgb.save(trainmodel, filename)
  }else{
    trainmodel <- xgb.load(filename)
  }
  
  data <- predictfun(trainmodel, data)
  
  if(test == 'Train'){
    data <- cohort(data)
  }else{
    age_cohort <- cut(data$c_Age, breaks =c(17.9, 30, 40, 50, 60, 86) ,right = T)
    data$cohort <- paste0(age_cohort,'*',data$c_Gender,'*',data$c_SMQ)
    cohortmean <- readRDS("cohortmean.RDS")
    data_c <- merge(data, cohortmean, by='cohort', all.x = T)
    data$Prediction_center <- (data_c$Prediction-data_c$mean)/data_c$sd
  }
  
  
  return(data)
}



xgboost_test <- function(data){
  features <- c("c_Age",
                "c_Gender",
                "c_SMQ",
                "c_BMI",
                "c_FHHtak",
                "c_FHAthma",
                "c_FHDibetes",
                "c_HighBP",
                "c_Insurance",
                "c_HeavyDrinker",
                "c_Anemia",
                "c_Asthma",
                "c_ConfusionMemory",
                "c_CHF",
                "c_CHD",
                "c_Htak",
                "c_Kidney",
                "c_WorkType",
                "c_TimeSinceHealcare", 
                "c_NumHopital",
                "c_Arthritis",
                "c_CancerType",
                "c_Diabetes")
  features_mort <- c(features, 'mortstat', 'ftime', 'wt')
  
  create_dmatrix <- function(data, subset){
    wt <-data[subset, ]$wt
    return(xgb.DMatrix(data = as.matrix(data[subset, features]),
                       weight = wt))
  }
  
  predictfun <- function(model, data){
    data_all <- create_dmatrix(data, subset = seq(1, nrow(data)))
    data[, 'Prediction'] <- predict(model, data_all)
    return(data)
  }
  
  trainmodel_filename <- 'model_coxph.model'
  trainmodel <- xgb.load(trainmodel_filename)
  
  data <- predictfun(trainmodel, data)
  
  age_cohort <- cut(data$c_Age, breaks =c(17.9, 30, 40, 50, 60, 86) ,right = T)
  data$cohort <- paste0(age_cohort,'*',data$c_Gender,'*',data$c_SMQ)
  cohortmean <- readRDS("cohortmean.RDS")
  data_c <- merge(data, cohortmean, by='cohort', all.x = T)
  data$Prediction_center <- (data_c$Prediction-data_c$mean)/data_c$sd
  
  return(data)
  
}
