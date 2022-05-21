library(tidyverse)
library(ggplot2)
library(tidymodels)
library(usemodels)
library(ggpubr)
library(vip)
library(patchwork)


## Import data
data <- read.csv("M1.csv")

# Splitting
set.seed(123)
data_split <- initial_split(data, prop = 0.70, strata = Group)
data_train <- training(data_split)
data_test <- testing(data_split)

## Re-partitioning training data: Bootstrap resamples
set.seed(243)
data_boot <- bootstraps(data_train, times = 100, strata = Group)
use_ranger(Group ~., data = data_train)     ## provides a recipe for models; copy recipe and edit

## Set model parameters
ranger_recipe <- 
	  recipe(formula = Group ~ ., data = data_train) %>% 
	  step_string2factor(one_of("Group")) 

ranger_spec <- 
	  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
	  set_mode("classification") %>% 
	  set_engine("ranger") 
  
ranger_workflow <- 
	  workflow() %>% 
	  add_recipe(ranger_recipe) %>% 
	  add_model(ranger_spec) 

doParallel::registerDoParallel()
set.seed(55513)
ranger_tune <-
  	tune_grid(ranger_workflow, resamples = data_boot, grid = 25)
  
## Explore results of model
autoplot(ranger_tune)

# finalize model
final_rf <- ranger_workflow %>% 
  	  finalize_workflow(select_best(ranger_tune))
  
##Fit model on test data
data_fit <- last_fit(final_rf, data_split)

##Evaluate
collect_metrics(data_fit)


##ROC curve
M24_roc_plot <- data_fit %>% 
	    collect_predictions() %>% 
	    roc_curve(Group, .pred_Condition) %>% 
	    ggplot(aes(1-specificity, sensitivity, colour = 'red')) +
	    geom_abline(lty = 2, color = "gray80", size = 1.5) +
	    geom_path(show.legend = F, alpha = 0.6, size = 1.2) + 
	    coord_equal() + theme_pubr()
  
  
## Feature importance
imp_spec <- ranger_spec %>% 
    finalize_model(select_best(ranger_tune)) %>% 
    set_engine("ranger", importance= "permutation")

malaria_imp <- workflow() %>% 
    add_recipe(ranger_recipe) %>% 
    add_model(imp_spec) %>% 
    fit(data_train) %>% 
    extract_fit_parsnip() %>% 
    vip(aesthetics = list(alpha = 0.8, fill = "midnightblue"), num_features = 30) + theme_pubr()

    
