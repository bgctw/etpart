# from https://www.r-bloggers.com/2021/04/quantile-regression-forests-for-prediction-intervals/
#install.packages(c("tidyverse","tidymodels","AmesHousing","gt"))
library(tidyverse)
library(tidymodels)
library(AmesHousing)
#library(gt)

ames <- make_ames() %>%
  mutate(Years_Old = Year_Sold - Year_Built,
         Years_Old = ifelse(Years_Old < 0, 0, Years_Old))
set.seed(4595)
data_split <- initial_split(ames, strata = "Sale_Price", prop = 0.75)
ames_train <- training(data_split)
ames_holdout  <- testing(data_split)

#RF models require comparably less pre-processing to linear models
rf_recipe <-
  recipe(
    Sale_Price ~ Lot_Area + Neighborhood  + Years_Old + Gr_Liv_Area +
      Overall_Qual + Total_Bsmt_SF + Garage_Area,
    data = ames_train
  ) %>%
  step_log(Sale_Price, base = 10) %>%
  step_other(Neighborhood, Overall_Qual, threshold = 50) %>%
  step_novel(Neighborhood, Overall_Qual) %>%
  step_dummy(Neighborhood, Overall_Qual)

rf_mod <- rand_forest() %>%
  set_engine("ranger", importance = "impurity", seed = 63233, quantreg = TRUE) %>%
  set_mode("regression")

set.seed(63233)
rf_wf <- workflows::workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(rf_recipe) %>%
  fit(ames_train)

# extract quantiles from fitting object
preds_bind <- function(data_fit, lower = 0.05, upper = 0.95){
  predict(
    rf_wf$fit$fit$fit,
    #workflows::pull_workflow_prepped_recipe(rf_wf) %>% bake(data_fit),
    workflows::extract_recipe(rf_wf) %>% bake(data_fit),
    type = "quantiles",
    quantiles = c(lower, upper, 0.50)
  ) %>%
    with(predictions) %>%
    as_tibble() %>%
    set_names(paste0(".pred", c("_lower", "_upper",  ""))) %>%
    mutate(across(contains(".pred"), ~10^.x)) %>%
    bind_cols(data_fit) %>%
    select(contains(".pred"), Sale_Price, Lot_Area, Neighborhood, Years_Old,
           Gr_Liv_Area, Overall_Qual, Total_Bsmt_SF, Garage_Area)
}
rf_preds_test <- preds_bind(ames_holdout)


pred_ranger_quantiles(rf_wf, ames_holdout)

# inspect prediction intervals
set.seed(1234)
rf_preds_test %>%
  mutate(pred_interval = ggplot2::cut_number(Sale_Price, 10)) %>%
  group_by(pred_interval) %>%
  #summarize(nrec = nrow(.))
  slice_sample(n = 2) %>%  #sample_n(2) %>%
  ggplot(aes(x = .pred))+
  geom_point(aes(y = .pred, color = "prediction interval"))+
  geom_errorbar(aes(ymin = .pred_lower, ymax = .pred_upper,
                    color = "prediction interval"))+
  geom_point(aes(y = Sale_Price, color = "actuals"))+
  scale_x_log10(labels = scales::dollar)+
  scale_y_log10(labels = scales::dollar)+
  labs(#title = "90% Prediction intervals on a holdout dataset",
       #subtitle = "Random Forest Model",
       y = "Sale_Price prediction")+
  theme_bw()+
  theme(legend.position = c(0.95,0.05), legend.justification = c(1,0)) +
  #coord_fixed() +
  theme()

