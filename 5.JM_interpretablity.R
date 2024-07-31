# cox模型可解释性

setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/interpretability")

library(survex)

# load data
load("JM_interpretability_pre.RData")

# 构建解释器
cox_exp <- explain(cox)


# global explanations(variable importance)
# permutation importance(排列重要性)：其他列不变，打乱某一列的随机排序后观察其预测结果的变化
variable_impor_auc <- model_parts(explainer = cox_exp,
                                  loss_function = loss_one_minus_cd_auc, # time-dependent auc
                                  type = "raw",
                                  output_type = "survival")
plot(variable_impor_auc)

variable_impor_intauc <- model_parts(explainer = cox_exp,
                                  loss_function = loss_one_minus_integrated_cd_auc, #based on integrated auc
                                  type="raw")
plot(variable_impor_intauc)

# local explanations(variable attributions, what variables contribute to the prediction)
# shap value
lmpre_HBS_BL$sex <- ifelse(lmpre_HBS_BL$sex == "M", "1", "2")
survshap <- predict_parts(explainer = cox_exp,
                          new_observation = lmpre_HBS_BL,
                          type = "survlime",
                          output_type = "survival")
plot(survshap)

# local explanations(variable attributions, what variables contribute to the prediction)
# shap value
survshap <- predict_parts(explainer = cox_exp,
                          new_observation = lmpre_HBS_BL,
                          type = "survshap",
                          output_type = "survival")
plot(survshap)
#ceteris paribus(how does a variable affect the prediction)
variable_profile <- predict_profile(explainer = cox_exp,
                                    new_observation = lmpre_HBS_BL[1,],
                                    type = "ceteris_paribus",
                                    output_type = "survival")

plot(variable_profile, numerical_plot_type = "contours")

# 所选变量是否会影响平均预测
cox_model_profile <- model_profile(explainer = cox_exp,
                                   type = "partial", # partial, accumulated or conditional
                                   output_type = "survival", # survival, chf (cumulative hazard function) or risk
                                   variables = c("moca","agediag","sex","educyrs","GBA.APOE4.PHS","pred"))

plot(cox_model_profile, numerical_plot_type = 'lines')

## model predictive performance
cox_performance <- model_performance(explainer = cox_exp, metrics = c('Integrated C/D AUC' = integrated_cd_auc, 'C/D AUC' = cd_auc))
plot(cox_performance)
