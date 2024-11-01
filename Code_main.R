## Task 1

gen_df = read.csv('C:/Users/85190/Desktop/stat/ass4/braindat.csv')
df = read.csv('C:/Users/85190/Desktop/stat/ass4/brain_sample_descriptions_PFC.csv')
df = df[,-1]
head(df)
head(gen_df)
library(dplyr)
library(tidyr)
# First, pivot gen_df to longer format
gen_df_long = gen_df %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression")

# Then, pivot wider to get genes as columns
gen_df_wide = gen_df_long %>%
  pivot_wider(names_from = GeneID, values_from = Expression)

# Check the result
head(gen_df_wide)

# Now, merge gen_df_wide with df on the sample IDs
# Ensure that the sample ID columns have the same name
df = df %>%
  rename(Sample = X.Sample_geo_accession)

# Merge the data frames
final_df <- df %>%
  inner_join(gen_df_wide, by = "Sample")

# View the final data frame
head(final_df)
final_df$Disease = ifelse(final_df$Disease == 'A',1,0)
final_df$Sex = ifelse(final_df$Sex == 'F',1,0)
rm(df,gen_df,gen_df_long,gen_df_wide)
# Normalization of gene data
pca_data = scale(final_df[,5:ncol(final_df)])
pca_result = prcomp(pca_data)
Loadings = pca_result$rotation
write.csv(Loadings,'D:/STATA/Loadings.csv')
explained_variance = (pca_result$sdev)^2 / sum((pca_result$sdev)^2)
plot(explained_variance[1:30],type = 'h',col='red',xlab = '',ylab = 'EV in %')
library(pheatmap)
pheatmap(pca_result$x[, 1:10], main = "Heatmap of First 10 Principal Components scores")

## Task 2

results = cbind(as.numeric(colnames(final_df)[5:ncol(final_df)]), 0, 0)

for (i in 5:ncol(final_df)) {
  model = lm(Disease ~ final_df[, i] + Age + Sex, data = final_df)
  results[i - 4, 2] = coef(model)[2]
  results[i - 4, 3] = summary(model)$coefficients[2, 4]
}


coefficients_df = data.frame(
  Gene_ID = results[, 1],
  Coefficient = as.numeric(results[, 2]),
  p_value = as.numeric(results[, 3])
)


write.csv(coefficients_df, 'D:/STATA/coefficients_df.csv', row.names = FALSE)
p_values_df = read.csv('D:/STATA/coefficients_df.csv')

adjusted_p_values = p.adjust(p_values_df$p_value,method = 'BH')

# Set FDR threshold
fdr_threshold = 0.1

# filt out significant
significant_genes = which(adjusted_p_values < fdr_threshold)
num_significant_genes = length(significant_genes)

# list first 30
if (num_significant_genes > 30) {
  top_genes = significant_genes[1:30]
} else {
  top_genes = significant_genes
}


significant_genes_df = p_values_df[top_genes,]
plot(log(1:length(p_values_df$p_value)), log(sort(p_values_df$p_value)), 
     xlab = "Log of Gene Index (ordered by p-value)", 
     ylab = "Log of Unadjusted p-value", 
     main = "P-value Plot with FDR Control Line")

abline(h = log(fdr_threshold), col = "red")
num_significant_genes
threshold_p_value = max(p_values_df$p_value[significant_genes])
threshold_p_value
top_genes_list = data.frame(Gene = p_values_df$Gene_ID[top_genes], Original_p_value = p_values_df$p_value[top_genes],coef = p_values_df$Coefficient[top_genes])
top_genes_list
library(knitr)

data = data.frame(
  Gene = c(10019481149, 10019687586, 10019799479, 10019927856, 10019948931, 
           10019975533, 10019977224, 10019977227, 10020014121, 10020044367,
           10020060841, 10020084138, 10020109481, 10020134015, 10020140475, 
           10020166212, 10020192767, 10020216227, 10020228202, 10020230502,
           10020234675, 10023804647, 10023804650, 10023804651, 10023804656, 
           10023804657, 10023804666, 10023804668, 10023804672, 10023804679),
  Original_p_value = c(0.029716265, 0.051135762, 0.024216107, 0.010350361, 
                       0.026922382, 1.39E-09, 3.01E-05, 0.000725046, 5.76E-12, 
                       9.61E-08, 1.18E-06, 9.22E-11, 0.01190894, 4.02E-07, 
                       2.40E-06, 1.65E-05, 0.000188067, 5.30E-13, 2.80E-07, 
                       3.81E-06, 0.03072534, 0.000350139, 9.03E-05, 0.037498767, 
                       2.07E-09, 1.38E-13, 6.58E-07, 0.004204457, 0.000466879, 
                       1.03E-07)
)


kable(data, caption = "Top 30 significant genes and their original p-values")

## Task 5  
# (The chosen classifier part is writen in python, 
# also included in '47616554_STAT3006_A4_supp.zip' named 'Code_task5_boosting.ipynb')

gen_df = read.csv('C:/Users/85190/Desktop/stat/ass4/braindat.csv')
df = read.csv('C:/Users/85190/Desktop/stat/ass4/brain_sample_descriptions_PFC.csv')
df = df[,-1]
gen_df_long = gen_df %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression")
gen_df_wide = gen_df_long %>%
  pivot_wider(names_from = GeneID, values_from = Expression)
df = df %>%
  rename(Sample = X.Sample_geo_accession)
final_df = df %>%
  inner_join(gen_df_wide, by = "Sample")
final_df$Disease = ifelse(final_df$Disease == 'A',1,0)
final_df$Sex = ifelse(final_df$Sex == 'F',1,0)
rownames(final_df) = final_df[,1]
final_df = final_df[,-1]
write.csv(final_df,'D:/STATA/Unormalized_df.csv')
library(glmnet)
?glmnet

X = as.matrix(final_df[, 5:ncol(final_df)])
y = as.factor(final_df$Disease)  
set.seed(47616554)  
cv_model = cv.glmnet(X, y, alpha = 1, family = "binomial")

best_lambda = cv_model$lambda.min
cat("Optimal lambda:", best_lambda, "\n")
plot(cv_model)
lasso_model = glmnet(X, y, alpha = 1, family = "binomial", lambda = best_lambda)

coefficients = as.matrix(coef(lasso_model))
coeff_df = data.frame(Gene = rownames(coefficients), Coefficient = coefficients[, 1])

important_genes = coeff_df[coeff_df$Coefficient != 0, ]

important_genes = important_genes[order(-abs(important_genes$Coefficient)), ]
print(important_genes)
set.seed(47616554)
overall_error = numeric(5)
class_1 = numeric(5)
class_2 = numeric(5)
folds = sample(rep(1:5, length.out = nrow(X)))

for (i in 1:5) {
  train_indices = which(folds!=i)
  test_indices = which(folds == i)
  model = glmnet(X[train_indices, ], y[train_indices], alpha = 1, family = "binomial", lambda = best_lambda)
  predictions_prob = as.numeric(predict(model, X[test_indices, ], type = "response")) 
  predicted_class = ifelse(predictions_prob>=0.5,1,0)
  overall_error[i] = sum(predicted_class!=y[test_indices])/length(test_indices)
  class_1[i] = sum(y[test_indices] == 1 & predicted_class!=y[test_indices])/length(test_indices)
  class_2[i] = sum(y[test_indices] == 0 & predicted_class!=y[test_indices])/length(test_indices)
}

mean(overall_error)
mean(class_2)
mean(class_1)

# Apparent error rate
model = glmnet(X, y, alpha = 1, family = "binomial", lambda = best_lambda)
predictions_prob = as.numeric(predict(model, X, type = "class")) 
sum(predictions_prob !=y)/230

## Task 6


loadings = pca_result$rotation
abs_loadings = abs(loadings)
abs_loadings_top3 = abs_loadings[, 1:3]

top_contributors = apply(abs_loadings_top3, 2, function(x) order(x, decreasing = TRUE)[1:10])

variable_names = rownames(loadings)

top_contributing_variables = apply(top_contributors, 2, function(indices) variable_names[indices])

top_contributing_variables



## Task 7

lasso_genes = c('10025907919','10025910280','10023828793','10023809844','10023814349','10025932435','10023815580','10025928439','10023827871','10023824734','10025908435','10023810924','10025905024','10023807849','10023825788','10025919391','10023821516','10031920093','10025909739','10025923114','10025907403','10031920282','10025930953','10025910042','10033668724','10025915075')
gens = cbind(final_df[,1],final_df[,lasso_genes,drop=F])
idx = which(final_df$Disease == 1)
yes_group = as.matrix(gens[idx,-1])
yes_group
no_group = as.matrix(gens[-idx,-1])
no_group
library(huge)
library(qgraph)
glasso.out = huge(no_group,method = 'glasso',lambda = seq(0.001,0.5,length.out=100))
glasso.sel = huge.select(glasso.out,criterion = 'ebic')
bestind = glasso.sel$opt.index
precm = glasso.sel$icov[[bestind]]
filtered_precm = precm
filtered_precm[abs(filtered_precm) < 0.1] = 0
qgraph(round(filtered_precm,2),graph = 'cor',edge.labels = T,filetype = "png", filename = "partial_correlation_graph2.png")
glasso.out = huge(yes_group,method = 'glasso',lambda = seq(0.001,0.5,length.out=100))
glasso.sel = huge.select(glasso.out,criterion = 'ebic')
bestind = glasso.sel$opt.index
precm = glasso.sel$icov[[bestind]]
filtered_precm = precm
filtered_precm[abs(filtered_precm) < 0.1] = 0
qgraph(round(filtered_precm,2),graph = 'cor',edge.labels = T,filetype = "png", filename = "partial_correlation_graph.png")


