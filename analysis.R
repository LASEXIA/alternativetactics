### Meta-analysis on alternative mating tactics ###
### Last update: May 2024 ###

# # # # # # # # # # # # # # # # # 
# Packages ----------------------
# # # # # # # # # # # # # # # # # 

library(metafor)
library(ape)
library(rotl)
library(scales)
library(dplyr)
library(ggplot2) 
library(orchaRd)
source("figure.R")

# # # # # # # # # # # # # # # # # 
# Data --------------------------
# # # # # # # # # # # # # # # # # 

# Importing data
tactics <- readxl::read_xlsx("data_tactics.xlsx", sheet = "data_tactics")

# Checking data
glimpse(tactics)
summary(tactics)

# Converting to factor
tactics <- tactics |> 
  mutate_if(is.character, as.factor)

# Importing variance-covariance matrix for study ID for all variables
diag1 <- readxl::read_xlsx("diag.xlsx", sheet = "diag", col_names = FALSE, na = "NA")

# Checking data
glimpse(diag1)

# Variance-covariance matrix as matrix
diag <- as.matrix(diag1)
colnames(diag) <- tactics$id
rownames(diag) <- tactics$id
isSymmetric.matrix(diag)

# # # # # # # # # # # # # # # # # 
# Phylogeny ---------------------
# # # # # # # # # # # # # # # # # 

# Converting it into character column sp2
tactics$sp2 <- as.character(tactics$species2)

# Creating a match between species in our database and species from tree web of life (TOL)
spp <- tnrs_match_names(unique(tactics$sp2, context_name = 'Animal'))
#View(spp)

# Creating the tree
my_tree <-  tol_induced_subtree(ott_ids=spp$ott_id)

# Removing TOL ID from tip label
otl_tips <- strip_ott_ids(my_tree$tip.label, remove_underscores=TRUE)

# Structuring
taxon_map <- structure(spp$search_string, names=spp$unique_name)

# Mapping taxon
my_tree$tip.label <- taxon_map[otl_tips] 

# No node label
my_tree$node.label <- NULL 

# Computing branch length
my_tree.ult <-  compute.brlen(my_tree, method = "Grafen")

# Lower case
my_tree.ult$tip.label <- stringr::str_to_lower(my_tree.ult$tip.label)

# Ploting tree
plot(my_tree.ult,no.margin = TRUE, cex=0.5) 

# Creating variance-covariance matrix
cov.matrix <-  vcv(my_tree.ult, corr = TRUE)


# # # # # # # # # # # # # # # # # 
# DESCRIPTION  ------------------
# # # # # # # # # # # # # # # # # 

## Overall data -----------------

# Checking number of studies (considering that each ID is a study)
unique(tactics$id) %>% length()

# Checking number of species
unique(tactics$species) %>% length()

# Checking order
unique(tactics$order) %>% length()

# Number of effect sizes per order
tactics %>% 
  group_by(order) %>% 
  count() %>% 
  arrange(desc(n))

# Number of the expression of tactics 
tactics %>% 
  group_by(type_tactics) %>% 
  count()


## Data without NA in type of tactic -----------

## Type of tactics: fixed or flexible

tactic_no_NA <- tactics |> 
  filter(tactics$type_tactics != "not_reported")

# Checking number of studies (considering that each ID is a study)
unique(tactic_no_NA$id) %>% length()

# Checking number of species
unique(tactic_no_NA$species) %>% length()

# Checking order
unique(tactic_no_NA$order) %>% length()

# Number of effect sizes per order
tactic_no_NA %>% 
  group_by(order) %>% 
  count() %>% 
  arrange(desc(n))

# Number of the type of tactics
tactic_no_NA %>% 
  group_by(type_tactics) %>% 
  count() %>% 
  arrange(desc(n))

# # # # # # # # # # # # # # # # # 
# TYPE OF MEASURE ---------------
# # # # # # # # # # # # # # # # # 

# Assumption test: there is no difference between the methods used to measure reproductive success (measures related to mating or to genetics) 

## Data -------------------------

## Genetic or mating
measure <- tactics |> 
  filter(measure1 != "") 

# Checking number of studies (considering that each ID is a study)
unique(measure$id) %>% length()

# Checking number of species
unique(measure$species) %>% length()

# Checking order
unique(measure$order) %>% length()

# Number of effect sizes per order
measure %>% 
  group_by(order) %>% 
  count() %>% 
  arrange(desc(n))

## Phylogeny -------------------------

# Dropping tips from the phylogeny by keeping only the tips from the measure dataset
## First: Creating an object with species names
measure_tips <- unique(as.character(measure$species2))

# Keeping tips from this data
measure_phylo <- keep.tip(phy = my_tree.ult, tip = measure_tips)

# Creating a new variance-covariance matrix
cov.matrix_measure <-  vcv(measure_phylo, corr = TRUE)

## Variance-covariance matrix for study ID ------------------------

# Importing variance-covariance matrix for study ID for measure data
measure_diag <- readxl::read_xlsx("diag.xlsx", sheet = "measure", col_names = FALSE, na = "NA")

# Checking data
glimpse(measure_diag)

# Variance-covariance matrix as matrix
measure_cvc <- as.matrix(measure_diag)
colnames(measure_cvc) <- measure$id
rownames(measure_cvc) <- measure$id
isSymmetric.matrix(measure_cvc)

## Model ----------------------------

m1 <- rma.mv(yi,vi, mods = ~measure1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_measure,
                      id = measure_cvc),
             control=list(optimizer="optim"),
             data = measure)

summary(m1)

## R2 ------------------------------

r2_ml(m1)

## Heterogeneity -------------------

i2_ml(m1) * 100

## Plot ----------------------------

# Creating a model without intercept to build the figure
m1_fig <- rma.mv(yi,vi, mods = ~measure1 -1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix,
                      id = measure_cvc),
             control=list(optimizer="optim"),
             data = measure)

# OrchaRd plot (modified version without prediction lines)
#jpeg("fig1_measure.jpg", width = 1400, height = 1000, res = 300)
tiff("fig1_measure_tif.tif", width = 1400, height = 1000, compression = "lzw", res = 300)
orchard(m1_fig, mod = "measure1", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  ylab("Measure") +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_y_discrete(labels = c("Mating", "Genetic"))+
  scale_size(range = c(0.5, 6), name = "Precision") + 
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()



# # # # # # # # # # # # # # # # # 
# OVERALL MODEL  -----------------
# # # # # # # # # # # # # # # # # 

# Hypothesis: males adopting main tactics will have higher reproductive success than males adopting secondary tactics

## Model ------------------------
m2 <- rma.mv(yi,vi, 
             random=list(~1|sp2, 
                         ~1|id, 
                         ~1|researcher, 
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs, 
                         ~1|species2), 
             R = list(sp2 = cov.matrix,
                      id = diag),
             data = tactics)

summary(m2)

## Heterogeneity --------------------------------------

i2_ml(m2, method = "ns") * 100

## Plot ---------------------------------------------

# Orchard plot (modified version without prediction lines)
#jpeg("fig2_overall.jpg", width = 1400, height = 1000, res = 300)
tiff("fig2_overall_tif.tif", width = 1400, height = 1000, compression = "lzw", res = 300)
orchard(m2, xlab= "Hedges' g") +
  labs(y = "Effect sizes") +
  scale_color_manual(values = "#808080") + 
  scale_fill_manual(values = "#0a0a0a") + 
  theme_classic() +
  theme(axis.title = element_text(size = 13), 
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_size(range = c(0.5, 6), name = "Precision") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()



# # # # # # # # # # # # # # # # # 
# EXPRESSION OF AMT  ------------
# # # # # # # # # # # # # # # # # 

# Hypothesis: the difference in the reproductive success should be lower for species that have fixed AMTs than for species that have flexible AMTs

## Phylogeny -------------------------

# Dropping tips from the phylogeny by keeping only the tips from the measure dataset
## First: Creating an object with species names
expression_tips <- unique(as.character(tactic_no_NA$species2))

# Keeping tips from this data
expression_phylo <- keep.tip(phy = my_tree.ult, tip = expression_tips)

# Creating a new variance-covariance matrix
cov.matrix_expression <-  vcv(expression_phylo, corr = TRUE)

## Variance-covariance matrix for study ID ------------------------

# Importing variance-covariance matrix for study ID for measure data
expression_diag <- readxl::read_xlsx("diag.xlsx", sheet = "expression", col_names = FALSE, na = "NA")

# Checking data
glimpse(expression_diag)

# Variance-covariance matrix as matrix
expression_cvc <- as.matrix(expression_diag)
colnames(expression_cvc) <- tactic_no_NA$id
rownames(expression_cvc) <- tactic_no_NA$id
isSymmetric.matrix(expression_cvc)

## Model ------------------------

m3 <- rma.mv(yi,vi, mods = ~type_tactics,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(m3)

## R2 ---------------------------------------------------

r2_ml(m3)


## Heterogeneity  ---------------------------------------

i2_ml(m3) * 100


## Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
m3_fig <- rma.mv(yi,vi, mods = ~type_tactics - 1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

tiff("fig3_AMTexpression_tif.tif", width = 1400, height = 1000, compression = "lzw", res = 300)
jpeg("fig3_AMTexpression.jpg", width = 1400, height = 1000, res = 300)
orchard(m3_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_y_discrete(labels = c("Fixed", "Flexible"))+
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()


# # # # # # # # # # # # # # # # # 
# BIAS  -------------------------
# # # # # # # # # # # # # # # # # 

# Funnel
funnel(tactics$yi, tactics$vi, yaxis="sei",
       xlab = "Effect size (Hedges' g)", digits = 2, las = 1) 

funnel(tactics$yi, tactics$vi, yaxis="seinv",
       #xlim = c(-3, 3),
       xlab = "Effect size (Hedges' g)",  digits = 2, las = 1)

# Calculating square root of variance
tactics$sei <- sqrt(tactics$vi)

## Study ---------------------------

# Bias of small study (equation 21 - Nakagawa et al, 2022)
pub_bias_study <- rma.mv(yi = yi, V = vi,
                           mod = ~1 + sei,
                           random=list(~1|sp2, 
                                       ~1|id, 
                                       ~1|researcher,
                                       ~1|experiment, 
                                       ~1|measure_timing, 
                                       ~1|comparison, 
                                       ~1|n_obs,
                                       ~1|species2), 
                           R = list(sp2 = cov.matrix,
                                    id = diag),
                           data=tactics)

summary(pub_bias_study)

# Creating predictions
study_pred <- predict(pub_bias_study)

# Organizing in a dataframe
dfstudy <- data.frame(sei = tactics$sei,
                     predi = study_pred$pred,
                     lower = study_pred$ci.lb,
                     upper = study_pred$ci.ub)

# Plot
study_fig <- ggplot(data = tactics, aes(x = sei, y = yi)) + 
  geom_point(aes(size = (1/sqrt(yi))), shape = 21, fill = "grey85", colour = "grey60", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) + 
  geom_line(data = dfstudy, aes(x = sei, y = predi), size = 1, colour = "#0a0a0a") + 
  geom_ribbon(data = dfstudy, aes(ymin = lower,  ymax = upper, y = 0), alpha = 0.3, fill = "#808080") + 
  labs(x = "Square root of variance (sei)", y = "Effect size (Hedges' g)", size = "Precision (1/SE)") + 
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        panel.grid = element_blank(),
        legend.text = element_text(size = 10), 
        legend.position = "top", 
        legend.title = element_text(size = 12), 
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))


## Time lag bias -------------------

# Time-lag bias (equation 23 - Nakagawa et al,  2022)
time_lag_overall <- rma.mv(yi = yi, V = vi,
                           mods= ~1 + year,
                           random=list(~1|sp2, 
                                       ~1|id, 
                                       ~1|researcher,
                                       ~1|experiment, 
                                       ~1|measure_timing, 
                                       ~1|comparison, 
                                       ~1|n_obs,
                                       ~1|species2),
                           R = list(sp2 = cov.matrix,
                                    id = diag),
                           data=tactics)

summary(time_lag_overall)

# Creating predictions
time_pred <- predict(time_lag_overall)

# Organizing a dataframe
dftime <- data.frame(year = tactics$year,
                     predi = time_pred$pred,
                     lower = time_pred$ci.lb,
                     upper = time_pred$ci.ub)

# Plot
time_fig <- ggplot(data = tactics, aes(x = year, y = yi)) + 
  geom_point(aes(size = (1/sqrt(yi))), shape = 21, fill = "grey85", colour = "grey60", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) + 
  geom_line(data = dftime, aes(x = year, y = predi), size = 1, colour = "#0a0a0a") + 
  geom_ribbon(data = dftime, aes(ymin = lower,  ymax = upper, y = 0), alpha = 0.3, fill = "#808080") + 
  labs(x = "Year of publication", y = "Effect size (Hedges' g)", size = "Precision (1/SE)") + 
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        panel.grid = element_blank(),
        legend.text = element_text(size = 10), 
        legend.position = "top", 
        legend.title = element_text(size = 12), 
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))

# Creating a figure with both biases
tiff("fig4_pubias_tif.tif", width = 1750, height = 1000, compression = "lzw", res = 300)
#jpeg("fig4_pubias.jpg", width = 1800, height = 1000, res = 300)
ggpubr::ggarrange(study_fig, time_fig, ncol = 2, common.legend = TRUE, legend = "top", labels = c("(a)", "(b)"))
dev.off()


# # # # # # # # # # # # # # # # # 
# Sensibility analysis ---------
# # # # # # # # # # # # # # # # # 

## Categories of AMT ------------

### Model ------------------------
s1 <- rma.mv(yi,vi, mods = ~comparison,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s1)

### R2 ---------------------------------------------------

r2_ml(s1)


### Heterogeneity  ---------------------------------------

i2_ml(s1) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s1_fig <- rma.mv(yi,vi, mods = ~comparison - 1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

orchard(s1_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c(rep("#808080", times = 11))) + 
  scale_fill_manual(values = c(rep("#0a0a0a", times = 11))) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))

## Categories of AMT: behavior x morphology ------------------------

### Model ------------------------

s4 <- rma.mv(yi,vi, mods = ~comparison2,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s4)

### R2 ---------------------------------------------------

r2_ml(s4)


### Heterogeneity  ---------------------------------------

i2_ml(m3) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s4_fig <- rma.mv(yi,vi, mods = ~comparison2 - 1,
                 random=list(~1|sp2, 
                             ~1|id,
                             ~1|researcher,
                             ~1|experiment, 
                             ~1|comparison, 
                             ~1|n_obs,
                             ~1|species2),
                 R = list(sp2 = cov.matrix_expression,
                          id = expression_diag),
                 control=list(optimizer="optim"), 
                 data = tactic_no_NA)


orchard(s4_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080", "#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a", "#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_y_discrete(labels = c("Fixed", "Flexible"))+
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()


## Taxonomic groups --------------------------------

### Model ------------------------

s2 <- rma.mv(yi,vi, mods = ~family,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s2)

### R2 ---------------------------------------------------

r2_ml(s2)


### Heterogeneity  ---------------------------------------

i2_ml(m3) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s2_fig <- rma.mv(yi,vi, mods = ~family - 1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)


orchard(s2_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c(rep("#808080", times = 35))) + 
  scale_fill_manual(values = c(rep("#0a0a0a", times = 35))) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()

## Experimental design --------------------------

### Model ------------------------

s3 <- rma.mv(yi,vi, mods = ~experiment,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s3)

### R2 ---------------------------------------------------

r2_ml(s3)


### Heterogeneity  ---------------------------------------

i2_ml(s3) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s3_fig <- rma.mv(yi,vi, mods = ~experiment - 1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|measure_timing, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

orchard(s3_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080", "#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a", "#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()

## Vertebrates x invertebrates  -----------------------------------

### Model ------------------------

s4 <- rma.mv(yi,vi, mods = ~vertebrate,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s4)

### R2 ---------------------------------------------------

r2_ml(s4)


### Heterogeneity  ---------------------------------------

i2_ml(m3) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s4_fig <- rma.mv(yi,vi, mods = ~measure_timing - 1,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = diag),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)


orchard(s4_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080", "#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a", "#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_y_discrete(labels = c("Fixed", "Flexible"))+
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()

## Timing -----------------------------------

### Model ----------------------------------

s4 <- rma.mv(yi,vi, mods = ~measure_timing,
             random=list(~1|sp2, 
                         ~1|id,
                         ~1|researcher,
                         ~1|experiment, 
                         ~1|comparison, 
                         ~1|n_obs,
                         ~1|species2),
             R = list(sp2 = cov.matrix_expression,
                      id = expression_cvc),
             control=list(optimizer="optim"), 
             data = tactic_no_NA)

summary(s4)

### R2 ---------------------------------------------------

r2_ml(s4)


### Heterogeneity  ---------------------------------------

i2_ml(m3) * 100


### Plot -------------------------------------------------

# Creating a model without intercept in order to build the figure
s4_fig <- rma.mv(yi,vi, mods = ~measure_timing - 1,
                 random=list(~1|sp2, 
                             ~1|id,
                             ~1|researcher,
                             ~1|experiment, 
                             ~1|comparison, 
                             ~1|n_obs,
                             ~1|species2),
                 R = list(sp2 = cov.matrix_expression,
                          id = diag),
                 control=list(optimizer="optim"), 
                 data = tactic_no_NA)


orchard(s4_fig, mod = "type_tactics", xlab = "Hedges' g") +
  scale_colour_manual(values = c("#808080", "#808080", "#808080", "#808080")) + 
  scale_fill_manual(values = c("#0a0a0a", "#0a0a0a", "#0a0a0a", "#0a0a0a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") +
  coord_flip() +
  scale_y_discrete(labels = c("Fixed", "Flexible"))+
  scale_size(range = c(0.5, 6), name = "Precision") +
  ylab("Expression of AMT") +
  guides(size = guide_legend(override.aes = list(colour = "#808080")))
dev.off()

# # # # # # # # # # # # # # # # # 
# Influential points ------------
# # # # # # # # # # # # # # # # # 

## Type of measure -----------------------

g1 <- rma.uni(yi,vi, data = measure)

# Testing influential points
leave1out <- leave1out(g1, progbar=TRUE)
hist(leave1out$estimate, breaks=40)
boxplot(leave1out$estimate)
boxplot(leave1out$ci.lb, ylim=c(0.1, 0.11))
boxplot(leave1out$ci.ub)

## Overall data ---------------------------
g2 <- rma.uni(yi,vi, data = tactics)

# Testing influential points
leave1out <- leave1out(g2, progbar=TRUE)
hist(leave1out$estimate, breaks=40)
boxplot(leave1out$estimate)
boxplot(leave1out$ci.lb, ylim=c(0.1, 0.11))
boxplot(leave1out$ci.ub)

## Expression of AMTs ---------------------------
g3 <- rma.uni(yi,vi, data = tactic_no_NA)

# Testing influential points
leave1out <- leave1out(g3, progbar=TRUE)
hist(leave1out$estimate, breaks=40)
boxplot(leave1out$estimate)
boxplot(leave1out$ci.lb, ylim=c(0.1, 0.11))
boxplot(leave1out$ci.ub)
