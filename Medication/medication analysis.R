library(mediation)
library(data.table)
library(tidyverse)
library(readxl)
library(tidyr)
library(ggplot2)
library(forcats)
library(dplyr)
library(openxlsx)
#load data
load("~/Data/ksads/covariat_all_18.Rdata")
load("~/Data//diag_all_name.Rdata")
load("~/Medication/hold_w24.Rdata")
load("~/Data/proteomics/proteomics_data_s.Rdata")
load("~/Data/proteomics/pc_scores.Rdata")
hold<-left_join(covariat_all_18%>%dplyr::select(1:4),proteomics_data_s,by="ABCNO")
diet<-read.csv("~/Data/New_Diet_PC2.csv")

#############################################################################################PC
pc<-cbind(hold%>%dplyr::select(1:4),pc_scores)


PM<-list(hold_w24,pc[,c(1,5:6)],proteomics_data_s, diet)%>%reduce(left_join,by="ABCNO")
PM[10,21:800] <- lapply(PM[10,21:800], function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
#X182 is Quinolinate 

#VEGFA########################################################3
model.M <- lm(X182 ~ VEGFA, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~VEGFA + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='VEGFA', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: VEGFA; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
#CCL3
model.M <- lm(X182 ~ CCL3, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~CCL3 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='CCL3', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: CCL3; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
#CCL3 
model.M <- lm(CCL3 ~ X182, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~CCL3 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='X182', mediator='CCL3',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
#####################################3CD5

model.M <- lm(X182 ~ CD5, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~CD5 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='CD5', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: CD5; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
######
model.M <- lm(CD5 ~ X182, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~CD5 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='X182', mediator='CD5',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)

#IL12B
model.M <- lm(X182 ~ IL_12B, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~IL_12B + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='IL_12B', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: IL-12B; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
##
model.M <- lm(IL_12B ~ X182, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~IL_12B + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='X182', mediator='IL_12B',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)

########MCP-1
model.M <- lm(X182 ~ MCP_1, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~ MCP_1 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='MCP_1', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
#FGF_23
model.M <- lm(X182 ~ FGF_23, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~ FGF_23 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='FGF_23', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)

#PC1
model.M <- lm(X182~PC1, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~PC1 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='PC1', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: Maternal inflammation PC1; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
##
model.M <- lm(PC1 ~ X182, na.action = na.omit,PM)
summary(model.M)

model.Y <- glm(adhd~PC1 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='X182', mediator='PC1',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
######################################BMI

model.M <- lm(X182 ~ mother_bmi, na.action = na.omit,PM)
summary(model.M)
model.X <- glm(adhd ~ mother_bmi, na.action = na.omit,family=binomial,PM)
summary(model.X)
model.Y <- glm(adhd~mother_bmi + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='mother_bmi', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

# Plot results
plot(results, main = "Treat: Mother BMI; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)

#####diabete

model.M <- lm(X182 ~ diabetes_all, na.action = na.omit,PM)
summary(model.M)
model.X <- glm(adhd ~ diabetes_all, na.action = na.omit,family=binomial,PM)
summary(model.X)
model.Y <- glm(adhd~diabetes_all + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='diabetes_all', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)



plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p  # Assuming ACME P-value is the same as d.avg.p

#########################diet PC2
model.M <- lm(X182 ~ Diet_PC2, na.action = na.omit,PM)
summary(model.M)
model.X <- glm(adhd ~ Diet_PC2, na.action = na.omit,family=binomial,PM)
summary(model.X)
model.Y <- glm(adhd~Diet_PC2 + X182, family=binomial,na.action = na.omit,PM)
summary(model.Y)

set.seed(234)
results <- mediate(model.M, model.Y, treat='Diet_PC2', mediator='X182',
                   boot=TRUE,sims = 1000)
summary(results)
plot(results)
prop_mediated <- results$n.avg * 100  # Convert to percentage
mediation_p_value <- results$n.avg.p
ACME_P_value <- results$d.avg.p 
# Plot results
plot(results, main = "Treat: Diet_PC2; Mediator: Quinolinate", 
     xlab = "Effect Sizes", ylab = "Mediation Pathways", 
     col = "darkblue", pch = 19, cex = 1.2)
grid(col = "gray")

# Dynamically adjust text position for the top-right corner
x_range <- par("usr")[1:2] # Get x-axis range
y_range <- par("usr")[3:4] # Get y-axis range

text_x <- max(x_range) - 0.18 * diff(x_range)  # Increase horizontal offset to move left
text_y_top <- max(y_range) - 0.05 * diff(y_range)  # Slightly below the top edge
text_y_middle <- max(y_range) - 0.1 * diff(y_range)  # Below the first text
text_y_bottom <- max(y_range) - 0.15 * diff(y_range)  # Below the second text

# Add annotations in the top-right corner
text(x = text_x, y = text_y_top, 
     labels = paste("Prop. Mediated (avg):", round(prop_mediated, 2), "%"), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_middle, 
     labels = paste("Prop. Mediated P-value:", signif(mediation_p_value, 3)), 
     col = "black", cex = 1.2)

text(x = text_x, y = text_y_bottom, 
     labels = paste("ACME P-value:", signif(ACME_P_value, 3)), 
     col = "black", cex = 1.2)
