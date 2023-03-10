Copyright 2023 Sophie Ostmeier

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# install all packages if necessary
install.packages(c("readxl", "dplyr", "ggpubr","gt", "tidyverse", "writexl", 
                    "gtExtras", "xtable", "boot", "rstatix","glue"))

# packages
rm(list=ls())
library(readxl)
library(dplyr)
library(ggpubr)
library(gt)
library(tidyverse)
library(writexl)
library(gtExtras)
library(xtable)
library(boot)
library(rstatix)
library(glue)

# variables, names
# full path to input folder name (...) with data in xlsx format
input <- "..."
# full path to output file (...) in .tex
output <- "..."

# names of input files in .xlsx format
list_data <- c(
  "interexpert_B",
  "modelexpert_B",
  "interexpert_C",
  "modelexpert_C"
)

columnnames <- c(
  "Categories",
  "Metric", 
  "Expert B to A", 
  "Expert B to Best Model",
  "ratio of variance 1",
  "CI 1",
  "p-value BA",
  "Expert C to A",
  "Expert C to Best Model",
  "ratio of variance 2",
  "CI 2",
  "p-value CA"
)

metricnames <- c(
  "VS",
  "AVD [ml]",
  "Dice",
  "Precision",
  "Recall",
  "HD 95 [mm]",
  "SDT 5mm"
)

metricnames_categories <- c(
  "Volume"," ",
  "Overlap"," "," ",
  "Distance"," "
)

# list of metrices you would like to include, they have to match the column name where the values for one metric is stored
metric_segmentation <- c(
  "('1', 'Volumetric Similarity')","('1', 'Volume Absolute Difference')",
  "('1', 'Dice')","('1', 'Precision')","('1', 'Recall')",
  "('1', 'Hausdorff Distance 95')",
  "('1', 'Surface Dice Variable 5')")

set.seed(12345)

# functions
# function for precision analysis
precision <- function(group1, group2, reps=1000) {
  results <- sapply(1:reps, function(r) {
    list <- sample(1:length(group1), size= length(group1), replace= TRUE)
    group1.resample <- (group1)[list]
    group2.resample <- (group2)[list]
    var(group1.resample, na.rm = TRUE)/var(group2.resample, na.rm = TRUE)
  })
  sd <- sd(results)
  print("model-expert")
  print(var(group1, na.rm = TRUE))
  print("inter-expert")
  print(var(group2, na.rm = TRUE))
  var_ratio <- var(group1, na.rm = TRUE)/var(group2, na.rm = TRUE)
  print("ratio")
  print(var_ratio)
  print("sd")
  print(sd)
  print("pnorm")
  print(pnorm(var_ratio,1.1,sd,lower.tail = TRUE))
  print(" ")
  
  data.frame(cbind(var(group2, na.rm = TRUE), var(group1, na.rm = TRUE), var_ratio, sd*1.96, pnorm(var_ratio,1.1,sd,lower.tail = TRUE)))
}


loop <- function(dataBA, dataBM, dataCA, dataCM, indices){
  BA <-dataBA[,indices]
  BModel <-dataBM[,indices]
  CA <-dataCA[,indices]
  CModel <-dataCM[,indices]
  
  # precision analysis
  print(metricnames[indices])
  B_pres <- precision(BModel, BA)
  C_pres <- precision(CModel, CA)
  
  data.frame(cbind(B_pres,C_pres))
}

#####################
# loading  ata
data_list_data <- data.frame(row.names = metricnames)

for (a in 1:length(list_data)){
  
  setwd(input)
  
  # create segmentation taskps= 1000, test, mups= 1000, test, mups= 1000, test, mu
  segmentation_data <- read_excel(glue("{list_data[a]}.xlsx"), sheet = "all")
  
  # data
  segmentation_needed <- select(segmentation_data, metric_segmentation)
  df <- data.frame(apply(segmentation_needed, 2, function(x) as.numeric(as.character(x))))
  
  assign(paste0(list_data[a]), df)
  
}

#####################
# Analysis
df_wilcox <- data.frame(row.names=metricnames)

for (i in 1:length(metricnames)){
  
  temp <- loop(interexpert_B,
                modelexpert_B,
                interexpert_C,
                modelexpert_C,
                i)
  
  df_loop <- round(rbind(df_wilcox, assign(paste0(metricnames[i]), temp)), digits=2)
  rm(temp)
}

results <- cbind(metricnames_categories,metricnames,df_loop)

colnames(results) <- columnnames
for (b in colnames(results[ , grepl('CI' , names(results))])){
  results[ , b] <-sub("^","?? ",results[ , b])
}

for (i in colnames(results[ , grepl('p-value' , names(results))])){
  results[ , i] <- p.adjust(results[ ,i], method = "holm", n = length(metricnames)*4)
  results[ , i] <- symnum(results[ , i], corr = FALSE, na = FALSE, 
                          cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), 
                          symbols = c("p<0.0001","p<0.001", "p<0.01", "p<0.05", "non-sig"))
}

# view table in RStudio
results %>%
  gt() %>%
  gt_theme_espn()

# saves tables as latex table
print(xtable(results, type = "latex"), 
      file = output, 
      include.rownames=FALSE
)

