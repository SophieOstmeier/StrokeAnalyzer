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

install.packages(c("readxl", "dplyr", "ggpubr","gt", "tidyverse", "writexl", 
                    "gtExtras", "xtable", "boot", "rstatix","glue"))

# non-inferior analyses
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

# names of input files folders in .xlsx format
list_data <- c(
  "interexpert_B",
  "modelexpert_B",
  "interexpert_C",
  "modelexpert_C"
)

columnnames <- c(
  "Categories",
  "Metric", 
  "Expert B to A", "CI 1", 
  "Expert B to Best Model","CI 2", 
  "p-value BA",
  "Expert C to A","CI 3", 
  "Expert C to Best Model","CI 4",
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

# boundary (or effect size) for metrics with from 0 to 1
mu_0to1 <- -0.1
test_0to1 <- "greater" # less

# boundary for AVD
mu_SI_volume <- 2
test_SI <- "less" # greater

# boundary for HD 95
mu_SI_distance <- 5
test_SI <- "less" # greater

# H0: x - y < -0.1, The difference is smaller than -0.1 (-0.1 to !-!1.0). (Agreement is more than 10% worse)
# H1: x - y > -0.1, Is the difference value !greater! than -0.1? (Agreement is less than -10% worse)

######################
# functions
foo <- function(data, indices){
  dt<-data[indices,]
  c(
    median(dt[,1], na.rm= TRUE), median(dt[,2], na.rm= TRUE), median(dt[,3], na.rm= TRUE),
    median(dt[,4], na.rm= TRUE), median(dt[,5], na.rm= TRUE), median(dt[,6], na.rm= TRUE),
    median(dt[,7], na.rm= TRUE))
}

wilcox <- function(dataBA, dataBM, dataCA, dataCM, indices, mu, test){
  BA <-dataBA[,indices]
  BModel <-dataBM[,indices]
  CA <-dataCA[,indices]
  CModel <-dataCM[,indices]
  
  # bias analysis
  B_wilcox <- wilcox.test(x= BModel, y= BA, mu= mu, paired= TRUE, alternative= test)$p.value 
  C_wilcox <- wilcox.test(x= CModel, y= CA, mu= mu, paired= TRUE, alternative= test)$p.value
  
  data.frame(cbind(B_wilcox,C_wilcox))
}

#####################
# loading  ata
data_list_data <- data.frame(row.names = metricnames)

for (a in 1:length(list_data)){
  
  setwd(input)
  
  # create segmentation taskps= 1000, test, mups= 1000, test, mups= 1000, test, mu
  segmentation_data <- read_excel(glue("{list_data[a]}.xlsx"), sheet = "all")
  
  # median
  segmentation_needed <- select(segmentation_data, metric_segmentation)
  segmentation_median <- data.frame(round(apply(segmentation_needed, 2, median, na.rm=TRUE),digits = 2))
  
  # bootstrap
  df <- data.frame(apply(segmentation_needed, 2, function(x) as.numeric(as.character(x))))
  myBootstrap <- boot(df, foo, R=1000)
  
  # CI
  SD_metric <- round(apply(myBootstrap$t,2,sd)*1.96, digits = 2) # 90% or 95% CI?
  
  # bind everything
  metric_all <- cbind(segmentation_median, SD_metric)
  colnames(metric_all) <- c(list_data[a], glue("CI {list_data[a]}"))
  rownames(metric_all) <- metricnames
  data_list_data <- cbind(data_list_data, metric_all)
  
  assign(paste0(list_data[a]), df)
  
}

#####################
# Analysis
df_wilcox <- data.frame(row.names=metricnames)

for (i in 1:length(metricnames)){
  mu <- mu_0to1
  test <- test_0to1
  if (i==2){ # for AVD
    mu <- mu_SI_volume
    test <- test_SI 
  }
  if (i==6){ # for HD
    mu <-  mu_SI_distance
    test <- test_SI 
  }
  temp <- wilcox(interexpert_B,
                 modelexpert_B,
                 interexpert_C,
                 modelexpert_C,
                 i, mu, test)
  
  df_wilcox <- rbind(df_wilcox, assign(paste0(metricnames[i]), temp))
  rm(temp)
}

results <- cbind(metricnames_categories,metricnames,
                 data_list_data[,c(1,2,3,4)],df_wilcox[,1],
                 data_list_data[,c(5,6,7,8)],df_wilcox[,2]
)

colnames(results) <- columnnames
for (b in colnames(results[ , grepl('CI' , names(results))])){
  results[ , b] <-sub("^","± ",results[ , b])
}

for (i in colnames(results[ , grepl('p-value' , names(results))])){
  results[ , i] <- p.adjust(results[ ,i], method = "holm", n = length(metricnames)*4)
  results[ , i] <- symnum(results[ , i], corr = FALSE, na = FALSE, 
                          cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), 
                          symbols = c("p<0.0001","p<0.001", "p<0.01", "p<0.05", "non-sig"))
}


# define title for both settings
title <- glue("One Sided")
subtitle <- glue("Boundary for relative metrics: {mu_0to1}, for volume metrics: {mu_SI_volume} and for distance metrics: {mu_SI_distance}")

# view table in RStudio
results %>%
  gt() %>%
  tab_header(
    title = title,
    subtitle =subtitle
  ) %>%
  gt_theme_espn()

# saves tables as latex table
print(xtable(results, type = "latex"), 
      file = output, 
      include.rownames=FALSE
)

