---
title: "breast-cancer-wisconsin-Dataset1 analysis"
author: "Jaah21"
date: "10 Oct 2017"
---
# breast-cancer-wisconsin-Dataset1
---
Analysis and ML tests on R
---

```{r}
library(formatR)
```

```{r}
## data loading
bc_data <- read.table("breast-cancer-wisconsin.data.txt", 
    header = FALSE, sep = ",")
colnames(bc_data) <- c("sample_code_number", "clump_thickness", 
    "uniformity_of_cell_size", "uniformity_of_cell_shape", 
    "marginal_adhesion", "single_epithelial_cell_size", 
    "bare_nuclei", "bland_chromatin", "normal_nucleoli", 
    "mitosis", "classes")
bc_data$classes <- ifelse(bc_data$classes == "2", "benign", 
    ifelse(bc_data$classes == "4", "malignant", NA))
bc_data[bc_data == "?"] <- NA
```

```{r}
x <- c('u','c','l','a')
x
```


### how many NAs are in the data
```{r, echo=FALSE,comment='##'}
length(which(is.na(bc_data)))
```


```{r, results='markup'}
x <- 1:10
y <- round(rnorm(10, x, 1), 2)
df <- data.frame(x, y)
df
```

