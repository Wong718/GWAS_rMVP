# Phenotype data processing

1. Input data : Use pandas.melt(id_vars='cul’) to generate the file

![Untitled](Phenotype%20data%20processing%20e0cdd574d8f94cb6a82c408e20590714/Untitled.png)

1. Data preprocessing

```r
> setwd("D:/GWAS_phe")
> raw_data <- read.csv("kernel20_filtered_melt.csv", header = T, check.names = F)
> dat=raw_data
> library('magrittr')
> head(dat)
  cul env kernelweight
1 Tzi10 06A           NA
2  K148 06A           NA
3 NC354 06A          3.3
4  Ki14 06A           NA
5  L317 06A          5.2
6  Mo46 06A           NA
> for(i in 1:2) dat[,i] = as.factor(dat[,i])  #前两列设置为factor
```

1. BLUE value calculating

```r
> library('lme4')
> m1 = lmer(kernelweight~ cul + (1|env), data = dat)
#固定因子：cul
#随机因子：env 
> as.data.frame(fixef(m1))
> library('lsmeans')
> re=lsmeans(m1,"cul") #结果中的lsmeans即为各品种该性状的BLUE值
```

1. Output results

```r
> blue.df<-data.frame(re)[,1:2]
> colnames(blue.df)<-c('Cul','Kernelweight')
> head(blue.df)
	Cul KernelWeight
1 33-16     4.497469
2  4226     4.132739
3  4722     2.932216
4  A188     5.214043
5  A239     5.775818
6  A272     4.021764
> write.table(blue.df, file = "KernelWeight_blue.txt", row.names = F, col.names = T, quote = F, sep = "\t")
```

1. If we like to process multiple phenotype …

```r
library('magrittr')
library('lme4')
library('lsmeans')

setwd("E:/Jupyter_notebook/杂种胚/GWAS/Phenotype/kernel_phe")
csv_files <- list.files(pattern = "\\.csv$")
combined_data <- data.frame()

for (file in csv_files) {
	phenotype_name <- sub("\\_filtered_melt.csv$", "", file)
  raw_data <- read.csv(file, header = TRUE, check.names = FALSE)
  for (i in 1:2) {
    raw_data[, i] = as.factor(raw_data[, i])
  }
  m1 <- lmer(as.formula(paste(phenotype_name, "~ cul + (1|env)")), data = raw_data)
  fixed_effects <- as.data.frame(fixef(m1))
  re <- lsmeans(m1, "cul")
  blue.df <- data.frame(re)[, 1:2]
  colnames(blue.df) <- c("Cul", phenotype_name )
  output_file <- gsub(".csv", "_blue.txt", file )
  write.table(blue.df, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

```