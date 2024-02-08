# GWAS with R::rMVP

> https://github.com/xiaolei-lab/rMVP
> 
1. Input files
- VCF files
    - Maize V5 mapping files: [https://datadryad.org/stash/dataset/doi:10.5061/dryad.bnzs7h4f1](https://datadryad.org/stash/dataset/doi:10.5061/dryad.bnzs7h4f1)
    - **Whole-genome resequencing data** from 1515 total samples were used in this analysis, including 1276 previously published samples (Brandenburg et al., [2017](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0002); Bukowski et al., [2018](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0005); Chen et al., [2022](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0006); Chia et al., [2012](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0008); Kistler et al., [2018](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0026); Qiu et al., [2021](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0044); Unterseer et al., [2014](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0055); Wang et al., [2017](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0059), [2020](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16123#tpj16123-bib-0058)) and 239 lines resequenced as part of this study.
    
    > Grzybowski, Marcin et al. (2023). A common resequencing‐based genetic marker dataset for global maize diversity [Dataset]. Dryad. [https://doi.org/10.5061/dryad.bnzs7h4f1](https://doi.org/10.5061/dryad.bnzs7h4f1)
    > 
- Phenotype data
1. Installation

```r
install.packages("rMVP")
library(rMVP)
```

1. Process input file
- Phenotype file

```
Cul     KernelWeight
33-16   4.49746883103842
4226    4.13273913519756
4722    2.9322164109919
A188    5.21404330937991
A239    5.77581815507938
```

- Process vcf file format

```r
# Full-featured function (Recommended)
MVP.Data(fileVCF="myVCF.vcf",
         filePhe="Phenotype.txt",
         fileKin=True,
         filePC=True,
         out="mvp.vcf")

```

1. Runing GWAS

```r
library(rMVP)

ph <- read.table("mvp.vcf.phe", header = TRUE)
geno <- attach.big.matrix("mvp.vcf.geno.desc")
map <- read.table("mvp.vcf.geno.map" , head = TRUE)
Kin <- attach.big.matrix("mvp.vcf.kin.desc")
pc <-  attach.big.matrix("mvp.vcf.pc.desc")[]

imMVP <- MVP(
  phe=ph,
  geno=geno,
  map=map,
	CV.GLM = pc
  CV.FarmCPU=pc,
  CV.MLM = pc,
  K = Kin,
  priority="speed",
	vc.method="BRENT",  ##only works for MLM
  ncpus=16,
  maxLoop=10,
  method.bin="FaST-LMM", ## "FaST-LMM", "static" (#only works for FarmCPU)
  method=c("GLM", "MLM", "FarmCPU"),
  file.output = c("pmap", "pmap.signal", "plot", "log"),
  permutation.threshold=TRUE,
  permutation.rep=100,
  threshold=0.05,
)

imMVP <- cbind(imMVP$map, imMVP$glm.results, imMVP$mlm.results, imMVP$farmcpu.results)
data.table::fwrite(imMVP, "./KernalWeight_All.csv.gz")
```

> **phe:** phenotype data
> 
> 
> **geno:** genotype data
> 
> **map:** map data
> 
> **K:** Kinship matrix
> 
> **CV.GLM:** Covariates added in GLM
> 
> **CV.MLM:** Covariates added in MLM
> 
> **CV.FarmCPU:** Covariates added in FarmCPU
> 
> **nPC.GLM:** number of first columns of Principal Components added in GLM
> 
> **nPC.MLM:** number of first columns of Principal Components added in MLM
> 
> **nPC.FarmCPU:** number of first columns of Principal Components added in FarmCPU
> 
> **priority:  "speed"** or **"memory"** when calculating the genomic relationship matrix
> 
> **ncpus:** number of CPUs used for parallel computation, If not set, all CPUs will be used by default
> 
> **vc.method:** methods of variance components analysis, three methods are avaiblable, "BRENT", "EMMA", and "HE"
> 
> **maxLoop:** a parameter for FarmCPU only, the maximum iterations allowed in FarmCPU
> 
> **method.bin:** a parameter for FarmCPU only, two options are available: 'static' or 'FaST-LMM'
> 
> **permutation.threshold:** if **TRUE**, a threshold of permutation will be used in manhattan plot. 
> 
> **permutation.rep:** number of permutation replicates, only used when **permutation.threshold** is
> 
> **TRUE**
> 
> **threshold** , 0.05/marker size, a cutoff line on manhattan plot
> 
> **method**, models for association tests, three models are available in MVP, **"GLM"**,**"MLM"**, and **"FarmCPU"**
> 
> **file.output**, a Boolean value or a string vector. If **TRUE**, output all types of files. If **FALSE**, no files are output. For string vectors, the available values are c("pmap", "pmap.signal", "plot", "log"). Among them, pmap represents all SNP P-Val files, pmap.signal represents significant SNP files, Plot represents visualization results, and log represents log files.
> 
- If have more than one phenotype

```r
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
	  phe=phenotype[, c(1, i)],
	  geno=geno,
	  map=map,
		CV.GLM = pc
	  CV.FarmCPU=pc,
	  CV.MLM = pc,
	  K = Kin,
	  priority="speed",
		vc.method="BRENT",  ##only works for MLM
	  ncpus=16,
	  maxLoop=10,
	  method.bin="FaST-LMM", ## "FaST-LMM", "static" (#only works for FarmCPU)
	  method=c("GLM", "MLM", "FarmCPU"),
	  file.output = c("pmap", "pmap.signal", "plot", "log"),
	  permutation.threshold=TRUE,
	  permutation.rep=100,
	  threshold=0.05,
  )
  gc()
}
```

1. Additional plots
- Parameters

```
plot.type, four options ("d", "c", "m", "q"); if "d", draw SNP-density plot
bin.size, the window size for counting SNP number
bin.max, maximum SNP number, for windows, which has more SNPs than bin.max, will be painted in same color
col, colors for separating windows with different SNP density
file.type, format of output figure
dpi, resolution of output figure
```

- PCA

```r
pca <- attach.big.matrix("mvp.pc.desc")[, 1:3]
#pca <- prcomp(t(as.matrix(genotype)))$x[, 1:3]
MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c("red", "green", "yellow"), file.type="jpg")
```

- Manhattan plot in Circular fashion

```r
> MVP.Report(pig60K,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
      threshold=c(1e-6,1e-4),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
      "blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
      bin.size=1e6,outward=FALSE,file.type="jpg",memo="",dpi=300)

#Note:
#1. if signal.line=NULL, the lines that crosse circles won't be added.
#2. if the length of parameter 'chr.den.col' is not equal to 1, SNP density that counts 
#   the number of SNP within given size('bin.size') will be plotted around the circle.
```

- Manhattan plot in Rectangular fashion for single trait or method

```r
> MVP.Report(pig60K, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
        col=c("grey60","grey30"), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
        chr.den.col=c("darkgreen", "yellow", "red"),bin.size=1e6,signal.col=c("red","green"),
        signal.cex=c(1,1),signal.pch=c(19,19),file.type="jpg",memo="",dpi=300)
        
#Note:
#if the length of parameter 'chr.den.col' is bigger than 1, SNP density that counts 
#   the number of SNP within given size('bin.size') will be plotted.
```

1. Extract high P-value loci from result file

```bash
ls *FarmCPU.csv | while read file;do
	phe=$(echo $file|cut -f1 -d'.')
	awk -F ',' '{print $1"\t"$NF}' $file | awk '$NF<0.00001{print $0}' > $phe.FarmCPU.filtered.csv
	awk 'BEGIN{OFS="\t"} {gsub(/"/, "", $1); split($1, arr, "_"); print arr[1], arr[2], arr[2]}' $phe.FarmCPU.filtered.csv | grep -v 'SNP' > $phe.FarmCPU.hq.snp.bed
	rm $phe.FarmCPU.filtered.csv
done

```

1. mapping SNPs to peak region

```bash
ls /data21/wongzj/Maize_embryo/GWAS/kernel_GWAS/hq_snp/*FarmCPU.hq.snp.bed | while read file;do
file_name=$(echo $file|cut -f8 -d'/')	
phe=$(echo $file_name|cut -f1 -d'.')
	for i in 7 9 11 13; do
		bedtools closest -a /data21/yjiang22/new_hybrid_BM_embryo/allele-specific-acr/ase_gene2peak/MB-B-type/MB.B-type.peak.DAP$i.bed -b /data21/wongzj/Maize_embryo/GWAS/kernel_GWAS/hq_snp/$file_name -d | awk '$NF!=-1{print $0}'| awk '$NF <5000{print $0}' > $phe.DAP$i.MB.B.closest
		bedtools closest -a /data21/yjiang22/new_hybrid_BM_embryo/allele-specific-acr/ase_gene2peak/MB-M-type/MB.M-type.peak.DAP$i.bed -b /data21/wongzj/Maize_embryo/GWAS/kernel_GWAS/hq_snp/$file_name -d | awk '$NF!=-1{print $0}'| awk '$NF <5000{print $0}' > $phe.DAP$i.MB.M.closest
		bedtools closest -a /data21/yjiang22/new_hybrid_BM_embryo/allele-specific-acr/ase_gene2peak/BM-B-type/BM.B-type.peak.DAP$i.bed -b /data21/wongzj/Maize_embryo/GWAS/kernel_GWAS/hq_snp/$file_name -d | awk '$NF!=-1{print $0}'| awk '$NF <5000{print $0}' > $phe.DAP$i.BM.B.closest
		bedtools closest -a /data21/yjiang22/new_hybrid_BM_embryo/allele-specific-acr/ase_gene2peak/BM-M-type/BM.M-type.peak.DAP$i.bed -b /data21/wongzj/Maize_embryo/GWAS/kernel_GWAS/hq_snp/$file_name -d | awk '$NF!=-1{print $0}'| awk '$NF <5000{print $0}' > $phe.DAP$i.BM.M.closest
	done
done
```

1. get enrich with cisDynet

```bash
library(cisDynet)
input_dir = "/public1/home/sc80041/wongzijie/Maize_embryo/cisDynet_analysis/input_file/"
euclideanNorm_file = paste(input_dir,"BM_delength_CPM_title.txt",sep = "")  

for (i in Sys.glob("./*.snp")){
    GWASEnrichment(euclideanNorm_file = euclideanNorm_file,
                   snp = i,trait = gsub(".snp","",basename(i)),
                   output_path = "./BM_gwas_enrich")
    }
pdf(file="./BM_gwas_enrich/heatmap.pdf")
getGER(enrichment_result_path = "./BM_gwas_enrich/", return_matrix = F)
dev.off()
pdf(file="./BM_gwas_enrich/GWAS.pdf")
plotGER(result = res)
dev.off()

```