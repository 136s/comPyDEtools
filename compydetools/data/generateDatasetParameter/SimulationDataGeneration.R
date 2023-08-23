# modified from `compareDEtools::generateDatasetParameter()`
# https://github.com/unistbig/compareDEtools/blob/7f7eecf301070a3b7c7c037da079b331b722622b/R/SimulationDataGeneration.R
# To modify the code minimally for outputting parameters to CSV for simulating KIRC and Bottomly.
# requires: edgeR, Biobase
# execution environment: Windows 10, R 4.1.3, edgeR 3.36.0, Biobase 2.54.0

library(edgeR)
library(Biobase)

#' Generate List containing estimated mean and dispersion parameters and filtered count from original count dataset.
#' @export
generateDatasetParameter = function(){
  #  Effect size 1.2/1.5~
  # Load TCGA KIRC RNA-seq data. 144 samples (72 cancer and matched normal, respectively)
  download.file("https://github.com/cran/SimSeq/raw/256a72d0b3cd24ff68541624ce35a6e9182cf450/data/kidney.rda", "kidney.rda")
  load("kidney.rda")
  k_count = kidney$counts # RNA-seq read count data
  index.cancer=(1:72)*2 # cancer sample index
  index.normal=index.cancer-1 # normal sample index

  k_count= k_count[,c(index.cancer, index.normal)] # Arrange samples for convinience

  # Get count mean and dispersion using edgeR package

  # Mean and dispersion values are obtained separately from cancer and normal samples when different dispersion is assummed between two sample types.
  # Mean and dispersion values from normal samples
  dge.normal=DGEList(counts=k_count[,73:144], group = factor(rep(2,72)))
  dge.normal=calcNormFactors(dge.normal)
  dge.normal=estimateCommonDisp(dge.normal)
  dge.normal=estimateTagwiseDisp(dge.normal)
  disp.normal = dge.normal$tagwise.dispersion # Dispersion
  mean.normal = apply(k_count[,73:144],1,mean)


  # Mean and dispersion values from cancer samples
  dge.cancer=DGEList(counts=k_count[,1:72], group=factor(rep(1,72)))
  dge.cancer=calcNormFactors(dge.cancer)
  dge.cancer=estimateCommonDisp(dge.cancer)
  dge.cancer=estimateTagwiseDisp(dge.cancer)
  disp.cancer = dge.cancer$tagwise.dispersion
  mean.cancer = apply(k_count[,1:72],1,mean)


  # Gene filtering: Genes having small read count (<10) are filtered
  k_mean.total = apply(k_count,1,mean)
  k_index.filter = which(k_mean.total < 10)
  k_mean.total = k_mean.total[-k_index.filter]
  disp.normal = disp.normal[-k_index.filter]
  disp.cancer = disp.cancer[-k_index.filter]
  mean.normal = mean.normal[-k_index.filter]
  mean.cancer = mean.cancer[-k_index.filter]


  # Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
  k_dge.total = DGEList(counts = k_count, group = factor(c(rep(1,72),rep(2,72))))
  k_dge.total = calcNormFactors(k_dge.total)
  k_dge.total = estimateCommonDisp(k_dge.total)
  k_dge.total = estimateTagwiseDisp(k_dge.total)
  k_disp.total = k_dge.total$tagwise.dispersion
  k_disp.total = k_disp.total[-k_index.filter]



  ###############################################
  #############bottomly###########################
  ################################################
  # Load Bottomly mouse RNA-seq data. 21 samples
  download.file("https://github.com/unistbig/compareDEtools/raw/7f7eecf301070a3b7c7c037da079b331b722622b/inst/extdata/bottomly_eset.RData", "bottomly_eset.RData")
  load("bottomly_eset.RData")
  b_count<-exprs(bottomly.eset) # RNA-seq read count data
  strain<-pData(bottomly.eset)[,'strain']
  index.C=which(strain=='C57BL/6J') # C57BL/6J sample index
  index.D=which(strain=='DBA/2J') # DBA/2J sample index


  b_count= b_count[,c(index.C, index.D)] # Arrange samples for convinience


  # Get count mean and dispersion using edgeR package


  # Mean and dispersion values are obtained separately from C57BL/6J and DBA/2J samples when different dispersion is assummed between two sample types.
  # Mean and dispersion values from C57BL/6J samples
  dge.C=DGEList(counts=b_count[,1:10], group=factor(rep(1,10)))
  dge.C=calcNormFactors(dge.C)
  dge.C=estimateCommonDisp(dge.C)
  dge.C=estimateTagwiseDisp(dge.C)
  disp.C = dge.C$tagwise.dispersion
  mean.C = apply(b_count[,1:10],1,mean)


  # Mean and dispersion values from DBA/2J samples
  dge.D=DGEList(counts=b_count[,11:21], group = factor(rep(2,11)))
  dge.D=calcNormFactors(dge.D)
  dge.D=estimateCommonDisp(dge.D)
  dge.D=estimateTagwiseDisp(dge.D)
  disp.D = dge.D$tagwise.dispersion # Dispersion
  mean.D = apply(b_count[,11:21],1,mean)


  # Gene filtering: Genes having small read count (<10) are filtered
  b_mean.total = apply(b_count,1,mean)
  b_index.filter = which(b_mean.total < 10)
  b_mean.total = b_mean.total[-b_index.filter]
  disp.C = disp.C[-b_index.filter]
  disp.D = disp.D[-b_index.filter]
  mean.C = mean.C[-b_index.filter]
  mean.D = mean.D[-b_index.filter]

  # Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
  b_dge.total = DGEList(counts = b_count, group = factor(c(rep(1,10),rep(2,11))))
  b_dge.total = calcNormFactors(b_dge.total)
  b_dge.total = estimateCommonDisp(b_dge.total)
  b_dge.total = estimateTagwiseDisp(b_dge.total)
  b_disp.total = b_dge.total$tagwise.dispersion
  b_disp.total = b_disp.total[-b_index.filter]


  ###############################################
  ########Modification for comPyDEtools##########
  ###############################################

  # Instead of returning the function's output as dataset.parameters, modify it using comPyDEtools to output as CSV.
  k_params = data.frame(
    k_total_mean = k_mean.total,
    k_total_disp = k_disp.total,
    k_cancer_mean = mean.cancer,
    k_cancer_disp = disp.cancer,
    k_normal_mean = mean.normal,
    k_normal_disp = disp.normal
  )
  write.csv(k_params, "../k_params.csv")
  b_params = data.frame(
    b_total_mean = b_mean.total,
    b_total_disp = b_disp.total,
    b_C_mean = mean.C,
    b_C_disp = disp.C,
    b_D_mean = mean.D,
    b_D_disp = disp.D
  )
  write.csv(b_params, "../b_params.csv")
  # As an extra, also output the counts of KIRC and Bottomly.
  write.csv(k_count[-k_index.filter,], "./k_count.csv")
  write.csv(b_count[-b_index.filter,], "./b_count.csv")
}

generateDatasetParameter()
