library(CMSclassifier)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

input = read.table(input_file,sep=',',header=TRUE,row.names=1,check.names=FALSE)

result = CMSclassifier::classifyCMS(input,method="SSP")

colnames(result$nearestCMS)<-c("nearestCMS")
colnames(result$predictedCMS)<-c("predictedCMS")

combined = cbind(result$predictedCMS,result$nearestCMS,result$SSP.details)

write.table(combined, output_file,sep='\t')