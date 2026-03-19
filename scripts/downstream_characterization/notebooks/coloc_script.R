library(coloc)
args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1])
bf1 <- setNames(df$logBF_1, df$snp)
bf2 <- setNames(df$logBF_2, df$snp)
my.res <- coloc.bf_bf(bf1 =  bf1 , bf2 =  bf2,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
  )
write.csv(my.res[[1]], args[2])