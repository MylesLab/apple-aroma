# Title     : Utility Functions
# Objective : This script contains all the functions used by this module
# Created by: tayabsoomro
# Created on: 2021-07-27

# R script to read the GRM binary file
ReadGRMBin=function(test, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(test,".grm.bin",sep="")
  NFileName=paste(test,".grm.N.bin",sep="")
  IDFileName=paste(test,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
