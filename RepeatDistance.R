
# Hamming distance between character arrays
hamming <- function(X, Y) {
if ( missing(Y) ) {
uniqs <- unique(as.vector(X))
U <- X == uniqs[1]
H <- t(U) %*% U
for ( uniq in uniqs[-1] ) {
U <- X == uniq
H <- H + t(U) %*% U
}
} else {
uniqs <- union(X, Y)
H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
for ( uniq in uniqs[-1] ) {
H <- H + t(X == uniq) %*% (Y == uniq)
}
}
nrow(X) - H
}
which_min <- function(x){
which(x == min(x, na.rm = TRUE))
}

# pairwise repeat alignment based on mutations, single repeat indels and single slippage
#
# Parameters:
# seq1: aligned unit sequences for first repeat
# seq2: aligned unit sequences for second repeat
# gendist: pre-computed genetic distance between all units in seq1 and seq2 (default: Hamming distance)
# wmut: weight for the genetic distance between units
# windel: weight for unit insertions/deletions
# wslippage: weight for single unit duplication/slippage
# output: could be the value of the final genetic distance between repeats ("distance"), or the pairwise alignment in terms of sequences ("sequence") or in terms of ordered units ("order"), or a list containing all these info ("all")

RepeatAlignment<-function(seq1, seq2, gendist=NULL, wmut, windel, wslippage, output=c("distance","order","sequence","all")){
seq1<-as.character(seq1)
seq2<-as.character(seq2)
m<-dim(seq1)[1]
n<-dim(seq2)[1]
l<-dim(seq1)[2]
stopifnot(l==dim(seq2)[2])
if(!is.null(gendist)){stopifnot(dim(gendist)[1]==m+n,dim(gendist)[2]==m+n)} else {gendist=hamming(t(rbind(seq1,seq2)))}
dm<-(gendist[1:m,(m+1):(m+n)]+t(gendist[(m+1):(m+n),1:m]))/2
ds<-gendist[1:m,1:m]
dt<-gendist[(m+1):(m+n),(m+1):(m+n)]
d<-matrix(0,nrow=m+1,ncol=n+1)
mfrom<-matrix(rep(list(),(m+1)*(n+1)),nrow=m+1,ncol=n+1)
d[2,1]<-windel
mfrom[2,1]<-list(c(3))
for(i in c(2:m)){
s_del<-d[i,1]+windel
s_ins<-s_del+1
s_delslip<-d[i,1]+wslippage+wmut*min(ds[i,i-1],ds[i,i+1-2*(i==m)])
s_insslip<-s_delslip+1
s_mut<-s_del+1
d[i+1,1]<-min(s_mut,s_ins,s_del,s_insslip,s_delslip)
mfrom[i+1,1]<-list(which_min(c(s_mut,s_ins,s_del,s_insslip,s_delslip)))
}
d[1,2]<-windel
mfrom[1,2]<-list(c(2))
for(j in c(2:n)){
s_ins<-d[1,j]+windel
s_del<-s_ins+1
s_insslip<-d[1,j]+wslippage+wmut*min(dt[j,j-1],dt[j,j+1-2*(j==n)])
s_delslip<-s_insslip+1
s_mut<-s_ins+1
d[1,j+1]<-min(s_mut,s_ins,s_del,s_insslip,s_delslip)
mfrom[1,j+1]<-list(which_min(c(s_mut,s_ins,s_del,s_insslip,s_delslip)))# CONVERT TO LIST
}
for (j in c(1:n)){
for (i in c(1:m)){
s_del <- d[i,j+1]+windel # deletion
s_ins <- d[i+1,j]+windel # an insertion
s_delslip <- d[i,j+1]+wslippage+wmut*ds[i,i-1] # slippage on the other
s_insslip <- d[i+1,j]+wslippage+wmut*dt[j,j-1] # slippage
s_mut <- d[i,j]+dm[i,j]*wmut  # a substitution
d[i+1,j+1]<-min(s_mut,s_ins,s_del,s_insslip,s_delslip)
mfrom[i+1,j+1]<-list(which_min(c(s_mut,s_ins,s_del,s_insslip,s_delslip)))
}
}
result1<-c()
result2<-c()
order1<-c()
order2<-c()
mdist<-d[m+1,n+1]
if(output=="distance"){return(mdist)}
#print(mdist)
#print(d)
i<-m+1
j<-n+1
while(i>1 | j>1){
#print(mfrom[i,j][[1]])
if(length(mfrom[i,j][[1]])==1){orig<-mfrom[i,j][[1]]} else {orig<-sample(mfrom[i,j][[1]],1)}
#print(orig)
if(orig==1){
result1<-c((seq1[i-1,]),result1)
result2<-c((seq2[j-1,]),result2)
order1<-c(i-1,order1)
order2<-c(j-1,order2)
i<-i-1
j<-j-1
} else {
if (orig==2 | orig==4){
result1<-c(rep("-",l),result1)
result2<-c((seq2[j-1,]),result2)
order1<-c(0,order1)
order2<-c(j-1,order2)
j<-j-1
} else {
result1<-c((seq1[i-1,]),result1)
result2<-c(rep("-",l),result2)
order1<-c(i-1,order1)
order2<-c(0,order2)
i<-i-1
}
}
}
names(result1)<-NULL
names(result2)<-NULL
names(order1)<-NULL
names(order2)<-NULL
#print(result1)
#print(result2)
result<-rbind(result1,result2)
rownames(result)<-c("seq1","seq2")
orders<-rbind(order1,order2)
rownames(orders)<-c("order1","order2")
if(output=="sequence"){return(result)}
if(output=="order"){return(orders)}
results<-list(distance=mdist,order=orders,sequence=result)
return(results)
}
