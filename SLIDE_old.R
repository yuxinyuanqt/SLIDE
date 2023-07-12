# This function is to estimate pairwise linkage disequilibrium with EM algorithm

ld.em<-function(genotype)
{
  
  n_site<-ncol(genotype)
  
  dat<-genotype
  n<-nrow(dat)
  
  LD0<-LD1<-LD2<-diag((colSums(dat)/2/n)*(1-(colSums(dat)/2/n)))
  for(i in 1:(n_site-1))
  {
    for(j in c(i:n_site)[which(c(i:n_site)!=i)])
    {
      p_A<-(2*n-sum(dat[,i]))/2/n
      p_B<-(2*n-sum(dat[,j]))/2/n
      p_AB<-p_A*p_B
      p_Ab<-p_A*(1-p_B)
      p_aB<-(1-p_A)*p_B
      p_ab<-(1-p_A)*(1-p_B)
      
      n_AABB<-length(which(dat[,i]==0 & dat[,j]==0))
      n_AABb<-length(which(dat[,i]==0 & dat[,j]==1))
      n_AaBB<-length(which(dat[,i]==1 & dat[,j]==0))
      n_AaBb<-length(which(dat[,i]==1 & dat[,j]==1))
      n_AAbb<-length(which(dat[,i]==0 & dat[,j]==2))
      n_Aabb<-length(which(dat[,i]==1 & dat[,j]==2))
      n_aabb<-length(which(dat[,i]==2 & dat[,j]==2))
      n_aaBb<-length(which(dat[,i]==2 & dat[,j]==1))
      n_aaBB<-length(which(dat[,i]==2 & dat[,j]==0))
      
      p_AB_i<-p_AB
      p_ab_i<-p_ab
      p_Ab_i<-p_Ab
      p_aB_i<-p_aB
      
      for(t in 1:100)
      {
        Full_prob <- ifelse((p_AB_i[t]*p_ab_i[t]+p_Ab_i[t]*p_aB_i[t])==0,1,p_AB_i[t]*p_ab_i[t]+p_Ab_i[t]*p_aB_i[t])
        a<-p_AB_i[t]*p_ab_i[t]/Full_prob
        
        b<-p_Ab_i[t]*p_aB_i[t]/Full_prob
        
        p_AB_i<-rbind(p_AB_i,(2*n_AABB+n_AABb+n_AaBB+n_AaBb*a)/2/n)
        
        p_ab_i<-rbind(p_ab_i,(2*n_aabb+n_aaBb+n_Aabb+n_AaBb*a)/2/n)
        
        p_Ab_i<-rbind(p_Ab_i,(2*n_AAbb+n_AABb+n_Aabb+n_AaBb*b)/2/n)
        
        p_aB_i<-rbind(p_aB_i,(2*n_aaBB+n_aaBb+n_AaBB+n_AaBb*b)/2/n)
        
        aa<-max(abs(c(p_AB_i[t+1]-p_AB_i[t],p_Ab_i[t+1]-p_Ab_i[t],p_aB_i[t+1]-p_aB_i[t],p_ab_i[t+1]-p_ab_i[t])))
        
        if(aa<=0.00001)   break
      }
      LD0[i,j]<-p_AB_i[t+1]*p_ab_i[t+1]-p_Ab_i[t+1]*p_aB_i[t+1]
      LD0[j,i]<-LD0[i,j]
      
      LD1[i,j]<-LD1[j,i]<-abs(LD0[i,j])/(min(p_A*(1-p_B),(1-p_A)*p_B)*sign(LD0[i,j]>0)+min((1-p_A)*(1-p_B), p_A*p_B)*sign(LD0[i,j]<0)+sign(LD0[i,j]==0))
      
      LD2[i,j]<- LD2[j,i]<-LD0[i,j]^2/(p_A*(1-p_B)*(1-p_A)*p_B)
      
    }
  }
  list(LD0,LD1,LD2)
}

## The permutation for case-control studies

permuC<-function(phenotype,genotype)
{
  n<-nrow(genotype)
  size1<-sum(phenotype==1)
  size2<-sum(phenotype==0)
  fa_per<-genotype[sample(c(1:n),n),]
  fa_per<-cbind(matrix(c(rep(1,size1),rep(0,size2)),ncol=1),fa_per)
  
  fa_per
}

# This function is to estimate pairwise linkage disequilibrium with EM algorithm

require(MASS)

SLIDE.old <- function(pheno_geno,method,LD=NULL,num_per=200)
{
  
  n_gv<-ncol(pheno_geno)-1
  
  fa_case<-pheno_geno[which(pheno_geno[,1]==1),-1]
  
  n_case<-nrow(fa_case)
  
  fa_con<-pheno_geno[which(pheno_geno[,1]==0),-1]
  n_con<-nrow(fa_con)
  
  U<-(colSums(fa_case)/n_case-colSums(fa_con)/n_con)/(1/n_case+1/n_con)
  
  
  if(method=="Cov") V1<-cov(fa_con)   else V1<-2*ld.em(fa_con)[[1]]
  
  LD_matrix<-list()
  
  LD_matrix[[1]]<-V1*n_case*n_con/(n_case+n_con)
  
  if(length(LD)>0) LD_matrix[[2]]<-2*LD*n_case*n_con/(n_case+n_con)
  
  
  LD_inverse<-list()
  
  TT<-NULL
  
  for(ii in 1:length(LD_matrix))
  {
    LD_inverse[[ii]]<-ginv(LD_matrix[[ii]])
    
    TT<-c(TT,U%*%LD_inverse[[ii]]%*%U)
    
  }
  
  for(B in 1:num_per)
  {
    
    fa_per<-permuC(pheno_geno[,1],pheno_geno[,-1])
    
    fa_case<-fa_per[which(fa_per[,1]==1),-1]
    
    fa_con<-fa_per[which(fa_per[,1]==0),-1]
    
    U<-(colSums(fa_case)/n_case-colSums(fa_con)/n_con)/(1/n_case+1/n_con)
    
    tt0<-NULL
    
    for(ii in 1:length(LD_inverse)) tt0<-c(tt0,U%*%LD_inverse[[ii]]%*%U)
    
    
    TT<-rbind(TT,tt0)
  }
  
  all<-(num_per+1-apply(TT,2,rank))/num_per
  
  if(length(LD)>0) p_asy_out_cov<-1-pchisq(TT[1,2],n_gv) else p_asy_out_cov<-NULL
  
  p_asy_cov<-1-pchisq(TT[1,1],n_gv)
  
  p.value<-c(all[1,],p_asy_cov,p_asy_out_cov)
  
  if(length(p.value)>2) {
    
    SLIDE.tem<-c(TT[1,1],p.value[c(1,3)]); names(SLIDE.tem)<-c("SLIDE","p.value.per","p.value.asy")
    SLIDE.pop<-c(TT[1,2],p.value[c(1,4)]); names(SLIDE.pop)<-c("SLIDE.pop","p.value.pop.per","p.value.pop.asy")
    
    statis<-list(SLIDE.tem,SLIDE.pop)
  } else {
    SLIDE.tem<-c(TT[1,1],p.value); names(SLIDE.tem)<-c("SLIDE","p.value.per","p.value.asy")
    statis<-list(SLIDE.tem)
  }
  
  
  return(statis)
}
