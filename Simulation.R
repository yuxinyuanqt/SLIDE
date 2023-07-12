##############################################################################

#ʱ��:2021��8��20��

#��״���������ݡ���������״��
#�ص㣺���Ǻ�������ͳ�����������ЧӦ
#Ŀ�ģ����� X Ⱦɫ���ϵĻ���������Ƿ��뼲�����������������������ܼ��г����������к������죻
#      ���Ƚϲ�ͬ��ֵ��ʽ��һ�����ͼ���Ч��

#��Ⱦɫ����ģ�⣬����SLIDE���µ�ģ��˼·

#################################���������R��################################

library(plyr)
library(tidyverse)
library(mvtnorm) #Q-Zmax����R��
library(expm) #sqrtm {expm}: ��ƽ��������ĺ���
#���м���
library(foreach)
library(doParallel)
library(snow)
library(reshape2)
library(RColorBrewer) #ȡɫ��
library(SKAT)
library(MASS)


###################################���庯��##################################

#����һ���õ����ַ�����һ�������߼���Ч�ܵĺ���
T1E_power <- function(seed,n_pool,n_ld_est,n,
                      n_rare,n_common,
                      beta_g,beta_0,
                      Nsim,rho=0,alpha=0.05){
  
  #seed: ������; 
  #n_pool: ���������������; #n_ld_est: ��������LD������ⲿ������������
  #n_case: ������; n_control: ��������
  #n: k*2�����ݿ�����������ֱ�Ϊ�������Ͷ��������ܵ�ȡֵ��
  
  #n_rare: ��������λ������n_common: ��������λ������
  #Ϊ����beta_g���Ӧ��n_rare��n_common�����Ƴɱ�����

  #beta_g: ��������ͳ�������λ������Ӧ�Ļع�ϵ������ɵ����ݿ�,ÿһ�д����ع�ϵ�����ܵ�һ��ȡֵ;
  #����ʵ��ģ���У���������ͳ�������Ļع�ϵ��ͨ������Գ��ֵ�,
  #�ʽ�beta_g����Ϊһ��k*(n_rare+n_common)�����ݿ�,
  #beta_gÿһ�е�ǰn_rare��Ԫ���Ǻ�������λ���Ӧ�Ļع�ϵ������n_common��Ԫ���ǳ�������λ���Ӧ�Ļع�ϵ����
  #��beta_g��ĳһ�о�Ϊ0��������ģ��һ����󣬷���ģ�����Ч�ܣ�
  
  #beta_0:�ؾ���Ļع�ϵ��;
  #Nsim: ģ�����; 
  #rho: AR(1)ģ������Ӧ�����ϵ��, ����ģ�����Ȥ������λ��֮���LD��Ϣ��
  #     ȡֵ��0~1֮��,Ĭ��Ϊ0;
  #alpha: ����һ�����. Ĭ��Ϊ0.05;
  
  #��n��n_rare��n_common��beta_rare��beta_common��alpha�⣬�����������Ϊ������������ÿһ��Ԫ�ش���һ�ο��ܵ�ȡֵ;
  
  
  if(is.vector(n)){
    n <- matrix(n,nrow=1)}
  colnames(n) <- c('n_case','n_control') #Ϊp����������ֱ�����Ϊn_case��n_control
  
  #���beta_g����������ֻ��һ�ֻع�ϵ��ȡֵ����beta_Gת��Ϊ1*(n_rare+n_common)�ľ���
  if(is.vector(beta_g)){
    beta_g <- matrix(beta_g,nrow=1)}

    
  
  #####����һЩ�Ӻ���
  
  #1.����һ���������ݿ�ĵѿ������ĺ����������� expand.grid �����������ö��������ݿ�
  expand.grid.df <- function(...) 
    Reduce(function(...) merge(..., by=NULL), list(...))
  
  
  
  #��������ģ����������
  situation <- expand.grid(seed=seed,n_pool=n_pool,
                           n_rare=n_rare,n_common=n_common,
                           beta_0=beta_0,Nsim=Nsim,rho=rho) %>%
    expand.grid.df(n,beta_g) %>%
    dplyr::mutate(
      id=1:n(),
      seed=paste(seed,id,sep = "")
    )
  
  #����ÿ�ֲ�������¸��ַ���ÿ��ģ���pֵ
  TestRes <- plyr::ddply(situation, .(id) ,function(sit){
    
    set.seed(sit$seed)
    
    #��������ͺ�������λ����
    n_loci <- sit$n_rare+sit$n_common
    
    #beta_g: λ��Ļع�ϵ������
    beta_g <- unlist(sit[paste('V',1:n_loci,sep = "")])
    
    #TF_causal: �Ƿ������λ�㣬����һ��������
    #           ����λ����������������Ӧ������Ԫ��ֵΪ1������ֵΪ0
    TF_causal <- ifelse(beta_g==0,0,1)
    
    #Cov_Z:Ԥ�ȸ�n_loci*n_loci����̬�ֲ��������Z�ķ���-Э�������������ڴ�ռ䣬֮���������Cov_Z���ɣ�
    #       ������ִ��ѭ���Ĺ����У����Ա�����ж���ڴ����
    #Z: ��̬�ֲ�����������������ɻ�����
    Cov_Z <- matrix(nrow = n_loci,ncol = n_loci)
    
    #ΪCov_Z�Խ�����Ԫ�ظ�ֵ
    diag(Cov_Z) <- 1
  
    #���ǵ�Cov_Z��һ���Գƾ���ֻ�迼�Ǿ���������ǲ��ּ���
    for (i in 1:(n_loci-1)) {
      #Indicative: ʾ�Ժ�����������i�����(i+1):n_lociλ���Ƿ��Ϊ����������Ϊ���������
      #            ���ǣ����ӦԪ��ֵΪ1������ֵΪ0
      Indicative <- ifelse(TF_causal[(i+1):n_loci]==TF_causal[i],1,0)
      #��λ��i��j��Ϊ������������Ϊ�����������Cov_Z[i,j]ֵΪrho^(abs(i-j))������ֵΪ0
      Cov_Z[(i+1):n_loci,i] <- Cov_Z[i,(i+1):n_loci] <- Indicative*sit$rho^c(1:(n_loci-i))
    }
    
    #threshold: ���������������ֵ
    threshold <- 1/sqrt(2*(sit$n_case+sit$n_control))
    #MAF: λ��Ĵε�λ����Ƶ�ʣ���������Ĵε�λ����Ƶ�ʴ�U(0.001,threshold)���������
    #     ��������Ĵε�λ����Ƶ�ʴ�U(threshold,0.5)���������
    #qnorm_MAF: MAF����Ӧ�ı�׼��̬�ֲ����·�λ��
    qnorm_MAF <- qnorm(c(runif(sit$n_rare,0.001,threshold),runif(sit$n_common,threshold,0.5)))
    
    #Genotype:����n_pool*n_loci�Ļ����;������������;������
    #ʹ�ùܵ����η�����ʡ�ڴ�
    Genotype <- t(MASS::mvrnorm(sit$n_pool,rep(0,n_loci),Cov_Z)   %>%
                      apply(1,function(x){x <= qnorm_MAF}))+
      t(MASS::mvrnorm(sit$n_pool,rep(0,n_loci),Cov_Z)   %>%
          apply(1,function(x){x <= qnorm_MAF}))
    
    #logit_p: logit P(yi=1)�����廼�����ʵ�logit�任��ȡֵ��-Inf~Inf֮��
    logit_p <- sit$beta_0+Genotype%*%beta_g
    
    #Phenotype: ���ɱ��ͣ�����ΪP(y=1)=exp(logit_p)/(1+exp(logit_p))
    Phenotype <- rbinom(n=sit$n_pool, size=1, prob=exp(logit_p)/(1+exp(logit_p)))
    
    #index_case��index_control: �����Ͷ��յ�����
    index_case <- which(Phenotype==1)
    index_control <- which(Phenotype==0)
    
    #�ӳ������г�ȡn_ld_est�ĸ���Ļ����;���������LD����
    #n_ld_est: ��Ϊ����
    #LD_est: �б����б���ÿһ����������ӳ������г�ȡn_ld_est�ĸ���Ļ����;���������LD����
    LD_est <- lapply(n_ld_est, function(n_ld_esti){LD.f(Genotype[sample(1:sit$n_pool,n_ld_esti),])$LD_f0})
    
    map_dfr(1:sit$Nsim,function(x){
      
      #index:�������������������ڴӳ������г���
      index <- c(sample(index_case,sit$n_case),sample(index_control,sit$n_control))
      
      #Genotype_s: �����Ļ����;���
      Genotype_s <- Genotype[index,]
      #Phenotype_s: �����ı�������
      Phenotype_s <- Phenotype[index]
      
      # #�����������ݹ��Ƴ���LD�����Ӧ��Pֵ
      # SLIDE_tem_p <- SILDE(Genotype_s,Phenotype_s,LD_estimate='sample')[2:3]
      # #���������Ķ������ݹ��Ƴ���LD�����Ӧ��Pֵ
      # SLIDE_con_p <- SILDE(Genotype_s,Phenotype_s,LD_estimate='control')[2:3]
      # 
      # #�����ⲿ����Ⱥ���ݹ��Ƴ���LD�����Ӧ��Pֵ
      # SLIDE_pop_p <- as.vector(sapply(LD_est,function(LD_esti){
      #   SILDE(Genotype_s,Phenotype_s,LD_estimate='population',
      #         LD_matrix=LD_esti)[2:3]
      # }))
      
      SLIDE_p <- as.vector(SILDE(Genotype_s,Phenotype_s,LD_estimate='all',
                     LD_matrix=LD_est)[,-1])
      
      
      #��Ͽ�ģ��
      obj <- SKAT::SKAT_Null_Model(Phenotype_s ~ 1, out_type="D",
                                   Adjustment=FALSE)
      #SKATO_p: SKAT-O��Pֵ
      SKATO_p <- SKAT::SKAT(Genotype_s, obj,method="optimal.adj")$p.value
      
      #SKAT_CR_p: ���Ǻ�������ͳ�������, Combined Sum Test of Rare- and Common-Variant Effects
      SKAT_CR_p <- SKAT::SKAT_CommonRare(Genotype_s, obj, method="C")$p.value
      Burden_CR_p <- SKAT::SKAT_CommonRare(Genotype_s, obj,r.corr.rare=1,
                                           r.corr.common=1, method="C")$p.value
      
      # #SKATO_R_p: ������������Ӧ��SKAT-O Pֵ
      # SKATO_R_p <- SKAT::SKAT(Genotype_s[,1:sit$n_rare], obj,method="optimal.adj")$p.value
      # #SKATO_C_p: ������������Ӧ��SKAT-O Pֵ
      # SKATO_C_p <- SKAT::SKAT(Genotype_s[,-(1:sit$n_rare)],  
      #                   weights.beta=c(0.5,0.5),obj,method="optimal.adj")$p.value
      # 
      # #��������ͳ��������Pֵ���п����ϲ�
      # SKATO_C <- 1/2-atan(sum(tan((0.5-c(SKATO_R_p,SKATO_C_p))*pi))/2)/pi
      
      # P_data <- c(SLIDE_p,SKATO_p,SKAT_CR_p,SKATO_R_p,SKATO_C_p,SKATO_C)
      # names(P_data) <- c('SLIDE_tem_asy','SLIDE_pop_asy',
      #                    'SLIDE_tem_per','SLIDE_pop_per',
      #                    'SKATO_p','SKAT_CR_p','SKATO_R_p',
      #                    'SKATO_C_p','SKATO_C')
      
      P_data <- c(SLIDE_p,SKATO_p,SKAT_CR_p,Burden_CR_p)
      names(P_data) <- c(paste(rep(c(paste('SLIDE_pop',n_ld_est,sep = ""),
                                     'SLIDE_tem','SLIDE_con'),2),
                               rep(c("asy","per"),each=length(n_ld_est)+2),
                               sep = "_"),
                         'SKATO','SKAT_CR','Burden_CR')
      
      return(P_data)
      
    })  %>%
      dplyr::bind_rows()
  })
  
  #������ֲ�������µ� Type I error ���߼���Ч��
  T1E_power <- plyr::ddply(TestRes,.(id),function(df){
    T1E_power <- apply(df[-1], 2, function(x) mean(x<alpha,na.rm=TRUE))
  }) %>%
    left_join(situation,by="id")
  return(T1E_power)

} 


####################���ú���������һ�������߼���Ч��#######################


seed <- 210825
n_pool <- 500000
n_ld_est <- 10000
n <- cbind(c(300,600),c(600,600))
n <- c(600,600)
n_rare <- 20
n_common <- 12
beta_g <- rbind(c(log(c(2,1/2,2,1/2,2,1/2,1.5,1.5)),rep(0,12),log(c(1.1,1/1.1,1.1,1.1)),rep(0,8)),
                rep(0,32),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2)),rep(0,12),log(c(1.1,1/1.1,1.1,1/1.1)),rep(0,8)))

beta_g <- rbind(c(log(c(2,1/2,2,1/2,2,1/2,1.5,1.5)),rep(0,12),log(c(1.1,1/1.1,1.1,1.1)),rep(0,8)),
                rep(0,32))
beta_0 <- log(0.1)
Nsim <- 1000
#rho <- c(0,0.3,0.5,0.7,0.9)
rho <- c(0,0.3,0.5,0.7,0.9)


##����������

cl <- makeCluster(detectCores()/2,type="SOCK") ## �����ض�cpu���ĵļ�Ⱥ
registerDoParallel(cl) ## ע�ἯȺ

# clusterExport(cl,)
clusterEvalQ(cl, {library(tidyverse);library(coin)}) ## ���ذ���cl��Ⱥ

#����һ���������Ч��
suppressWarnings(myT1E_power1 <- T1E_power(seed,n_pool,n_ld_est,n,
                                           n_rare,n_common,
                                           beta_g,beta_0,
                                           Nsim,rho,alpha=0.05))

#����
newmyT1E_power1 <- myT1E_power1[order(myT1E_power1$theta,
                                      myT1E_power1$beta_g,
                                      myT1E_power1$ratio,
                                      myT1E_power1$pm,
                                      myT1E_power1$pf,
                                      myT1E_power1$gamma),]

##ֹͣ��Ⱥ���ͷ��ڴ�  
suppressWarnings({
  stopImplicitCluster()
  stopCluster(cl)
})


seed <- 210827
n_pool <- 500000
n_ld_est <- 20000
n <- cbind(c(300,600),c(600,600))
beta_0 <- log(0.1)
Nsim <- 1000
rho <- c(0,0.3,0.5,0.7,0.9)

################################Scenario 1-3

n_rare <- 32
n_common <- 0
beta_g <- rbind(c(log(c(3.5,1/3,3.5,1/3)),rep(0,28)),
                c(log(c(3.5,1/3.5,3.5,1/3.5)),rep(0,28)),
                c(log(c(3,1/3,3,1/3,2,1/2,2,2)),rep(0,24)),
                c(log(c(3,1/3,3,1/3,2.5,1/2.5,2,1/2)),rep(0,24)),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2,2,1/2,2,2)),rep(0,20)),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2,2,1/2,2,1/2)),rep(0,20)),
                rep(0,32))

system.time(myT1E_power1_3 <- T1E_power(seed,n_pool,n_ld_est,n,
                          n_rare,n_common,
                          beta_g,beta_0,
                          Nsim,rho,alpha=0.05))

################################Scenario 4-5

n_rare <- 22
n_common <- 10
beta_g <- rbind(c(log(c(2,1/2,2,1/2,2,1/2,1.5,1.5,2,1/2,2,1/2)),rep(0,20)),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2,2,1/2,2,1/2)),rep(0,20)),
                c(log(c(2,1/2,2,1/2,1.5,1/1.5,1.5,1/1.5,1.5,1.5)),rep(0,12),log(c(1/1.1,1.15)),rep(0,8)),
                c(log(c(2.5,1/2.5,2,1/2,2,1/2,1.5,1/1.5,1.5,1/1.5)),rep(0,12),log(c(1/1.1,1.1)),rep(0,8)),
                rep(0,32))


system.time(myT1E_power4_5 <- T1E_power(seed,n_pool,n_ld_est,n,
                                      n_rare,n_common,
                                      beta_g,beta_0,
                                      Nsim,rho,alpha=0.05))

################################Scenario 6


n_rare <- 18
n_common <- 14
beta_g <- rbind(c(log(c(2,1/2,2,1/2,1.5,1/1.5,1.5,1/1.5,1.5,1.5)),rep(0,8),log(c(1/1.1,1.15)),rep(0,12)),
                c(log(c(2.5,1/2.5,2,1/2,2,1/2,1.5,1/1.5,1.5,1/1.5)),rep(0,8),log(c(1/1.1,1.1)),rep(0,12)),
                rep(0,32))


system.time(myT1E_power6 <- T1E_power(seed,n_pool,n_ld_est,n,
                                        n_rare,n_common,
                                        beta_g,beta_0,
                                        Nsim,rho,alpha=0.05))

################################Scenario 7

n_rare <- 20
n_common <- 12
beta_g <- rbind(c(log(c(2,1/2,2,1/2,2,1/2,1.5,1.5)),rep(0,12),log(c(1.1,1/1.1,1.1,1.1)),rep(0,8)),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2)),rep(0,12),log(c(1.1,1/1.1,1.1,1/1.1)),rep(0,8)),
                rep(0,32))

system.time(myT1E_power7 <- T1E_power(seed,n_pool,n_ld_est,n,
                                      n_rare,n_common,
                                      beta_g,beta_0,
                                      Nsim,rho,alpha=0.05))

################################Scenario 8

n_rare <- 16
n_common <- 16
beta_g <- rbind(c(log(c(2,1/2,2,1/2,2,1/2,1.5,1.5)),rep(0,8),log(c(1.1,1/1.1,1.1,1.1)),rep(0,12)),
                c(log(c(2,1/2,2,1/2,2,1/2,2,1/2)),rep(0,8),log(c(1.1,1/1.1,1.1,1/1.1)),rep(0,12)),
                rep(0,32))


system.time(myT1E_power8 <- T1E_power(seed,n_pool,n_ld_est,n,
                                      n_rare,n_common,
                                      beta_g,beta_0,
                                      Nsim,rho,alpha=0.05))

################################Scenario 9

n_rare <- 0
n_common <- 32
beta_g <- rbind(c(log(c(1.3,1/1.3,1.3,1/1.2)),rep(0,28)),
                c(log(c(1.3,1/1.3,1.3,1/1.3)),rep(0,28)),
                rep(0,32))

system.time(myT1E_power9 <- T1E_power(seed,n_pool,n_ld_est,n,
                                      n_rare,n_common,
                                      beta_g,beta_0,
                                      Nsim,rho,alpha=0.05))