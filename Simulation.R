##############################################################################

#时间:2021年8月20日

#性状：独立数据、二分类性状；
#特点：考虑罕见变异和常见变异联合效应
#目的：检验 X 染色体上的基因或区域是否与疾病相关联，这个基因或区域可能既有常见变异又有罕见变异；
#      并比较不同赋值方式的一类错误和检验效能

#常染色体结果模拟，基于SLIDE文章的模拟思路

#################################加载所需的R包################################

library(plyr)
library(tidyverse)
library(mvtnorm) #Q-Zmax所需R包
library(expm) #sqrtm {expm}: 求平方根矩阵的函数
#并行计算
library(foreach)
library(doParallel)
library(snow)
library(reshape2)
library(RColorBrewer) #取色板
library(SKAT)
library(MASS)


###################################定义函数##################################

#定义一个得到各种方法的一类错误或者检验效能的函数
T1E_power <- function(seed,n_pool,n_ld_est,n,
                      n_rare,n_common,
                      beta_g,beta_0,
                      Nsim,rho=0,alpha=0.05){
  
  #seed: 种子数; 
  #n_pool: 抽样池里的样本量; #n_ld_est: 用来估计LD矩阵的外部数据样本量；
  #n_case: 病例数; n_control: 对照数；
  #n: k*2的数据框，其两个组件分别为病例数和对照数可能的取值；
  
  #n_rare: 罕见变异位点数；n_common: 常见变异位点数；
  #为了与beta_g相对应，n_rare和n_common均限制成标量；

  #beta_g: 罕见变异和常见变异位点所对应的回归系数所组成的数据框,每一行代表回归系数可能的一次取值;
  #由于实际模拟中，罕见变异和常见变异的回归系数通常是配对出现的,
  #故将beta_g定义为一个k*(n_rare+n_common)的数据框,
  #beta_g每一行的前n_rare个元素是罕见变异位点对应的回归系数，后n_common个元素是常见变异位点对应的回归系数，
  #若beta_g的某一行均为0向量，则模拟一类错误，否则模拟检验效能；
  
  #beta_0:截距项的回归系数;
  #Nsim: 模拟次数; 
  #rho: AR(1)模型所对应的相关系数, 用于模拟感兴趣区域中位点之间的LD信息，
  #     取值在0~1之间,默认为0;
  #alpha: 名义一类错误. 默认为0.05;
  
  #除n、n_rare、n_common、beta_rare、beta_common和alpha外，其余变量均可为向量，向量的每一个元素代表一次可能的取值;
  
  
  if(is.vector(n)){
    n <- matrix(n,nrow=1)}
  colnames(n) <- c('n_case','n_control') #为p的两个组件分别命名为n_case和n_control
  
  #如果beta_g是向量，即只有一种回归系数取值，则将beta_G转化为1*(n_rare+n_common)的矩阵
  if(is.vector(beta_g)){
    beta_g <- matrix(beta_g,nrow=1)}

    
  
  #####定义一些子函数
  
  #1.定义一个生成数据框的笛卡尔积的函数，类似于 expand.grid 函数，但作用对象是数据框
  expand.grid.df <- function(...) 
    Reduce(function(...) merge(..., by=NULL), list(...))
  
  
  
  #生成所有模拟参数的组合
  situation <- expand.grid(seed=seed,n_pool=n_pool,
                           n_rare=n_rare,n_common=n_common,
                           beta_0=beta_0,Nsim=Nsim,rho=rho) %>%
    expand.grid.df(n,beta_g) %>%
    dplyr::mutate(
      id=1:n(),
      seed=paste(seed,id,sep = "")
    )
  
  #计算每种参数组合下各种方法每次模拟的p值
  TestRes <- plyr::ddply(situation, .(id) ,function(sit){
    
    set.seed(sit$seed)
    
    #常见变异和罕见变异位点数
    n_loci <- sit$n_rare+sit$n_common
    
    #beta_g: 位点的回归系数向量
    beta_g <- unlist(sit[paste('V',1:n_loci,sep = "")])
    
    #TF_causal: 是否是因果位点，返回一个向量，
    #           若该位点是因果变量，其对应的向量元素值为1，否则值为0
    TF_causal <- ifelse(beta_g==0,0,1)
    
    #Cov_Z:预先给n_loci*n_loci的正态分布随机向量Z的方差-协方差矩阵分配了内存空间，之后依次填充Cov_Z即可，
    #       这样在执行循环的过程中，可以避免进行多次内存分配
    #Z: 正态分布随机向量，用于生成基因型
    Cov_Z <- matrix(nrow = n_loci,ncol = n_loci)
    
    #为Cov_Z对角线上元素赋值
    diag(Cov_Z) <- 1
  
    #考虑到Cov_Z是一个对称矩阵，只需考虑矩阵的上三角部分即可
    for (i in 1:(n_loci-1)) {
      #Indicative: 示性函数向量，第i个与第(i+1):n_loci位点是否均为因果变量或均为非因果变量
      #            若是，则对应元素值为1，否则值为0
      Indicative <- ifelse(TF_causal[(i+1):n_loci]==TF_causal[i],1,0)
      #若位点i和j均为因果变量，或均为非因果变量，Cov_Z[i,j]值为rho^(abs(i-j))，否则值为0
      Cov_Z[(i+1):n_loci,i] <- Cov_Z[i,(i+1):n_loci] <- Indicative*sit$rho^c(1:(n_loci-i))
    }
    
    #threshold: 罕见变异的理论阈值
    threshold <- 1/sqrt(2*(sit$n_case+sit$n_control))
    #MAF: 位点的次等位基因频率，罕见变异的次等位基因频率从U(0.001,threshold)中随机产生
    #     常见变异的次等位基因频率从U(threshold,0.5)中随机产生
    #qnorm_MAF: MAF所对应的标准正态分布的下分位数
    qnorm_MAF <- qnorm(c(runif(sit$n_rare,0.001,threshold),runif(sit$n_common,threshold,0.5)))
    
    #Genotype:生成n_pool*n_loci的基因型矩阵，两条单倍型矩阵相加
    #使用管道修饰符，节省内存
    Genotype <- t(MASS::mvrnorm(sit$n_pool,rep(0,n_loci),Cov_Z)   %>%
                      apply(1,function(x){x <= qnorm_MAF}))+
      t(MASS::mvrnorm(sit$n_pool,rep(0,n_loci),Cov_Z)   %>%
          apply(1,function(x){x <= qnorm_MAF}))
    
    #logit_p: logit P(yi=1)，个体患病概率的logit变换，取值在-Inf~Inf之间
    logit_p <- sit$beta_0+Genotype%*%beta_g
    
    #Phenotype: 生成表型，概率为P(y=1)=exp(logit_p)/(1+exp(logit_p))
    Phenotype <- rbinom(n=sit$n_pool, size=1, prob=exp(logit_p)/(1+exp(logit_p)))
    
    #index_case和index_control: 病例和对照的索引
    index_case <- which(Phenotype==1)
    index_control <- which(Phenotype==0)
    
    #从抽样池中抽取n_ld_est的个体的基因型矩阵来估计LD矩阵
    #n_ld_est: 可为向量
    #LD_est: 列表，列表的每一个组件代表从抽样池中抽取n_ld_est的个体的基因型矩阵来估计LD矩阵
    LD_est <- lapply(n_ld_est, function(n_ld_esti){LD.f(Genotype[sample(1:sit$n_pool,n_ld_esti),])$LD_f0})
    
    map_dfr(1:sit$Nsim,function(x){
      
      #index:样本的索引向量，用于从抽样池中抽样
      index <- c(sample(index_case,sit$n_case),sample(index_control,sit$n_control))
      
      #Genotype_s: 样本的基因型矩阵
      Genotype_s <- Genotype[index,]
      #Phenotype_s: 样本的表型向量
      Phenotype_s <- Phenotype[index]
      
      # #基于样本数据估计出的LD矩阵对应的P值
      # SLIDE_tem_p <- SILDE(Genotype_s,Phenotype_s,LD_estimate='sample')[2:3]
      # #基于样本的对照数据估计出的LD矩阵对应的P值
      # SLIDE_con_p <- SILDE(Genotype_s,Phenotype_s,LD_estimate='control')[2:3]
      # 
      # #基于外部的人群数据估计出的LD矩阵对应的P值
      # SLIDE_pop_p <- as.vector(sapply(LD_est,function(LD_esti){
      #   SILDE(Genotype_s,Phenotype_s,LD_estimate='population',
      #         LD_matrix=LD_esti)[2:3]
      # }))
      
      SLIDE_p <- as.vector(SILDE(Genotype_s,Phenotype_s,LD_estimate='all',
                     LD_matrix=LD_est)[,-1])
      
      
      #拟合空模型
      obj <- SKAT::SKAT_Null_Model(Phenotype_s ~ 1, out_type="D",
                                   Adjustment=FALSE)
      #SKATO_p: SKAT-O的P值
      SKATO_p <- SKAT::SKAT(Genotype_s, obj,method="optimal.adj")$p.value
      
      #SKAT_CR_p: 考虑罕见变异和常见变异, Combined Sum Test of Rare- and Common-Variant Effects
      SKAT_CR_p <- SKAT::SKAT_CommonRare(Genotype_s, obj, method="C")$p.value
      Burden_CR_p <- SKAT::SKAT_CommonRare(Genotype_s, obj,r.corr.rare=1,
                                           r.corr.common=1, method="C")$p.value
      
      # #SKATO_R_p: 罕见变异所对应的SKAT-O P值
      # SKATO_R_p <- SKAT::SKAT(Genotype_s[,1:sit$n_rare], obj,method="optimal.adj")$p.value
      # #SKATO_C_p: 常见变异所对应的SKAT-O P值
      # SKATO_C_p <- SKAT::SKAT(Genotype_s[,-(1:sit$n_rare)],  
      #                   weights.beta=c(0.5,0.5),obj,method="optimal.adj")$p.value
      # 
      # #罕见变异和常见变异的P值进行柯西合并
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
  
  #计算各种参数组合下的 Type I error 或者检验效能
  T1E_power <- plyr::ddply(TestRes,.(id),function(df){
    T1E_power <- apply(df[-1], 2, function(x) mean(x<alpha,na.rm=TRUE))
  }) %>%
    left_join(situation,by="id")
  return(T1E_power)

} 


####################调用函数，计算一类错误或者检验效能#######################


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


##并行运算框架

cl <- makeCluster(detectCores()/2,type="SOCK") ## 创建特定cpu核心的集群
registerDoParallel(cl) ## 注册集群

# clusterExport(cl,)
clusterEvalQ(cl, {library(tidyverse);library(coin)}) ## 加载包到cl集群

#计算一类错误或检验效能
suppressWarnings(myT1E_power1 <- T1E_power(seed,n_pool,n_ld_est,n,
                                           n_rare,n_common,
                                           beta_g,beta_0,
                                           Nsim,rho,alpha=0.05))

#排序
newmyT1E_power1 <- myT1E_power1[order(myT1E_power1$theta,
                                      myT1E_power1$beta_g,
                                      myT1E_power1$ratio,
                                      myT1E_power1$pm,
                                      myT1E_power1$pf,
                                      myT1E_power1$gamma),]

##停止集群，释放内存  
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