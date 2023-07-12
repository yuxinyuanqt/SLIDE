##############################################################################

#时间:2021年8月11日

#性状：独立数据、二分类性状；
#特点：考虑罕见变异和常见变异联合效应
#目的：检验 X 染色体上的基因或区域是否与疾病相关联，这个基因或区域可能既有常见变异又有罕见变异；
#      并比较不同赋值方式的一类错误和检验效能




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

#LD函数: 估计出女性和男性基因型向量Gf和Gm所对应的LD矩阵
#LD_f0和LD_m0: 对角线上元素pi(1-pi),pi是等位基因频率
#              非对角线上元素是加权平均的LD值

LD <- function(Genotype,sex,weight=2,
               em_max=100,tolerance=0.00001){
  
  #Genotype: n*m的基因型矩阵, n是样本量, m是感兴趣位点个数
  #          女性基因型赋值为0、1和2，男性基因型赋值为0和1
  #sex:性别向量，用于拆分基因型矩阵，男性赋值为0, 女性赋值为1
  #weight: 用于对女性数据和男性数据估计出的LD值进行加权平均得到人群估计的LD值，
  #        女性和男性LD值对应的权重是weight*nf/(weight*nf+nm)和1-weight*nf/(weight*nf+nm),默认是2，
  #        即女性LD对应的权重是2*nf/(2*nf+nm)
  #em_max: 容许的最大的EM迭代次数, 默认值是100
  #tolerance: 前一次迭代和后一次迭代最大绝对误差的容许值，默认值为0.00001
  
  #nm、nf：分别代表男性和女性个体数, m代表感兴趣基因或区域的位点数
  nm <- sum(sex==0)
  nf <- sum(sex==1)
  m <- ncol(Genotype)
  
  #对男性数据和女性数据分别调用LD.m和LD.f函数，得到两个列表LD_m和LD_f
  LD_m <- LD.m(Genotype[sex==0,])
  LD_f <- LD.f(Genotype[sex==1,],em_max,tolerance)
  
  #LD_m0、LD_f0: 预先给男性和女性LD矩阵分配了内存空间
  LD_m0 <- LD_f0 <- matrix(nrow = m,ncol = m)
  
  #lower_tri、upper_tri: 逻辑矩阵，是否为下三角或上三角元素
  lower_tri <- lower.tri(LD_m0)
  upper_tri <- upper.tri(LD_m0)
  
  #weight_ld: 女性和男性估计的ld值进行加权合并的权重
  weight_ld <- weight*nf/(weight*nf+nm)
  
  #为LD_m0和LD_f0非对角线区域赋值
  LD_m0[lower_tri] <- LD_f0[lower_tri] <- (1-weight_ld)*LD_m$LD_m0[lower_tri]+weight_ld*LD_f$LD_f0[lower_tri]
  LD_m0[upper_tri] <- LD_f0[upper_tri] <- t(LD_m0)[upper_tri]
  
  #为LD_m0和LD_f0对角线区域赋值
  diag(LD_m0) <- diag(LD_m$LD_m0)
  diag(LD_f0) <- diag(LD_f$LD_f0)
  
  #返回一个列表，包括4个组件,LD_m0和LD_f0分别是男性和女性的LD矩阵,
  #pm和pf则是男性和女性感兴趣基因各位点估计出的等位基因频率向量
  return(list(LD_m0=LD_m0,LD_f0=LD_f0,pm=LD_m$pm,pf=LD_f$pf))
}

#估计出女性基因型向量Gf所对应的方差-协方差矩阵
Cov.Gf <- function(LD_f,pf,XCI='XCI',gamma=1){
  
  #LD_f: 女性的LD矩阵，对角线上元素是pf*(1-pf)，非对角线上元素是位点i和位点j的LD系数ldij
  #pf: 等位基因频率向量
  #XCI:失活模式,可为逃逸失活("XCI-E")或非逃逸失活("XCI")，默认是XCI
  #gamma: 若XCI为非逃逸失活, 则需指定偏移失活系数gamma的值, gamma取值在[0,2]，默认是1
  
  #m代表感兴趣基因或区域的位点数
  m <- ncol(LD_f)
  
  if(XCI=='XCI-E' | (XCI=='XCI' & gamma==1)){
    Cov_Gf <- 2*LD_f
  } else{
    
    #E_Gf:女性基因型向量Gf对应的期望
    E_Gf <- 2*pf*(gamma*(1-pf)+pf)
    
    
    #Cov_Gf: 预先给女性基因型向量Gf对应的方差-协方差矩阵分配了内存空间
    Cov_Gf <- matrix(nrow = m,ncol = m)
    
    #为Cov_Gf对角线上元素赋值，即基因型向量各分量对应的方差
    diag(Cov_Gf) <- 2*pf*(1-pf)*gamma^2+4*pf^2-E_Gf^2
    
    #考虑到Cov_Gf是一个对称矩阵，只需考虑矩阵的上三角部分即可
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        
        pf_AB <- pf[i]*pf[j]+LD_f[i,j]
        pf_Ab <- pf[i]*(1-pf[j])-LD_f[i,j]
        pf_aB <- (1-pf[i])*pf[j]-LD_f[i,j]
        pf_ab <- (1-pf[i])*(1-pf[j])+LD_f[i,j]
        
        #Cov_Gf[i,j]=Cov_Gf[j,i]为位点i和位点j的协方差
        Cov_Gf[i,j] <- Cov_Gf[j,i] <- 2*gamma^2*(pf_ab*pf_AB+pf_Ab*pf_aB)+
          4*gamma*pf_AB*(pf_Ab+pf_aB)+4*pf_AB^2-E_Gf[i]*E_Gf[j]
      }
    }
  }
  return(Cov_Gf)
}

#使用EM算法估计出女性的基因型向量Gf所对应的LD矩阵，
#LD_f0-2第i个对角线上元素是pi(1-pi), 非对角线元素是位点i和位点j的ld系数

LD.f <- function(Genotype_f,em_max=100,tolerance=0.00001){
  
  #Genotype_f: nf*m维的基因型矩阵
  #em_max: 容许的最大的EM迭代次数, 默认值是100
  #tolerance: 前一次迭代和后一次迭代最大绝对误差的容许值，默认值为0.00001
  
  #nf代表女性个体数,m代表感兴趣基因或区域的位点数
  nf <- nrow(Genotype_f)
  m <- ncol(Genotype_f)
  
  #dat_Gf: 将基因型矩阵转换为列表，并将每一列均转化为水平为(0、1、2)的因子,
  #        保证后续使用table函数能得出想要的结果
  dat_Gf <- lapply(as.data.frame(Genotype_f),factor,levels=0:2)
  
  #LD_f0-2:预先给m*m的Gf的方差-协方差矩阵分配了内存空间，之后依次填充Cov_Gf即可，
  #       这样在执行循环的过程中，可以避免进行多次内存分配
  #LD_f0: 非对角线上元素是位点i和位点j的LD系数
  #LD_f1: 非对角线上元素是Lewontin（1964）提出的 LD 标准化参数，即相对连锁不平衡系数
  #LD_f2: 非对角线上元素r^2
  LD_f0 <- LD_f1 <- LD_f2 <- matrix(nrow = m,ncol = m)
  
  #pf:利用Genotype_f估计出来的女性感兴趣位点的次等位基因频率
  pf <- colSums(Genotype_f)/(2*nf)
  
  #为LD_f0-2对角线上元素赋值
  diag(LD_f0) <- diag(LD_f1) <- diag(LD_f2) <- pf*(1-pf)
  
  #考虑到LD_f是一个对称矩阵，只需考虑矩阵的上三角部分即可
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      
      #预先给向量分配内存空间，避免进行多次内存分配
      #pf_AB:pf_ab:每次迭代估计出的单倍型频率所组成的向量
      pf_AB <- pf_Ab <- pf_aB <- pf_ab <- vector(length = (em_max+1))
      
      #nf_geno:对AABB(频数为nf_geno[3,3])和AaBb(频数为nf_geno[2,2])等可能的9种基因型进行计数
      nf_geno <- table(dat_Gf[[i]],dat_Gf[[j]])
      
      #pf_AB_b:每次迭代中固定的部分(base)，代表各单倍型中相已知的部分
      #        （未知部分在AaBb中），避免重复计算
      pf_AB_b <- 2*nf_geno[3,3]+nf_geno[3,2]+nf_geno[2,3]
      pf_ab_b <- 2*nf_geno[1,1]+nf_geno[1,2]+nf_geno[2,1]
      pf_Ab_b <- 2*nf_geno[3,1]+nf_geno[3,2]+nf_geno[2,1]
      pf_aB_b <- 2*nf_geno[1,3]+nf_geno[1,2]+nf_geno[2,3]
      
      #单倍型频率初始值，即在连锁平衡的条件下的单倍型频率
      pf_AB[1] <- pf[i]*pf[j]
      pf_Ab[1] <- pf[i]*(1-pf[j])
      pf_aB[1] <- (1-pf[i])*pf[j]
      pf_ab[1] <- (1-pf[i])*(1-pf[j])
      
      for(t in 1:em_max){
        
        #a和b: E步，已知基因型为AaBb，女性两条单倍型分别为AB和ab、Ab和aB的频率
        
        #Full_prob: 保证分母不为0
        Full_prob <- ifelse((pf_AB[t]*pf_ab[t]+pf_Ab[t]*pf_aB[t])==0,1,pf_AB[t]*pf_ab[t]+pf_Ab[t]*pf_aB[t])
        a <- pf_AB[t]*pf_ab[t]/Full_prob
        b <- pf_Ab[t]*pf_aB[t]/Full_prob
        
        #M步，该次迭代单倍型频率的更新值
        pf_AB[t+1] <- (pf_AB_b+a*nf_geno[2,2])/(2*nf)
        pf_ab[t+1] <- (pf_ab_b+a*nf_geno[2,2])/(2*nf)
        pf_Ab[t+1] <- (pf_Ab_b+b*nf_geno[2,2])/(2*nf)
        pf_aB[t+1] <- (pf_aB_b+b*nf_geno[2,2])/(2*nf)

        #tol: 前一次迭代和后一次迭代的最大绝对误差
        tol <- max(abs(c(pf_AB[t+1]-pf_AB[t],pf_Ab[t+1]-pf_Ab[t],
                           pf_aB[t+1]-pf_aB[t],pf_ab[t+1]-pf_ab[t])))
        
        #当tol <= tolerance，迭代终止
        if(tol <= tolerance)   break
      }
      #LD_f0[i,j]=LD_f0[j,i]=ldij,ldij即位点i和位点j的LD系数
      LD_f0[i,j] <- LD_f0[j,i] <- pf_AB[t+1]*pf_ab[t+1]-pf_Ab[t+1]*pf_aB[t+1]
      
      #LD_f1: 非对角线上元素是Lewontin（1964）提出的 LD 标准化参数，即相对连锁不平衡系数
      LD_f1[i,j] <- LD_f1[j,i] <- abs(LD_f0[i,j])/(min(pf[i]*(1-pf[j]),pf[j]*(1-pf[i]))*sign(LD_f0[i,j]>0)+
                                                     min(pf[i]*pf[j],(1-pf[i])*(1-pf[j]))*sign(LD_f0[i,j]<0)+
                                                     sign(LD_f0[i,j]==0))
      
      #LD_f2: 非对角线上元素r^2
      LD_f2[i,j] <- LD_f2[j,i] <- LD_f0[i,j]^2/(pf[i]*pf[j]*(1-pf[i])*(1-pf[j]))
      }
  }
  return(list(LD_f0=LD_f0,LD_f1=LD_f1,LD_f2=LD_f2,pf=pf))
}


#估计出男性的基因型向量Gm所对应的LD矩阵，
#LD_m0-2第i个对角线上元素是pi(1-pi), 非对角线元素是位点i和位点j的ld系数

LD.m <- function(Genotype_m){
  
  #Genotype_m: nm*m维的基因型矩阵
  
  #nm代表男性个体数,m代表感兴趣基因或区域的位点数
  nm <- nrow(Genotype_m)
  m <- ncol(Genotype_m)
  
  #dat_Gm: 将基因型矩阵转换为数据框，并将每一列均转化为水平为(0、1)的因子,
  #        保证后续使用table函数能得出想要的结果
  dat_Gm <- as.data.frame(lapply(as.data.frame(Genotype_m),factor,levels=0:1))
  
  #Cov_Gm0-2:预先给m*m的Gf的方差-协方差矩阵分配了内存空间，之后依次填充Cov_Gm即可，
  #       这样在执行循环的过程中，可以避免进行多次内存分配
  #LD_m0: 非对角线上元素是位点i和位点j的LD系数
  #LD_m1: 非对角线上元素是Lewontin（1964）提出的 LD 标准化参数，即相对连锁不平衡系数
  #LD_m2: 非对角线上元素r^2
  LD_m0 <- LD_m1 <- LD_m2 <- matrix(nrow = m,ncol = m)
  
  
  #pm:利用Genotype_m估计出来的男性感兴趣位点的次等位基因频率
  pm <- colSums(Genotype_m)/(nm)
  
  #为LD_f0-2对角线上元素赋值
  diag(LD_m0) <- diag(LD_m1) <- diag(LD_m2) <- pm*(1-pm)
  
  
  #考虑到LD_m是一个对称矩阵，只需考虑矩阵的上三角部分即可
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      
      #nm_geno:对AB(频数为nm_geno[2,2])和Ab(频数为nm_geno[2,1])等可能的4种基因型进行计数
      nm_geno <- table(dat_Gm[,i],dat_Gm[,j])
      
      #LD_m0[i,j]=LD_m0[j,i]=ldij,ldij即位点i和位点j的LD系数
      LD_m0[i,j] <- LD_m0[j,i] <- (nm_geno[2,2]*nm_geno[1,1]-nm_geno[1,2]*nm_geno[2,1])/(nm^2)
      
      #LD_m1: 非对角线上元素是Lewontin（1964）提出的 LD 标准化参数，即相对连锁不平衡系数
      LD_m1[i,j] <- LD_m1[j,i] <- abs(LD_m0[i,j])/(min(pm[i]*(1-pm[j]),pm[j]*(1-pm[i]))*sign(LD_m0[i,j]>0)+
                                                     min(pm[i]*pm[j],(1-pm[i])*(1-pm[j]))*sign(LD_m0[i,j]<0)+
                                                     sign(LD_m0[i,j]==0))
      
      #LD_m2: 非对角线上元素r^2
      LD_m2[i,j] <- LD_m2[j,i] <- LD_m0[i,j]^2/(pm[i]*pm[j]*(1-pm[i])*(1-pm[j]))
    }
  }
  return(list(LD_m0=LD_m0,LD_m1=LD_m1,LD_m2=LD_m2,pm=pm))
}


SILDE <- function(Genotype,Phenotype,LD_estimate='sample',
                  LD_matrix,n_per=1000,
                  em_max=100,tolerance=0.00001){
  
  #Genotype: 基因型矩阵
  #Phenotype:表型向量, 患病个体赋值为1，未患病个体赋值为0
  #LD_estimate: 估计出基因型向量G所对应的方差-协方差矩阵的方法，可为c('control','population','both)
  #             'sample': 使用样本的数据来估计LD值
  #             'control': 使用样本的对照数据来估计LD值
  #             'population': 使用外部的人群数据来估计LD值，若选择该选项，则需给LD_matrix赋值
  #             'all': 提供'sample'、'control'和'population'所对应的结果
  #LD_matrix: LD_estimate='population'时需提供, 根据外部的人群数据估计出来的LD矩阵
  #p_value: 计算p值的方法，可为c('Asymptotically','permutation','both)，
  #         'Asymptotically': 使用卡方分布得到近似的p值
  #         'permutation': 使用permutation来计算p值
  #         'both': 提供卡方分布得到近似的p值和permutation计算的p值
  #n_per: permutation过程的次数，默认是1000
  #em_max: 容许的最大的EM迭代次数, 默认值是100
  #tolerance: 前一次迭代和后一次迭代最大绝对误差的容许值，默认值为0.00001
  
  
  #n_phe:病例-对照数目，对照数为n_phe[1]，病例数为n_phe[2]
  n_phe <- table(Phenotype)
  
  #n:样本量; prod_n_phe:病例数与对照数相乘; m: 位点数
  n <- sum(n_phe)
  prod_n_phe <- prod(n_phe)
  m <- ncol(Genotype)
  
  #U统计量，即病例组和对照组之间多位点平均基因型得分的差异
  U <- prod_n_phe*(colSums(Genotype[Phenotype==1,])/n_phe[2]-colSums(Genotype[Phenotype==0,])/n_phe[1])/n
  
  # #U统计量的另一种表示
  # G_I <- t(Genotype)%*%(diag(n)-(rep(1,n)%*%t(rep(1,n)))/n)
  # #U <- G_I%*%Phenotype
  # U <- rowSums(G_I[,Phenotype])
  
  if(LD_estimate=='sample'){
    
    #使用样本的数据来估计LD值
    #ginv_Cov_U_tem: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_tem <- n*ginv(LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_tem: SLIDE 统计量
    T_tem <- as.vector(U%*%ginv_Cov_U_tem%*%U)
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_tem_asy <- pchisq(T_tem,m,lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_tem_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_tem_per <- U_per%*%ginv_Cov_U_tem%*%U_per
    })
    
    #使用permutation过程来计算p值
    p_value_tem_per <- sum(T_tem_per>=T_tem)/n_per
    
    SLIDE_tem <- c(T_tem,p_value_tem_asy,p_value_tem_per)
    names(SLIDE_tem) <- c('T_tem','p_value_tem_asy','p_value_tem_per')
    return(SLIDE_tem)
    
  } else if(LD_estimate=='control'){
    
    #使用样本的对照数据来估计LD值
    #ginv_Cov_U_tem: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_con <- n*ginv(LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_con: SLIDE 统计量
    T_con <- as.vector(U%*%ginv_Cov_U_con%*%U)
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_con_asy <- pchisq(T_con,m,lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_con_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_con_per <- U_per%*%ginv_Cov_U_con%*%U_per
    })
    
    #使用permutation过程来计算p值
    p_value_con_per <- sum(T_con_per>=T_con)/n_per
    
    SLIDE_con <- c(T_con,p_value_con_asy,p_value_con_per)
    names(SLIDE_con) <- c('T_con','p_value_con_asy','p_value_con_per')
    return(SLIDE_con)
    
  } else if(LD_estimate=='population'){
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    
    
    #使用外部的人群数据来估计LD值
    #ginv_Cov_U_pop: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_pop <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    #T_pop: SLIDE 统计量
    T_pop <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U%*%ginv_Cov_U_popi%*%U})
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_pop_asy <- pchisq(T_pop,m, lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_pop_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_pop_per <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U_per%*%ginv_Cov_U_popi%*%U_per})
    })
    
    #使用permutation过程来计算p值
    p_value_pop_per <- rowSums(T_pop_per>=T_pop)/n_per
    
    SLIDE_pop<- cbind(T_pop,p_value_pop_asy,p_value_pop_per)
    rownames(SLIDE_pop) <- paste('SLIDE_pop',1:length(LD_matrix),sep = "")
    return(SLIDE_pop)
    
  }else if(LD_estimate=='all'){
    
    #提供'sample'、'control'和'population'所对应的结果
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    LD_matrix$ld_sap <- LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0
    LD_matrix$ld_con <- LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0
    
    ginv_Cov_U <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    SLIDE_T <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U%*%ginv_Cov_Ui%*%U})
    
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_asy <- pchisq(SLIDE_T,m, lower.tail = FALSE)
    
    #返回一个3*n_per的矩阵
    #第一行是基于样本数据估计出LD对应的n_per次的permutation统计量
    #第二行是基于样本的对照组数据估计出LD对应的n_per次的permutation统计量
    #第三行是基于人群数据估计出LD对应的n_per次的permutation统计量
    T_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      
      T_per <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U_per%*%ginv_Cov_Ui%*%U_per})
      return(T_per)
    })
    
    #使用permutation过程来计算p值,返回一个长度为3的向量
    #向量的元素分别是T_tem、T_com和T_pop对应的p值
    p_value_per <- rowSums(T_per>=SLIDE_T)/n_per
    
    SLIDE <- cbind(SLIDE_T,p_value_asy,p_value_per)
    rownames(SLIDE) <- c(paste('SLIDE_pop',1:(length(LD_matrix)-2),sep = ""),
                         'SLIDE_tem','SLIDE_con')
    
    return(SLIDE)
  }
  
}
  

SILDE.X <- function(Genotype,Phenotype,sex,LD_estimate='sample',
                    Cov_G=FALSE,Cov_Gf,Cov_Gm,LD_f,LD_m,
                    gamma=seq(0,2,length=5),n_per=1000,
                    em_max=100,tolerance=0.00001){
  
  #Genotype: 基因型矩阵，女性赋值为0/1/2、男性赋值为0/1
  #Phenotype:表型向量, 患病个体赋值为1，未患病个体赋值为0
  #sex: 性别向量，用于拆分基因型矩阵，男性赋值为0, 女性赋值为1
  #LD_estimate: 估计出基因型向量G所对应的方差-协方差矩阵的方法，可为c('control','population','both)
  #             'sample': 使用样本的数据来估计LD值
  #             'control': 使用样本的对照数据来估计LD值
  #             'population': 使用外部的人群数据来估计LD值，若选择该选项，则需给LD_matrix赋值
  #             'all': 提供'sample'、'control'和'population'所对应的结果
  #LD_matrix: LD_estimate='population'时需提供, 根据外部的人群数据估计出来的LD矩阵
  
  #gamma: 偏倚失活系数可能的取值，默认是0.0、0.5、1.0、1.5、2.0
  #使用网格搜索法找到最大的统计量
  
  #p_value: 计算p值的方法，可为c('Asymptotically','permutation','both)，
  #         'Asymptotically': 使用卡方分布得到近似的p值
  #         'permutation': 使用permutation来计算p值
  #         'both': 提供卡方分布得到近似的p值和permutation计算的p值
  #n_per: permutation过程的次数，默认是1000
  #em_max: 容许的最大的EM迭代次数, 默认值是100
  #tolerance: 前一次迭代和后一次迭代最大绝对误差的容许值，默认值为0.00001
  
  
  #n_phe:病例-对照数目，对照数为n_phe[1]，病例数为n_phe[2]
  n_phe <- table(Phenotype)
  
  #n:样本量; prod_n_phe:病例数与对照数相乘; m: 位点数
  n <- sum(n_phe)
  prod_n_phe <- prod(n_phe)
  m <- ncol(Genotype)
  
  #U统计量，即病例组和对照组之间多位点平均基因型得分的差异
  U <- prod_n_phe*(colSums(Genotype[Phenotype==1,])/n_phe[2]-colSums(Genotype[Phenotype==0,])/n_phe[1])/n
  
  # #U统计量的另一种表示
  # G_I <- t(Genotype)%*%(diag(n)-(rep(1,n)%*%t(rep(1,n)))/n)
  # #U <- G_I%*%Phenotype
  # U <- rowSums(G_I[,Phenotype])
  
  if(LD_estimate=='sample'){
    
    #使用样本的数据来估计LD值
    #ginv_Cov_U_tem: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_tem <- n*ginv(LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_tem: SLIDE 统计量
    T_tem <- as.vector(U%*%ginv_Cov_U_tem%*%U)
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_tem_asy <- pchisq(T_tem,m,lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_tem_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_tem_per <- U_per%*%ginv_Cov_U_tem%*%U_per
    })
    
    #使用permutation过程来计算p值
    p_value_tem_per <- sum(T_tem_per>=T_tem)/n_per
    
    SLIDE_tem <- c(T_tem,p_value_tem_asy,p_value_tem_per)
    names(SLIDE_tem) <- c('T_tem','p_value_tem_asy','p_value_tem_per')
    return(SLIDE_tem)
    
  } else if(LD_estimate=='control'){
    
    #使用样本的对照数据来估计LD值
    #ginv_Cov_U_tem: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_con <- n*ginv(LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_con: SLIDE 统计量
    T_con <- as.vector(U%*%ginv_Cov_U_con%*%U)
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_con_asy <- pchisq(T_con,m,lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_con_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_con_per <- U_per%*%ginv_Cov_U_con%*%U_per
    })
    
    #使用permutation过程来计算p值
    p_value_con_per <- sum(T_con_per>=T_con)/n_per
    
    SLIDE_con <- c(T_con,p_value_con_asy,p_value_con_per)
    names(SLIDE_con) <- c('T_con','p_value_con_asy','p_value_con_per')
    return(SLIDE_con)
    
  } else if(LD_estimate=='population'){
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    
    
    #使用外部的人群数据来估计LD值
    #ginv_Cov_U_pop: U统计量的方差-协方差矩阵的逆矩阵
    ginv_Cov_U_pop <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    #T_pop: SLIDE 统计量
    T_pop <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U%*%ginv_Cov_U_popi%*%U})
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_pop_asy <- pchisq(T_pop,m, lower.tail = FALSE)
    
    #返回一个长度为n_per的向量, 即n_per次的permutation统计量
    T_pop_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_pop_per <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U_per%*%ginv_Cov_U_popi%*%U_per})
    })
    
    #使用permutation过程来计算p值
    p_value_pop_per <- rowSums(T_pop_per>=T_pop)/n_per
    
    SLIDE_pop<- cbind(T_pop,p_value_pop_asy,p_value_pop_per)
    rownames(SLIDE_pop) <- paste('SLIDE_pop',1:length(LD_matrix),sep = "")
    return(SLIDE_pop)
    
  }else if(LD_estimate=='all'){
    
    #提供'sample'、'control'和'population'所对应的结果
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    LD_matrix$ld_sap <- LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0
    LD_matrix$ld_con <- LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0
    
    ginv_Cov_U <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    SLIDE_T <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U%*%ginv_Cov_Ui%*%U})
    
    
    #使用自由度为m的卡方分布得到渐近的p值
    p_value_asy <- pchisq(SLIDE_T,m, lower.tail = FALSE)
    
    #返回一个3*n_per的矩阵
    #第一行是基于样本数据估计出LD对应的n_per次的permutation统计量
    #第二行是基于样本的对照组数据估计出LD对应的n_per次的permutation统计量
    #第三行是基于人群数据估计出LD对应的n_per次的permutation统计量
    T_per <- sapply(1:n_per,function(i){
      #permutation统计量
      #对表型重抽样
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      
      T_per <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U_per%*%ginv_Cov_Ui%*%U_per})
      return(T_per)
    })
    
    #使用permutation过程来计算p值,返回一个长度为3的向量
    #向量的元素分别是T_tem、T_com和T_pop对应的p值
    p_value_per <- rowSums(T_per>=SLIDE_T)/n_per
    
    SLIDE <- cbind(SLIDE_T,p_value_asy,p_value_per)
    rownames(SLIDE) <- c(paste('SLIDE_pop',1:(length(LD_matrix)-2),sep = ""),
                         'SLIDE_tem','SLIDE_con')
    
    return(SLIDE)
  }
  
}