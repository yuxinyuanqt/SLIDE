##############################################################################

#ʱ��:2021��8��11��

#��״���������ݡ���������״��
#�ص㣺���Ǻ�������ͳ�����������ЧӦ
#Ŀ�ģ����� X Ⱦɫ���ϵĻ���������Ƿ��뼲�����������������������ܼ��г����������к������죻
#      ���Ƚϲ�ͬ��ֵ��ʽ��һ�����ͼ���Ч��




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

#LD����: ���Ƴ�Ů�Ժ����Ի���������Gf��Gm����Ӧ��LD����
#LD_f0��LD_m0: �Խ�����Ԫ��pi(1-pi),pi�ǵ�λ����Ƶ��
#              �ǶԽ�����Ԫ���Ǽ�Ȩƽ����LDֵ

LD <- function(Genotype,sex,weight=2,
               em_max=100,tolerance=0.00001){
  
  #Genotype: n*m�Ļ����;���, n��������, m�Ǹ���Ȥλ�����
  #          Ů�Ի����͸�ֵΪ0��1��2�����Ի����͸�ֵΪ0��1
  #sex:�Ա����������ڲ�ֻ����;������Ը�ֵΪ0, Ů�Ը�ֵΪ1
  #weight: ���ڶ�Ů�����ݺ��������ݹ��Ƴ���LDֵ���м�Ȩƽ���õ���Ⱥ���Ƶ�LDֵ��
  #        Ů�Ժ�����LDֵ��Ӧ��Ȩ����weight*nf/(weight*nf+nm)��1-weight*nf/(weight*nf+nm),Ĭ����2��
  #        ��Ů��LD��Ӧ��Ȩ����2*nf/(2*nf+nm)
  #em_max: ����������EM��������, Ĭ��ֵ��100
  #tolerance: ǰһ�ε����ͺ�һ�ε�����������������ֵ��Ĭ��ֵΪ0.00001
  
  #nm��nf���ֱ�������Ժ�Ů�Ը�����, m��������Ȥ����������λ����
  nm <- sum(sex==0)
  nf <- sum(sex==1)
  m <- ncol(Genotype)
  
  #���������ݺ�Ů�����ݷֱ����LD.m��LD.f�������õ������б�LD_m��LD_f
  LD_m <- LD.m(Genotype[sex==0,])
  LD_f <- LD.f(Genotype[sex==1,],em_max,tolerance)
  
  #LD_m0��LD_f0: Ԥ�ȸ����Ժ�Ů��LD����������ڴ�ռ�
  LD_m0 <- LD_f0 <- matrix(nrow = m,ncol = m)
  
  #lower_tri��upper_tri: �߼������Ƿ�Ϊ�����ǻ�������Ԫ��
  lower_tri <- lower.tri(LD_m0)
  upper_tri <- upper.tri(LD_m0)
  
  #weight_ld: Ů�Ժ����Թ��Ƶ�ldֵ���м�Ȩ�ϲ���Ȩ��
  weight_ld <- weight*nf/(weight*nf+nm)
  
  #ΪLD_m0��LD_f0�ǶԽ�������ֵ
  LD_m0[lower_tri] <- LD_f0[lower_tri] <- (1-weight_ld)*LD_m$LD_m0[lower_tri]+weight_ld*LD_f$LD_f0[lower_tri]
  LD_m0[upper_tri] <- LD_f0[upper_tri] <- t(LD_m0)[upper_tri]
  
  #ΪLD_m0��LD_f0�Խ�������ֵ
  diag(LD_m0) <- diag(LD_m$LD_m0)
  diag(LD_f0) <- diag(LD_f$LD_f0)
  
  #����һ���б�������4�����,LD_m0��LD_f0�ֱ������Ժ�Ů�Ե�LD����,
  #pm��pf�������Ժ�Ů�Ը���Ȥ�����λ����Ƴ��ĵ�λ����Ƶ������
  return(list(LD_m0=LD_m0,LD_f0=LD_f0,pm=LD_m$pm,pf=LD_f$pf))
}

#���Ƴ�Ů�Ի���������Gf����Ӧ�ķ���-Э�������
Cov.Gf <- function(LD_f,pf,XCI='XCI',gamma=1){
  
  #LD_f: Ů�Ե�LD���󣬶Խ�����Ԫ����pf*(1-pf)���ǶԽ�����Ԫ����λ��i��λ��j��LDϵ��ldij
  #pf: ��λ����Ƶ������
  #XCI:ʧ��ģʽ,��Ϊ����ʧ��("XCI-E")�������ʧ��("XCI")��Ĭ����XCI
  #gamma: ��XCIΪ������ʧ��, ����ָ��ƫ��ʧ��ϵ��gamma��ֵ, gammaȡֵ��[0,2]��Ĭ����1
  
  #m��������Ȥ����������λ����
  m <- ncol(LD_f)
  
  if(XCI=='XCI-E' | (XCI=='XCI' & gamma==1)){
    Cov_Gf <- 2*LD_f
  } else{
    
    #E_Gf:Ů�Ի���������Gf��Ӧ������
    E_Gf <- 2*pf*(gamma*(1-pf)+pf)
    
    
    #Cov_Gf: Ԥ�ȸ�Ů�Ի���������Gf��Ӧ�ķ���-Э�������������ڴ�ռ�
    Cov_Gf <- matrix(nrow = m,ncol = m)
    
    #ΪCov_Gf�Խ�����Ԫ�ظ�ֵ����������������������Ӧ�ķ���
    diag(Cov_Gf) <- 2*pf*(1-pf)*gamma^2+4*pf^2-E_Gf^2
    
    #���ǵ�Cov_Gf��һ���Գƾ���ֻ�迼�Ǿ���������ǲ��ּ���
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        
        pf_AB <- pf[i]*pf[j]+LD_f[i,j]
        pf_Ab <- pf[i]*(1-pf[j])-LD_f[i,j]
        pf_aB <- (1-pf[i])*pf[j]-LD_f[i,j]
        pf_ab <- (1-pf[i])*(1-pf[j])+LD_f[i,j]
        
        #Cov_Gf[i,j]=Cov_Gf[j,i]Ϊλ��i��λ��j��Э����
        Cov_Gf[i,j] <- Cov_Gf[j,i] <- 2*gamma^2*(pf_ab*pf_AB+pf_Ab*pf_aB)+
          4*gamma*pf_AB*(pf_Ab+pf_aB)+4*pf_AB^2-E_Gf[i]*E_Gf[j]
      }
    }
  }
  return(Cov_Gf)
}

#ʹ��EM�㷨���Ƴ�Ů�ԵĻ���������Gf����Ӧ��LD����
#LD_f0-2��i���Խ�����Ԫ����pi(1-pi), �ǶԽ���Ԫ����λ��i��λ��j��ldϵ��

LD.f <- function(Genotype_f,em_max=100,tolerance=0.00001){
  
  #Genotype_f: nf*mά�Ļ����;���
  #em_max: ����������EM��������, Ĭ��ֵ��100
  #tolerance: ǰһ�ε����ͺ�һ�ε�����������������ֵ��Ĭ��ֵΪ0.00001
  
  #nf����Ů�Ը�����,m��������Ȥ����������λ����
  nf <- nrow(Genotype_f)
  m <- ncol(Genotype_f)
  
  #dat_Gf: �������;���ת��Ϊ�б�������ÿһ�о�ת��ΪˮƽΪ(0��1��2)������,
  #        ��֤����ʹ��table�����ܵó���Ҫ�Ľ��
  dat_Gf <- lapply(as.data.frame(Genotype_f),factor,levels=0:2)
  
  #LD_f0-2:Ԥ�ȸ�m*m��Gf�ķ���-Э�������������ڴ�ռ䣬֮���������Cov_Gf���ɣ�
  #       ������ִ��ѭ���Ĺ����У����Ա�����ж���ڴ����
  #LD_f0: �ǶԽ�����Ԫ����λ��i��λ��j��LDϵ��
  #LD_f1: �ǶԽ�����Ԫ����Lewontin��1964������� LD ��׼�������������������ƽ��ϵ��
  #LD_f2: �ǶԽ�����Ԫ��r^2
  LD_f0 <- LD_f1 <- LD_f2 <- matrix(nrow = m,ncol = m)
  
  #pf:����Genotype_f���Ƴ�����Ů�Ը���Ȥλ��Ĵε�λ����Ƶ��
  pf <- colSums(Genotype_f)/(2*nf)
  
  #ΪLD_f0-2�Խ�����Ԫ�ظ�ֵ
  diag(LD_f0) <- diag(LD_f1) <- diag(LD_f2) <- pf*(1-pf)
  
  #���ǵ�LD_f��һ���Գƾ���ֻ�迼�Ǿ���������ǲ��ּ���
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      
      #Ԥ�ȸ����������ڴ�ռ䣬������ж���ڴ����
      #pf_AB:pf_ab:ÿ�ε������Ƴ��ĵ�����Ƶ������ɵ�����
      pf_AB <- pf_Ab <- pf_aB <- pf_ab <- vector(length = (em_max+1))
      
      #nf_geno:��AABB(Ƶ��Ϊnf_geno[3,3])��AaBb(Ƶ��Ϊnf_geno[2,2])�ȿ��ܵ�9�ֻ����ͽ��м���
      nf_geno <- table(dat_Gf[[i]],dat_Gf[[j]])
      
      #pf_AB_b:ÿ�ε����й̶��Ĳ���(base)��������������������֪�Ĳ���
      #        ��δ֪������AaBb�У��������ظ�����
      pf_AB_b <- 2*nf_geno[3,3]+nf_geno[3,2]+nf_geno[2,3]
      pf_ab_b <- 2*nf_geno[1,1]+nf_geno[1,2]+nf_geno[2,1]
      pf_Ab_b <- 2*nf_geno[3,1]+nf_geno[3,2]+nf_geno[2,1]
      pf_aB_b <- 2*nf_geno[1,3]+nf_geno[1,2]+nf_geno[2,3]
      
      #������Ƶ�ʳ�ʼֵ����������ƽ��������µĵ�����Ƶ��
      pf_AB[1] <- pf[i]*pf[j]
      pf_Ab[1] <- pf[i]*(1-pf[j])
      pf_aB[1] <- (1-pf[i])*pf[j]
      pf_ab[1] <- (1-pf[i])*(1-pf[j])
      
      for(t in 1:em_max){
        
        #a��b: E������֪������ΪAaBb��Ů�����������ͷֱ�ΪAB��ab��Ab��aB��Ƶ��
        
        #Full_prob: ��֤��ĸ��Ϊ0
        Full_prob <- ifelse((pf_AB[t]*pf_ab[t]+pf_Ab[t]*pf_aB[t])==0,1,pf_AB[t]*pf_ab[t]+pf_Ab[t]*pf_aB[t])
        a <- pf_AB[t]*pf_ab[t]/Full_prob
        b <- pf_Ab[t]*pf_aB[t]/Full_prob
        
        #M�����ôε���������Ƶ�ʵĸ���ֵ
        pf_AB[t+1] <- (pf_AB_b+a*nf_geno[2,2])/(2*nf)
        pf_ab[t+1] <- (pf_ab_b+a*nf_geno[2,2])/(2*nf)
        pf_Ab[t+1] <- (pf_Ab_b+b*nf_geno[2,2])/(2*nf)
        pf_aB[t+1] <- (pf_aB_b+b*nf_geno[2,2])/(2*nf)

        #tol: ǰһ�ε����ͺ�һ�ε��������������
        tol <- max(abs(c(pf_AB[t+1]-pf_AB[t],pf_Ab[t+1]-pf_Ab[t],
                           pf_aB[t+1]-pf_aB[t],pf_ab[t+1]-pf_ab[t])))
        
        #��tol <= tolerance��������ֹ
        if(tol <= tolerance)   break
      }
      #LD_f0[i,j]=LD_f0[j,i]=ldij,ldij��λ��i��λ��j��LDϵ��
      LD_f0[i,j] <- LD_f0[j,i] <- pf_AB[t+1]*pf_ab[t+1]-pf_Ab[t+1]*pf_aB[t+1]
      
      #LD_f1: �ǶԽ�����Ԫ����Lewontin��1964������� LD ��׼�������������������ƽ��ϵ��
      LD_f1[i,j] <- LD_f1[j,i] <- abs(LD_f0[i,j])/(min(pf[i]*(1-pf[j]),pf[j]*(1-pf[i]))*sign(LD_f0[i,j]>0)+
                                                     min(pf[i]*pf[j],(1-pf[i])*(1-pf[j]))*sign(LD_f0[i,j]<0)+
                                                     sign(LD_f0[i,j]==0))
      
      #LD_f2: �ǶԽ�����Ԫ��r^2
      LD_f2[i,j] <- LD_f2[j,i] <- LD_f0[i,j]^2/(pf[i]*pf[j]*(1-pf[i])*(1-pf[j]))
      }
  }
  return(list(LD_f0=LD_f0,LD_f1=LD_f1,LD_f2=LD_f2,pf=pf))
}


#���Ƴ����ԵĻ���������Gm����Ӧ��LD����
#LD_m0-2��i���Խ�����Ԫ����pi(1-pi), �ǶԽ���Ԫ����λ��i��λ��j��ldϵ��

LD.m <- function(Genotype_m){
  
  #Genotype_m: nm*mά�Ļ����;���
  
  #nm�������Ը�����,m��������Ȥ����������λ����
  nm <- nrow(Genotype_m)
  m <- ncol(Genotype_m)
  
  #dat_Gm: �������;���ת��Ϊ���ݿ򣬲���ÿһ�о�ת��ΪˮƽΪ(0��1)������,
  #        ��֤����ʹ��table�����ܵó���Ҫ�Ľ��
  dat_Gm <- as.data.frame(lapply(as.data.frame(Genotype_m),factor,levels=0:1))
  
  #Cov_Gm0-2:Ԥ�ȸ�m*m��Gf�ķ���-Э�������������ڴ�ռ䣬֮���������Cov_Gm���ɣ�
  #       ������ִ��ѭ���Ĺ����У����Ա�����ж���ڴ����
  #LD_m0: �ǶԽ�����Ԫ����λ��i��λ��j��LDϵ��
  #LD_m1: �ǶԽ�����Ԫ����Lewontin��1964������� LD ��׼�������������������ƽ��ϵ��
  #LD_m2: �ǶԽ�����Ԫ��r^2
  LD_m0 <- LD_m1 <- LD_m2 <- matrix(nrow = m,ncol = m)
  
  
  #pm:����Genotype_m���Ƴ��������Ը���Ȥλ��Ĵε�λ����Ƶ��
  pm <- colSums(Genotype_m)/(nm)
  
  #ΪLD_f0-2�Խ�����Ԫ�ظ�ֵ
  diag(LD_m0) <- diag(LD_m1) <- diag(LD_m2) <- pm*(1-pm)
  
  
  #���ǵ�LD_m��һ���Գƾ���ֻ�迼�Ǿ���������ǲ��ּ���
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      
      #nm_geno:��AB(Ƶ��Ϊnm_geno[2,2])��Ab(Ƶ��Ϊnm_geno[2,1])�ȿ��ܵ�4�ֻ����ͽ��м���
      nm_geno <- table(dat_Gm[,i],dat_Gm[,j])
      
      #LD_m0[i,j]=LD_m0[j,i]=ldij,ldij��λ��i��λ��j��LDϵ��
      LD_m0[i,j] <- LD_m0[j,i] <- (nm_geno[2,2]*nm_geno[1,1]-nm_geno[1,2]*nm_geno[2,1])/(nm^2)
      
      #LD_m1: �ǶԽ�����Ԫ����Lewontin��1964������� LD ��׼�������������������ƽ��ϵ��
      LD_m1[i,j] <- LD_m1[j,i] <- abs(LD_m0[i,j])/(min(pm[i]*(1-pm[j]),pm[j]*(1-pm[i]))*sign(LD_m0[i,j]>0)+
                                                     min(pm[i]*pm[j],(1-pm[i])*(1-pm[j]))*sign(LD_m0[i,j]<0)+
                                                     sign(LD_m0[i,j]==0))
      
      #LD_m2: �ǶԽ�����Ԫ��r^2
      LD_m2[i,j] <- LD_m2[j,i] <- LD_m0[i,j]^2/(pm[i]*pm[j]*(1-pm[i])*(1-pm[j]))
    }
  }
  return(list(LD_m0=LD_m0,LD_m1=LD_m1,LD_m2=LD_m2,pm=pm))
}


SILDE <- function(Genotype,Phenotype,LD_estimate='sample',
                  LD_matrix,n_per=1000,
                  em_max=100,tolerance=0.00001){
  
  #Genotype: �����;���
  #Phenotype:��������, �������帳ֵΪ1��δ�������帳ֵΪ0
  #LD_estimate: ���Ƴ�����������G����Ӧ�ķ���-Э�������ķ�������Ϊc('control','population','both)
  #             'sample': ʹ������������������LDֵ
  #             'control': ʹ�������Ķ�������������LDֵ
  #             'population': ʹ���ⲿ����Ⱥ����������LDֵ����ѡ���ѡ������LD_matrix��ֵ
  #             'all': �ṩ'sample'��'control'��'population'����Ӧ�Ľ��
  #LD_matrix: LD_estimate='population'ʱ���ṩ, �����ⲿ����Ⱥ���ݹ��Ƴ�����LD����
  #p_value: ����pֵ�ķ�������Ϊc('Asymptotically','permutation','both)��
  #         'Asymptotically': ʹ�ÿ����ֲ��õ����Ƶ�pֵ
  #         'permutation': ʹ��permutation������pֵ
  #         'both': �ṩ�����ֲ��õ����Ƶ�pֵ��permutation�����pֵ
  #n_per: permutation���̵Ĵ�����Ĭ����1000
  #em_max: ����������EM��������, Ĭ��ֵ��100
  #tolerance: ǰһ�ε����ͺ�һ�ε�����������������ֵ��Ĭ��ֵΪ0.00001
  
  
  #n_phe:����-������Ŀ��������Ϊn_phe[1]��������Ϊn_phe[2]
  n_phe <- table(Phenotype)
  
  #n:������; prod_n_phe:����������������; m: λ����
  n <- sum(n_phe)
  prod_n_phe <- prod(n_phe)
  m <- ncol(Genotype)
  
  #Uͳ��������������Ͷ�����֮���λ��ƽ�������͵÷ֵĲ���
  U <- prod_n_phe*(colSums(Genotype[Phenotype==1,])/n_phe[2]-colSums(Genotype[Phenotype==0,])/n_phe[1])/n
  
  # #Uͳ��������һ�ֱ�ʾ
  # G_I <- t(Genotype)%*%(diag(n)-(rep(1,n)%*%t(rep(1,n)))/n)
  # #U <- G_I%*%Phenotype
  # U <- rowSums(G_I[,Phenotype])
  
  if(LD_estimate=='sample'){
    
    #ʹ������������������LDֵ
    #ginv_Cov_U_tem: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_tem <- n*ginv(LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_tem: SLIDE ͳ����
    T_tem <- as.vector(U%*%ginv_Cov_U_tem%*%U)
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_tem_asy <- pchisq(T_tem,m,lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_tem_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_tem_per <- U_per%*%ginv_Cov_U_tem%*%U_per
    })
    
    #ʹ��permutation����������pֵ
    p_value_tem_per <- sum(T_tem_per>=T_tem)/n_per
    
    SLIDE_tem <- c(T_tem,p_value_tem_asy,p_value_tem_per)
    names(SLIDE_tem) <- c('T_tem','p_value_tem_asy','p_value_tem_per')
    return(SLIDE_tem)
    
  } else if(LD_estimate=='control'){
    
    #ʹ�������Ķ�������������LDֵ
    #ginv_Cov_U_tem: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_con <- n*ginv(LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_con: SLIDE ͳ����
    T_con <- as.vector(U%*%ginv_Cov_U_con%*%U)
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_con_asy <- pchisq(T_con,m,lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_con_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_con_per <- U_per%*%ginv_Cov_U_con%*%U_per
    })
    
    #ʹ��permutation����������pֵ
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
    
    
    
    #ʹ���ⲿ����Ⱥ����������LDֵ
    #ginv_Cov_U_pop: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_pop <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    #T_pop: SLIDE ͳ����
    T_pop <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U%*%ginv_Cov_U_popi%*%U})
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_pop_asy <- pchisq(T_pop,m, lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_pop_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_pop_per <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U_per%*%ginv_Cov_U_popi%*%U_per})
    })
    
    #ʹ��permutation����������pֵ
    p_value_pop_per <- rowSums(T_pop_per>=T_pop)/n_per
    
    SLIDE_pop<- cbind(T_pop,p_value_pop_asy,p_value_pop_per)
    rownames(SLIDE_pop) <- paste('SLIDE_pop',1:length(LD_matrix),sep = "")
    return(SLIDE_pop)
    
  }else if(LD_estimate=='all'){
    
    #�ṩ'sample'��'control'��'population'����Ӧ�Ľ��
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    LD_matrix$ld_sap <- LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0
    LD_matrix$ld_con <- LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0
    
    ginv_Cov_U <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    SLIDE_T <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U%*%ginv_Cov_Ui%*%U})
    
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_asy <- pchisq(SLIDE_T,m, lower.tail = FALSE)
    
    #����һ��3*n_per�ľ���
    #��һ���ǻ����������ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    #�ڶ����ǻ��������Ķ��������ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    #�������ǻ�����Ⱥ���ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    T_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      
      T_per <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U_per%*%ginv_Cov_Ui%*%U_per})
      return(T_per)
    })
    
    #ʹ��permutation����������pֵ,����һ������Ϊ3������
    #������Ԫ�طֱ���T_tem��T_com��T_pop��Ӧ��pֵ
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
  
  #Genotype: �����;���Ů�Ը�ֵΪ0/1/2�����Ը�ֵΪ0/1
  #Phenotype:��������, �������帳ֵΪ1��δ�������帳ֵΪ0
  #sex: �Ա����������ڲ�ֻ����;������Ը�ֵΪ0, Ů�Ը�ֵΪ1
  #LD_estimate: ���Ƴ�����������G����Ӧ�ķ���-Э�������ķ�������Ϊc('control','population','both)
  #             'sample': ʹ������������������LDֵ
  #             'control': ʹ�������Ķ�������������LDֵ
  #             'population': ʹ���ⲿ����Ⱥ����������LDֵ����ѡ���ѡ������LD_matrix��ֵ
  #             'all': �ṩ'sample'��'control'��'population'����Ӧ�Ľ��
  #LD_matrix: LD_estimate='population'ʱ���ṩ, �����ⲿ����Ⱥ���ݹ��Ƴ�����LD����
  
  #gamma: ƫ��ʧ��ϵ�����ܵ�ȡֵ��Ĭ����0.0��0.5��1.0��1.5��2.0
  #ʹ�������������ҵ�����ͳ����
  
  #p_value: ����pֵ�ķ�������Ϊc('Asymptotically','permutation','both)��
  #         'Asymptotically': ʹ�ÿ����ֲ��õ����Ƶ�pֵ
  #         'permutation': ʹ��permutation������pֵ
  #         'both': �ṩ�����ֲ��õ����Ƶ�pֵ��permutation�����pֵ
  #n_per: permutation���̵Ĵ�����Ĭ����1000
  #em_max: ����������EM��������, Ĭ��ֵ��100
  #tolerance: ǰһ�ε����ͺ�һ�ε�����������������ֵ��Ĭ��ֵΪ0.00001
  
  
  #n_phe:����-������Ŀ��������Ϊn_phe[1]��������Ϊn_phe[2]
  n_phe <- table(Phenotype)
  
  #n:������; prod_n_phe:����������������; m: λ����
  n <- sum(n_phe)
  prod_n_phe <- prod(n_phe)
  m <- ncol(Genotype)
  
  #Uͳ��������������Ͷ�����֮���λ��ƽ�������͵÷ֵĲ���
  U <- prod_n_phe*(colSums(Genotype[Phenotype==1,])/n_phe[2]-colSums(Genotype[Phenotype==0,])/n_phe[1])/n
  
  # #Uͳ��������һ�ֱ�ʾ
  # G_I <- t(Genotype)%*%(diag(n)-(rep(1,n)%*%t(rep(1,n)))/n)
  # #U <- G_I%*%Phenotype
  # U <- rowSums(G_I[,Phenotype])
  
  if(LD_estimate=='sample'){
    
    #ʹ������������������LDֵ
    #ginv_Cov_U_tem: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_tem <- n*ginv(LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_tem: SLIDE ͳ����
    T_tem <- as.vector(U%*%ginv_Cov_U_tem%*%U)
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_tem_asy <- pchisq(T_tem,m,lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_tem_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_tem_per <- U_per%*%ginv_Cov_U_tem%*%U_per
    })
    
    #ʹ��permutation����������pֵ
    p_value_tem_per <- sum(T_tem_per>=T_tem)/n_per
    
    SLIDE_tem <- c(T_tem,p_value_tem_asy,p_value_tem_per)
    names(SLIDE_tem) <- c('T_tem','p_value_tem_asy','p_value_tem_per')
    return(SLIDE_tem)
    
  } else if(LD_estimate=='control'){
    
    #ʹ�������Ķ�������������LDֵ
    #ginv_Cov_U_tem: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_con <- n*ginv(LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0)/(2*prod_n_phe)
    #T_con: SLIDE ͳ����
    T_con <- as.vector(U%*%ginv_Cov_U_con%*%U)
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_con_asy <- pchisq(T_con,m,lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_con_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_con_per <- U_per%*%ginv_Cov_U_con%*%U_per
    })
    
    #ʹ��permutation����������pֵ
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
    
    
    
    #ʹ���ⲿ����Ⱥ����������LDֵ
    #ginv_Cov_U_pop: Uͳ�����ķ���-Э�������������
    ginv_Cov_U_pop <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    #T_pop: SLIDE ͳ����
    T_pop <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U%*%ginv_Cov_U_popi%*%U})
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_pop_asy <- pchisq(T_pop,m, lower.tail = FALSE)
    
    #����һ������Ϊn_per������, ��n_per�ε�permutationͳ����
    T_pop_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      T_pop_per <- sapply(ginv_Cov_U_pop,function(ginv_Cov_U_popi){U_per%*%ginv_Cov_U_popi%*%U_per})
    })
    
    #ʹ��permutation����������pֵ
    p_value_pop_per <- rowSums(T_pop_per>=T_pop)/n_per
    
    SLIDE_pop<- cbind(T_pop,p_value_pop_asy,p_value_pop_per)
    rownames(SLIDE_pop) <- paste('SLIDE_pop',1:length(LD_matrix),sep = "")
    return(SLIDE_pop)
    
  }else if(LD_estimate=='all'){
    
    #�ṩ'sample'��'control'��'population'����Ӧ�Ľ��
    
    if(is.matrix(LD_matrix)){
      LD_matrix <- list(LD_matrix)
    }else if(!is.list(LD_matrix)){
      stop("'LD_matrix' must be a matrix or list")
    }
    
    LD_matrix$ld_sap <- LD.f(Genotype,em_max=em_max,tolerance=tolerance)$LD_f0
    LD_matrix$ld_con <- LD.f(Genotype[Phenotype==0,],em_max=em_max,tolerance=tolerance)$LD_f0
    
    ginv_Cov_U <- lapply(LD_matrix,function(LD_matrixi){n*ginv(LD_matrixi)/(2*prod_n_phe)})
    SLIDE_T <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U%*%ginv_Cov_Ui%*%U})
    
    
    #ʹ�����ɶ�Ϊm�Ŀ����ֲ��õ�������pֵ
    p_value_asy <- pchisq(SLIDE_T,m, lower.tail = FALSE)
    
    #����һ��3*n_per�ľ���
    #��һ���ǻ����������ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    #�ڶ����ǻ��������Ķ��������ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    #�������ǻ�����Ⱥ���ݹ��Ƴ�LD��Ӧ��n_per�ε�permutationͳ����
    T_per <- sapply(1:n_per,function(i){
      #permutationͳ����
      #�Ա����س���
      samp_Phe <- sample(Phenotype)
      U_per <- prod_n_phe*(colSums(Genotype[samp_Phe==1,])/n_phe[2]-colSums(Genotype[samp_Phe==0,])/n_phe[1])/n
      
      T_per <- sapply(ginv_Cov_U,function(ginv_Cov_Ui){U_per%*%ginv_Cov_Ui%*%U_per})
      return(T_per)
    })
    
    #ʹ��permutation����������pֵ,����һ������Ϊ3������
    #������Ԫ�طֱ���T_tem��T_com��T_pop��Ӧ��pֵ
    p_value_per <- rowSums(T_per>=SLIDE_T)/n_per
    
    SLIDE <- cbind(SLIDE_T,p_value_asy,p_value_per)
    rownames(SLIDE) <- c(paste('SLIDE_pop',1:(length(LD_matrix)-2),sep = ""),
                         'SLIDE_tem','SLIDE_con')
    
    return(SLIDE)
  }
  
}