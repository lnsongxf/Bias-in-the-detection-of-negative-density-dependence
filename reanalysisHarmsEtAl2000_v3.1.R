# PART definitions for litter files (for reference) 
# 0  reproductive buds
# 1  mature fruit
# 2  seeds
# 3  capsules
# 4  fragments
# 5  immature fruit
# 6  flowers (female or perfect)
# 7  fruit with insect emergence holes / insect damaged
# 8  aborted fruit
# 9  male flowers
#10  fruto comido por animal (ave, mono, ardilla y otros).
#11  categor?a para hojas 
#12  "desconocido"
#13  pedicelo de hojas
#14  pedazos de ramas, less than 2 cm in diameter
#15  basura fina
#16  pedicelo de flor
#17  pedicelo de fruto
#18  caca de mono

## Prepare and load data
load("pairedseedtrapsApril2018.Rdata")
SeedEq = read.csv("SeedEquivalent.csv")
harmsdat <- read.table("Nature2000_DataRecreated_20180404.txt",header=TRUE)

require(MASS)
source('rbinom2.r')

# subset and pair traps
pairtraps=subset(trapdata,
                 (trap<=299|trap>=400) & year>=2011 & fruit==1 & (year==2011&month>=3) & (part==1|part==2|part==7|part==10) )


Atraps=subset(pairtraps,trap<=200)
Btraps=subset(pairtraps,trap>=201)

Atraps$trap <- as.factor(Atraps$trap)
Btraps$trap <- as.factor(Btraps$trap)


spp <- table(Atraps$sp)
N <- spp[spp>10]
spp <- names(spp)[spp>10]

SE=numeric(length(spp))+1
for(i in 1:length(spp)){
  
  use=SeedEq[,1]==spp[i]
  if (sum(use)==1)
    SE[i] = SeedEq[use,2]
  
}



results <- data.frame(error=as.numeric(N),
                      N=as.numeric(N),
                      rho=as.numeric(N),
                      spp=spp)

for(i in 1:length(spp)){
  
  sp1A <- subset(Atraps,sp==spp[i])
  sp1B <- subset(Btraps,sp==spp[i])
  
  if (is.na(SE[i])==FALSE)
  {
    sp1A$se <- ifelse(sp1A$part==1,sp1A$quan*SE[i],sp1A$quan) 
    sp1B$se <- ifelse(sp1B$part==1,sp1B$quan*SE[i],sp1B$quan)
  } else {
    sp1A$se <- sp1A$quan
    sp1B$se <- sp1B$quan
  }
  
  
  sptrapnA<-as.numeric(xtabs(sp1A$se~sp1A$trap))
  sptrapnB<-as.numeric(xtabs(sp1B$se~sp1B$trap))
  
  results$rho[i] <- cor(sptrapnA,sptrapnB)
  results$spp[i] <- spp[i] 
  results$N[i]   <- mean(sptrapnA+sptrapnB)/2
  
  results$error[i] <- mean(abs(log(sptrapnA+1)-log(sptrapnB+1))/sqrt(2))
  

}


################################################################################
## Re-analysis of the Harms et al 2000 Nature paper
################################################################################


spcodes <- as.character(unique(harmsdat$sp))
Rep <- 1000
mod.rnd0=mod.rnd1=mod.rnd2=mod <- numeric(Rep)
H0=H1=H2 <- data.frame(b1=numeric(length(spcodes)),
                      lb1=numeric(length(spcodes)),
                      ub1=numeric(length(spcodes)))
                      
r0=r=p0=p=rho=b1=f=N=bhat=numeric(length(spcodes))
rho.avg=median(results$rho, na.rm=TRUE)
Z=numeric(Rep)


for(i in 1:length(spcodes)) {
  temp <- subset(harmsdat,sp==spcodes[i])
  R=as.numeric(temp$R*4*3)    # total number of recruits measured per site
  S=as.numeric(temp$S*4/2)    # total number of equilvalent seeds measured per trap
  
  N[i]=mean(temp$S)
  use=results$spp==spcodes[i]
  if (sum(use)>0){
  rho[i]=abs(results$rho[use])
  
  } else {
    rho[i]=NA
  }

  f[i]=sum(R)/sum(6*S)
  
  # parameter estimation of the negative binomial distribution using method of moments
  m = mean(S)
  V = var(S)
  r = m^2/(V-m)
  p = m/V
  
    if (is.na(rho[i])){
      alfa = rho.avg*r

  } else {
    
    alfa=rho[i]*r
    }
      
      
  
  # fit model and get confidence intervals
  S.0 = S*2/4     # per m-2 yr-1
  R.0 = R/3/4     # per m-2 yr-1
  
  # use0 = R.0>0|S.0>0
  use0 = S.0>=0
  fit = lm(log(R.0[use0]+1)~log(S.0[use0]+1))
  H0[i,] <- c(coef(fit)[2], confint(fit,'log(S.0 + 1)',level=0.95))

  
    for(j in 1:Rep){
      cat("\r",i,":",Rep-j,"\r")
    
    # null model 1
    R.r = rbinom2(200,S*6,f[i])
    R.0 = R.r/3/4   # per m-2 yr-1
    mod.rnd1[j] <- coef(lm(log(R.0[use0]+1)~log(S.0[use0]+1)))[2]
    
    # null model 2 (collocation error)
    S.r=0
    A <- rbeta(200,alfa,r-alfa)
    z=numeric(200)
    for (k in 1:200){
      if (S[k]>0 & A[k]<1)
        z[k] <- rbinom2(1,6*S[k],A[k])
    }
    
    S.r <- z+rnbinom(200,6*(r-alfa),p)
    R.r <- rbinom(200,S.r,f[i])
    R.0 = R.r/3/4   # per m-2 yr-1
    mod.rnd2[j] <- coef(lm(log(R.0[use0]+1)~log(S.0[use0]+1)))[2]
    
    }
  

  H1[i,] <- c(mean(mod.rnd1,na.rm=TRUE),quantile(mod.rnd1,prob=c(0.025,0.975),na.rm=TRUE))
  H2[i,] <- c(mean(mod.rnd2,na.rm=TRUE),quantile(mod.rnd2,prob=c(0.025,0.975),na.rm=TRUE))

}



## Make Harms et al re-analysis Figure
b1=H0[,1]
pdf("HarmsFigure.pdf",width=8,height=8)

## null models for all species
par(bty="n",mar=c(5,6,2,2),cex.lab=1.2)

layoutmatrix <- matrix(c(
    1,1,1,1,1,2,2,2,2,
    1,1,1,1,1,2,2,2,2,
    1,1,1,1,1,3,3,3,3,
    1,1,1,1,1,3,3,3,3,
    1,1,1,1,1,4,4,4,4,
    1,1,1,1,1,4,4,4,4),
    ncol=9,byrow=TRUE)
layout(layoutmatrix)
    
spN <- 1:length(spcodes)
xl=c(-1.2,.5)
plot(b1-1,spN,pch=16,cex=1,xlim=c(xl[1],xl[2]),yaxt="n",ylab="",xlab=expression(italic(b)[OLS] - italic(b)['NULL']),main='A')

text(y=seq(1,52),  -1.1, 
     labels = spcodes,srt = 0, pos = 2,cex=.7)


points(b1-H1$b1,spN-0.1,pch=16,cex=1,xlim=xl,yaxt="n",col='red')

points(b1-H2$b1,spN+0.1,pch=16,cex=1,xlim=xl,yaxt="n",col='blue')
# legend(x= "topleft",c('H0','H1','H2'))

segments(H0$ub1-1,spN,H0$lb1-1,spN)
segments(b1-H1$ub1,spN-0.1,b1-H1$lb1,spN-0.1,col='red')
segments(b1-H2$ub1,spN+0.1,b1-H2$lb1,spN+0.1,col='blue')
abline(v=0); abline(v=-1.1)


# plot with rho and germination probability
hist(results$rho[as.character(results$spp)%in%spcodes], xlab='paired trap correlation',
     ylab='number of species',main='B',
     breaks=10,freq=TRUE)
# lines(density(na.omit(results$rho)),col="grey",lwd=4)

plot(rho,b1-H1[,1],main='C',xlab='paired trap correlation',ylab=expression(italic(b)[OLS] - H['1']))
abline(lm(b1-H1[,1]~rho))
# abline(h=0,lty=2,col="grey")
m0 <- (lm(b1-H1[,1]~rho))
R2 <- signif(summary(m0)$r.sq,2)
legend("bottomright",legend=bquote(R^2==.(R2)*"***"),bty="n")


kk=(N+1)/(N+1/f)*rho
# plot(kk,b1,main='D',ylab='b',xlab=bquote(frac(1+bar(S)^-1,1+bar(R)^-1)*" "*rho))
plot(kk,b1,main='D',ylab=expression(italic(b)[OLS]),xlab='analytical model')
mod0=lm(b1~kk)
# abline(mod0)
abline(0,1)
summary(mod0)
R2 <- signif(summary(mod0)$r.sq,2)
legend("bottomright",legend=bquote(R^2==.(R2)*"***"),bty="n")

dev.off()
