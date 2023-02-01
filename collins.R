library(pracma)
mutualinformation = function(state,action,alpha){
  uS = unique(state);uA = unique(action)
  N = matrix(0,length(uS),length(uA))
  for(x in 1:length(uS)){
    for(y in 1:length(uA)){
      N[x,y] = alpha+sum(state==x & action==y)
    }
  }
  n = sum(N)
  nA = colSums(N)
  nS = rowSums(N)
  psi1 = sweep(digamma(N+1),2,digamma(nA+1))
  psi2 = sweep(psi1,2,digamma(nS+1))
  P =  sweep(psi2,2,digamma(n+1),FUN="+")
  I = sum(colSums(N*P))/n
}
logsumexp = function(x,N){
  y = x[order(x,decreasing=T)[1:N]]
  xx = x - y
  s = y + log(rowSums(exp(xx)))
  return(s)
}
beta = seq(0.1,15,length.out=50)
blahut_arimoto = function(Ps,Q,b){
  A = ncol(Q)
  niter = 50
  R = V = rep(0,length(b)); Pa = matrix(0,length(b),A)
  
  for(j in 1:length(b)){
    FF = b[j]*Q
    v0 = mean(colMeans(Q))
    q = rep(1,A) / A
    for(g in 1:niter){
      logP = sweep(FF,2,log(q),FUN="+")
      Z = logsumexp(logP,3)
      Psa = exp(logP-Z)
      q = Ps%*%Psa
      v = sum(colSums(Ps*(Psa*Q)))
      v0=v
    }
    Pa[j,] = q
    V[j] = v
    R[j] = b[j]*v - Ps%*%Z
  }
  return(list(V=V,
              R=R))
}

data = read.csv("~/Dropbox (University of Oregon)/Collins18_data.csv")
subs = unique(data$ID)

V=R = list()

V_dat=R_dat=matrix(0,length(subs),2)
Varr=Rarr=array(0,dim=c(length(subs),50,2))
for(sub in 1:length(subs)){
  s1 = data[data$ID==subs[sub],]
  B = unique(s1$learningblock)
  cond = rep(0,length(B))
  R_data = rep(0,length(B))
  V_data = rep(0,length(B))
  V=R=matrix(0,14,50)
  for(bb in 1:length(B)){
    ix = s1[s1$learningblock==B[bb] & s1$phase==0,]
    ix = ix[!ix$choice==-1,]
    state = ix$stim
    cc = ix$corchoice
    action = ix$choice
    R_data[bb] = mutualinformation(state,action,0.1)
    V_data[bb] = mean(ix$cor,na.rm=T)
    
    S = unique(state)
    Q = matrix(0,length(S),3)
    Ps = rep(0,length(S))
    for(i in 1:length(S)){
      ii = state==S[i]
      Ps[i] = mean(ii,na.rm=T)
      a = cc[ii]; a = a[1]
      Q[i,a] = 1
    }
    ba = blahut_arimoto(Ps,Q,beta)
    V[bb,] = ba$V; R[bb,] = ba$R
    if(length(S)==3)cond[bb] = 1
    if(length(S)!=3)cond[bb] = 2
  }
  
  for(ccc in 1:2){
    Varr[sub,,ccc] = colMeans(V[cond==ccc,],na.rm=T)
    Rarr[sub,,ccc] = colMeans(R[cond==ccc,],na.rm=T)
    V_dat[sub,ccc] = mean(V_data[cond==ccc],na.rm=T)
    R_dat[sub,ccc] = mean(R_data[cond==ccc],na.rm=T)
  }
}

pracma::interp1(Rarr[1,,1],Varr[1,,1],R_dat[,1])


plot(Rarr[,,1],Varr[,,1],ylim=c(0.25,1.2),xlim=c(0,0.9),type="l")
lines(Rarr[,,2],Varr[,,2])
points(R_dat[,1],V_dat[,1])
points(R_dat[,2],V_dat[,2])
