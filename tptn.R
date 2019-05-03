library(corrplot)
mat1=(abs(as.matrix(read.table('COVGGMn100-s200-c3-d6-nwc0-corr60vary.txt')))[1:18,1:18])
mat2=(abs(as.matrix(read.table('COVGGMn100-s200-c3-d5-nwc0-corr60vary.txt')))[1:18,1:18])
mat3=(abs(as.matrix(read.table('COVGGMn100-s200-c3-d4-nwc0-corr60vary.txt')))[1:18,1:18])
colnames(mat1)=c(1:18);colnames(mat2)=c(1:18);colnames(mat3)=c(1:18)
par(mfrow=c(1,3))
corrplot(mat1)
corrplot(mat2)
corrplot(mat3)

tptn<-function(X,covar, neigh,thresh){
X=abs(X)

tempX= abs(X[1:(covar*neigh),1:(covar*neigh)]); tempX=(tempX>thresh)+0
wantX= matrix(0,covar*neigh, covar*neigh); for(iii in 1:covar){wantX[(1+(iii-1)*neigh):(1+(iii-1)*neigh+neigh-1) , (1+(iii-1)*neigh):(1+(iii-1)*neigh+neigh-1) ]=1}
#print(wantX)
diag(tempX)=100
diag(wantX)=100 
resmat=wantX*2+tempX
tp=sum(resmat==3)/2 
tn=sum(resmat==0)/2
fp=sum(resmat==1)/2
fn=sum(resmat==2)/2
tpr=tp/(tp+fn); tnr=tn/(tn+fp)
return(list(tpr,tnr,tp,tn,fp,fn))

}

tprtnror <- function(covar,neigh, ori, thre1, sol, thre2){
  ori=abs(ori)
  sol=abs(sol)
  ori[1:(covar*neigh),1:(covar*neigh)]=0 
  sol[1:(covar*neigh),1:(covar*neigh)]=0 
  ori[upper.tri(ori)] <- t(ori)[upper.tri(ori)]
  sol[upper.tri(sol)] <- t(sol)[upper.tri(sol)]
  ori=(ori>thre1)
  sol=(sol>thre2)
  diag(ori)=100
  diag(sol)=100
  #print(dim(ori))
  #print(dim(sol))
  resmat=ori*2+sol
  tp=sum(resmat==3)/2 
  tn=sum(resmat==0)/2-(covar*neigh)*(covar*neigh-1)/2
  fp=sum(resmat==1)/2
  fn=sum(resmat==2)/2
  tpr= tp/(tp+fn)
  tnr=tn/(tn+fp)
  return(list(tpr,tnr,tp,tn,fp,fn))
}

tprtnrmm <- function(covar,neigh,x,thresh){
  x=abs(x)
  tp=0;tn=0;fp=0;fn=0
  for(i in 1:covar){
      tp=tp+sum(x[i, ((1+(i-1)*neigh):(neigh+(i-1)*neigh))]>thresh)
      fn=fn+neigh-sum(x[i, ((1+(i-1)*neigh):(neigh+(i-1)*neigh))]>thresh)
      tn=sum(x[i,-((1+(i-1)*neigh):(neigh+(i-1)*neigh))]<thresh)+tn
      fp=fp+dim(x)[2]-neigh-sum(x[i,-((1+(i-1)*neigh):(neigh+(i-1)*neigh))]<thresh)
  }
  tpr=tp/(tp+fn)
  tnr=tn/(tn+fp)
  return(list(tpr,tnr,tp,tn,fp,fn))
}



##################################
###########

for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 6; cv=3; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  justggm=abs(as.matrix(read.table(paste('ONLYGGMn',f,'-s',n,'-c',0,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.001)
  covbin1=tptn(covbin,cv,d,0.001)
  nondic1=tptn(noncovdichbin,cv,d,0.001)
  covdic1=tptn(covdichbin,cv,d,0.001)
  nonggm1=tptn(noncovggm,cv,d,0.001)
  covggm1=tptn(covggm,cv,d,0.001)
  justggm1=tptn(justggm,cv,d,0.001)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.001)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.001)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.001)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.001)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.001)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.001)
  justggm2=tprtnror(cv,d,ggmoriginal,0.09,justggm,0.001)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  #print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  #print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  #print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  #print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  #print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  #print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
  
  print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]],justggm1[[1]]+justggm1[[2]])/2)
  print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]],justggm2[[1]]+justggm2[[2]])/2)
  print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)
  
}
######
#######
for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 5; cv=3; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  justggm=abs(as.matrix(read.table(paste('ONLYGGMn',f,'-s',n,'-c',0,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.001)
  covbin1=tptn(covbin,cv,d,0.001)
  nondic1=tptn(noncovdichbin,cv,d,0.001)
  covdic1=tptn(covdichbin,cv,d,0.001)
  nonggm1=tptn(noncovggm,cv,d,0.001)
  covggm1=tptn(covggm,cv,d,0.001)
  justggm1=tptn(justggm,cv,d,0.001)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.001)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.001)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.001)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.001)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.001)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.001)
  justggm2=tprtnror(cv,d,ggmoriginal,0.09,justggm,0.001)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  #print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  #print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  #print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  #print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  #print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  #print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
  
  print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]],justggm1[[1]]+justggm1[[2]])/2)
  print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]],justggm2[[1]]+justggm2[[2]])/2)
  print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)
  
}




#########
##########
for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 4; cv=3; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  justggm=abs(as.matrix(read.table(paste('ONLYGGMn',f,'-s',n,'-c',0,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.001)
  covbin1=tptn(covbin,cv,d,0.001)
  nondic1=tptn(noncovdichbin,cv,d,0.001)
  covdic1=tptn(covdichbin,cv,d,0.001)
  nonggm1=tptn(noncovggm,cv,d,0.001)
  covggm1=tptn(covggm,cv,d,0.001)
  justggm1=tptn(justggm,cv,d,0.001)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.001)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.001)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.001)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.001)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.001)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.001)
  justggm2=tprtnror(cv,d,ggmoriginal,0.09,justggm,0.001)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  #print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  #print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  #print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  #print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  #print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  #print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
  
  print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]],justggm1[[1]]+justggm1[[2]])/2)
  print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]],justggm2[[1]]+justggm2[[2]])/2)
  print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)
  
}

###

f=100; n= 200; d= 10; cv=1; nwc=10; corrnwc=20#d*nwc/5

binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}

noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}

ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
justggm=abs(as.matrix(read.table(paste('ONLYGGMn',f,'-s',n,'-c',0,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]


#1) CA
nonbin1=tptn(noncovbin,cv,d,0.01)
covbin1=tptn(covbin,cv,d,0.01)
nondic1=tptn(noncovdichbin,cv,d,0.01)
covdic1=tptn(covdichbin,cv,d,0.01)
nonggm1=tptn(noncovggm,cv,d,0.01)
covggm1=tptn(covggm,cv,d,0.01)
justggm1=tptn(justggm,cv,d,0.01)

#2) original 

nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
justggm2=tprtnror(cv,d,ggmoriginal,0.09,justggm,0.01)

#3) metal-metabolite

nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
covbin3=tprtnrmm(cv,d,covbinextra,0.09)
nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
covggm3=tprtnrmm(cv,d,covggmextra,0.09)

print(paste('sample',idk))
#print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
#print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
#print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
#print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
#print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
#print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))

print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]],justggm1[[1]]+justggm1[[2]])/2)
print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]],justggm2[[1]]+justggm2[[2]])/2)
print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)












####################33



###########

for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 5; cv=3; nwc=1; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.011)
  covbin1=tptn(covbin,cv,d,0.011)
  nondic1=tptn(noncovdichbin,cv,d,0.011)
  covdic1=tptn(covdichbin,cv,d,0.011)
  nonggm1=tptn(noncovggm,cv,d,0.011)
  covggm1=tptn(covggm,cv,d,0.011)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  #print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  #print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  #print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  #print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  #print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  #print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
  
  print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]])/2)
  print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]])/2)
  print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)
  
}


#################
#############################


for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 5; cv=3; nwc=2; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.01)
  covbin1=tptn(covbin,cv,d,0.01)
  nondic1=tptn(noncovdichbin,cv,d,0.01)
  covdic1=tptn(covdichbin,cv,d,0.01)
  nonggm1=tptn(noncovggm,cv,d,0.01)
  covggm1=tptn(covggm,cv,d,0.01)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  #print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  #print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  ##print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  #print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  ##print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  #print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
  
  print(c(nonbin1[[1]]+nonbin1[[2]],covbin1[[1]]+covbin1[[2]],nondic1[[1]]+nondic1[[2]],covdic1[[1]]+covdic1[[2]],nonggm1[[1]]+nonggm1[[2]],covggm1[[1]]+covggm1[[2]])/2)
  print(c(nonbin2[[1]]+nonbin2[[2]],covbin2[[1]]+covbin2[[2]],nondic2[[1]]+nondic2[[2]],covdic2[[1]]+covdic2[[2]],nonggm2[[1]]+nonggm2[[2]],covggm2[[1]]+covggm2[[2]])/2)
  print(c(nonbin3[[1]]+nonbin3[[2]],covbin3[[1]]+covbin3[[2]],nondic3[[1]]+nondic3[[2]],covdic3[[1]]+covdic3[[2]],nonggm3[[1]]+nonggm3[[2]],covggm3[[1]]+covggm3[[2]])/2)
  
}

############################## nwc=4


for (idk in c(110,140,170,200,230,260)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= idk; d= 5; cv=3; nwc=4; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovbinextra)){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covbinextra)){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovdichbinextra)){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covdichbinextra)){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(noncovggmextra)){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(covggmextra)){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.01)
  covbin1=tptn(covbin,cv,d,0.01)
  nondic1=tptn(noncovdichbin,cv,d,0.01)
  covdic1=tptn(covdichbin,cv,d,0.01)
  nonggm1=tptn(noncovggm,cv,d,0.01)
  covggm1=tptn(covggm,cv,d,0.01)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('sample',idk))
  print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
}


##################### 


for (idk in c(3:15)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= 200; d= idk; cv=1; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovbinextra))){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covbinextra))){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovdichbinextra))){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covdichbinextra))){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovggmextra))){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covggmextra))){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.01)
  covbin1=tptn(covbin,cv,d,0.01)
  nondic1=tptn(noncovdichbin,cv,d,0.01)
  covdic1=tptn(covdichbin,cv,d,0.01)
  nonggm1=tptn(noncovggm,cv,d,0.01)
  covggm1=tptn(covggm,cv,d,0.01)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('fea',idk))
  print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
}

ising=c(0,0,0,0,0,0,0,0,0,0)
liising= c(0.33,0,0,0,0,0,0,0,0,0)
ggm=c(1,1,0.5,0,0,0,0,0,0,0)
liggm=c(1,1,1,1,0.8,0.4,0.25,0.11,0.054,0.03)

plot(ising,ylim=c(0,1), x=c(3:12) ,type = "o", col = "red", xlab = "Degree", ylab = "True Positive Rate of Detecting Latent Interactions",main = "TPR vs Degree")
  lines(ising, x=c(3:12), type="b",lty=1, lwd=3.5, col='red') 
  lines(liising, x=c(3:12), type="b",lty=2, lwd=3.5, col='black') 
  lines(ggm, x=c(3:12), type="b", lwd=3.5, lty=3, col='blue')
  lines(liggm, x=c(3:12), type="b", lwd=3.5, lty=4, col='green')
  legend(c(10,12), c(0.7,1),c('Ising', 'Integrated \n Ising','GGM','IGGM'),lty=c(1,2,3,4), cex=1, col=c('red','black','blue','green'), title="Models")

### feat= 100 120 140 160 180 

for (idk in c(100,120,140,160,180)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=idk; n= 200; d= 5; cv=3; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovbinextra))){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covbinextra))){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovdichbinextra))){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covdichbinextra))){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovggmextra))){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covggmextra))){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.01)
  covbin1=tptn(covbin,cv,d,0.01)
  nondic1=tptn(noncovdichbin,cv,d,0.01)
  covdic1=tptn(covdichbin,cv,d,0.01)
  nonggm1=tptn(noncovggm,cv,d,0.01)
  covggm1=tptn(covggm,cv,d,0.01)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('fea',idk))
  print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
}




##################### cov= 1-5


for (idk in c(1:5)){
  # varied coefficients gaucorr-n100s200c3d5 nwc=0 
  f=100; n= 200; d= 5; cv=idk; nwc=0; corrnwc=d*nwc/5
  
  binoriginal=abs(as.matrix(read.table(paste('bincorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));binoriginal=binoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovbin= abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovbinextra=abs(as.matrix(read.table(paste('NONCOVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovbinextra))){noncovbinextra= t(as.matrix(noncovbinextra))}
  covbin=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covbinextra=abs(as.matrix(read.table(paste('COVBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covbinextra))){covbinextra= t(as.matrix(covbinextra))}
  
  noncovdichbin= abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovdichbinextra=abs(as.matrix(read.table(paste('NONCOVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovdichbinextra))){noncovdichbinextra= t(as.matrix(noncovdichbinextra))}
  covdichbin=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covdichbinextra=abs(as.matrix(read.table(paste('COVDICHBINn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covdichbinextra))){covdichbinextra= t(as.matrix(covdichbinextra))}
  
  ggmoriginal=abs(as.matrix(read.table(paste('gaucorr-n',f,'s',n,'c',cv,'d',d,'cor0.65','nwc',corrnwc,'.txt',sep='') )));ggmoriginal=ggmoriginal[(cv+1):(f+cv),(cv+1):(f+cv)]
  noncovggm= abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  noncovggmextra=abs(as.matrix(read.table(paste('NONCOVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(noncovggmextra))){noncovggmextra= t(as.matrix(noncovggmextra))}
  covggm=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[1:f,1:f]
  covggmextra=abs(as.matrix(read.table(paste('COVGGMn',f,'-s',n,'-c',cv,'-d',d,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))))[(f+1):(f+cv),1:f]
  if(is.null(dim(covggmextra))){covggmextra= t(as.matrix(covggmextra))}
  
  #1) CA
  nonbin1=tptn(noncovbin,cv,d,0.01)
  covbin1=tptn(covbin,cv,d,0.01)
  nondic1=tptn(noncovdichbin,cv,d,0.01)
  covdic1=tptn(covdichbin,cv,d,0.01)
  nonggm1=tptn(noncovggm,cv,d,0.01)
  covggm1=tptn(covggm,cv,d,0.01)
  
  #2) original 
  
  nonbin2=tprtnror(cv,d,binoriginal,0.09,noncovbin,0.01)
  covbin2=tprtnror(cv,d,binoriginal,0.09,covbin,0.01)
  nondic2=tprtnror(cv,d,binoriginal,0.09,noncovdichbin,0.01)
  covdic2=tprtnror(cv,d,binoriginal,0.09,covdichbin,0.01)
  nonggm2=tprtnror(cv,d,ggmoriginal,0.09,noncovggm,0.01)
  covggm2=tprtnror(cv,d,ggmoriginal,0.09,covggm,0.01)
  
  #3) metal-metabolite
  
  nonbin3=tprtnrmm(cv,d,noncovbinextra,0.09)
  covbin3=tprtnrmm(cv,d,covbinextra,0.09)
  nondic3=tprtnrmm(cv,d,noncovdichbinextra,0.09)
  covdic3=tprtnrmm(cv,d,covdichbinextra,0.09)
  nonggm3=tprtnrmm(cv,d,noncovggmextra,0.09)
  covggm3=tprtnrmm(cv,d,covggmextra,0.09)
  
  print(paste('fea',idk))
  print(c(nonbin1[[1]],covbin1[[1]],nondic1[[1]],covdic1[[1]],nonggm1[[1]],covggm1[[1]]))
  print(c(nonbin1[[2]],covbin1[[2]],nondic1[[2]],covdic1[[2]],nonggm1[[2]],covggm1[[2]]))
  print(c(nonbin2[[1]],covbin2[[1]],nondic2[[1]],covdic2[[1]],nonggm2[[1]],covggm2[[1]]))
  print(c(nonbin2[[2]],covbin2[[2]],nondic2[[2]],covdic2[[2]],nonggm2[[2]],covggm2[[2]]))
  print(c(nonbin3[[1]],covbin3[[1]],nondic3[[1]],covdic3[[1]],nonggm3[[1]],covggm3[[1]]))
  print(c(nonbin3[[2]],covbin3[[2]],nondic3[[2]],covdic3[[2]],nonggm3[[2]],covggm3[[2]]))
}




