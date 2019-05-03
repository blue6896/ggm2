#set.seed(100)
library(glmnet)
library(SparseM)

#set.seed(16)

simu100AIC <- function ( X, selectedNodes,numNodes,numSamples,alphaval,validtype,covadj){
##########
hi=c()
bye=c()
lst=read.table('lambdaset.txt')

 for(i in 1:dim(X)[2]){
#if (sum(X[,i])>1 | sum(X[,i])<(dim(X)[1]-1)){ hi=c(hi,i)}
#if (sum(X[,i])<=1 | sum(X[,i])>=(dim(X)[1]-1)){ bye=c(bye,i)}
}
if (!is.null(bye)){
#X=X[,-bye]
}
#print(dim(X))
#print(colnames(X))
#print(length(colnames(X)))

#print(dim(X))
#print(rowSums(X))
#print(length(microname))

#selectedNodes<-0
numNodes<-dim(X)[2]
numSamples<-dim(X)[1]

  # truegraph: true interaction graph 
  # X: generated samples for true graph 
  # selectedNodes: selected markers 
  # numNodes: number of nodes
  # numSamples: number of Samples 
  
  selectedindex<-1           #<-head(order(colSums(truegraph!=0) ), selectedNodes)
if(selectedNodes>0){
selectedindex<- c(1:selectedNodes)
}else if (selectedNodes==0){
selectedindex=c()
}
mindeviance=c()

  Xold<- X
  X<- X[,(selectedNodes+1):numNodes]
  
  graphL1red<-c()
  devdev<-c()
if(covadj){
  p.fac= c(rep(1,(numNodes-selectedNodes-1)),rep(0,selectedNodes)  )
}
else{
 p.fac= c(rep(1,(numNodes-selectedNodes-1)),rep(1,selectedNodes)  )
}

for(lamlam in 1:99){
print(lamlam)
  graphL1red<-c()
  devdev<-c()

  for (s in 1: (numNodes-selectedNodes) ){
    Vbs = setdiff(c(1:numNodes), selectedindex)-selectedNodes
    Vbs = setdiff(Vbs,s)

    # next four lines are commendted 
    if(length(unique(X[,s]))<=1){
      graphL1red<-cbind(graphL1red,rep(0, (numNodes)))
      next
    }
newX= cbind(X[,Vbs], Xold[,selectedindex])
#print(colSums(newX))
#print(sum(X[,s]))   
      fit <- glmnet( newX,(X[,s]),family="binomial",alpha=alphaval, penalty.factor=p.fac, lambda=lst[lamlam])
#print(length(Vbs))
small.lambda.betas <- fit$beta
#print(small.lambda.betas)
#print(X[,Vbs])
#print(length(small.lambda.betas))

#print(dim(X[,Vbs]))
#print(class(small.lambda.betas))

#print(class(X[,Vbs]))
#k= sum((fit$beta)!=0)

tLL <- fit$nulldev - deviance(fit)
k <- fit$df
n <- fit$nobs

if (validtype=='BIC'){
#print(Vbs)
#print(X[,Vbs])
#print(small.lambda.betas)
devdev=c(devdev, log(n)*k-tLL)
#devdev=c(devdev, -2*(sum(X[,s]* ( newX%*%small.lambda.betas)   -log(1+exp(  newX%*%small.lambda.betas))))   +   log(numSamples)* k)
}
else if (validtype=='AIC'){
#devdev=c(devdev, -2*(sum(X[,s]* (  newX%*%small.lambda.betas)   -log(1+exp(  newX%*%small.lambda.betas))))   +   2*k)
devdev=c(devdev,-tLL+2*k)
#devdev=c(devdev, numSamples * log( sum(( X[,s]- (  newX%*%small.lambda.betas))^2)/numSamples )    +   2*k)
}
else if (validtype=='AICC'){
devdev= c(devdev,-tLL+2*k+2*k*(k+1)/(n-k-1))
#devdev=c(devdev, -2*(sum(X[,s]* ( newX%*%small.lambda.betas)   -log(1+exp(  newX%*%small.lambda.betas))))   +   2*k + 2*k*(k+1)/(numNodes-k-1)  )
}
         
    colnames(small.lambda.betas)<-NULL
    hall<- small.lambda.betas# gt$beta
    hall2<-c()
    hall<-as.matrix(hall) #
   # print('checkin')
    hall2<-c()
    if (s==1){
      hall2<-rbind(0,hall)
      graphL1red<-as.vector(hall2)
  #    print(graphL1red)
      
    }
    else {
      hall2<-c(hall[1:(s-1)],0,hall[s:length(hall)])
      graphL1red<-cbind(graphL1red,(hall2))
    }
    
  }

mindeviance= c(mindeviance,  mean(devdev) ) 
#print(mindeviance)
}



  for (s in 1: (numNodes-selectedNodes) ){
    Vbs = setdiff(c(1:numNodes), selectedindex)-selectedNodes
    Vbs = setdiff(Vbs,s)

    # next four lines are commendted
    if(length(unique(X[,s]))<=1){
      graphL1red<-cbind(graphL1red,rep(0, (numNodes-selectedNodes)))
      next
    }
newX= cbind(X[,Vbs], Xold[,selectedindex])
#      fit <- glmnet( newX,(X[,s]),family="binomial",alpha=alphaval, penalty.factor=p.fac, lambda=0.01*minind)
 fit <- glmnet( newX,(X[,s]),family="binomial",alpha=alphaval, penalty.factor=p.fac, lambda=lst[(which(mindeviance== min(mindeviance))[1])])
small.lambda.betas <- fit$beta

    colnames(small.lambda.betas)<-NULL
    hall<- small.lambda.betas# gt$beta
    hall2<-c()
    hall<-as.matrix(hall) #
   # print('checkin')
    hall2<-c()
    if (s==1){
      hall2<-rbind(0,hall)
      graphL1red<-as.vector(hall2)
  #    print(graphL1red)

    }
    else {
      hall2<-c(hall[1:(s-1)],0,hall[s:length(hall)])
      graphL1red<-cbind(graphL1red,(hall2))
    }

  }
# just add last columns (covariate) by using cbind 
print(dim(graphL1red))
print(numNodes)
print(selectedNodes)
#print(cbind(graphL1red, matrix(1,(numNodes+selectedNodes), selectedNodes)))
print('dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd')

return(as.matrix(graphL1red))
                                                                            #print(graphL1red[1:30,1:30])
#################  return(graphL1red[1:(numNodes-selectedNodes),1:(numNodes-selectedNodes)])
}


###results=graphL1redgenerator3( X, selectedNodes,numNodes,numSamples,1,'AIC')
#print(results)
#print(dim(results))
#print(numNodes)
#print(selectedNodes)
###write.table(results,'results100AIC.txt',col.names=hi[1:dim(X)[2]])
###mainmat= abs(results)
###mainmat= (mainmat>0.01)+0

###print(mainmat)