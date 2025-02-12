cat_custom_delta<-function(ZZod,Z,Z_y,Z_list,zm,Q,nvar,method,Qs){

  .x = NULL
  a = NULL
  b = NULL
  blocks = NULL
  id = NULL
  # Qs=ncol(Z)
  n=nrow(Z)
  Qlist=split(Q,1:length(Q))
  blocks = bdiag(map(.x=Qlist,~matrix(1,nrow=.x,ncol=.x)))

  Qlisty=split(c(Q,ncol(Z_y)),1:length(c(Q,ncol(Z_y))))
  blocksy = bdiag(map(.x=Qlisty,~matrix(1,nrow=.x,ncol=.x)))


  create_delta <- function(Z,nvar,p,chi2=FALSE){
    if(chi2==FALSE){
      if(nvar !=1)
        Dsel<-as.matrix(dist(Z,method="minkowski", p=p))/((3-p)*(nvar-1)) # divide by 2 if ln=1
      else
        ### CHECK IF THIS IS OK 30-Oct-2024
        Dsel<-as.matrix(dist(Z,method="minkowski", p=p))/((3-p)*(nvar)) # divide by 2 if ln=1
    }else{
      if(nvar !=1)
        Dsel<-as.matrix(dist(Z,method="minkowski", p=p)^2)/((3-p)*(nvar-1)) # divide by 2 if ln=1
      else
        ### CHECK IF THIS IS OK 30-Oct-2024
        Dsel<-as.matrix(dist(Z,method="minkowski", p=p))/((3-p)*(nvar)) # divide by 2 if ln=1
    }
    return(Dsel)
  }



  if(method=="tot_var_dist"){
    p=1

    full_delta = blocks * create_delta(ZZod, nvar,p,chi2=FALSE)


  }
  # if(method=="euclid_prof"){
  #   #print("in euclid_prof")
  #   p=2
  #   full_delta = blocks * create_delta(ZZod, nvar,p,chi2=FALSE)
  #
  # }
  if(method=="gifi_chi2" ){
    #print("in chi2")
    #ZZod2
    p=2
   # full_delta = blocks *  create_delta(t(t(ZZod)/sqrt(zm)), nvar,p,chi2=TRUE)
    full_delta = blocks *  create_delta(t(t(ZZod)/sqrt(zm/n)), nvar,p,chi2=TRUE) # New: Corrected probability zm/n. Takes average over nvar-1
    
    
  }
  # if(method=="w_chi2"){
  #   #print("in w_chi2")
  #   #ZZod4
  #   p=2
  #   delta_input=t(t(ZZod)/sqrt(zm)) * zm^(.5)
  #   full_delta = blocks * create_delta(delta_input,nvar,p,chi2=FALSE)
  #
  # }
  # if(method=="std_res"){
  #   #print("in std_res")
  #   #ZZod3
  #   p=2
  #   full_delta = blocks * create_delta(t(t(ZZod)/sqrt(zm)) / zm^(.5), nvar,p,chi2=FALSE)
  #
  # }
  if(method=="supervised"){
    #  print("in suvervised")
    #ZZyod
    full_delta = blocks * create_delta((t(Z)%*%Z_y)/zm,3/2,1,chi2=FALSE) # so that ((3-1)*(nv-1))= 1

  }

  if(method=="supervised_full"){
    # print("in suvervised")
    #ZZyod
    ZZy <- cbind(Z,Z_y)
    zmy <- c(zm,colSums(t(Z_y)%*%Z_y))
    BBy <- (t(ZZy)%*%ZZy)/zmy
    full_delta = blocksy * create_delta(BBy,3/2,1,chi2=FALSE) # so that ((3-1)*(nv-1))= 1
    full_delta = full_delta[1:ncol(Z),1:ncol(Z)]
  }

  if(method=="matching"){
    #print("in matching")
    full_delta = blocks ### /length(Q) WARNING 3 OCT SH1TSHOW

    diag(full_delta) = 0

  }
  if(method=="eskin"){

    # print("in eskin")

    full_delta = bdiag(map(.x=Z_list,.f=function(x=.x){
      Qi=ncol(x)
      sk=Qi^2/(Qi^2+2)
      return(
        #matrix((1/sk -1),nrow=Qi,ncol=Qi)/nvar
        matrix((1/sk -1),nrow=Qi,ncol=Qi)
      )
    })
    )

    diag(full_delta)=0


  }
  if(method=="goodall_3"){
    #print("in good3")
    full_delta = (blocks* (rep(1,Qs,Qs)-
                             diag(1 - (zm *  (zm- 1)) / (n*(n-1)),
                                  nrow = Qs,ncol = Qs)))/nvar

  }
  if(method=="goodall_4"){
    #print("in good4")
    full_delta = (blocks* (rep(1,Qs,Qs)-
                             diag(((zm *  (zm- 1)) / (n*(n-1))),
                                  nrow = Qs,ncol = Qs)))/nvar

  }
  if(method=="iof"){
    #print("in iof")
    full_delta = 1/(1+(log(zm)%*%t(log(zm))))
    diag(full_delta) = 1
    full_delta = (blocks*(1/full_delta -1))/nvar

  }
  if(method=="of"){
    #print("in of")
    full_delta = 1/(1+(log(n/zm)%*%t(log(n/zm))))
    diag(full_delta) = 1
    full_delta = (blocks*(1/full_delta -1))/nvar
    #  if (sum(is.nan(full_delta))>0)
    #    print('NaNs!!')
  }

  if(method=="lin"){
    prop = as.matrix(zm/n)

    pv<-prop
    pr<-pv %*% rep(1,Qs)
    pc<-rep(1,Qs) %*% t(pv)
    pp<-pr+pc
    pplog<-log(pr)+log(pc)
    diag(pp)<-prop
  #  linsim2<- (pplog - 2*log(pp))/2*log(pp)
    linsim2<- (pplog - 2*log(pp))/(2*log(pp)) # NEW CODE, ADDED BRACKETS FOR DIVISION AND ADDED EPS 27-12-2024
    # print(blocks)
    full_delta<-as.matrix(blocks*linsim2)
    #full_delta<-as.matrix(blocks*linsim2)/nvar
    # For all 2x2 blocks corresponding to q_j=2, set full_delta equal to 0
    start_index <- 1
    for (j in seq_along(Q)) {
      if (Q[j] == 2) { # Check if q_j = 2
        end_index <- start_index + 1
        # full_delta[start_index:end_index, start_index:end_index] <- 0 # Set 2x2 block to 0
        # # Create simple matching dissimilarity matrix for this block
        simple_matching_block <- matrix(1, nrow = 2, ncol = 2) - diag(2)
        full_delta[start_index:end_index, start_index:end_index] <- simple_matching_block
      }
      start_index <- start_index + Q[j] # Move to the next block
    }
  }

  if(method=="var_entropy"){
    prop = as.matrix(zm/n)
    full_delta = matrix(0,Qs,Qs)
    pos=0
    for(i in 1:nvar){
      sk= -1/log(Q[i])
      plogp<-prop[(pos+1):(pos+Q[i])] %*% log(prop[(pos+1):(pos+Q[i])])
      full_delta[(pos+1):(pos+Q[i]),(pos+1):(pos+Q[i])]<-1 - sk*diag(rep(plogp,Q[i]))
      pos<-pos+Q[i]
    }
    full_delta=full_delta/nvar
  }

  if(method=="var_mutability"){
    prop = as.matrix(zm/n)
    full_delta = matrix(0,Qs,Qs)
    pos=0
    for(i in 1:nvar){
      sk2<-Q[i]/(Q[i]-1)
      pp<-1-prop[(pos+1):(pos+Q[i])] %*% prop[(pos+1):(pos+Q[i])]
      full_delta[(pos+1):(pos+Q[i]),(pos+1):(pos+Q[i])]<- 1 - sk2*diag(rep(pp,Q[i]))
      pos<-pos+Q[i]
    }
    full_delta=full_delta/nvar
  }
  return(full_delta)
}
