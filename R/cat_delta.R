cat_delta <- function(x, y = NULL, method = NULL, method_cat="tot_var_dist", mkw_p = 1){
  
  a = NULL
  b = NULL
  blocks = NULL
  .x = NULL
  id = NULL
  
  # Define is_scalar_vector function
  is_scalar_vector <- function(x) {
    is.atomic(x) && length(x) == 1
  }
  
  # if(method=="association_based"){
  #   ab_method=method
  method=method_cat
  # }else{
  #   ab_method=method
  # }
  catdiss = method
  # print(map(x,class))
  # print(map_dbl(x,nlevels))
  x = map_df(x,fct_drop)
  # print(map_dbl(x,nlevels))
  
  
  
  if(is.null(dim(x))){
    Q=nlevels(x)
  }else{
    Q=map_dbl(x,nlevels)
  }
  
  
  n=nrow(x)
  nvar=length(Q)
  
  z_prep = z_preproc(x=x,y=y,Q=Q)
  #print("preprocessed")
  
  Z_names= colnames(z_prep$Z)
  
  Z = z_prep$Z %>% data.matrix()
  
  ZZod = z_prep$ZZod
  zm = z_prep$zm
  
  Z_list = z_prep$Z_list
  # level_pos = z_prep$level_pos
  level_pos = data.table(start=c(1,cumsum(Q)[-length(Q)]+1),stop=cumsum(Q))
  if(!is.null(y)){
    Z_y = z_prep$Z_y
  }else{Z_y=NULL}
  
  
  
  Qs=ncol(Z)
  
  
  if(method %in% c("tot_var_dist", "gifi_chi2",
                   "supervised","supervised_full","matching","eskin",
                   "goodall_3","goodall_4","iof","of","lin","var_entropy","var_mutability")){
    
    full_delta = cat_custom_delta(ZZod=ZZod,Z=Z,Z_y=Z_y,Z_list=Z_list,
                                  zm=zm,Q=Q,nvar=nvar,method=method,Qs=Qs)
    
    
  }else{
    
    
    crs = crossing(a=1:(nvar),b=1:nvar) %>% filter(a!=b)
    blocks_id_a = level_pos[crs$a,]
    blocks_id_b = level_pos[crs$b,] %>% rename(`end_start`=start,`end_stop`=stop)
    block_ids=cbind(blocks_id_a,blocks_id_b)
    
    pull_block <-function(start=1,stop=1,end_start=1,end_stop=1,squared=TRUE){
      if(squared==TRUE){
        return(ZZod[start:stop,end_start:end_stop])
      }else{
       return(Z[,start:stop])
      }
  
    }
    
    pull_conditional_block <- function(start, stop, end_start, end_stop) {
      # Extract the joint probability block between variables j and k
      joint_block <- ZZod[start:stop, end_start:end_stop]
      
      # Convert to conditional probabilities P(X_j | X_k)
      # Each column represents P(X_j = x_j | X_k = x_k)
      conditional_block <- sweep(joint_block, 2, colSums(joint_block), FUN = "/")
      
      # Handle cases where colSums might be zero (avoid division by zero)
      conditional_block[is.nan(conditional_block)] <- 0
      conditional_block[is.infinite(conditional_block)] <- 0
    
      return(conditional_block)
    }
    # distance_blocks = tibble(row_ind = crs$a,col_ind=crs$b,
    #                          blocks = pmap(block_ids, ~pull_block(start=..1,stop=..2,
    #                                                               squared=FALSE)
    #                          )
    # )
    # }else{
    distance_blocks = tibble(row_ind = crs$a,col_ind=crs$b,
                             blocks = pmap(block_ids, ~pull_block(start=..1,stop=..2,
                                                                  end_start=..3,end_stop=..4)
                             )
    )
    # }
    #print("distance_blocks$blocks")
    #print(distance_blocks$blocks)
    
    if(catdiss=="le_and_ho"){
      
      # Uses pull_conditional_block
      distance_blocks = tibble(row_ind = crs$a,col_ind=crs$b,
                               blocks = pmap(block_ids, ~pull_conditional_block(start=..1,stop=..2,
                                                                                end_start=..3,end_stop=..4)
                               ))
      
      distance_blocks = distance_blocks %>%
        mutate(block_dist=map(
          .x=blocks,.f = function(x=.x){
            x[is.na(x)]=0
            # x=R
            # ORIGINAL CODE
            dist_long = crossing(a=1:nrow(x), b=1:nrow(x)) |> filter(a != b) |> # ORIGINAL CODE
              mutate(kl = map2_dbl(.x=a,.y=b,~philentropy::kullback_leibler_distance(P =x[.x,],Q=x[.y,],unit = "log2",testNA=FALSE,epsilon=0.00000001))) # ORIGINAL CODE
            
            # NEW CODE: IF YOU WANT TO INCLUDE R^{j|j}
            # dist_long = crossing(a=1:nrow(x), b=1:nrow(x)) |> # NEW CODE
            #   mutate(kl = map2_dbl(.x=a,.y=b,~philentropy::kullback_leibler_distance(P =x[.x,],Q=x[.y,],unit = "log2",testNA=FALSE,epsilon=0.00000001))) # NEW CODE
            
            phil_dist = matrix(0,nrow(x),nrow(x))
            for(i in 1:nrow(dist_long)){
              phil_dist[dist_long$a[i],dist_long$b[i]]=dist_long$kl[i]
              # phil_dist[dist_long$b[i],dist_long$a[i]]= phil_dist[dist_long$b[i],dist_long$a[i]] + dist_long$kl[i]
            }
            phil_dist = phil_dist + t(phil_dist)
            
            
            if(is_scalar_vector(phil_dist)){
              # print(is_scalar_vector(phil_dist))
              phil_dist=matrix(phil_dist,2,2);
              diag(phil_dist)=0
            }
            
            return(phil_dist)
          }
          
        )
        )
    }else{
      distance_blocks = distance_blocks %>%
        mutate(block_dist=map(
          .x=blocks,.f = function(x=.x){
            x[is.na(x)]=0
            phil_dist = philentropy::distance(x = x,method=catdiss,
                                              mute.message = TRUE,p=mkw_p);
            
            if(is_scalar_vector(phil_dist)){
              # print(is_scalar_vector(phil_dist))
              phil_dist=matrix(phil_dist,2,2);
              diag(phil_dist)=0}
            
            return(phil_dist)
          }
        )
        )
      
    }
    
    # print("distance_blocks$block_dist")
    # print(distance_blocks$block_dist)
    #
    
    ####################################################################
    ####################################################################
    ### THE WEIGHTS SHOULD GO HERE #####################################
    ####################################################################
    ####################################################################
    
    delta_blocks = tibble(id=as.list(1:nvar)) %>%
      mutate(diag_delta = map(.x=id,
                              .f=~Reduce("+", distance_blocks %>%
                                           filter(row_ind==.x) %>%
                                           pull(block_dist))
      )
      )
    
    ####################################################################
    ####################################################################
    
    #print("delta_blocks$diag_delta")
    #print(delta_blocks$diag_delta)
    
    full_delta = bdiag(delta_blocks$diag_delta) %>% as.matrix
    full_delta = full_delta/(nvar-1)
  }
  
  # if(ab_method=="association_based"){
  #   method=ab_method
  # }
  
  out=list()
  out$delta_names = Z_names
  # full_delta[is.na(full_delta)]=0
  out[[method_cat]] = full_delta %>% as.matrix
  # out$delta_blocks = delta_blocks
  out$Z = Z
  
  return(out)
}
