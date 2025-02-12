cdist<-function(x,response=NULL,validate_x=NULL,method="tot_var_dist", commensurable = FALSE, weights=1){
  
  full_delta <- NULL
  level_stop <- NULL
  level_start <- NULL
  # response <- NULL
  delta_tmp <- NULL
  all_of <- NULL
  map2 <- NULL
  .x <- NULL
  a <- NULL
  b <- NULL
  blocks <- NULL
  id <- NULL
  delta_ms <- NULL
  
  
  y=response
  
  x=as_tibble(x) |> mutate(across(everything(),fct_drop))
  
  if ((length(method) > 1 ) & (any(method %in% c("none","st_dev","HL","cat_dis","HLeucl","mca")))) {
    warning("Specifying a different method for each variable is not currently possible for 'st_dev','HL','cat_dis','HLeucl','mca'.")
    return()
  }
  
  if(length(method)==1){
    
    if(method %in% c("none", "st_dev","HL","cat_dis","HLeucl","mca")) {
      
      outindbased <- indicator_based(x=x,validate_x=validate_x,commensurable=commensurable, scaling=method, weights=1)
      distance_mat <- outindbased$distance_mat
      delta <- outindbased$delta
      delta_ms <- outindbased$delta_ms
      delta_names <- outindbased$delta_names
      # colnames(delta) = delta_names
      # rownames(delta) = delta_names
        
    } else {
      
      out_delta = cat_delta(x=x,y=y,method_cat=method)
      delta = out_delta[[method]] |> data.matrix()
      delta_names = out_delta$delta_names
    
      colnames(delta) = delta_names
      rownames(delta) = delta_names
      
      
      ########################################################################
      ########################################################################
      ########################################################################
      # delta[is.na(delta)]=0
      ########################################################################
      ########################################################################
      ########################################################################
      
      Z = out_delta$Z |> data.matrix()
      
      if(is.null(dim(x))){
        Q=nlevels(x)
      }else{
        Q=map_dbl(x,nlevels)
      }
      
      if(commensurable==TRUE){
        
        distance_mat <- commensurability_for_cat(x=x,y=y,validate_x=validate_x,delta=delta)
        
      }else{ 
        # print("non-commensurable")
        #not commensurable
        if(is.null(dim(weights))){
          if(length(weights) == 1){
            # print("length weights =1")
            if(is.null(validate_x)){
              distance_mat = Z  %*% delta %*%  t(Z)
            }else{
              # tib_x=
              val_Z_tib =  recipe(~.,data = x)|>
                step_dummy(all_nominal_predictors(),one_hot = TRUE) |>
                prep(training = x) |>
                bake(new_data=validate_x) 
              
              val_Z = val_Z_tib |> as.matrix()
              
              
              distance_mat = val_Z  %*% delta %*%  t(Z)
            }
            
          }else{
            weightsexp =NULL
            for (i in 1:ncol(x)) {
              weightsexp = c(weightsexp,rep(weights[i],Q[i]))
            }
            W=diag(weightsexp,nrow=length(weightsexp),ncol=length(weightsexp))
            
            if(is.null(validate_x)){
              distance_mat = Z %*% W %*% delta %*% W %*% t(Z)
            }else{
              val_Z = validate_x |> recipe(~.)|>
                step_dummy(all_predictors(),one_hot = TRUE) |>
                prep(training = as_tibble(x)) |>
                bake(new_data=validate_x) |> as.matrix()
              
              distance_mat = val_Z  %*% W %*% delta %*% W %*% t(Z)
            }  
          }
          
        }
      }
    }
    
  }else{ #differnt method for each variable
    
    
    method_vec = method
    
    if(is.null(dim(x))){
      Q=nlevels(x)
    }else{
      Q=map_dbl(x,nlevels)
    }
    nvar=length(Q)
    level_pos = data.table(start=c(1,cumsum(Q)[-length(Q)]+1),stop=cumsum(Q))
    
    delta_structure = tibble(method = method_vec) |> mutate(
      x=map(.x=method_vec,~as_tibble(x)),
      delta_tmp=map2(.x = x,.y = method,.f=~cat_delta(x=.x,method=.y)),
      delta_names = map(.x=delta_tmp,.f=~.x[[1]]),
      full_delta = map(.x=delta_tmp,.f=~.x[[2]]),
      level_start = level_pos$start,
      level_stop = level_pos$stop,
      delta_block=pmap(.l=list(a = full_delta,b=level_start,c=level_stop),
                       .f=function(a,b,c)
                       {return(a[b:c,b:c])}
      )
    )
    
    Z = delta_structure$delta_tmp[[1]]$Z |> as.matrix()
    
    delta = bdiag(delta_structure$delta_block) |> as.matrix()
    delta_names = delta_structure$delta_names |> unlist()
    
    if(commensurable==TRUE){
      #   weight_vec=commensurable_weight(cat_data=cat_data,delta_cat=delta)
      #    W=diag(weight_vec,nrow=length(weight_vec),ncol=length(weight_vec))
      #    distance_mat = Z %*% W %*% delta %*% W %*% t(Z)
      distance_mat <- commensurability_for_cat(x=x,y=y,validate_x=validate_x,delta=delta)
      
    }else{
      if(is.null(dim(weights))){
        if(length(weights) == 1){
          if(is.null(validate_x)){
            distance_mat = Z  %*% delta %*%  t(Z)
          }else{
            val_Z = validate_x |> recipe(~.)|>
              step_dummy(all_predictors(),one_hot = TRUE) |>
              prep(training = as_tibble(x)) |>
              bake(new_data=validate_x) |> as.matrix()
            
            distance_mat = val_Z  %*% delta %*%  t(Z)
          }
        }else{
          
          weightsexp =NULL
          for (i in 1:ncol(x)) {
            weightsexp = c(weightsexp,rep(weights[i],Q[i]))
          }
          W=diag(weightsexp,nrow=length(weightsexp),ncol=length(weightsexp))
          if(is.null(validate_x)){
            distance_mat = Z %*% W %*% delta %*% W %*% t(Z)
          }else{
            val_Z = validate_x |> recipe(~.)|>
              step_dummy(all_predictors(),one_hot = TRUE) |>
              prep(training = as_tibble(x)) |>
              bake(new_data=validate_x) |> as.matrix()
            
            distance_mat = val_Z %*% W %*% delta %*% W %*% t(Z)
          }
          
        }
      }else{ # weights is a matrix
        weights = diag(weights)
        weightsexp =NULL
        for (i in 1:ncol(x)) {
          weightsexp = c(weightsexp,rep(weights[i],Q[i]))
        }
        W=diag(weightsexp)
        if(is.null(validate_x)){
          distance_mat = Z %*% W %*% delta %*% W %*% t(Z)
        }else{
          val_Z = validate_x |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(x)) |>
            bake(new_data=validate_x) |> as.matrix()
          
          distance_mat = val_Z %*% W %*% delta %*% W %*% t(Z)
        }
        
      }
      
    }
  }
  # print(dim(distance_mat))
  # if(is.null(out_delta)){
  #    return(distance_mat)
  #  }else{
  out_catdist = list()
  out_catdist$distance_mat = distance_mat
  out_catdist$delta = delta
  out_catdist$delta_ms = delta_ms
  out_catdist$delta_names = delta_names
  return(out_catdist)
  # }
}


