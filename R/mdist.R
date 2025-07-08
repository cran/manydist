mdist <- function(x,validate_x=NULL,response=NULL, distance_cont="manhattan", distance_cat="tot_var_dist",
                  commensurable = FALSE,scaling_cont="none",
                  ncomp=ncol(x), threshold = NULL,preset = "custom"){#,prop_nn=0.1, alpha=.5){
  
  # source("gower_recipe.R")
  # source("ndist.R")
  # source("cdist.R")
  
  .x = NULL
  a <- NULL
  b <- NULL
  gowdist <- NULL
  cat_data  = x %>% dplyr::select(where(is.factor))
  cont_data = x %>% dplyr::select(where(is.numeric))
  if (ncol(cat_data) == 0)
    cat_data = NULL
  
  if (ncol(cont_data) == 0)
    cont_data = NULL
  
  
  # Check if tot_var_dist is specified but only one categorical variable exists
  if (!is.null(cat_data) && distance_cat == "tot_var_dist" && ncol(cat_data) == 1) {
    warning("'tot_var_dist' requires more than one categorical variable. Switching to 'matching' distance.")
    distance_cat <- "matching"
  }
  
  # Check if pc_scores scaling_cont is specified but only one continuous variable exists
  if (!is.null(cont_data) && scaling_cont == "pc_scores" && ncol(cont_data) == 1) {
    warning("With only one variable, PCA produces a single component identical to standardization. Consider using scaling_cont = \"std\" instead.")
    
  }

  if(!is.null(validate_x)){
    cat_data_val  = validate_x %>% dplyr::select(where(is.factor))
    cont_data_val = validate_x %>% dplyr::select(where(is.numeric))
    if (ncol(cat_data_val) == 0)
      cat_data_val = NULL
    
    if (ncol(cont_data_val) == 0)
      cont_data_val = NULL
  }
  
  
  #### ACTUAL MIXED DATA
  if(!is.null(cont_data) & !is.null(cat_data)){
    if(preset == "gower"){
      if(is.null(validate_x)){
        if (commensurable == FALSE){
          distance_mat <-  as.matrix(daisy(x, metric = "gower"))
        } else {
          
          gowerlist = x %>% map(~daisy(as_tibble(.x),metric="gower") %>% as.matrix())
          gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT : CHECK IF DISTA with gower works
        gow_prep <- gower_recipe(data=x) |> prep(training=x)
        x = gow_prep |> bake(new_data = NULL)
        validate_x = gow_prep |> bake(new_data = validate_x)
        
        if (commensurable == FALSE){
          # distance_mat <-  as.matrix(dist(x, method = "manhattan"))[1:5,1:5]
          
          distance_mat <-  Rfast::dista(xnew = validate_x, 
                                        x = x,
                                        type = "manhattan") |> as.matrix()
          
          
        }else{
          gowerlist = map2(.x=x,.y=validate_x,
                           ~Rfast::dista(xnew = .y,x = .x,
                                         type = "manhattan") |> as.matrix()
          )
          
          gowerlist = tibble(gowdist = gowerlist) |>  
            mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }
      
    }else if(preset == "catdissim"){
      distance_cont = "manhattan"
      distance_cat = "matching"
      commensurable = TRUE
      # scaling_cont="none"
      scaling_cont="std"
      
      if(is.null(validate_x)){
        cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x = cont_data,validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,validate_x=cat_data_val,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }
    }else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      distance_cat = "tot_var_dist"
      commensurable = TRUE
      cont_scaling_cont = "pc_scores"
      # cont_scaling_cont="none"
      # cat_scaling_cont="none"
     # scaling_cont="std"
      
      # Additional check for this preset
    #  if (!is.null(cat_data) && ncol(cat_data) == 1) {
    #    warning("'unbiased_dependent' preset with 'tot_var_dist' requires more than one categorical variable. Switching to 'matching' distance.")
    #    distance_cat <- "matching"
    #  }
    
      if(is.null(validate_x)){
        cont_dist_mat = ndist(x = cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x = cont_data, validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,validate_x= cat_data_val,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }
    } else if(preset == "euclidean_onehot"){
      
      distance_cont = "euclidean"
      commensurable = FALSE
      scaling_cont="std"
      
      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal(),one_hot = TRUE)
      
      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)
      
      if(is.null(validate_x)){
        cont_dist_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = ndist(x=cat_data_dummy, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x=cont_data,validate_x=cont_data_val, 
                              method = distance_cont,
                              commensurable = commensurable,
                              scaling=scaling_cont)  |>  as.matrix()
        
        cat_data_val_dummy = dummy_recipe |>
          prep(training = cat_data) |>
          bake(new_data=cat_data_val)
        
        cat_dist_mat = ndist(x=cat_data_dummy,
                             validate_x=cat_data_val_dummy, 
                             method = distance_cont,
                             commensurable = commensurable,
                             scaling=scaling_cont)  |>  as.matrix()
        
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        
      }
      }else if(preset=="entropy_based"){
        # if(weight_cont != "commensurable"){weight_cont = weight_cont}
        #  if(weight_cat != "commensurable"){weight_cat = weight_cat}
        # 
        # n_cont=ncol(cont_data)
        # x=cbind(cont_data,cat_data)
        # distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
        
        #   }else if(preset=="indicator_based"){
        #      if(distance_cont!="manhattan"){distance_cont = distance_cont}
        #      if(distance_cat!="tot_var_dist"){distance_cat = distance_cat}
        #  if(weight_cont != "commensurable"){weight_cont = weight_cont}
        #  if(weight_cat != "commensurable"){weight_cat = weight_cat}
        #     if(cont_scaling_cont!="none"){cont_scaling_cont=cont_scaling_cont}
        
        #      cat_dist_mat <- indicator_based(x,commensurable = commensurable, scaling_cont=cat_scaling_cont, weights=1)$distance_mat
        #     cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling_cont=cont_scaling_cont)  |>  as.matrix()
        # #  print(cont_dist_mat[1:5,1:5])
        #    if ((distance_cont == "euclidean") | (distance_cat == "euclidean"))
        #      distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        #    else
        #      distance_mat = cat_dist_mat + cont_dist_mat
        
        
      }else if(preset=="custom"){
        #    if(weight_cont != "commensurable"){weight_cont = weight_cont}
        #    if(weight_cat != "commensurable"){weight_cat = weight_cat}
        if(is.null(validate_x)){
          cont_dist_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
          cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
          distance_mat = cat_dist_mat + cont_dist_mat
          if ((distance_cont == "euclidean") & (distance_cat=="HLeucl")){
            distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
          }
        }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT    
          cont_dist_mat = ndist(x=cont_data,validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
          cat_dist_mat = cdist(x=cat_data,validate_x=cat_data_val, response=response,method=distance_cat,commensurable = commensurable)$distance_mat
          distance_mat = cat_dist_mat + cont_dist_mat
          if ((distance_cont == "euclidean") & (distance_cat=="HLeucl"))
            distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        }
        
        #### 
      }
      
    } else if(is.null(cont_data) & !is.null(cat_data)){
      
      if(preset == "gower"){
        
        distance_mat=ncol(cat_data)*daisy(cat_data,metric = "gower") %>% as.matrix()
        
        
      }else if(preset == "catdissim"){
        
        distance_cat = "matching"
        commensurable = TRUE
        cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable=commensurable)$distance_mat
        distance_mat = cat_dist_mat
        
      }else if(preset == "unbiased_dependent"){
        
        distance_cat = "tot_var_dist"
        weight_cat = "commensurable"
        cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable=commensurable)$distance_mat
        distance_mat = cat_dist_mat
        
      }else if(preset == "euclidean_onehot"){
        
        dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal())
        cat_data_dummy = dummy_recipe |>
          prep(training = cat_data) |>
          bake(new_data=NULL)
        
        cat_dist_mat = ndist(cat_data_dummy, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        
      }else if(preset=="entropy_based"){
        
        n_cont=0
        x = cat_data
        # distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
        
      }
      else if(preset=="custom"){
        
        distance_mat = cdist(cat_data, method = distance_cat,commensurable=commensurable)$distance_mat  |>  as.matrix()
        
      }
      # Continuous only
      
      
    }else if(!is.null(cont_data) & is.null(cat_data)){
      
      if(preset == "gower"){
        
        distance_mat=ncol(cont_data)*daisy(cont_data,metric= "gower") %>% as.matrix()
        
      }else if(preset == "catdissim"){
        distance_cont = "manhattan"
        commensurable=TRUE
        # scaling_cont="none"
        scaling_cont="std"
        cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
        distance_mat = cont_dist_mat
        
      }else if(preset == "unbiased_dependent"){
        distance_cont = "manhattan"
        commensurable=TRUE
        scaling_cont="std"
        cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
        distance_mat = cont_dist_mat
      }else if(preset == "euclidean_onehot"){
        distance_cont = "euclidean"
        commensurable=FALSE
        scaling_cont="std"
        distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
      }else if(preset=="entropy_based"){
        
        n_cont=ncol(cont_data)
        x=cont_data
        #    distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
      }else if(preset=="custom"){
        distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
        
      }
      
    }
    
    return(distance_mat)
  }
  
  
