indicator_based<-function(x,validate_x,commensurable=FALSE, scaling="none", weights=1){
  factor_name = NULL
  delta = NULL
  by_var_dist = NULL
  mean_by_var_dist = NULL
  .x = NULL
  eta = NULL
  scaled_Zs = NULL
  comm_dist = NULL
  scaled_val_Zs = NULL
  
  
  #  if (weight_sys != "commensurable") {
  delta_output = cat_delta(x=x,method_cat="matching")
  Z = delta_output$Z
  
  
  
  
  delta_m = delta_output$matching
  delta_ms <- NULL
  delta_names = delta_output$delta_names
  p <- ncol(x)
  
  if(scaling=="none"){
    delta_m = 2 * delta_m#sqrt(2) * delta_m
    #  cat_dist_mat = Z %*%  delta_m  %*% t(Z)
    
    if(!commensurable) {
      
      if(is.null(validate_x)){
        cat_dist_mat = Z %*%  delta_m  %*% t(Z)
      }else{
        
        validate_Z = validate_x |> recipe(~.)|>
          step_dummy(all_predictors(),one_hot = TRUE) |>
          prep(training = as_tibble(x)) |>
          bake(new_data=validate_x) |> as.matrix()
        
        
        
        cat_dist_mat = validate_Z %*% delta_m %*% t(Z)
      }
    } else {
      #### Commensurability
      prep_list = x |> map(
        ~as_tibble(.x) |> recipe(~.)|>
          step_dummy(all_predictors(),one_hot = TRUE) |>
          prep(training = as_tibble(.x)) 
      )
      
      Z_list =  prep_list |> 
        map(~.x |> bake(new_data=NULL))
      
      
      
      Q=map_dbl(x,nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      
      #delta_out = cat_delta(cats,method = method)
      
      if(is.null(validate_x)){
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_m[levels_identifier==.x,
                                      levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
      }else{
        
        
        validate_Z_list =  map2(.x=validate_x, .y=prep_list,
                                ~ .y |> bake(new_data=as_tibble(.x)))
        
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(
            delta = map(.x=factor_name,
                        ~delta_m[levels_identifier==.x,
                                 levels_identifier==.x] |> as.matrix()
            ),
            Zs=Z_list,
            val_Zs=validate_Z_list,
            by_var_dist = pmap(.l=list(a=Z_list,b=delta,c=validate_Z_list),
                               .f=function(a,b,c){
                                 return(as.matrix(c) %*% b %*% t(as.matrix(a)))}),
            mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
            comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                             ~.x /.y)
          )
        # commensurable_dist_structure |> print()
      }
      cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
      
    }
    
  }else if(scaling=="st_dev")
  {
    
    # Create the qs_vec
    qs_vec = x |>
      as_tibble() |>
      map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) |>
      unlist()
    
    # Create the qs_diag
    #  qs_diag = diag(1 / unlist(qs_vec))
    qs_diag = diag(1 / sqrt(unlist(qs_vec))) # ADDED SQRT 27-12-2024
    # Dummy encode and center categorical features
    prep_c_Z = x |>
      as_tibble() |>
      recipe(~.) |>
      step_dummy(all_predictors(), one_hot = TRUE) |>
      step_center(all_predictors())  |> 
      prep(training = x) 
    
    c_Z= prep_c_Z |> 
      bake(new_data = NULL) |> 
      as.matrix()
    
    
    
    # Calculate the covariance matrix S
    S = (1 / nrow(c_Z)) * t(c_Z) %*% c_Z  # EQ 17
    
    # Calculate the inverse square root of the diagonal of S
    inv_sq_S_d = diag(1 / sqrt(diag(S)))
    
    # Dummy encode without centering
    prep_Z = x |>
      recipe(~.) |>
      step_dummy(all_nominal(), one_hot = TRUE) |>
      prep(training = x)
    
    Z = prep_Z |>
      bake(new_data = NULL) |>
      as.matrix()
    
    
    # Convert Z to numeric
    Z = apply(Z, 2, as.numeric)
    
    #delta_ms <- ncol(x)*inv_sq_S_d %*% qs_diag %*% delta_m %*% inv_sq_S_d
    delta_ms <- inv_sq_S_d %*% qs_diag %*% delta_m %*% inv_sq_S_d # NEW
    # Compute the category dissimilarity scaled distance matrix
    #   cat_dist_mat = Z %*% delta_sd %*% t(Z)
    #      print(cat_dist_mat[1:4,1:5])
    if(!commensurable) {
      # Compute the category dissimilarity scaled distance matrix
      if(is.null(validate_x)){
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      }else{
        val_Z = prep_Z |> 
          bake(new_data=validate_x) |> as.matrix()
        
        val_Z = apply(val_Z, 2, as.numeric)
        
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
      #    print(cat_dist_mat[1:4,1:5])
    } else {
      #### Commensurability
      Z_prep_list = x |> map(
        ~as_tibble(.x) |> recipe(~.)|>
          step_dummy(all_predictors(),one_hot = TRUE) |>
          prep(training = as_tibble(.x))
      )
      
      Z_list = Z_prep_list |> 
        map(~.x |> bake(new_data=NULL))
      
      Q=map_dbl(x,nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      #delta_out = cat_delta(cats,method = method)
      
      if(is.null(validate_x)){
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
        
        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
        
      }else{
        
        val_Z_list  = map2(.x = Z_prep_list,.y = validate_x,
                           ~.x |>
                             bake(new_data=as_tibble(.y)))
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          val_Zs=val_Z_list,
          by_var_dist = pmap(.l=list(a=Z_list,b=delta,c=val_Z_list),
                             function(a,b,c) as.matrix(c) %*% b %*% t(as.matrix(a))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
        
        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
        
      }
    }
    
  }else if(scaling=="HL"){
    
    eta_vec = x |>as_tibble() |>
      map(function(x=.x){
        as.vector(rep(fpc::distancefactor(
          cat=nlevels(x),catsizes=table(x)
        ),nlevels(x))
        )
      }
      )
    
    eta_diag = diag(unlist(eta_vec))
    delta_ms <- 2*delta_m %*% eta_diag
    
    
    if(!commensurable) {
      if(is.null(validate_x)){
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      }else{
        
        prep_Z = x |>
          recipe(~.) |>
          step_dummy(all_nominal(), one_hot = TRUE) |>
          prep(training = x)
        
        val_Z = prep_Z |> 
          bake(new_data=validate_x) |> as.matrix()
        
        val_Z = apply(val_Z, 2, as.numeric)
        
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
    } else {
      #### Commensurability
      if(is.null(validate_x)){
        
        prep_Z = x |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)))
        
        #Fix:2 July 2025
        Z_list = prep_Z |> 
          map(~.x |> bake(new_data = NULL))
        
        Q=map_dbl(x,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
        
      }else{
        
        prep_Z = x |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)))
        
        
        Z_list = prep_Z  |> 
          map(~.x |> bake(new_data=NULL)
          )
        
        
        val_Z_list  = map2(.x = prep_Z,.y = validate_x,
                           ~.x |>
                             bake(new_data=as_tibble(.y)))
        
        Q=map_dbl(x,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          val_Zs=val_Z_list,
          by_var_dist = pmap(.l=list(a=Z_list,b=delta,c=val_Z_list),
                             function(a,b,c) as.matrix(c) %*% b %*% t(as.matrix(a))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
      }
      cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
    }
    
  }else if(scaling=="cat_dis"){
    
    # Create the qs_vec
    qs_vec = x |>
      as_tibble() |>
      map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) |>
      unlist()
    
    # Create the qs_diag
    qs_diag = diag(1 / unlist(qs_vec))
    
    # Dummy encode and center categorical features
    c_Z = x |>
      as_tibble() |>
      recipe(~.) |>
      step_dummy(all_predictors(), one_hot = TRUE) |>
      step_center(all_predictors()) |>
      prep(training = x) |>
      bake(new_data = NULL) |>
      as.matrix()
    
    # Calculate the covariance matrix S
    S = (1 / nrow(c_Z)) * t(c_Z) %*% c_Z  # EQ 17
    
    # Calculate the inverse square root of the diagonal of S
    inv_sq_S_d = diag(1 / sqrt(diag(S)))
    
    # Dummy encode without centering
    prep_Z = x  |> 
      recipe(~.) |> 
      step_dummy(all_nominal(), one_hot = TRUE) |>
      prep(training = x)
    
    Z = prep_Z   |> 
      bake(new_data = NULL) |> 
      as.matrix()
    
    # Convert Z to numeric
    Z = apply(Z, 2, as.numeric)
    
    Zs = Z %*% inv_sq_S_d %*% diag(1/sqrt(unlist(qs_vec))) # EQ 18
    
    #delta_ms = sqrt(1/qs_vec)*(inv_sq_S_d %*% delta_m %*% inv_sq_S_d)
    delta_ms = (qs_diag %*% inv_sq_S_d %*% delta_m %*% inv_sq_S_d) # Adjusted 27-12-2024
    # Compute the category dissimilarity scaled distance matrix
    # cat_dist_mat = Z %*% inv_sq_S_d %*% delta_m %*% inv_sq_S_d %*% t(Z)
    #      print(cat_dist_mat[1:4,1:4])
    if(!commensurable) {
      if(is.null(validate_x)){
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      }else{
        
        prep_Z = x |>
          recipe(~.) |>
          step_dummy(all_nominal(), one_hot = TRUE) |>
          prep(training = x)
        
        val_Z = prep_Z |> 
          bake(new_data=validate_x) |> as.matrix()
        
        val_Z = apply(val_Z, 2, as.numeric)
        
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
    } else {
      #### Commensurability
      if(is.null(validate_x)){
        
        prep_Z = x |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)))
        
        #Fix:2 July 2025
        Z_list = prep_Z  |> 
          map(~.x |> bake(new_data=NULL)
          )
        
        Q=map_dbl(x,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
        
      }else{
        
        prep_Z = x |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)))
        
        Z_list = prep_Z  |> 
          map(~.x |> bake(new_data=NULL)
          )
        
        val_Z_list  = map2(.x = prep_Z,.y = validate_x,
                           ~.x |>
                             bake(new_data=as_tibble(.y)))
        
        Q=map_dbl(x,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          val_Zs=val_Z_list,
          by_var_dist = pmap(.l=list(a=Z_list,b=delta,c=val_Z_list),
                             function(a,b,c) as.matrix(c) %*% b %*% t(as.matrix(a))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
      }
      cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
      
    }
    # }
  }else if(scaling=="HLeucl"){
    
    delta_ms <- delta_m
    
    
    Z <- x |>
      recipe(~.) |>
      step_select(all_nominal()) |>
      step_dummy(all_nominal(), one_hot = TRUE) |>
      prep() |>
      bake(new_data = NULL)
    
    ### PREVIOUS eta_vec 7-Oct-24
    #  eta_vec = x |>as_tibble() |>
    #    map(function(x=.x){
    #      as.vector(rep(distancefactor(
    #        cat=nlevels(x),catsizes = table(x)
    #      ),nlevels(x))
    #      )
    #    }
    #    )
    
    eta_vec = x |>as_tibble() |>
      map(function(x=.x){
        as.vector(rep(fpc::distancefactor(
          cat=nlevels(x),catsizes = table(x)
        ),nlevels(x))
        )
      }
      )
    # Apply Hennig and Liao scaling to categorical data
    Zs <- (as.matrix(Z) %*% diag(unlist(eta_vec)))
    
    ### PREVIOUS eta_vec 7-Oct-24
    # Zs <- (as.matrix(Z) %*% diag(unlist(eta_vec))/2)
    
    if(!commensurable) {
      #   if(weight_cat != "commensurable") {
      if(is.null(validate_x)){
         cat_dist_mat = daisy(Zs,metric = "euclidean",warnType=FALSE) |> as.matrix()
      }else{
        prep_Z = x |>
          recipe(~.) |>
          step_select(all_nominal()) |>
          step_dummy(all_nominal(), one_hot = TRUE) |>
          prep(training = x)
        
        val_Z = prep_Z |> 
          bake(new_data=validate_x) |> as.matrix()
        
        val_Zs=(as.matrix(val_Z) %*% diag(unlist(eta_vec)))
        
        
        cat_dist_mat = Rfast::dista(xnew=val_Zs,Zs,type = "euclidean") |> as.matrix()
      }
      #    }
    } else {
      ####### NEW CODE
      if(is.null(validate_x)){
        # Dummy encode categorical variables
        Z_list <- x |>
          map(~ as_tibble(.x) |>
                recipe(~.) |>
                step_dummy(all_nominal(), one_hot = TRUE) |>
                prep() |>
                bake(new_data = NULL)
          )
        
        # Apply Hennig and Liao scaling to each categorical variable
        eta_vec <- x |>
          as_tibble() |>
          map(function(x) {
            as.vector(rep(fpc::distancefactor(
              cat = nlevels(x),
              catsizes = table(x)
            ), nlevels(x)))
          })
        
        # Create a distance matrix for each categorical variable
        commensurable_dist_structure <- tibble(factor_name = names(x)) |>
          mutate(
            Zs = Z_list,
            eta = eta_vec,
            #   scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y) / 2), # Hennig-Liao scaling
            scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y)), # Hennig-Liao scaling
            by_var_dist = map(scaled_Zs, ~ daisy(.x, metric = "euclidean",warnType=FALSE) |> as.matrix()), # Compute distance matrix
            mean_by_var_dist = map_dbl(by_var_dist, ~ mean(.x)), # Compute mean distance
            comm_dist = map2(by_var_dist, mean_by_var_dist, ~ .x / .y) # Normalize for commensurability
          )
      }else{
        
        prep_Z <- x |>
          map(~ as_tibble(.x) |>
                recipe(~.) |>
                step_dummy(all_nominal(), one_hot = TRUE) |>
                prep(training = .x)
          )
        
        Z_list = prep_Z |> 
          map(~.x |> bake(new_data = NULL))
        
        val_Z_list = map2(.x= prep_Z,.y=validate_x,~.x |> bake(new_data = as_tibble(.y)))
        
        # Apply Hennig and Liao scaling to each categorical variable
        eta_vec <- x |>
          as_tibble() |>
          map(function(x) {
            as.vector(rep(fpc::distancefactor(
              cat = nlevels(x),
              catsizes = table(x)
            ), nlevels(x)))
          })
        
        # Create a distance matrix for each categorical variable
        commensurable_dist_structure <- tibble(factor_name = names(x)) |>
          mutate(
            Zs = Z_list,
            val_Zs = val_Z_list,
            eta = eta_vec,
            #   scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y) / 2), # Hennig-Liao scaling
            scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y)), # Hennig-Liao scaling
            scaled_val_Zs = map2(val_Zs, eta, ~ as.matrix(.x) %*% diag(.y)), # Hennig-Liao scaling
            by_var_dist = map2(.x=scaled_val_Zs,.y=scaled_Zs, ~ Rfast::dista(xnew=.x, x=.y, type = "euclidean") |> as.matrix()), # Compute distance matrix
            mean_by_var_dist = map_dbl(by_var_dist, ~ mean(.x)), # Compute mean distance
            comm_dist = map2(by_var_dist, mean_by_var_dist, ~ .x / .y) # Normalize for commensurability
          )
        
        
      }
      # Aggregate the commensurable distance matrices by summing them up
      cat_dist_mat <- Reduce(`+`, commensurable_dist_structure |> pull(comm_dist))
    }
    ######
    
    
  }else if(scaling=="mca"){
    
    # Create the qs_vec
    qs_vec <- x |>
      dplyr::select(where(is.factor)) |>
      as_tibble() |>
      map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) |>
      unlist()
    
    # Dummy encode categorical features without centering
    prep_Z <- x |>
      recipe(~.) |>
      step_select(all_nominal()) |>
      step_dummy(all_nominal(), one_hot = TRUE) |>
      prep(training = x) 
    
    
    Z = prep_Z |>
      bake(new_data = NULL) |>
      as.matrix()
    
    # Convert the matrix to numeric to avoid any issues with factors/characters
    Z <- apply(Z, 2, as.numeric)
    
    # Calculate the diagonal matrix of category proportions
    D_p <- diag(colSums(Z) / nrow(Z))
    
    # Calculate the inverse square root of D_p
    inv_sq_D_p <- diag(1 / sqrt(diag(D_p)))
    #      print('INSIDE')
    # Compute the category dissimilarity scaled distance matrix
    #cat_dist_mat <- Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
    
    # Compute delta_ms 
    delta_ms <- inv_sq_D_p %*% delta_m %*% inv_sq_D_p # NOTE: ORIGINAL IMPLEMENTATION DOES NOT UPDATE DELTA_MS!
    
    
    if(!commensurable) {
      if(is.null(validate_x)){
        cat_dist_mat = Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
      }else{
        val_Z = bake(prep_Z, new_data = validate_x) |> as.matrix()
        val_Z <- apply(val_Z, 2, as.numeric)
        cat_dist_mat = val_Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
      }
    } else {
      #### Commensurability
      
      prep_Z = x |> map(
        ~as_tibble(.x) |> recipe(~.)|>
          step_dummy(all_predictors(),one_hot = TRUE) |>
          prep(training = as_tibble(.x))
      ) 
      
      Z_list = prep_Z |>
        map(~.x |> bake(new_data=NULL))
      
      Q=map_dbl(x,nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      #delta_out = cat_delta(cats,method = method)
      
      if(is.null(validate_x)){
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_m[levels_identifier==.x,
                                      levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
      }else{
        
        val_Z_list = map2(
          .x=prep_Z, .y=validate_x,
          ~.x |> bake(new_data=as_tibble(.y))
        )
        
        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_m[levels_identifier==.x,
                                      levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          val_Zs = val_Z_list,
          by_var_dist = pmap(.l=list(a=Z_list,b=delta,c=val_Z_list),
                             .f=function(a,b,c){return(as.matrix(c) %*% b %*% t(as.matrix(a)))}),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )
      }
      cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
      
      
      
      
      #       print(cat_dist_mat[1:4,1:4])
    }
  }
  # }
  #    else {
  
  #   }
  
  #}
  
  out_catdist = list()
  out_catdist$distance_mat = cat_dist_mat
  out_catdist$delta = delta_m
  out_catdist$delta_ms = delta_ms
  out_catdist$delta_names = delta_names
  return(out_catdist)
}

