## Add a warning message when weight_sys is not user-defined and weights is differnt that all ones
## Add a check when the length of weights vector is different from nvars

ndist <- function(x,validate_x=NULL, commensurable=FALSE, method = "manhattan",sig=NULL,
                  scaling = "none" ,ncomp=ncol(x),threshold=NULL,weights=rep(1,ncol(x))){
  
 # source("manydist/R/mahalanobis_like.R")
  
  .x <- NULL
  y <- NULL
  by_var_dist_w <- NULL
  
  # x = iris |> select(where(is.numeric))
  # x_split=initial_split(x,prop=0.7)
  # x = training(x_split)
  # validate_x = testing(x_split)
  
  
  rec_x = recipe(~., data = x)
  tr_x=x
  
  if(scaling == "std"){
    
    x = rec_x |> step_scale(all_predictors()) |> prep(training = tr_x) |> bake(new_data=NULL)
    if(!is.null(validate_x)){
      validate_x = rec_x |> step_scale(all_predictors()) |> prep(training = tr_x) |> bake(new_data=validate_x)
    }
  }else if(scaling == "range"){
    x = rec_x |> step_range(all_predictors()) |> prep(training = tr_x) |> bake(new_data=NULL)
    if(!is.null(validate_x)){
      validate_x = rec_x  |> step_range(all_predictors())  |> prep(training = tr_x) |> bake(new_data=validate_x)
    }
  }else if(scaling == "pc_scores" ){
    
    if(is.null(threshold)){
      x = rec_x |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
        step_pca(all_predictors(),num_comp=ncomp) |> prep(training = tr_x) |> bake(new_data=NULL)
    
      if(!is.null(validate_x)){
        validate_x = rec_x  |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
          step_pca(all_predictors(),num_comp=ncomp)  |> prep(training = tr_x) |> bake(new_data=validate_x)
      } 
      
    }else{
      x = rec_x |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
        step_pca(all_predictors(),threshold=threshold) |> prep(training = tr_x) |> bake(new_data=NULL)
      
      if(!is.null(validate_x)){
        
        validate_x = rec_x |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
          step_pca(all_predictors(),threshold=threshold) |> prep(training = tr_x) |> bake(new_data=validate_x)
      }
      
    }
  }else if(scaling=="robust"){
    # old_x = x  |>  map(.f=function(x=.x){
    #   x=(x-median(x))/IQR(x)
    #   return(x)
    # }) |> dplyr::bind_cols()
    
    
    x = rec_x |> step_mutate(across(everything(), ~ (. - median(.))/IQR(.))) |> 
      prep(training = tr_x) |> bake(new_data=NULL)
    
    if(!is.null(validate_x)){
      pre_validate_x = rec_x |> 
        step_mutate(across(everything(), ~ (. - median(.))/IQR(.))) |>  
        prep(training = tr_x) |> 
        bake(new_data=validate_x)
    }
    
  }else if(scaling=="none"){
    # print("scaling=none does nothing")
    
  }
  
  if(method=="mahalanobis"){
    
      distance=mahalanobis_like(x=x,validate_x=validate_x,sig=sig)
  }else{
    
    if(commensurable == TRUE){
      if(is.null(validate_x)){
        by_var_dist = map(.x=as_tibble(x),.f=function(x=.x){
          b_v_d = daisy(data.frame(x),metric=method, warnBin = FALSE) %>% as.matrix()
          b_v_d = b_v_d/mean(b_v_d)
          return(b_v_d)
        })}else{
          by_var_dist = map2(.x=as_tibble(x),.y=as_tibble(validate_x),.f=function(x=.x,xnew=y){
            b_v_d = Rfast::dista(xnew=xnew,x=x,type=method) %>% as.matrix()
            b_v_d = b_v_d/mean(b_v_d)
            return(b_v_d)
          })
        }
    }else if (commensurable==FALSE){
    
    if(is.null(validate_x)){
      by_var_dist = map(.x=as_tibble(x),.f=function(x=.x){
        b_v_d = daisy(data.frame(x),metric=method, warnBin = FALSE)
        return(b_v_d)
      })
    }else{
      by_var_dist = map2(.x=as_tibble(x),.y=as_tibble(validate_x),.f=function(x=.x,xnew=y){
        b_v_d = Rfast::dista(xnew=xnew,x=x,type=method)  |>  as.matrix()
        return(b_v_d)
      })
    }
    }
    if (method == "euclidean") {
      by_var_structure = tibble(by_var_dist=by_var_dist,weights=weights) |>
        #  mutate(by_var_dist_w= map(.x=by_var_dist,~.x^2))
        mutate(by_var_dist_w= map2(.x=by_var_dist,.y=weights,~(.x^2)*.y))
      
      distance = as.matrix(Reduce(`+`,by_var_structure |> pull(by_var_dist_w)))
      distance = sqrt(distance)
      
    } else {
      
      by_var_structure = tibble(by_var_dist=by_var_dist,weights=weights) |>
        mutate(by_var_dist_w= map2(.x=by_var_dist,.y=weights,~.x*.y))
      
      distance = as.matrix(Reduce(`+`,by_var_structure |> pull(by_var_dist_w)))
      
    }
    
  }

return(distance)

}

