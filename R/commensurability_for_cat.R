commensurability_for_cat<-function(x,y=NULL,validate_x=NULL, delta) {

#  library(tidymodels)
#  library(tidyverse)
#  library(data.table)
#  library(Matrix)
  .x = NULL
  by_var_dist = NULL
  mean_by_var_dist = NULL
  factor_name = NULL
  by_var_dist = NULL
  comm_dist = NULL
  weight = NULL


  cats = x

  prep_Z = cats |> map(
    ~as_tibble(.x) |> recipe(~.)|>
      step_dummy(all_predictors(),one_hot = TRUE) |>
      prep(training = as_tibble(.x))) 
    
    Z_list = prep_Z |>
      map(~bake(.x,new_data=NULL))
  

  Q=map_dbl(cats,nlevels)
  levels_identifier = rep(names(Q), times = as.vector(Q))
  #delta_out = cat_delta(cats,method = method)

if(is.null(validate_x)){
    commensurable_dist_structure = tibble(factor_name = names(Q)) |>
    mutate(delta = map(.x=factor_name,
                       ~delta[levels_identifier==.x,
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
  
  validate_Z_list =  map2(.x=validate_x, .y=prep_Z,
                          ~ .y |> bake(new_data=as_tibble(.x)))
  
  
  commensurable_dist_structure = tibble(factor_name = names(Q)) |>
    mutate(
      delta = map(.x=factor_name,
                  ~delta[levels_identifier==.x,
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
  
}
  distance = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
  
  
  return(distance)

}
