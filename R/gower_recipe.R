gower_recipe <- function(data){
  
  .x = NULL
  dname = NULL
  n_cat= data |> select(where(is.factor)) |> ncol()    
  n_con= data |> select(where(is.numeric)) |> ncol()    
  
  if(n_cat>0){
    dummy_names = recipe(~.,data=data) |> 
      step_dummy(all_nominal_predictors(),one_hot = TRUE,id="dummy") |> prep(training=data) |> 
      tidy(id="dummy") |> unite(col="dname",c(terms,columns)) |> pull(dname)
    
    gow_recipe = recipe(~.,data=data) |> 
      step_range(all_numeric_predictors(), min = 0, max = 1) |>
      step_dummy(all_nominal_predictors(),one_hot = TRUE,id="dummy") |> 
      step_mutate(across(all_of(dummy_names), 
                         function(x=.x){
                           return((x/2))
                         }
      )
      )
    
  }else{
    gow_recipe = recipe(~.,data=data) |> 
      step_range(all_numeric_predictors(), min = 0, max = 1) 
  }
  
  
  # n_levs = data |> select(where(is.factor)) |> map(~rep(nlevels(.x)-1,nlevels(.x))) |> unlist()
  # names(n_levs)=dummy_names
  
  return(gow_recipe)
}
