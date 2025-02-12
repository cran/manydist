commensurable_weight<-function(cat_data,delta_cat){

  delta_list = NULL
  margins = NULL
  weight = NULL

  Q=map_dbl(cat_data,nlevels)

  # delta_cat= matrix(1,sum(Q),sum(Q))
  #
  # diag(delta_cat)=0

  nvar=length(Q)
  delta_structure = data.table(start=c(1,cumsum(Q)[-length(Q)]+1),stop=cumsum(Q)) |> as_tibble() |>
    mutate(
      delta_list = map2(.x=start,.y=stop,~delta_cat[.x:.y,.x:.y]),
      # # Lin quick fix: remove if not needed
      # delta_list = map(.x=delta_list,function(x=.x){
      #   if(sum(x)==0){
      #   x=matrix(.5,nrow(x),ncol(x))
      #   diag(x)=0
      # }
      # return(x)
      #   }),
      #### end of Lin quick fix
  margins= map(.x=cat_data,~fct_count(.x,prop=TRUE)|> pull(p)),
  weight = 1/map2_dbl(.x=delta_list,.y=margins,.f = ~t(.y)%*% .x %*% .y) #|> sqrt()
  )

  w = rep(delta_structure |> pull(weight),Q)


  return(w)

}

