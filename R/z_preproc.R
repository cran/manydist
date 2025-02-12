z_preproc = function(x,y=NULL,Q){
  out=list()

  .x = NULL

  Z_tib = dummy_cols(x, remove_selected_columns = TRUE) %>% as_tibble()
  Z = dummy_cols(x, remove_selected_columns = TRUE) %>% data.matrix()
  Z_names=colnames(Z)

  if(!is.null(y)){
    Z_y = as_tibble(factor(y)) %>% dummy_cols(remove_selected_columns = TRUE) %>% data.matrix()
    out$Z_y = Z_y
  }


  level_pos = split(data.table(start=c(1,cumsum(Q)[-length(Q)]+1),stop=cumsum(Q)),
                         f = 1:length(Q))
  Z_list = map(.x=level_pos,.f=function(x=.x) Z[,x$start:x$stop])


  ZZ<-t(Z)%*%Z
  zm = colSums(Z)

  ZZod = ZZ
  diag(ZZod)=0
  ZZod = ZZod / zm

  out$ZZod=ZZod
  out$zm=zm
  out$Z= Z_tib
  # out$Z = Z
  out$Z_list = Z_list

  out$level_pos = level_pos
  return(out)
}
