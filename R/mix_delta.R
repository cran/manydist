mix_delta <- function(x,y=NULL,prop_nn=.10, number_of_components=2, 
                      method_cat = "tot_var_dist", 
                      con_on_cat = NULL,
                      mahalanobize=TRUE){# 
  
 # source("manydist/R/cont_delta_single_attr.R")
#  source("manydist/R/cat_delta.R")
  

  .x=NULL
  y = NULL
  
  mixed_dataset  = x

  
  if(con_on_cat=="neighbors_based"){
    n_neighb = (nrow(mixed_dataset)*prop_nn) |> round(digits=0)
    
    mixed_dataset = mixed_dataset |> 
      mutate(across(where(is.factor),.fns = function(x=.x,n_ne = n_neighb){
        low_freq = min(table(x))
        nn_thresh = n_ne/2
        if(low_freq <  nn_thresh){
          x = fct_lump_lowfreq(x)
        }else{
          x = x
        }
        return(x)
      }
      )
      )
    
    
    
    cont_data_only = mixed_dataset |> dplyr::select(where(is.numeric))
    
    
    ### amarkos: 7.2.2025. Check if method should be NULL
    cat_delta_from_cats = cat_delta(x = mixed_dataset |> dplyr::select(where(is.factor)),y=y,
                                    method = NULL, method_cat = method_cat)
    
    
    
    delta_cat = cat_delta_from_cats[[method_cat]] |>  data.matrix()
    
    cat_delta_names = cat_delta_from_cats$delta_names
    Z_cat = cat_delta_from_cats$Z
    
    mh_dist = cont_data_only |>  dist() |>  as.matrix() 
    
    delta_con_list = mixed_dataset  |>  dplyr::select(where(is.factor)) |>
      map(.f = function(x=.x){
        out_delta = cont_delta_single_attr(x = x, y = mh_dist, neighbors = n_neighb)$delta
        
      })
    
    
    
    if(mahalanobize){
      #Mahalanobis Euclidean
      mh_dist = cont_data_only |>  as.matrix() |>  distances(normalize = "mahalanobize") |> as.matrix() 
    }
    
    delta_con = bdiag(delta_con_list)  |>  as.matrix()
    
    
    delta_mix = delta_cat + delta_con
    # delta_mix = delta_cat * delta_con
    
    ####################################
    ####################################  
    # delta_mix = delta_cat
    ####################################
    ####################################
  }else if(con_on_cat=="relevance_based"){
    warning("no GUDMM-based relevance index implementation is available yet")
  } else{
    
    cont_data_only = mixed_dataset |> dplyr::select(where(is.numeric))
    cat_delta_from_cats = cat_delta(x = mixed_dataset |> dplyr::select(where(is.factor)),y=y,
                                    method = NULL, method_cat = method_cat)
    
    
    
    delta_cat = cat_delta_from_cats[[method_cat]] |> data.matrix()
    
    cat_delta_names = cat_delta_from_cats$delta_names
    Z_cat = cat_delta_from_cats$Z
    
    mh_dist = cont_data_only |> dist() |> as.matrix() 
    
    if(mahalanobize){
      mh_dist = cont_data_only |>  as.matrix() |>  distances(normalize = "mahalanobize") |> as.matrix() 
    }
    
    delta_mix = delta_cat
  }
  
  out=list()
  out$preproc_train=mixed_dataset
  out$delta = delta_mix
  out$delta_cat = delta_cat
  out$cat_delta_names = cat_delta_names
  out$cont_distance = mh_dist
  out$Z = Z_cat 
  
  return(out)
}
