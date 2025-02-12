  cont_delta_single_attr <- function(x, y = NULL,neighbors = NULL, method = NULL){
#  library(data.table)
#  library(tidyverse)
 # source("R/cont_on_cat_pair.R")
  # x is the categorical attribute
  # y is the distance matrix (mahalanobis) of the continuous
  # obtained computing the Euclidean on the principal components
  
  ####################
  ####################
   #  y = dist(cont_data) |> as.matrix()
   #  x = cat_data$a
   # neighbors=10
  # ####################
  ####################
  
a <- NULL
b <- NULL
n_a <- NULL
n_b <- NULL
sub_dist_mats <- NULL
prop_well_classified <- NULL
phi_weight <-  NULL

  k=neighbors
  # print(neighbors)
  
  x = as.factor(x)
  nlevs=nlevels(x)
  levels(x) = 1:nlevs
  
  single_delta = matrix(0,nlevs,nlevs)
  
  att_freqs = x %>% fct_count()  
  
  tmp_structure  = crossing(a = levels(x) %>% parse_number(), b = levels(x) %>% parse_number()) %>%
    filter(a!=b) %>% 
    mutate(n_a= map_dbl(.x=a,~att_freqs %>% filter(f==.x) %>% pull(n)),
           n_b= map_dbl(.x=b,~att_freqs %>% filter(f==.x) %>% pull(n)),
           prop_a= map2_dbl(.x = n_a,.y = n_b, ~.x/(.x+.y)),
           attr = map2(.x=a,.y=b,~x[x==.x|x==.y]),
           # attr=map(a,~x),
           sub_dist_mats = map2(.x=a,.y=b,~y[x==.x|x==.y,x==.x]),# %>% as.data.table()),
           neighbors=k,
           prop_well_classified = pmap(.l = list(..1 = sub_dist_mats, ..2 = attr, ..3 = neighbors, ..4 = a),
                                       .f = ~cont_on_cat_pair(distance_matrix = ..1, attribute = ..2, neighbors = ..3, target_cat = as.character(..4))$well_class),
           phi_weight = map_dbl(.x = prop_well_classified,~mean(.x>.5)-.5)
           # phi_weight = map_dbl(.x = prop_well_classified,~mean(.x>.5))
    )
  
  selected_structure = tmp_structure %>% dplyr::select(a, b,phi_weight)
  
  for(i in 1:nrow(selected_structure)){
    #single_delta[selected_structure$a[i], selected_structure$b[i]] = single_delta[selected_structure$a[i], selected_structure$b[i]]+ max(0,selected_structure$phi_weight[i])
    #single_delta[selected_structure$b[i], selected_structure$a[i]] = single_delta[selected_structure$b[i], selected_structure$a[i]]+ max(0,selected_structure$phi_weight[i])
    # single_delta[selected_structure$a[i], selected_structure$b[i]] = single_delta[selected_structure$a[i], selected_structure$b[i]]+ selected_structure$phi_weight[i]
    single_delta[selected_structure$a[i], selected_structure$b[i]] = selected_structure$phi_weight[i]
    # single_delta[selected_structure$b[i], selected_structure$a[i]] = single_delta[selected_structure$a[i], selected_structure$b[i]]+ selected_structure$phi_weight[i]
  }
  
  single_delta = (single_delta + t(single_delta))
  single_delta[single_delta<0]=0
  out=list()
  out$structure = tmp_structure
  out$delta=single_delta
  out$single_weight = single_delta[2,1]
  # out$structure=tmp_structure
  return(out)
}




# tmp_structure  = crossing(a = levels(x) %>% parse_number(), b = levels(x) %>% parse_number()) %>%
#   filter(a!=b) %>% 
#   mutate(n_a= map_dbl(.x=a,~att_freqs %>% filter(f==.x) %>% pull(n)),
#          n_b= map_dbl(.x=b,~att_freqs %>% filter(f==.x) %>% pull(n)),
#          prop_a= map2_dbl(.x = n_a,.y = n_b, ~.x/(.x+.y)),
#          attr=map(a,~x),
#          sub_dist_mats = map2(.x=a,.y=b,~rbind(y[x==.x,x==.x],y[x==.y,x==.x])),# %>% as.data.table()),
#          # neighbours=k,
#          neighbours_by_cat = map_dbl(.x=prop_a,~round(k*.x, digits=0)),
#          a_names =  map(.x=sub_dist_mats,~rownames(.x)),
#          b_names =  map(.x=sub_dist_mats,~colnames(.x)),
#          preds_pos = map2(.x = sub_dist_mats,.y= neighbours_by_cat,
#                           .f= function(x = .x,k=.y){
#                             nbs = x %>% as.data.table() %>% map_dfc(~order(.x)[2:(k+1)])
#                           }),
#          preds_names = map2(.x = preds_pos,.y = a_names,
#                             .f= function(x=.x, y = .y){
#                               preds = map_dfc(.x = x, ~y[as.matrix(.x)] %>% parse_number())
#                             }
#          ),
#          preds = map2(.x = preds_names,.y = attr,
#                       .f= function(x=.x, y = .y){
#                         preds = map_dfc(.x = x, ~y[.x] )
#                       }
#          ),
#          prop_hat_a = map2(.x = preds,.y = a,
#                            .f= function(x = .x, y = .y){
#                              
#                              preds = map_dfc(.x = x, ~mean(.x==y))
#                              
#                            }
#          ),
#          phi_weight = map2_dbl(.x=prop_hat_a,.y=prop_a,~mean(.x>=.y)),
#          phi_weight = map_dbl(.x = phi_weight,~ifelse((.x-.5)>0,(.x-.5),0))
#          # phi_weight = map_dbl(.x=prop_hat_a,~mean(.x>=.5))
#   )
