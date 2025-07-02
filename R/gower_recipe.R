gower_recipe <- function(data) {
  .x = NULL
  n_cat <- data |> dplyr::select(where(is.factor)) |> ncol()
  n_con <- data |> dplyr::select(where(is.numeric)) |> ncol()
  
  if (n_cat > 0) {
    gow_recipe <- recipe(~., data = data) |> 
      step_range(all_numeric_predictors(), min = 0, max = 1) |> 
      step_dummy(all_nominal_predictors(), one_hot = TRUE, id = "dummy") |> 
      step_mutate(
        dplyr::across(matches("^f.*_"), ~ .x / 2)
      )
  } else {
    gow_recipe <- recipe(~., data = data) |> 
      step_range(all_numeric_predictors(), min = 0, max = 1)
  }
  
  return(gow_recipe)
}
