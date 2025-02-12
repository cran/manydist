\name{cdist}
\alias{cdist}
\title{Calculation of Pairwise Distances for Categorical Data}
\description{
Computes a distance matrix for categorical variables with support for validation data, multiple distance metrics, and variable weighting. The function implements various distance calculation approaches as described in van de Velden et al. (2024), including commensurable distances and supervised options when response variable is provided.
}
\usage{
cdist(x, response = NULL, validate_x = NULL, method = "tot_var_dist",
      commensurable = FALSE, weights = 1)
}
\arguments{
  \item{x}{A data frame or matrix of categorical variables (factors).}
  \item{response}{Optional response variable for supervised distance calculations. Default is \code{NULL}.}
  \item{validate_x}{Optional validation data frame or matrix. If provided, distances are computed between observations in \code{validate_x} and \code{x}. Default is \code{NULL}.}
  \item{method}{Character string or vector specifying the distance metric(s). Options include:
    \itemize{
      \item \code{"tot_var_dist"}: Total variation distance (default)
      \item \code{"HL"}, \code{"HLeucl"}: Hennig-Liao distance
      \item \code{"cat_dis"}: Category-based dissimilarity
      \item \code{"mca"}: Multiple correspondence analysis based
      \item \code{"st_dev"}: Standard deviation based
      \item \code{"matching"}, \code{"eskin"}, \code{"iof"}, \code{"of"}: Various coefficients
      \item \code{"goodall_3"}, \code{"goodall_4"}: Goodall-based distances
      \item \code{"gifi_chi2"}: Gifi chi-square distance
      \item \code{"lin"}: Lin's similarity measure
      \item \code{"var_entropy"}, \code{"var_mutability"}: Variability-based measures
      \item \code{"supervised"}, \code{"supervised_full"}: Response-guided distances
      \item \code{"le_and_ho"}: Le and Ho distance
      \item Additional methods from philentropy package
    }
    Can be a single string or vector for different methods per variable.}
  \item{commensurable}{Logical. If \code{TRUE}, standardizes each variable's distance matrix by dividing by its mean distance. Default is \code{FALSE}.}
  \item{weights}{Numeric vector or matrix of weights. If vector, must have length equal to number of variables. If matrix, must match the dimension of level-wise distances. Default is 1 (equal weighting).}
}
\details{
The \code{cdist} function provides a comprehensive framework for categorical distance calculations:

\itemize{
  \item Supports multiple distance calculation methods that can be specified globally or per variable
  \item Handles validation data through \code{validate_x} parameter
  \item Implements supervised distances when response variable is provided
  \item Supports commensurable distances for better comparability across variables
  \item Provides flexible weighting schemes at variable and level granularity
}

\strong{Important notes:}
\itemize{
  \item Input variables are automatically converted to factors with dropped unused levels
  \item Different methods per variable is not supported for \code{"none"}, \code{"st_dev"}, \code{"HL"}, \code{"cat_dis"}, \code{"HLeucl"}, \code{"mca"}
  \item Weight vector length must match the number of variables when specified as a vector
  \item Variables should be factors; numeric variables will cause errors
}
}
\value{
A list containing:
  \item{distance_mat}{Matrix of pairwise distances. If \code{validate_x} is provided, rows correspond to validation observations and columns to training observations.}
  \item{delta}{Matrix or list of matrices containing level-wise distances for each variable.}
  \item{delta_names}{Vector of level names used in the delta matrices.}
}
\examples{
library(palmerpenguins)
library(rsample)

# Prepare data with complete cases for both categorical variables and response
complete_vars <- c("species", "island", "sex", "body_mass_g")
penguins_complete <- penguins[complete.cases(penguins[, complete_vars]), ]
penguins_cat <- penguins_complete[, c("species", "island", "sex")]
response <- penguins_complete$body_mass_g

# Create training-test split
set.seed(123)
penguins_split <- initial_split(penguins_cat, prop = 0.8)
tr_penguins <- training(penguins_split)
ts_penguins <- testing(penguins_split)
response_tr <- response[penguins_split$in_id]
response_ts <- response[-penguins_split$in_id]

# Basic usage
result <- cdist(tr_penguins)

# With validation data
val_result <- cdist(x = tr_penguins, 
                   validate_x = ts_penguins,
                   method = "tot_var_dist")
                   
# ...and commensurability
val_result_COMM <- cdist(x = tr_penguins, 
                   validate_x = ts_penguins,
                   method = "tot_var_dist",
                   commensurable = TRUE)

# Supervised distance with response variable
sup_result <- cdist(x = tr_penguins, 
                   response = response_tr,
                   method = "supervised")

# Supervised with validation data
sup_val_result <- cdist(x = tr_penguins,
                       validate_x = ts_penguins,
                       response = response_tr,
                       method = "supervised")

# Commensurable distances with custom weights
comm_result <- cdist(tr_penguins,
                    commensurable = TRUE,
                    weights = c(2, 1, 1))

# Different methods per variable
multi_method <- cdist(tr_penguins,
                     method = c("matching", "goodall_3", "tot_var_dist"))

}
\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}
\seealso{
\code{\link{mdist}} for mixed-type data distances, \code{\link{ndist}} for continuous data distances
}
