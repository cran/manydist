\name{ndist}
\alias{ndist}
\title{Calculation of Pairwise Distances for Continuous Data}
\description{
Computes a distance matrix for continuous data with support for multiple distance metrics, scaling methods, dimensionality reduction, and validation data. The function implements various distance calculation approaches as described in van de Velden et al. (2024), including options for commensurable distances and variable weighting.
}
\usage{
ndist(x, validate_x = NULL, commensurable = FALSE, method = "manhattan",
      sig = NULL, scaling = "none", ncomp = ncol(x), threshold = NULL,
      weights = rep(1, ncol(x)))
}
\arguments{
  \item{x}{A data frame or matrix of continuous input variables.}
  \item{validate_x}{Optional data frame or matrix for validation data. If provided, distances are computed between observations in \code{validate_x} and \code{x}. Default is \code{NULL}.}
  \item{commensurable}{Logical. If \code{TRUE}, standardizes each variable's distance matrix by dividing by its mean distance, making distances comparable across variables. Default is \code{FALSE}.}
  \item{method}{Character string specifying the distance metric. Options include \code{"manhattan"}, \code{"euclidean"}, and \code{"mahalanobis"}. Default is \code{"manhattan"}.}
  \item{sig}{Covariance matrix to be used when \code{method = "mahalanobis"}. If \code{NULL}, computed from the data. Default is \code{NULL}.}
  \item{scaling}{Character string specifying the scaling method. Options:
    \itemize{
      \item \code{"none"}: No scaling
      \item \code{"std"}: Standardization (zero mean, unit variance)
      \item \code{"range"}: Min-max scaling to [0,1]
      \item \code{"pc_scores"}: PCA-based dimensionality reduction
      \item \code{"robust"}: Robust scaling using median and IQR
    }
    Default is \code{"none"}.}
  \item{ncomp}{Number of principal components to retain when \code{scaling = "pc_scores"}. Default is the number of columns in \code{x}.}
  \item{threshold}{Proportion of variance to retain when \code{scaling = "pc_scores"}. If specified, overrides \code{ncomp}. Default is \code{NULL}.}
  \item{weights}{Numeric vector of weights for each variable. Must have length equal to the number of variables in \code{x}. Default is a vector of ones.}
}
\details{
  The \code{ndist} function provides a comprehensive framework for distance calculations in continuous data:
  
  \itemize{
    \item When \code{validate_x} is provided, computes distances between observations in \code{validate_x} and \code{x}.
    \item Supports multiple scaling methods that can be applied before distance calculation.
    \item PCA-based dimensionality reduction can be controlled either by number of components or variance threshold.
    \item For Mahalanobis distance, handles singular covariance matrices with appropriate error messages.
    \item Implements commensurable distances for better comparability across variables.
  }
  
  \strong{Warning:} The function validates:
  \itemize{
    \item Weight vector length must match the number of variables
    \item Covariance matrix singularity for Mahalanobis distance
    \item Compatibility of \code{x} and \code{validate_x} dimensions
  }
}
\value{
A distance matrix where element [i,j] represents the distance between:
\itemize{
  \item observation i and j of \code{x} if \code{validate_x} is \code{NULL}
  \item observation i of \code{validate_x} and observation j of \code{x} if \code{validate_x} is provided
}
}
\examples{
library(palmerpenguins)
library(rsample)

penguins_cont <- palmerpenguins::penguins[, c("bill_length_mm",
"bill_depth_mm", "flipper_length_mm", "body_mass_g")]
penguins_cont <- penguins_cont[complete.cases(penguins_cont), ]

# Basic usage
dist_matrix <- ndist(penguins_cont)

# Commensurable distances with standardization
dist_matrix <- ndist(penguins_cont, 
                    commensurable = TRUE, 
                    scaling = "std")

# PCA-based dimensionality reduction
dist_matrix <- ndist(penguins_cont, 
                    scaling = "pc_scores", 
                    threshold = 0.95)

# Mahalanobis distance
dist_matrix <- ndist(penguins_cont, 
                    method = "mahalanobis")

# Weighted Euclidean distance
dist_matrix <- ndist(penguins_cont, 
                    method = "euclidean",
                    weights = c(1, 0.5, 2, 1))
                    
# Training-test split example with validation data
set.seed(123)
# Create training-test split using rsample
penguins_split <- initial_split(penguins_cont, prop = 0.8)
tr_penguins <- training(penguins_split)
ts_penguins <- testing(penguins_split)

# Basic usage with training data only
dist_matrix <- ndist(tr_penguins)

# Computing distances between test and training sets
val_dist_matrix <- ndist(x = tr_penguins, 
                        validate_x = ts_penguins,
                        method = "euclidean")

# Using validation data with standardization
val_dist_matrix_std <- ndist(x = tr_penguins,
                            validate_x = ts_penguins,
                            scaling = "std",
                            method = "manhattan")

# Validation with PCA and commensurability
val_dist_matrix_pca <- ndist(x = tr_penguins,
                            validate_x = ts_penguins,
                            scaling = "pc_scores",
                            ncomp = 2,
                            commensurable = TRUE)

# Validation with robust scaling and custom weights
val_dist_matrix_robust <- ndist(x = tr_penguins,
                               validate_x = ts_penguins,
                               scaling = "robust",
                               weights = c(1, 0.5, 2, 1))

# Mahalanobis distance with validation data
val_dist_matrix_mahal <- ndist(x = tr_penguins,
                              validate_x = ts_penguins,
                              method = "mahalanobis")

}
\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}
\seealso{
\code{\link{mdist}} for mixed-type data distances, \code{\link{cdist}} for categorical data distances.
}
