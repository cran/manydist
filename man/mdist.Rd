\name{mdist}
\alias{mdist}
\title{Calculation of Pairwise Distances for Mixed-Type Data}
\description{
Computes pairwise distances between observations described by numeric and/or categorical attributes, with support for validation data. The function provides options for computing independent, dependent, and practice-based distances as defined in van de Velden et al. (2024), with support for various continuous and categorical distance metrics, scaling, and commensurability adjustments.
}
\usage{
mdist(x, validate_x = NULL, response = NULL, distance_cont = "manhattan", 
      distance_cat = "tot_var_dist", commensurable = FALSE, scaling_cont = "none",
      ncomp = ncol(x), threshold = NULL, preset = "custom")
}
\arguments{
  \item{x}{A dataframe or tibble containing continuous (coded as numeric), categorical (coded as factors), or mixed-type variables.}
  \item{validate_x}{Optional validation data with the same structure as \code{x}. If provided, distances are computed between observations in \code{validate_x} and \code{x}. Default is \code{NULL}.}
  \item{response}{An optional factor for supervised distance calculation in categorical variables, applied only if \code{distance_cat = "supervised"}. Default is \code{NULL}.}
  \item{distance_cont}{Character string specifying the distance metric for continuous variables. Options include \code{"manhattan"} (default) and \code{"euclidean"}.}
  \item{distance_cat}{Character string specifying the distance metric for categorical variables. Options include \code{"tot_var_dist"} (default), \code{"HL"}, \code{"HLeucl"}, \code{cat_dis}, \code{mca}, \code{st_dev},  \code{"matching"}, \code{"eskin"}, \code{"iof"}, \code{"of"}, \code{"goodall_3"}, \code{"goodall_4"}, \code{"gifi_chi2"}, \code{"lin"}, \code{"var_entropy"}, \code{"var_mutability"}, \code{"supervised"}, \code{"supervised_full"}, \code{"le_and_ho"} and all the options in the package philentropy.}
  \item{commensurable}{Logical. If \code{TRUE}, the function adjusts each variable's contribution to ensure equal average influence in the overall distance. Default is \code{FALSE}.}
  \item{scaling_cont}{Character string specifying the scaling method for continuous variables. Options include \code{"none"} (default), \code{"std"}, \code{"range"}, \code{"pc_scores"}, and \code{"robust"}.}
  \item{ncomp}{Integer specifying the number of components to retain when \code{scaling_cont = "pc_scores"}. Default is \code{ncol(x)}.}
  \item{threshold}{Numeric value specifying the percentage of variance explained by retained components when \code{scaling_cont = "pc_scores"}. Overrides \code{ncomp} if specified. Default is \code{NULL}.}
  \item{preset}{Character string specifying pre-defined combinations of arguments. Options include:
    \itemize{
      \item{\code{"custom"} (default): Use specified distance metrics and parameters}
      \item{\code{"gower"}: Gower's distance for mixed data}
      \item{\code{"unbiased_dependent"}: Total variation distance for categorical and Manhattan for standardized continuous}
      \item{\code{"euclidean_onehot"}: Euclidean distance on one-hot encoded categorical and standardized continuous}
      \item{\code{"catdissim"}: Matching distance for categorical and Manhattan for standardized continuous}
    }
  }
}
\value{
A matrix of pairwise distances. If \code{validate_x} is provided, rows correspond to validation observations and columns to training observations.
}

\examples{
library(palmerpenguins)
library(rsample)

# Prepare complete data
pengmix <- palmerpenguins::penguins[complete.cases(palmerpenguins::penguins), ]

# Create training-test split
set.seed(123)
pengmix_split <- initial_split(pengmix, prop = 0.8)
tr_pengmix <- training(pengmix_split)
ts_pengmix <- testing(pengmix_split)

# Example 1: Basic usage with validation data
dist_matrix <- mdist(x = tr_pengmix, 
                    validate_x = ts_pengmix)

# Example 2: Gower preset with validation
dist_gower <- mdist(x = tr_pengmix, 
                   validate_x = ts_pengmix,
                   preset = "gower", 
                   commensurable = TRUE)

# Example 3: Euclidean one-hot preset with validation
dist_onehot <- mdist(x = tr_pengmix, 
                    validate_x = ts_pengmix,
                    preset = "euclidean_onehot")

# Example 4: Custom preset with standardization
dist_custom <- mdist(x = tr_pengmix,
                    validate_x = ts_pengmix,
                    preset = "custom",
                    distance_cont = "manhattan",
                    distance_cat = "matching",
                    commensurable = TRUE,
                    scaling_cont = "std")

# Example 5: PCA-based scaling with threshold
dist_pca <- mdist(x = tr_pengmix,
                 validate_x = ts_pengmix,
                 distance_cont = "euclidean",
                 scaling_cont = "pc_scores",
                 threshold = 0.85)

# Example 6: Categorical variables only
cat_vars <- c("species", "island", "sex")
dist_cat <- mdist(tr_pengmix[, cat_vars],
                 validate_x = ts_pengmix[, cat_vars],
                 distance_cat = "tot_var_dist")

# Example 7: Continuous variables only
num_vars <- c("bill_length_mm", "bill_depth_mm", 
              "flipper_length_mm", "body_mass_g")
dist_cont <- mdist(tr_pengmix[, num_vars],
                  validate_x = ts_pengmix[, num_vars],
                  distance_cont = "manhattan",
                  scaling_cont = "std")

# Example 8: Supervised distance with response
response_tr <- tr_pengmix$body_mass_g
dist_sup <- mdist(tr_pengmix,
                 validate_x = ts_pengmix,
                 response = response_tr,
                 distance_cat = "supervised")


}

\references{
van de Velden, M., Iodice D'Enza, A., Markos, A., Cavicchia, C. (2024). (Un)biased distances for mixed-type data. \emph{arXiv preprint}. Retrieved from \url{https://arxiv.org/abs/2411.00429}.
}

\seealso{
\code{\link{cdist}} for categorical-only distances, \code{\link{ndist}} for continuous-only distances
}
