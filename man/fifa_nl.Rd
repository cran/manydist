\name{fifa_nl}
\alias{fifa_nl}
\docType{data}
\title{FIFA 21 Player Data - Dutch League}
\description{
The \code{fifa_nl} dataset contains information on players in the Dutch League from the FIFA 21 video game. This dataset includes various attributes of players, such as demographics, club details, skill ratings, and physical characteristics.
}
\usage{data("fifa_nl")}
\format{
   A data frame with observations on various attributes describing the players.
  \describe{
    \item{\code{player_positions}}{Primary playing positions of the player.}
    \item{\code{nationality}}{The country the player represents.}
    \item{\code{team_position}}{Player's assigned position within their club.}
    \item{\code{club_name}}{Name of the club the player is part of.}
    \item{\code{work_rate}}{The player's work rate, describing defensive and attacking intensity.}
    \item{\code{weak_foot}}{Skill rating for the player's non-dominant foot, ranging from 1 to 5.}
    \item{\code{skill_moves}}{Skill moves rating, indicating technical skill and ability to perform complex moves, on a scale of 1 to 5.}
    \item{\code{international_reputation}}{Player's reputation on an international scale, from 1=local to 3=global star.}
    \item{\code{body_type}}{Body type of the player ( Lean, Normal, Stocky.}
    \item{\code{preferred_foot}}{Dominant foot of the player, either Left or Right.}
    \item{\code{age}}{Age of the player in years.}
    \item{\code{height_cm}}{Height of the player in centimeters.}
    \item{\code{weight_kg}}{Weight of the player in kilograms.}
    \item{\code{overall}}{Overall skill rating of the player out of 100.}
    \item{\code{potential}}{Potential skill rating the player may achieve in the future.}
    \item{\code{value_eur}}{Estimated market value of the player in Euros.}
    \item{\code{wage_eur}}{Player's weekly wage in Euros.}
    \item{\code{release_clause_eur}}{Release clause value in Euros, which other clubs must pay to buy out the player's contract.}
    \item{\code{pace}}{Speed rating of the player out of 100.}
    \item{\code{shooting}}{Shooting skill rating out of 100.}
    \item{\code{passing}}{Passing skill rating out of 100.}
    \item{\code{dribbling}}{Dribbling skill rating out of 100.}
    \item{\code{defending}}{Defending skill rating out of 100.}
    \item{\code{physic}}{Physicality rating out of 100.}
  }
}
\details{
This dataset provides a snapshot of player attributes and performance indicators as represented in FIFA 21 for players in the Dutch League. It can be used to analyze player characteristics, compare skills across players, and explore potential relationships among variables such as age, position, and various skill ratings.
}

\references{
Stefano Leone. (2021). \emph{FIFA 21 Complete Player Dataset}. Retrieved from \url{https://www.kaggle.com/datasets/stefanoleone992/fifa-21-complete-player-dataset}.
}

\examples{
data(fifa_nl)
summary(fifa_nl)
}
\keyword{datasets}
