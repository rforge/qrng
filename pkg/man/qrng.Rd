\name{qrng}
\alias{korobov}
\alias{ghalton}
\title{Compute Quasi-Random Sequences}
\description{
  Computing the following quasi-random sequences: Korobov, generalize Halton.
}
\usage{
korobov(n, d, generator)
ghalton(n, d, method=c("generalized", "halton"))
}
\arguments{
  \item{n}{Number \eqn{n} of points to be generated \eqn{\ge 2}{>= 2}.}
  \item{d}{Dimension \eqn{d}.}
  \item{generator}{A \code{\link{numeric}} of length \eqn{d} or length 1
    (in which case it is appropriately extended to length \eqn{d}). All numbers
    must be in \eqn{\{1,\dots,n\}}{{1,...,n}} and must be (coercible to)
    integers.}
  \item{method}{A \code{\link{character}} string indicating which
    sequence is generated, generalized Halton or (plain) Halton.}
}
\value{
  \code{korobov()} and \code{ghalton()} return an
  \eqn{(n,d)}-\code{\link{matrix}}; for \eqn{d=1} an \eqn{n}-vector
  is returned.
}
\details{
  Note that these procedures call fast C code.

  The choice of parameters for \code{korobov()} is crucial for the quality of
  this quasi-random sequence (only basic sanity checks are
  conducted). For more details, see l'Ecuyer and Lemieux (2000).

  For \eqn{d\ge360}{d>=360}, the generalized Halton sequences uses the
  identity to TODO.
}
\author{Marius Hofert}
\references{
  Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
  A Comparative Study. \emph{ACM-TOMACS} \bold{19}(4), Article 15.

  l'Ecuyer, P., Lemieux, C. (2000). Variance Reduction via Lattice
  Rules. \emph{Stochastic Models and Simulation}, 1214--1235.

  Lemieux, C., Cieslak, M., Luttmer, K. (2004). RandQMC User's guide.
  See https://www.math.uwaterloo.ca/~clemieux/randqmc/guide.pdf
}
%% \seealso{
%% }
\examples{
n <- 2^10
d <- 4
generator <- 17
u <- korobov(n, d=d, generator=generator)
pairs(u, gap=0, pch=".", labels=as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

u <- ghalton(n, d)
pairs(u, gap=0, pch=".", labels=as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))
}
\keyword{distribution}