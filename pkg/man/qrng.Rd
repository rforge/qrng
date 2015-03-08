\name{qrng}
\alias{korobov}
\alias{ghalton}
\title{Compute Quasi-Random Sequences}
\description{
  Computing the following quasi-random sequences: Korobov, generalize Halton.
}
\usage{
korobov(n, d, generator)
ghalton(n, d, n.random)
}
\arguments{
  \item{n}{Number \eqn{n} of points to be generated \eqn{\ge 2}{>= 2}.}
  \item{d}{Dimension \eqn{d}.}
  \item{generator}{A \code{\link{numeric}} of length \eqn{d} or length 1
    (in which case it is appropriately extended to length \eqn{d}). All numbers
    must be in \eqn{\{1,\dots,n\}}{{1,...,n}} and must be (coercible to)
    integers.}
  \item{n.random}{An \code{\link{integer}} (or coercible to one) giving
    the number of randomization points.}
}
\value{
  \code{korobov()} and \code{ghalton()} return an
  \eqn{(n,d)}-\code{\link{matrix}}; for \eqn{d=1} an \eqn{n}-vector
  is returned.
}
\details{
  The choice of parameters can be crucial for the quality of
  these quasi-random sequences (only basic sanity checks are
  conducted). Also note that these procedures call fast C code.
}
\author{Marius Hofert}
\references{
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
}
\keyword{distribution}