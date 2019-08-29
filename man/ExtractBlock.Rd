\name{extractBlock}
\alias{extractBlock}

\title{
Extracting a block from a SCM
}
\description{
This function allows to extract the bottom-left block in the sample covariance matrices constructed with the Super Trial.
}
\usage{
extractBlock(scm_st_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scm_st_list}{a list of sample covariance matrices constructed with the Super Trial.
}
}

\value{
It returns a list containing the block of every sample covariance matrix in input list.
}
\references{
Barachant, Alexandre e Marco Congedo (2014). "A Plug&Play P300 BCI Using Information Geometry". In: url: http://arxiv.org/abs/1409.0107.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\seealso{
\code{\link{SuperTrial}}, \code{\link{SCM}}
}

\keyword{ ~block }

