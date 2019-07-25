\name{perm.fdr}
\alias{perm.fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single variable FDR correction from Joaquín Sanz}
\description{Input empirical P-values }
\usage{
perm.fdr(input_df,perm_df,Pvals_col_name,name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_df}{
  Input_df is N genes by X Variables (can include other betas). Column "Pvals_col_name" contains the vector of empirical p-values for the permutation.
  }
  \item{perm_df}{Perm_df is (N genes x # of permutations; e.g. 100 genes *10 permutations, 1000 rows in total) by 1 (Permutation variable of interest)}
  \item{Pvals_col_name}{Column name in "input_df" that contains the empirical p-values of interest}
  \item{name}{Name of output DF with empirical FDR corrections -- FDR_"name"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
