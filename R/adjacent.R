#R

## Tobias.Kockmann <ethz.ch>
## 20170315

# $HeadURL: svn+ssh://r21/home/tobiasko/__SVNREPO/posprot/posprot/R/adjacent.R $
# $Id: adjacent.R 5 2017-06-26 13:42:53Z tobiasko $
# $Date: 2017-06-26 15:42:53 +0200 (Mon, 26 Jun 2017) $

.adjacent <- function(x){

## helper function for split+apply+combine

  s <- match(end(x)+1, start(x))
  q <- 1:length(x)
  y <- Hits(from = q[!is.na(s)], to = s[!is.na(s)], nLnode = length(x), nRnode = length(x))
  return(y)

}
