x <- rnorm(8, 10, 1)
y <- rnorm(8, 10, 1)
x.scaled <- x/x[1]
y.scaled <- y/y[1]
t.test(x = x.scaled[c(1,3,5,7,9)], y = y.scaled[c(2,4,6,8,10)])

library("plyr")
A <- array(rnorm(4000, 10, 2), dim = c(1000, 4))
B <- array(rnorm(4000, 10, 2), dim = c(1000, 4))
A.scaled <- A/A[,1]
B.scaled <- B/B[,1]
#AB <- cbind(A,B)
AB <- cbind(A.scaled,B.scaled)
head(AB)
tt <- function(x){
  return(t.test(x = x[c(1,3,5,7)], y = x[c(2,4,6,8)]))
}
l <- alply(AB, .margins = 1, .fun = tt)
get.p <- function(x){
  return(x$p.value)
}
df <- ldply(l, .fun = get.p)
hist(df$V1)
summary(df$V1 < 0.05)
