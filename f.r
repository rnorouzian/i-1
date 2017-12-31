
foo <- function(x)
{
  UseMethod("foo")
}

foo.d <- function(x){
  
  return(c(x + 1, x - 1))
}