
foo <- function(x)
{
  UseMethod("foo")
}

foo.default <- function(x){
  
  return(c(x + 1, x - 1))
}
