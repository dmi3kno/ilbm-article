# initial values for the g-and-h model
make_initials_gnh <- function(){
  set.seed(42)
  matrix(c(runif(4,2,7), runif(4, 0.1,0.5)),
                ncol = 2, byrow = FALSE, dimnames = list(NULL, c("g","h")))
}
