# Nested stage measurements must not force a full collection per fit.
slcn_time <- function(expr) base::system.time(expr, gcFirst = FALSE)
