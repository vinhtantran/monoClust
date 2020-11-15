anneal <- function(val_vec,
                   max_tries,
                   init_temp,
                   alpha,
                   stop_temp,
                   init_val,
                   verbose = c("all", "some", "none"),
                   datamems, cutsmems, dist, mems) {

  verbose <- match.arg(verbose)

  current_temp <- init_temp
  current_val <- init_val

  if (verbose == "all" || verbose == "some") {
    print("ORI: temp = ", current_temp,
          ";\tvalue = ", current_val,
          ";\tinertia = ", inertia(val_vec, current_val))
  }

  while (current_temp >= stop_temp) {
    for (i in 1:max_tries) {
      next_val <- pertub(current_val)

      # Objective function is to minimize the inertia
      delta_inertia =
        inertia(val_vec, next_val) -
        inertia(val_vec, current_val)

      # If inertia decrease, or if it increases but the temp is high
      if ((delta_inertia < 0) ||
          exp(-delta_inertia / current_temp) >= runif(0, 10)) {
        current_val <- next_val
      }
    }

    current_temp <- alpha * current_temp

    if (verbose == "all") {
      print("RUN: temp = ", current_temp,
            ";\tvalue = ", current_val,
            ";\tinertia = ", inertia(val_vec, current_val))
    }

  }

  if (verbose == "all" || verbose == "some") {
    print("END: temp = ", current_temp,
          ";\tvalue = ", current_val,
          ";\tinertia = ", inertia(val_vec, current_val))
  }
}

pertub <- function(current_val, val_vec) {
  return(sample(val_vec, 1))
}

mult_inertia <- function(current_val, datamems, cutsmems, dist, mems) {
  data_col <- dplyr::pull(datamems, current_val)
  cuts_col <- dplyr::pull(cutsmems, current_val)

  new_inertia <- purrr::map_dbl(cuts_col, function(x) {
    mems_a <- mems[which(data_col < x)]
    mems_b <- setdiff(mems, mems_a)
    ifelse(length(mems_a) * length(mems_b) == 0L,
           NA,
           inertia_calc(dist[mems_a, mems_a]) +
             inertia_calc(dist[mems_b, mems_b]))
  })

  return(new_inertia)
}

inertia_calc <- function(X) {

  if (!is.numeric(X) && !is.matrix(X))
    stop("X has to be a numerical value or matrix.")

  # If singleton cluster, inertia is 0
  inertia_value <- ifelse(length(X) == 1 && is.numeric(X),
                          0,
                          sum(X^2) / (dim(X)[1] * 2))
  return(inertia_value)
}



library(cluster)
data(ruspini)
val_vec <- ruspini$x
max_tries <- round(length(val_vec) / 4)

avg_cost_increase <- 200
num_cost_increases <- 100
acc_ratio <- 0.75
prob_e <- 0.00000000001
beta <- 0.125
num_temp <- 200

init_temp <-
  avg_cost_increase /
  log(num_cost_increases /
        ((num_cost_increases * acc_ratio) -
           (1 - acc_ratio) *
           (max_tries - num_cost_increases)))

stop_temp <- -beta * avg_cost_increase / log(prob_e)

alpha <- (stop_temp / init_temp)^(1 / num_temp)

anneal(val_vec,
       max_tries,
       init_temp,
       alpha,
       stop_temp,
       init_val = val_vec[1],
       "all")
