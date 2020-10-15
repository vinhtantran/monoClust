#' Parallel Coordinates Plot with Circular Variables
#'
#' Making a parallel coordinates plot with the circular variables are plotted
#' as ellipses. The function currently works well with data with one circular
#' variable.
#'
#' @param data Data set.
#' @param circ.var Circular variable(s) in the data set, indicated by names
#'   or index in the data set.
#' @param order.appear The order of appearance of the variables, listed by a
#'   vector of names or indices.
#' @param linetype Line type. Default is solid line. See details in
#'   `vignette("ggplot2-specs")`.
#' @param size Size of a line is its width in mm. Default is 0.5. See details in
#'   `vignette("ggplot2-specs")`.
#' @param object A `MonoClust` object which contains the cluster results
#'   and memberships.
#' @param cluster.col Color of clusters, indicating by a vector. The length of this
#'   vector must be equal to the number of clusters in the object object.
#' @param is.degree Whether the unit of the circular variables is degree or not
#'   (radian). Default is `TRUE`.
#' @param shift The shift (offset, rotation) of the circular variable, in
#'   radians. Default is 0 (no rotation).
#' @param alpha The transparency of the lines. Default is 0.1.
#' @param show.medoids Whether to highlight the median lines or not. Default is
#'   `TRUE`.
#' @param north What value of the circular variable is labeled North. Default is
#'   0 radian.
#' @param cw Which direction of the circular variable is considered increasing
#'   in value, clockwise (`TRUE`) or counter-clockwise (`FALSE`). Default is
#'   `TRUE`.
#' @param labelsize The size of labels on the plot. Default is 4.
#' @param xlab Labels for x-axis.
#' @param ylab Labels for y-axis.
#' @param legend.cluster Labels for group membership. Implemented by setting
#'   label for ggplot `color` aesthetics.
#'
#' @return A ggplot2 object.
#' @import dplyr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' # Set color constant
#' COLOR4 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#' # Reduce the size of the data for for sake of example speed
#' set.seed(12345)
#' wind_reduced <- wind_sensit_2007[sample.int(nrow(wind_sensit_2007), 50), ]
#'
#' sol42007 <- MonoClust(wind_reduced, cir.var = 3, nclusters = 4)
#'
#' library(ggplot2)
#' pcp.gg(data = wind_reduced,
#'        circ.var = "WDIR",
#'        # To improve aesthetics
#'        shift = pi*3/4-0.3,
#'        order.appear = c("WDIR", "has.sensit", "WS"),
#'        alpha = 0.5,
#'        object = sol42007,
#'        cluster.col = COLOR4) +
#'   theme(panel.background = element_rect(color = "white"),
#'         panel.border = element_rect(color = "white", fill = NA),
#'         panel.grid.major = element_line(color = "#f0f0f0"),
#'         panel.grid.minor = element_blank(),
#'         axis.line = element_line(color = "black"),
#'         legend.key = element_rect(color = NA),
#'         legend.position = "bottom",
#'         legend.direction = "horizontal",
#'         legend.title = element_text(face = "italic"),
#'         legend.justification = "center")
pcp.gg <- function(data,
                   circ.var = NULL, is.degree = TRUE, shift = 0, north = 0,
                   cw = FALSE,
                   order.appear = NULL, linetype = 1, size = 0.5, alpha = 0.1,
                   object, cluster.col = NULL, show.medoids = TRUE,
                   labelsize = 4,
                   xlab = "Variables", ylab = NULL, legend.cluster = "groups") {

  if (!is.data.frame(data))
    stop("\"data\" must be a data frame.")
  if (!all(purrr::map_lgl(data, is.numeric)))
    stop("Function only supports numerical variables.")
  if (!inherits(object, "MonoClust"))
    stop("\"object\" must be a MonoClust object.")

  # Commonly used variables
  all_var <- colnames(data)

  if (!is.null(circ.var) && !(circ.var %in% all_var))
    stop("If exist, \"circ.var\" must be columns in \"data\".")

  # Commonly used variables
  numeric_var <- setdiff(colnames(data), circ.var)

  if (!is.null(order.appear)) {
    # Check if order.appear includes variables in the data
    if (!all(order.appear %in% all_var)) {
      stop("The \"order.appear\" has to be a subset of data variables.")
    }
  } else {
    order.appear <- all_var
  }

  # Commonly used variables
  cluster_mem <- object$Membership

  if (is.null(cluster.col)) {
    cluster.col <- seq_len(length(unique(cluster_mem)))
  } else {
    if (length(cluster.col) != length(unique(cluster_mem))) {
      message("The number of colors does not match the number of clusters.
            Default values will be used.")
      cluster.col <- seq_len(length(unique(cluster_mem)))
    }

    if (is.null(names(cluster.col)) |
        !identical(sort(names(cluster.col)),
                   sort(as.character(unique(cluster_mem))))) {
      if (!is.null(names(cluster.col)))
        message("Named \"cluster.col\" does not match cluster names. Default
              values will be used.")
      names(cluster.col) <- sort(unique(cluster_mem))
    }
  }

  # Circular variable transformation
  # Transform angle to negative if the positive direction is clockwise
  # Because we use sin/cos, the R default is counter-clockwise
  if (cw) {
    data <-
      data %>%
      mutate(across(all_of(circ.var), ~ .x * (-1)))
  }
  # Make sure circ.var is in radian and positive
  if (is.degree) {
    data <-
      data %>%
      mutate(across(all_of(circ.var), ~ torad(.x)))
  } else {
    data <-
      data %>%
      mutate(across(all_of(circ.var), ~ .x %cr+% 0))
  }

  # Transform variables: circular to ellipse, numeric between 0 and 1
  data_transformed <-
    data %>%
    mutate(across(all_of(circ.var), ~ sin(.x %cr+% shift) / 2 + .5),
           across(all_of(numeric_var), ~ (.x - min(.x))/max(.x)),
           id = row_number(),
           member = cluster_mem)

  # Make a longer data form
  pcp_data <-
    data_transformed %>%
    tidyr::pivot_longer(all_of(all_var),
                        names_to = "var_name",
                        values_to = "value") %>%
    mutate(var_x = as.double(match(.data$var_name, order.appear)))

  # Override x values for circular variables
  pcp_data_2x <-
    data %>%
    select(all_of(circ.var)) %>%
    mutate(across(everything(), ~ cos(.x %cr+% shift) / 4))

  for (i in circ.var) {
    pcp_data[pcp_data$var_name == i, "var_x"] <-
      which(order.appear == i) + pcp_data_2x[[i]]
  }

  # Create values for annotation on the plot
  mins <- data %>% summarize(across(everything(), min))
  maxs <- data %>% summarize(across(everything(), max))

  # Create an ellipse shape
  ellipse <- seq(1, 360)
  ellipse_x <- cos(torad(ellipse)) / 4 + which(order.appear == circ.var)
  ellipse_y <- sin(torad(ellipse)) / 2 + 0.5
  ellipse_dat <- tibble::tibble(var_x = ellipse_x, value = ellipse_y)

  pos <- seq_len(length(order.appear))

  # Create vectors of position for NEWS label on circular variable
  put_degrees <- c(0, pi / 2, pi, 3 * pi / 2) * ifelse(cw, -1, 1) + shift
  news_x <- cos(put_degrees) / 4 + which(order.appear == circ.var)
  news_y <- 0.05 + sin(put_degrees) / 2 + 0.5

  # If medians plot are TRUE, show them in bold and no transparent
  medoid_size <- if_else(show.medoids, 1, size)
  medoid_alpha <- if_else(show.medoids, 1, alpha)

  medoids <- object$medoids

  # Make the plot
  p <-
    ggplot2::ggplot(pcp_data, ggplot2::aes(x = .data$var_x,
                                           y = .data$value,
                                           group = .data$id,
                                           col = as.character(.data$member))) +
    ggplot2::geom_line(alpha = alpha, size = size, linetype = linetype) +
    ggplot2::geom_line(data = filter(pcp_data, .data$id %in% medoids),
                       alpha = medoid_alpha,
                       size = medoid_size) +
    ggplot2::scale_x_continuous(breaks = seq(length(order.appear)),
                                labels = order.appear) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 1.1)) +
    ggplot2::scale_color_manual(values = cluster.col) +
    #scale_color_discrete(#labels = as.character(unique(pcp_data$groups)),
    #  name = "Cluster") +
    ggplot2::annotate("text", x = pos[-which(order.appear == circ.var)],
                      y = -0.05,
                      label = as.character(
                        mins[order.appear[-which(order.appear == circ.var)]]),
                      size = labelsize) +
    ggplot2::annotate("text", x = pos[-which(order.appear == circ.var)],
                      y = 1.05,
                      label = as.character(
                        maxs[order.appear[-which(order.appear == circ.var)]]),
                      size = labelsize) +
    ggplot2::annotate("text", x = news_x, news_y, label = c("N", "E", "S", "W"),
                      size = labelsize) +
    # make an ellipse shape
    ggplot2::geom_path(data = ellipse_dat,
                       ggplot2::aes(x = .data$var_x,
                                    y = .data$value,
                                    group = NULL,
                                    col = NULL),
                       linetype = 2) +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  color = legend.cluster) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())

  return(p)
}

#' Transform Degree to Radian
#'
#' This function transforms a circular angle from degree to radian.
#'
#' @param x A degree value if `torad` or radian value if `todeg`.
#'
#' @return A radian value if `torad` or degree value if `todeg`.
#' @name to_deg_rad
#'
#' @examples
#' torad(90)
#'
#' torad(-45)
#'
#' todeg(pi/2)
NULL

#' @export
#' @rdname to_deg_rad
torad <- function(x) {
  return(((x / 180) * pi) %cr+% 0)
}

#' @export
#' @rdname to_deg_rad
todeg <- function(x) {
  return(((x / pi) * 180) %cd+% 0)
}
