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
#'   vector of names or index. If set, length has to be equal to the number of
#'   variables in the data set.
#' @param linetype Line type. Default is solid line. See details in
#'   `vignette("ggplot2-specs")`.
#' @param size Size of a line is its width in mm. Default is 0.5. See details in
#'   `vignette("ggplot2-specs")`.
#' @param cluster.col Color of clusters, indicating by a vector. If set, the
#'   length of this vector must be equal to the number of clusters in
#'   `clustering`.
#' @param is.degree Whether the unit of the circular variables is degree or not
#'   (radian). Default is `TRUE`.
#' @param rotate The rotate (offset, shift) of the circular variable, in
#'   radians. Default is 0 (no rotation).
#' @param alpha The transparency of the lines. Default is 0.1.
#' @param show.medoids Whether to highlight the median lines or not. Default is
#'   `FALSE`.
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
#' @param clustering Cluster membership.
#' @param medoids Vector of medoid observations of cluster. Only required when
#' `show.medoids = TRUE`.
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
#' ggpcp(data = wind_reduced,
#'       circ.var = "WDIR",
#'       # To improve aesthetics
#'       rotate = pi*3/4-0.3,
#'       order.appear = c("WDIR", "has.sensit", "WS"),
#'       alpha = 0.5,
#'       clustering = sol42007$membership,
#'       medoids = sol42007$medoids,
#'       cluster.col = COLOR4,
#'       show.medoids = TRUE) +
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
ggpcp <- function(data, circ.var = NULL, is.degree = TRUE, rotate = 0,
                   north = 0, cw = FALSE, order.appear = NULL, linetype = 1,
                   size = 0.5, alpha = 0.5, clustering, medoids = NULL,
                   cluster.col = NULL, show.medoids = FALSE, labelsize = 4,
                   xlab = "Variables", ylab = NULL, legend.cluster = "groups") {

  if (missing(data) || !is.data.frame(data))
    stop("\"data\" must be a data frame.")
  if (!all(purrr::map_lgl(data, is.numeric)))
    stop("Function only supports numerical variables.")

  if (missing(clustering))
    stop("\"clustering\" is required.")
  if (length(clustering) != nrow(data))
    stop("\"clustering\" must have the same length as \"data\" observations.")

  if (show.medoids && is.null(medoids))
    stop("\"medoids\" must be set when \"show.medoids = TRUE\"")

  # clustering
  clustering <- as.character(clustering)

  # Commonly used variables
  all_var <- colnames(data)

  # circ.var
  if (!is.null(circ.var)) {
    # Check if circ.var is an index
    if (is.numeric(circ.var))
      circ.var <- all_var[circ.var]
    # Check if circ.var is a valid column
    if (!all(circ.var %in% all_var))
      stop("\"circ.var\" must be columns in \"data\".")
  }
  # rotate
  if (length(rotate) == 1) {
    rotate <- rep(rotate, length(circ.var))
  } else if (length(rotate) > length(circ.var))
    rotate <- rotate[seq_len(circ.var)]
  # Commonly used variables
  numeric_var <- setdiff(all_var, circ.var)
  # order.appear
  if (!is.null(order.appear)) {
    # Check if order.appear is a vector of index
    if (is.numeric(order.appear))
      order.appear <- all_var[order.appear]
    # Check if order.appear includes variables in the data
    if (!all(order.appear %in% all_var))
      stop("\"order.appear\" has to be a subset of data variables.")
    order.appear <- c(order.appear, setdiff(all_var, order.appear))
  } else {
    order.appear <- all_var
  }
  # cluster.col
  if (is.null(cluster.col)) {
    cluster.col <- seq_len(length(unique(clustering)))
  } else {
    if (length(cluster.col) != length(unique(clustering))) {
      message("The number of colors does not match the number of clusters.
            Default values will be used.")
      cluster.col <- seq_len(length(unique(clustering)))
    }

    if (is.null(names(cluster.col)) ||
        !identical(sort(names(cluster.col)),
                   sort(as.character(unique(clustering))))) {
      if (!is.null(names(cluster.col)))
        message("Named \"cluster.col\" does not match cluster names. Default
              values will be used.")
      names(cluster.col) <- sort(unique(clustering))
    }
  }

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

  circ_var_tbl <-
    tibble(circ.var,
           rotate)

  pcp_data <-
    data %>%
    # Transform numericals between 0 and 1
    mutate(across(all_of(numeric_var), ~ (.x - min(.x))),
           across(all_of(numeric_var), ~ (.x / max(.x))),
           # Add ID for each observation
           id = row_number(),
           # Add cluster membership
           member = clustering) %>%
    # Longer version
    tidyr::pivot_longer(all_of(all_var),
                        names_to = "var_name",
                        values_to = "value") %>%
    # Join to match rotate, assign 0 if not circular
    left_join(circ_var_tbl, by = c("var_name" = "circ.var")) %>%
    tidyr::replace_na(list(rotate = 0)) %>%
    # Rotate circular and transform to ellipse
    # Okay to use for all variable because numericals are now between 0 and 1
    mutate(value_y = if_else(.data$var_name %in% circ.var,
                             sin(.data$value %cr+% .data$rotate) / 2 + .5,
                             .data$value),
           value_x = if_else(.data$var_name %in% circ.var,
                             cos(.data$value %cr+% .data$rotate) / 4 +
                               as.double(match(.data$var_name, order.appear)),
                             as.double(match(.data$var_name, order.appear))))

  # Create values for annotation on the plot
  mins <- data %>% summarize(across(everything(), min))
  maxs <- data %>% summarize(across(everything(), max))

  # Create an ellipse shape
  ellipse <- seq(1, 360)
  ellipse_x <- cos(torad(ellipse)) / 4
  ellipse_y <- sin(torad(ellipse)) / 2 + 0.5
  ellipse_dat <-
    tibble::tibble(value = rep(ellipse_x, times = length(circ.var)),
                   value_y = rep(ellipse_y, times = length(circ.var)),
                   var_name = rep(circ.var, each = length(ellipse))) %>%
    mutate(value_x = .data$value + match(.data$var_name, order.appear))

  pos <- seq_len(length(order.appear))

  # Create vectors of position for NEWS label on circular variable
  put_degrees <-
    rep(c(0, pi / 2, pi, 3 * pi / 2) * ifelse(cw, -1, 1),
        times = length(circ.var)) %cr+% rep(rotate, each = 4)
  news_x <- cos(put_degrees) / 4 + match(rep(circ.var, each = 4), order.appear)
  news_y <- 0.05 + sin(put_degrees) / 2 + 0.5

  # If medians plot are TRUE, show them in bold and no transparent
  medoid_size <- if_else(show.medoids, 1, size)
  medoid_alpha <- if_else(show.medoids, 1, alpha)

  # Make the plot
  p <-
    ggplot2::ggplot(pcp_data, ggplot2::aes(x = .data$value_x,
                                           y = .data$value_y,
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
                        mins[order.appear[-match(circ.var, order.appear)]]),
                      size = labelsize) +
    ggplot2::annotate("text", x = pos[-match(circ.var, order.appear)],
                      y = 1.05,
                      label = as.character(
                        maxs[order.appear[-which(order.appear == circ.var)]]),
                      size = labelsize) +
    ggplot2::annotate("text", x = news_x, y = news_y,
                      label = rep(c("N", "E", "S", "W"),
                                  times = length(circ.var)),
                      size = labelsize) +
    # make an ellipse shape
    ggplot2::geom_path(data = ellipse_dat,
                       ggplot2::aes(x = .data$value_x,
                                    y = .data$value_y,
                                    group = .data$var_name,
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
