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
#' @param cluster.sol A `MonoClust` object which contains the cluster results
#'   and memberships.
#' @param cols Color of clusters, indicating by a vector. The length of this
#'   vector must be equal to the number of clusters in the cluster.sol object.
#' @param is.degree Whether the unit of the circular variables is degree or not
#'   (radian). Default is `TRUE`.
#' @param shift The shift (offset, rotation) of the circular variable, in
#'   radians. Default is 0 (no rotation).
#' @param alpha The transparency of the lines. Default to be 0.1.
#' @param show.medoids Whether to highlight the median lines or not. Default is
#'   `TRUE`.
#' @param north What value of the circular variable is labeled North. Default is
#'   0 radian.
#' @param cw Which direction of the circular variable is considered increasing
#'   in value, clockwise (`TRUE`) or counter-clockwise (`FALSE`). Default is
#'   `TRUE`.
#' @param labelsize The size of labels on the plot, default is 4.
#'
#' @return A ggplot2 object
#' @export
#' @importFrom dplyr `%>%`
#' @importFrom rlang .data
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
#' pcp.gg(wind_reduced,
#'        "WDIR",
#'        c("WDIR", "has.sensit", "WS"),
#'        sol42007,
#'        COLOR4,
#'        shift = pi*3/4-0.3,
#'        alpha = 0.1,
#'        labelsize = 8) +
#'   theme(panel.background = element_rect(color = NA),
#'         plot.background = element_rect(color = NA),
#'         panel.border = element_rect(color = NA),
#'         panel.grid.major = element_line(color = "#f0f0f0"),
#'         panel.grid.minor = element_blank(),
#'         axis.line = element_line(color = "black"),
#'         legend.key = element_rect(color = NA),
#'         legend.position = "bottom",
#'         legend.direction = "horizontal",
#'         legend.title = element_text(face = "italic"),
#'         legend.justification = "center",
#'         text = element_text(size = 24))
pcp.gg <- function(data, circ.var, order.appear, cluster.sol,
                   cols, is.degree = TRUE, shift = 0, alpha = 0.1,
                   show.medoids = TRUE, north = 0, cw = FALSE,
                   labelsize = 4) {
  # Constants, may change into parameters later
  LT <- 1.5
  REG.SIZE <- 0.5
  # If medians plot are TRUE, show them in bold and no transparent
  MED.SIZE <- ifelse(show.medoids, 1, REG.SIZE)
  MED.ALPHA <- ifelse(show.medoids, 1, alpha)

  medoids <- cluster.sol$medoids

  # Get variables from data set
  variables <- colnames(data)

  # Default values for parameters
  if (is.null(order.appear))
    order.appear <- variables

  # Check if order.appear includes variables in the data
  if (!all(order.appear %in% variables))
    stop("The order.appear has to be a subset of data variables.")

  # Take cluster memberships out from cluster solution
  cluster_mem <- factor(cluster.sol$Membership)

  if (length(cols) != length(unique(cluster_mem))) {
    warning("The number of colors does not match the number of clusters.
            Default values will be used.")
    cols <- seq_len(length(unique(cluster_mem)))
  }

  # Transform angle to negative if the positive direction is clockwise
  # Because we use sin/cos, the R default is counter-clockwise
  if (cw) data[, circ.var] <- data[, circ.var] * (-1)

  # Create a data set with standardized variables
  pcp_data_1 <- data %>%
    dplyr::mutate_all(list(~ (.data - min(.data)))) %>%
    dplyr::mutate_all(list(~ (.data / max(.data)))) %>%
    dplyr::mutate_all(id = seq_len(nrow(.data)), groups = cluster_mem)

  # Create y values for circular variables override y values
  if (!is.null(circ.var)) {
    if (is.degree)
      data[, circ.var] <- torad(data[, circ.var])
    pcp_data_2y <- data[, circ.var] %>%
      dplyr::mutate_all(list(~sin(.data %cr+% shift) / 2 + .5))

    # Combine the value into the plotting data set
    pcp_data_1[, circ.var] <- pcp_data_2y
  }

  pcp_data <- pcp_data_1 %>%
    tidyr::gather(key = .data$var.name, value = .data$value, !! order.appear,
                  factor_key = TRUE) %>%
    dplyr::mutate_all(var.x = match(.data$var.name, order.appear)) %>%
    dplyr::mutate_all(colors = cols[.data$groups])
  summary(pcp_data)

  # Override x values for circular variables
  if (!is.null(circ.var)) {
    pcp_data_2x <- data[, circ.var] %>%
      dplyr::mutate_all(list(~cos(.data %cr+% shift) / 4))

    # Combine the value into the plotting data set
    pcp_data[pcp_data$var.name == circ.var, "var.x"] <-
      pcp_data[pcp_data$var.name == circ.var, "var.x"] +
      pcp_data_2x
  }

  # Create values for annotation on the plot
  mins <- data %>% dplyr::summarize_all(min)
  maxs <- data %>% dplyr::summarize_all(max)

  # Create an ellipse shape
  ellipse <- seq(1, 360)
  ellipse_x <- cos(torad(ellipse)) / 4 + which(order.appear == circ.var)
  ellipse_y <- sin(torad(ellipse)) / 2 + 0.5
  ellipse_dat <- tibble::tibble(var.x = ellipse_x, value = ellipse_y)

  pos <- seq_len(length(order.appear))

  # Create vectors of position for NEWS label on circular variable
  put_degrees <- c(0, pi / 2, pi, 3 * pi / 2) * ifelse(cw, -1, 1) + shift
  news_x <- cos(put_degrees) / 4 + which(order.appear == circ.var)
  news_y <- 0.05 + sin(put_degrees) / 2 + 0.5

  # Make the plot
  ggplot2::ggplot(pcp_data, ggplot2::aes(x = .data$var.x,
                                         y = .data$value,
                                         group = .data$id,
                                         col = .data$groups)) +
    ggplot2::geom_line(alpha = alpha, size = REG.SIZE, linetype = LT) +
    ggplot2::geom_line(data = dplyr::filter(pcp_data, .data$id %in% medoids),
                       alpha = MED.ALPHA,
                       size = MED.SIZE) +
    ggplot2::scale_x_continuous(breaks = seq(length(order.appear)),
                                labels = order.appear,
                                name = "Variables") +
    ggplot2::scale_y_continuous(name = "",
                                limits = c(-0.1, 1.1)) +
    ggplot2::scale_color_manual(values = cols) +
    #scale_color_discrete(#labels = as.character(unique(pcp_data$groups)),
    #  name = "Cluster") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank()) +
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
                       ggplot2::aes(x = .data$var.x,
                                    y = .data$value,
                                    group = NULL,
                                    col = NULL),
                       linetype = 2)
}

#' Transform Degree to Radian
#'
#' @param x A degree value.
#'
#' @return A radian value.
#' @export
#'
#' @examples
#' torad(90)
#'
#' torad(45)
torad <- function(x) {
  (x / 180) * pi
}
