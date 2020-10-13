#' Plot a parallel coordinate plot with one circular variable
#'
#' @param data The full data set
#' @param cir.var The index or name of the circular variable in the data set
#' @param membership A vector of cluster membership for all observations in the
#'   data set
#' @param order The order of appearance of the variables in the plot
#' @param col List of color for the membership. It should match the number of
#'   groups in membership parameter
#' @param zero The position of the zero degree value in the circular variable.
#'   Can be chosen in c('north', 'east', 'west', 'south'). Default is 'east'
#'   side of the circle (usual 0 position in calculus circle)
#' @param nslabel Whether to print the North and South label on the circular
#'   variable. Default is TRUE
#' @param rotation The positive direction of the circular variable. Default is
#'   'counter'-clockwise. Can be changed to 'clock'.
#' @param meds A list of median values of the cluster. The length of this vector
#'   should match the unique values in membership.
#'
#' @return A PCP plot with an ellipse for circular variable
#' @export
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics axis points segments text
#'
#' @examples
#' library(cluster)
pcp.circ.plot <- function(data, cir.var = NULL, membership = NULL, order = NULL,
                          col = NULL, zero = "east", nslabel = TRUE,
                          rotation = "counter", meds = NULL) {

    direction <- c("east", "north", "west", "south")
    shift <- c(0, 90, 180, 270)
    if (rotation == "clock") {
        data[, cir.var] <- data[, cir.var] * (-1)
    }

    if (!(zero %in% direction)) {
        cat("Zero value needs to go with either \"east\", \"north\", \"west\", or \"south\". Default is \"east\".")
        return(0)
    } else {
        data[, cir.var] <- data[, cir.var] %c+% shift[which(direction == zero)]
    }

    # Circular variable
    if (!is.null(cir.var)) {
        cirvar <- torad(data[, cir.var])
    }

    linearvar_y <- apply(data, 2, function(x) (x - min(x))/max(x - min(x)))

    if (is.null(order)) {
        order <- 1:ncol(data)
    }

    if (!is.null(col)) {
        col1 <- col[membership]
    }

    xprime <- seq(0, 7, length.out = length(order))
    x <- xprime
    x[order] <- x
    # data[,cir.var] <- seg_dat <- data.frame(cbind((cos(aspect)/2), (sin(aspect)/2 + .5), elev_x,
    # elev_y)) seg_dat1 <- data.frame(cbind(elev_x, elev_y, density_x, density_y))

    ########################################### work on circular plot set up data for making plot
    theta <- seq(0, 360, 0.01)
    x1 <- cos(pi * theta/180)
    y1 <- sin(pi * theta/180)
    y1_stand <- y1/2 + 0.5
    x1_stand <- x1/2 + x[cir.var]

    plot(y1_stand ~ x1_stand, type = "l", xlim = c(-0.5, 7), ylim = c(0, 1), xaxt = "none", xlab = "",
        yaxt = "none", ylab = "", main = "PCP plot - circular")
    labs <- colnames(data)
    axis(1, at = x, labels = labs, las = 1)

    points(cos(cirvar)/2 + x[cir.var], sin(cirvar)/2 + 0.5, cex = 0.75, col = membership - 3)

    linear.l <- 1:ncol(data)
    linear.l <- linear.l[-cir.var]
    for (i in linear.l) {
        points(rep(x[i], length(data[, i])), linearvar_y[, i], col = membership - 3)
    }

    for (i in 1:(ncol(data) - 1)) {
        if (order[i] == cir.var) {
            segments(cos(cirvar)/2 + x[cir.var], sin(cirvar)/2 + 0.5, rep(xprime[i + 1], length(data[,
                i])), linearvar_y[, order[i + 1]], col = adjustcolor(col = membership - 3, alpha.f = 0.2))
            if (!is.null(meds)) {
                segments(cos(cirvar[meds])/2 + x[cir.var], sin(cirvar[meds])/2 + 0.5, rep(xprime[i +
                  1], length(data[, i])), linearvar_y[meds, order[i + 1]], col = as.numeric(names(meds)) -
                  3, lwd = 2)
            }
        } else if (order[i + 1] == cir.var) {
            segments(rep(xprime[i], length(data[, i])), linearvar_y[, order[i]], cos(cirvar)/2 +
                x[cir.var], sin(cirvar)/2 + 0.5, col = adjustcolor(col = membership - 3, alpha.f = 0.2))
            if (!is.null(meds)) {
                segments(rep(xprime[i], length(data[, i])), linearvar_y[meds, order[i]], cos(cirvar[meds])/2 +
                  x[cir.var], sin(cirvar[meds])/2 + 0.5, col = as.numeric(names(meds)) - 3, lwd = 2)
            }
        } else {
            segments(rep(xprime[i], length(data[, i])), linearvar_y[, order[i]], rep(xprime[i +
                1], length(data[, i])), linearvar_y[, order[i + 1]], col = adjustcolor(col = membership -
                3, alpha.f = 0.2))
            if (!is.null(meds)) {
                segments(rep(xprime[i], length(data[, i])), linearvar_y[meds, order[i]], rep(xprime[i +
                  1], length(data[, i])), linearvar_y[meds, order[i + 1]], col = as.numeric(names(meds)) -
                  3, lwd = 2)
            }
        }
    }
    if (nslabel) {
        text(x[cir.var] - 0.45, 0, "S")
        text(x[cir.var] - 0.45, 1, "N")
    }
}
