orchard <- function (object, mod = "Int", xlab, N = "none", alpha = 0.5, 
                     angle = 90, cb = TRUE, k = TRUE, transfm = c("none", "tanh")) 
{
  transfm <- match.arg(transfm)
  if (any(class(object) %in% c("rma.mv", "rma"))) {
    if (mod != "Int") {
      object <- mod_results(object, mod)
    }
    else {
      object <- mod_results(object, mod = "Int")
    }
  }
  mod_table <- object$mod_table
  data <- object$data
  data$moderator <- factor(data$moderator, levels = mod_table$name, 
                           labels = mod_table$name)
  data$scale <- (1/sqrt(data[, "vi"]))
  legend <- "Precision (1/SE)"
  if (any(N != "none")) {
    data$scale <- N
    legend <- "Sample Size (N)"
  }
  if (transfm == "tanh") {
    cols <- sapply(mod_table, is.numeric)
    mod_table[, cols] <- Zr_to_r(mod_table[, cols])
    data$yi <- Zr_to_r(data$yi)
    label <- xlab
  }
  else {
    label <- xlab
  }
  mod_table$K <- as.vector(by(data, data[, "moderator"], function(x) length(x[, 
                                                                              "yi"])))
  group_no <- nrow(mod_table)
  cbpl <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
            "#CC79A7", "#56B4E9", "#999999")
  plot <- ggplot2::ggplot(data = mod_table, aes(x = estimate, 
                                                y = name)) + ggbeeswarm::geom_quasirandom(data = data, 
                                                                                          aes(x = yi, y = moderator, size = scale, colour = moderator), 
                                                                                          groupOnX = FALSE, alpha = alpha) + 
    #ggplot2::geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR), height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) + 
    ggplot2::geom_errorbarh(aes(xmin = lowerCL, 
                                                                                                                                                                                                xmax = upperCL), height = 0, show.legend = FALSE, size = 1.2) + 
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", 
                        alpha = alpha) + ggplot2::geom_point(aes(fill = name), 
                                                             size = 3, shape = 21) + ggplot2::theme_bw() + ggplot2::guides(fill = "none", 
                                                                                                                           colour = "none") + ggplot2::theme(legend.position = c(1, 
                                                                                                                                                                                 0), legend.justification = c(1, 0)) + ggplot2::theme(legend.title = element_text(size = 9)) + 
    ggplot2::theme(legend.direction = "horizontal") + ggplot2::theme(legend.background = element_blank()) + 
    ggplot2::labs(x = label, y = "", size = legend) + ggplot2::theme(axis.text.y = element_text(size = 10, 
                                                                                                colour = "black", hjust = 0.5, angle = angle))
  if (cb == TRUE) {
    plot <- plot + scale_fill_manual(values = cbpl) + scale_colour_manual(values = cbpl)
  }
  if (k == TRUE) {
    plot <- plot + ggplot2::annotate("text", x = (max(data$yi) + 
                                                    (max(data$yi) * 0.1)), y = (seq(1, group_no, 1)), label = paste("italic(k)==", mod_table$K), 
                                     parse = TRUE, vjust = "right", size = 3.5)
  }
  return(plot)
}
