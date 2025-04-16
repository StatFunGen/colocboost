#' @rdname colocboost_plot
#'
#' @title Plot visualization plot from a ColocBoost output.
#'
#' @description `colocboost_plot` generates visualization plots for colocalization events from a ColocBoost analysis.
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param y Specifies the y-axis values, default is "log10p" for -log10 transformed marginal association p-values.
#' @param pos Optional plotting range of x-axis to zoom in to a specific region.
#' @param plot_target_only Logical, if TRUE only plots colocalization with target outcome, default is FALSE.
#' @param plot_cos_idx Optional indices of CoS to plot
#' @param outcome_idx Optional indices of outcomes to include in the plot. \code{outcome_idx=NULL} to plot only the outcomes having colocalization.
#' @param points_color Background color for non-colocalized variables, default is "grey80".
#' @param cos_color Optional custom colors for CoS.
#' @param add_vertical Logical, if TRUE adds vertical lines at specified positions, default is FALSE
#' @param add_vertical_idx Optional indices for vertical lines.
#' @param outcome_names Optional vector of outcomes names for the subtitle of each figure. \code{outcome_names=NULL} for the outcome name shown in \code{data_info}.
#' @param plot_cols Number of columns in the plot grid, default is 2. If you have many colocalization. please consider increasing this.
#' @param variant_coord Logical, if TRUE uses variant coordinates on x-axis, default is FALSE. This is required the variable names including position information.
#' @param show_hits Logical, if TRUE shows top variables for each CoS, default is FALSE
#' @param show_cos_to_uncoloc Logical, if TRUE shows colocalization to uncolocalized outcomes to diagnose, default is FALSE
#' @param show_cos_to_uncoloc_idx Optional indices for showing CoS to all uncolocalized outcomes
#' @param show_cos_to_uncoloc_outcome Optional outcomes for showing CoS to uncolocalized outcomes
#' @param gene_name Optional gene name to display in plot title
#' @param ylim_each Logical, if TRUE uses separate y-axis limits for each plot, default is TRUE
#' @param outcome_legend_pos Position for outcome legend, default is "top"
#' @param outcome_legend_size Size for outcome legend text, default is 1.2
#' @param cos_legend_pos Position for colocalization set legend, default is "bottomleft"
#' @param show_variable Logical, if TRUE displays variant IDs, default is FALSE
#' @param lab_style Vector of two numbers for label style (size, boldness), default is c(2, 1)
#' @param axis_style Vector of two numbers for axis style (size, boldness), default is c(2, 1)
#' @param title_style Vector of two numbers for title style (size, boldness), default is c(2.5, 2)
#' @param ... Additional parameters passed to `plot` functions
#'
#' @return Visualization plot for each colcoalization event.
#'
#' @importFrom utils head tail
#' @importFrom graphics abline axis legend mtext par points text
#' @importFrom grDevices adjustcolor
#'
#' @keywords cb_plot
#' @export
colocboost_plot <- function(cb_output, y = "log10p",
                            pos = NULL,
                            plot_target_only = FALSE,
                            plot_cos_idx = NULL,
                            outcome_idx = NULL,
                            points_color = "grey80",
                            cos_color = NULL,
                            add_vertical = FALSE,
                            add_vertical_idx = NULL,
                            outcome_names = NULL,
                            plot_cols = 2,
                            variant_coord = FALSE,
                            show_hits = FALSE,
                            show_cos_to_uncoloc = FALSE,
                            show_cos_to_uncoloc_idx = NULL,
                            show_cos_to_uncoloc_outcome = NULL,
                            gene_name = NULL,
                            ylim_each = TRUE,
                            outcome_legend_pos = "top",
                            outcome_legend_size = 1.2,
                            cos_legend_pos = "bottomleft",
                            show_variable = FALSE,
                            lab_style = c(2, 1),
                            axis_style = c(2, 1),
                            title_style = c(2.5, 2),
                            ...) {
  if (!inherits(cb_output, "colocboost")) {
    stop("Input of colocboost_plot must be a 'colocboost' object!")
  }

  # get cb_plot_input data from colocboost results
  cb_plot_input <- get_input_plot(cb_output,
    plot_cos_idx = plot_cos_idx,
    plot_target_only = plot_target_only,
    variant_coord = variant_coord,
    outcome_names = outcome_names,
    show_cos_to_uncoloc = show_cos_to_uncoloc,
    show_cos_to_uncoloc_idx = show_cos_to_uncoloc_idx,
    show_cos_to_uncoloc_outcome = show_cos_to_uncoloc_outcome
  )
  # get initial set up of plot
  cb_plot_init <- plot_initial(cb_plot_input,
    y = y, points_color = points_color, cos_color = cos_color,
    ylim_each = ylim_each, gene_name = gene_name,
    outcome_legend_pos = outcome_legend_pos, outcome_legend_size = outcome_legend_size,
    cos_legend_pos = cos_legend_pos,
    show_variable = show_variable, lab_style = lab_style, axis_style = axis_style,
    title_style = title_style, ...
  )

  colocboost_plot_basic <- function(cb_plot_input, cb_plot_init,
                                    outcome_idx = NULL, pos = NULL,
                                    plot_cols = 2,
                                    add_vertical = FALSE, add_vertical_idx = NULL,
                                    show_hits = TRUE,
                                    ...) {
    args <- list(...)
    args <- c(args, cb_plot_init[c("xlab", "ylab")])
    args$col <- cb_plot_init$bg
    if (is.null(pos)) {
      args$x <- cb_plot_init$x
      y <- cb_plot_init$y
    } else {
      args$x <- cb_plot_init$x[pos]
      y <- lapply(cb_plot_init$y, function(yy) yy[pos])
    }
    args$pch <- cb_plot_init$pch
    args$cex.axis <- cb_plot_init$axis_size
    args$cex.lab <- cb_plot_init$lab_size
    args$font.lab <- cb_plot_init$lab_face
    # - change position
    cb_plot_init$outcome_legend_pos <- switch(cb_plot_init$outcome_legend_pos,
      "right" = 4,
      "left" = 2,
      "top" = 3,
      "bottom" = 1
    )

    # - begin plotting
    coloc_cos <- cb_plot_input$cos
    outcomes <- cb_plot_input$outcomes
    if (is.null(outcome_idx)) {
      if (is.null(coloc_cos)) {
        # - no colocalized effects, draw all outcomes in this region
        if (length(cb_plot_input$outcomes) == 1) {
          message("There is no fine-mapped causal effect in this region!. Showing margianl for this outcome!")
        } else {
          message("There is no colocalization in this region!. Showing margianl for all outcomes!")
        }
        outcome_idx <- 1:length(y)
      } else {
        n.coloc <- length(coloc_cos)
        coloc_index <- cb_plot_input$coloc_index
        outcome_idx <- Reduce(union, coloc_index)
      }
      if (!is.null(cb_plot_input$target_outcome)) {
        p_target <- grep(cb_plot_input$target_outcome, outcomes)
        include_target <- sapply(cb_plot_input$coloc_index, function(ci) {
          p_target %in% ci
        })
        if (any(include_target)) {
          coloc_index <- cb_plot_input$coloc_index[order(include_target == "FALSE")]
          coloc_index <- Reduce(union, coloc_index)
          outcome_idx <- c(p_target, setdiff(coloc_index, p_target))
        }
      }
    }
    if (length(outcome_idx) == 1) {
      plot_cols <- 1
    }
    nrow <- ceiling(length(outcome_idx) / plot_cols)
    if (!is.null(cb_plot_init$xtext)) {
      bottom <- 6
    } else {
      bottom <- 2
    }
    if (!is.null(cb_plot_init$title)) {
      par(mfrow = c(nrow, plot_cols), mar = c(bottom, 5, 2, 1), oma = c(0, 0, 3, 0))
    } else {
      par(mfrow = c(nrow, plot_cols), mar = c(bottom, 5, 2, 1), oma = c(0, 0, 1, 0))
    }
    for (iy in outcome_idx) {
      args$y <- y[[iy]]
      args$ylim <- c(cb_plot_init$ymin[iy], cb_plot_init$ymax[iy])
      if (!is.null(cb_plot_init$xtext)) {
        args$xaxt <- "n"
        do.call(plot, args)
        axis(1, at = args$x, labels = FALSE)
        text(x = args$x, y = par("usr")[3] - 0.1, labels = cb_plot_init$xtext, srt = 45, adj = 1, xpd = TRUE)
      } else {
        do.call(plot, args)
      }
      mtext(outcomes[iy],
        side = cb_plot_init$outcome_legend_pos, line = 0.2, adj = 0.5,
        cex = cb_plot_init$outcome_legend_size, font = 1
      )
      if (add_vertical) {
        for (iii in 1:length(add_vertical_idx)) {
          abline(v = add_vertical_idx[iii], col = "#E31A1C", lwd = 1.5, lty = "dashed")
        }
      }

      # mark variables in CoS to colocalized outcomes
      if (!is.null(coloc_cos)) {
        n.coloc <- length(coloc_cos)
        coloc_index <- cb_plot_input$coloc_index
        legend_text <- list(col = vector())
        legend_text$col <- head(cb_plot_init$col, n.coloc)

        # check which coloc set for this outcome
        p.coloc <- sapply(coloc_index, function(idx) sum(idx == iy) != 0)
        p.coloc <- which(p.coloc)
        coloc_cos.idx <- coloc_cos[p.coloc]
        for (i.cs in p.coloc) {
          # add the points with specific color
          cs <- as.numeric(coloc_cos[[i.cs]])
          x0 <- intersect(args$x, cs)
          y1 <- args$y[match(x0, args$x)]
          points(x0, y1, pch = 21, bg = legend_text$col[i.cs], col = NA, cex = 2.5, lwd = 2.5)
          if (show_hits) {
            # add the hits points with "red"
            cs_hits <- as.numeric(cb_plot_input$cos_hits[[i.cs]])
            x_hits <- intersect(args$x, cs_hits)
            y_hits <- args$y[match(x_hits, args$x)]
            points(x_hits, y_hits, pch = 21, bg = legend_text$col[i.cs], col = "#E31A1C", cex = 3, lwd = 3)
          }
        }

        # mark variables in CoS to uncolocalized outcomes
        uncoloc <- cb_plot_input$uncoloc
        if (!is.null(uncoloc)) {
          p.uncoloc <- sapply(uncoloc$outcome_to_uncoloc, function(idx) sum(idx == iy) != 0)
          p.uncoloc <- which(p.uncoloc)
          texts <- shape_col <- texts_col <- c()
          for (i.uncoloc in p.uncoloc) {
            uncoloc_outcome <- uncoloc$outcome_to_uncoloc[[i.uncoloc]]
            if (iy %in% uncoloc_outcome) {
              # add the points with specific color
              cs <- as.numeric(uncoloc$cos_to_uncoloc[[i.uncoloc]])
              x0 <- intersect(args$x, cs)
              y1 <- args$y[match(x0, args$x)]
              points(x0, y1,
                pch = 4, col = adjustcolor(legend_text$col[i.uncoloc], alpha.f = 0.3),
                cex = 1.5, lwd = 1.5
              )
              texts <- c(texts, uncoloc$cos_uncoloc_texts[i.cs])
              shape_col <- c(shape_col, adjustcolor(legend_text$col[i.uncoloc], alpha.f = 1))
              texts_col <- c(texts_col, adjustcolor(legend_text$col[i.uncoloc], alpha.f = 0.8))
            }
          }
          if (length(texts) == 0) {
            next
          }
          legend(cb_plot_init$cos_legend_pos, texts,
            bty = "n", col = shape_col, text.col = texts_col,
            cex = 1.5, pt.cex = 1.5, pch = 4, x.intersp = 0.1, y.intersp = 0.3
          )
        }
      }
    }
    if (!is.null(cb_plot_init$title)) {
      mtext(cb_plot_init$title,
        side = 3, line = 0, outer = TRUE,
        cex = cb_plot_init$title_size, font = cb_plot_init$title_face
      )
    }
    return(invisible())
  }

  colocboost_plot_basic(cb_plot_input, cb_plot_init,
    pos = pos,
    outcome_idx = outcome_idx, plot_cols = plot_cols,
    add_vertical = add_vertical, add_vertical_idx = add_vertical_idx,
    show_hits = show_hits,
    ...
  )
}



# get input data for cb_plot
get_input_plot <- function(cb_output, plot_cos_idx = NULL,
                           plot_target_only = FALSE,
                           variant_coord = FALSE,
                           outcome_names = NULL,
                           show_cos_to_uncoloc = FALSE,
                           show_cos_to_uncoloc_idx = NULL,
                           show_cos_to_uncoloc_outcome = NULL) {
  # redefined outcome names
  if (!is.null(outcome_names)) {
    cb_output$data_info$outcome_info$outcome_names <- outcome_names
    names(cb_output$data_info$z) <- outcome_names
  }

  # extract results from colocboost
  analysis_outcome <- cb_output$data_info$outcome_info$outcome_names
  target_outcome_idx <- which(cb_output$data_info$outcome_info$is_target)
  if (length(target_outcome_idx) != 0) {
    target_outcome <- analysis_outcome[target_outcome_idx]
  } else {
    target_outcome <- NULL
  }
  # check if target cos
  target_cos <- cb_output$cos_summary$cos_id[cb_output$cos_summary$target_outcome != FALSE]
  if_target <- !is.na(match(names(cb_output$cos_details$cos$cos_variables), target_cos))
  # extract z-scores
  variables <- cb_output$data_info$variables
  Z <- cb_output$data_info$z
  coef <- cb_output$data_info$coef

  # if finemapping
  if (cb_output$data_info$n_outcomes == 1) {
    cb_output$cos_details$cos$cos_variables <- cb_output$ucos_details$ucos$ucos_variables
    cb_output$cos_details$cos$cos_index <- cb_output$ucos_details$ucos$ucos_index
    cb_output$cos_details$cos_top_variables <- cb_output$ucos_details$ucos_top_variables
    cb_output$cos_details$cos_outcomes <- cb_output$ucos_details$ucos_outcomes
    if (!is.null(cb_output$cos_details$cos$cos_variables)) {
      cb_output$cos_details$cos_vcp <- cb_output$ucos_details$ucos_weight
    }
    if_target <- rep(TRUE, length(cb_output$cos_details$cos$cos_variables))
  }
  # extract coloc_cos
  coloc_variables <- cb_output$cos_details$cos$cos_variables
  coloc_cos <- cb_output$cos_details$cos$cos_index
  coloc_index <- cb_output$cos_details$cos_outcomes$outcome_index
  # top_variables
  coloc_hits <- lapply(names(coloc_cos), function(cn) {
    p <- grep(cn, rownames(cb_output$cos_details$cos_top_variables))
    cb_output$cos_details$cos_top_variables$top_index[p]
  })
  names(coloc_hits) <- names(coloc_cos)
  if (!is.null(coloc_cos)) {
    # extract vcp each outcome
    vcp <- lapply(1:length(analysis_outcome), function(iy) {
      pos <- which(sapply(coloc_index, function(idx) iy %in% idx))
      if (length(pos) != 0) {
        w <- do.call(cbind, cb_output$cos_details$cos_vcp[pos])
        return(1 - apply(1 - w, 1, prod))
      } else {
        return(rep(0, length(variables)))
      }
    })
    ncos <- length(cb_output$cos_details$cos$cos_index)

    select_cs <- 1:ncos
    if (plot_target_only) {
      if (sum(if_target) == 0) {
        message("No target CoS, draw all CoS.")
      } else {
        select_cs <- which(if_target)
      }
    } else {
      if (!is.null(plot_cos_idx)) {
        if (length(setdiff(plot_cos_idx, c(1:ncos))) != 0) {
          stop("please check plot_cos_idx!")
        }
        select_cs <- plot_cos_idx
      }
    }
    coloc_variables <- coloc_variables[select_cs]
    coloc_cos <- coloc_cos[select_cs]
    coloc_index <- coloc_index[select_cs]
    coloc_hits <- coloc_hits[select_cs]
  } else {
    if (cb_output$data_info$n_outcomes == 1) {
      warnings("No fine-mapped causal effects in this region!")
    } else {
      warnings("No colocalized effects in this region!")
    }
    coloc_index <- NULL
    vcp <- rep(0, cb_output$data_info$n_variables)
  }
  # extract x axis
  if (variant_coord) {
    variable.info <- do.call(rbind, lapply(variables, function(variable) strsplit(variable, ":")[[1]]))
    chrom <- as.numeric(sapply(variable.info[, 1], function(ss) strsplit(ss, "chr")[[1]][2]))
    x <- data.frame(
      "chrom" = chrom,
      "pos" = as.numeric(variable.info[, 2])
    )
    coloc_cos <- lapply(coloc_cos, function(cos) x$pos[cos])
    coloc_hits <- lapply(coloc_hits, function(cos) x$pos[cos])
  } else {
    x <- data.frame(
      "chrom" = NA,
      "pos" = c(1:length(variables))
    )
  }
  plot_input <- list(
    "outcomes" = analysis_outcome,
    "target_outcome" = target_outcome,
    "variables" = variables,
    "x" = x,
    "Zscores" = Z,
    "vcp" = vcp,
    "coef" = coef,
    "cos" = coloc_cos,
    "cos_hits" = coloc_hits,
    "coloc_index" = coloc_index
  )

  # check if plot cos to uncolocalized outcome
  if (show_cos_to_uncoloc & !is.null(coloc_cos)) {
    if (is.null(show_cos_to_uncoloc_idx)) {
      cos_to_uncoloc <- coloc_cos
      cos_idx_to_uncoloc <- 1:length(coloc_index)
      if (is.null(show_cos_to_uncoloc_outcome)) {
        warning("Show all CoSs to uncolocalized outcomes.")
        outcome_to_uncoloc <- lapply(coloc_index, function(cidx) {
          setdiff(1:length(analysis_outcome), cidx)
        })
      } else {
        warning("Show all CoSs to uncolocalized outcomes ", paste(show_cos_to_uncoloc_outcome, collapse = ","))
        outcome_to_uncoloc <- lapply(coloc_index, function(cidx) {
          setdiff(show_cos_to_uncoloc_outcome, cidx)
        })
      }
    } else {
      if (show_cos_to_uncoloc_idx > length(coloc_cos)) {
        warning(
          "There are only ", length(coloc_cos), " CoS in this region. ",
          "Cannot show the ordered ", paste(show_cos_to_uncoloc_idx, collapse = ","),
          " CoS, which does not exist"
        )
        outcome_to_uncoloc <- cos_to_uncoloc <- cos_idx_to_uncoloc <- NULL
      } else {
        cos_idx_to_uncoloc <- show_cos_to_uncoloc_idx
        cos_to_uncoloc <- coloc_cos[show_cos_to_uncoloc_idx]
        if (is.null(show_cos_to_uncoloc_outcome)) {
          warning(
            "Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","),
            " CoS for all uncolocalized outcomes."
          )
          outcome_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx) {
            l <- list(setdiff(1:length(analysis_outcome), coloc_index[[idx]]))
            names(l) <- names(coloc_index[idx])
            return(l)
          })
        } else {
          warning(
            "Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","),
            " CoS for outcomes ", paste(show_cos_to_uncoloc_outcome, collapse = ",")
          )
          outcome_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx) {
            l <- list(setdiff(show_cos_to_uncoloc_outcome, coloc_index[[idx]]))
            names(l) <- names(coloc_index[idx])
            return(l)
          })
        }
      }
    }
    if (!is.null(outcome_to_uncoloc)) {
      cos_uncoloc_texts <- rep("Uncolocalized effect", length(outcome_to_uncoloc))
    } else {
      cos_uncoloc_texts <- NULL
    }
    draw_uncoloc <- sapply(outcome_to_uncoloc, length) != 0
    if (any(draw_uncoloc)) {
      pos <- which(draw_uncoloc)
      outcome_to_uncoloc <- outcome_to_uncoloc[pos]
      cos_to_uncoloc <- cos_to_uncoloc[pos]
      cos_idx_to_uncoloc <- cos_idx_to_uncoloc[pos]
      cos_uncoloc_texts <- cos_uncoloc_texts[pos]
    }
    uncoloc <- list(
      "outcome_to_uncoloc" = outcome_to_uncoloc,
      "cos_to_uncoloc" = cos_to_uncoloc,
      "cos_idx_to_uncoloc" = cos_idx_to_uncoloc,
      "cos_uncoloc_texts" = cos_uncoloc_texts
    )
    plot_input$uncoloc <- uncoloc
  }
  class(plot_input) <- "colocboost"
  return(plot_input)
}


#' @importFrom stats pnorm
plot_initial <- function(cb_plot_input, y = "log10p",
                         points_color = "grey80", cos_color = NULL,
                         ylim_each = TRUE, gene_name = NULL,
                         outcome_legend_size = 1.5,
                         outcome_legend_pos = "right",
                         cos_legend_pos = "bottomleft",
                         show_variable = FALSE,
                         lab_style = c(2, 1),
                         axis_style = c(1, 1),
                         title_style = c(2.5, 2),
                         ...) {
  args <- list(...)
  if (!exists("pch", args)) args$pch <- 16

  # - set background point color and cos color pools
  args$bg <- points_color
  if (is.null(cos_color)) {
    # cos_color <- c("dodgerblue2", "#6A3D9A", "#FF7F00", "#FB9A99", "#33A02C",
    #                "#A6CEE3",  "gold1", "#01665E","#FDBF6F", "#CAB2D6", "#B2DF8A",
    #                "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#01665E",
    #                "#35978F", "#80CDC1", "#C7EAE5", "#003C30")
    # cos_color <- c("#1B9E77", "#D95F02", "#7570B3", "#1F78B4", "#66A61E",
    #                "#E6AB02", "#A6761D", "#666666", "#E7298A", "#B2182B",
    #                "#D73027", "#4575B4", "#313695", "#542788", "#74ADD1",
    #                "#F46D43", "#4DAF4A", "#984EA3", "#FF7F00", "#A50F15")
    cos_color <- c(
      "#377EB8", "#E69F00", "#33A02C", "#984EA3", "#F46D43",
      "#A65628", "#1F4E79", "#B2182B", "#D73027", "#F781BF",
      "#4DAF4A", "#E41A1C", "#FF7F00", "#6A3D9A", "#1B7837",
      "#E6AB02", "#542788", "#74ADD1", "#A50F15", "#01665E"
    )


    # cos_color <- c("#1F70A9", "#33A02C", "#CAB2D6", "#EA7827")
  }
  args$col <- cos_color

  # - set data and x-lab and y-lab
  if (y == "log10p") {
    plot_data <- lapply(cb_plot_input$Zscores, function(z) {
      -log10(2 * pnorm(-abs(z)))
    })
    ylab <- "-log10(p)"
  } else if (y == "z_original") {
    plot_data <- cb_plot_input$Zscores
    ylab <- "Z score"
  } else if (y == "vcp") {
    plot_data <- cb_plot_input$vcp
    ylab <- "VCP"
    if (length(cb_plot_input$outcomes) == 1) {
      ylab <- "PIP"
    }
    args$ylim <- c(0, 1)
  } else if (y == "coef") {
    plot_data <- cb_plot_input$coef
    ylab <- "Coefficients"
  } else {
    stop("Invalid y value! Choose from 'z' or 'z_original'")
  }
  if (!exists("xlab", args)) args$xlab <- "variables"
  if (!exists("ylab", args)) args$ylab <- ylab
  args$lab_size <- as.numeric(lab_style[1])
  args$lab_face <- lab_style[2]

  # - set title format
  args$title <- gene_name
  args$title_size <- as.numeric(title_style[1])
  args$title_face <- title_style[2]

  # - set x-axis and y-axis
  args$x <- cb_plot_input$x$pos
  args$y <- plot_data
  if (show_variable) {
    args$xaxt <- "n"
    args$xtext <- cb_plot_input$variables
  }
  args$axis_size <- as.numeric(axis_style[1])
  args$axis_face <- axis_style[2]

  # - set ylim for each subfigure
  if (exists("ylim", args)) {
    ymax <- rep(args$ylim[2], length(args$y))
    ymin <- rep(args$ylim[1], length(args$y))
  } else {
    ymax <- NULL
    ymin <- rep(0, length(args$y))
  }
  if (ylim_each & is.null(ymax)) {
    ymax <- sapply(plot_data, function(p) {
      valid_values <- p[is.finite(p)]
      if (length(valid_values) == 0) {
        ymax <- 330 # Default if no valid values
      } else {
        ymax <- max(valid_values) * 1.05
      }
      return(ymax)
    })
    if (y == "coef") {
      ymin <- sapply(plot_data, function(p) min(p) * 1.05)
    }
  }
  args$ymax <- ymax
  args$ymin <- ymin

  # - set legend text position and format
  args$outcome_legend_pos <- outcome_legend_pos
  args$outcome_legend_size <- outcome_legend_size
  if (outcome_legend_pos == "right") {
    args$outcome_legend_angle <- 90
  } else if (outcome_legend_pos == "left") {
    args$outcome_legend_angle <- 270
  } else {
    args$outcome_legend_angle <- 0
  }

  if (!(cos_legend_pos %in% c(
    "bottomright", "bottom", "bottomleft", "left",
    "topleft", "top", "topright", "right", "center"
  ))) {
    cos_legend_pos <- "bottomleft"
  }
  args$cos_legend_pos <- cos_legend_pos

  return(args)
}
