library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

volcano_plot <- function(diff_obj, mode = "mu", evidence_thresh = 0.9,
                               change_thresh = log(1.5), size = 1.6, xlab = "Log odds ratio (A/B) change",
                               ylab = "Posterior evidence probability", title = NULL, nlabel = 0) {
  
  dt <- diff_obj %>% as.data.table
  if (mode == "mu") {
    dt %>% setnames(c("mean_LOR","mean_diff_tail_prob","mean_diff_test"),c("x","y","test"))
    x_text <- 5
  } else if (mode == "gamma") {
    dt %>% setnames(c("disp_LOR","disp_diff_tail_prob","disp_diff_test"),c("x","y","test"))
    dt[test=="ExcludedFromTesting",test:="NoDiff"]
    x_text <- 3
  } else if (mode == "epsilon") {
    dt %>% setnames(c("res_disp_LOR","res_disp_diff_tail_prob","res_disp_diff_test"),c("x","y","test"))
    x_text <- 1.5
  }
  dt %>% .[,abs_x:=abs(x)] %>% setorder(-y,-x) %>% .[,abs_x:=NULL]
  
  gg <- ggplot(dt, aes(x = x, y = y, fill = test)) +
    geom_jitter(aes(size = test, alpha = test), shape=21, stroke=0.15, height=0.01, width=0)
  
  # Add N= labels
  foo <- unique(dt$test); foo <- foo[foo!="NoDiff"]
  A_hits <- dt[test==foo[1],.N]
  B_hits <- dt[test==foo[2],.N]
  gg <- gg + annotate("text", x=-x_text, y=1.1, size=4.5, label=sprintf("N=%d (%.2f%%)",A_hits,100*(A_hits/dt[,.N])))
  gg <- gg + annotate("text", x=x_text, y=1.1, size=4.5, label=sprintf("N=%d (%.2f%%)",B_hits,100*(B_hits/dt[,.N])))
  
  gg <- gg +
    geom_hline(yintercept = evidence_thresh, color = "orange", linetype = "dashed") +
    geom_vline(xintercept = c(-change_thresh, change_thresh), color = "dodgerblue4", linetype = "dashed") +
    coord_cartesian(ylim=c(0,1.15)) +
    xlab(xlab) + ylab(ylab) + 
    scale_fill_manual(name = "Hits", values = c("#595959", "lightpink3", "darkolivegreen3","blue")) +
    scale_size_manual(values = c(0.6,1.75,1.75)) +
    scale_alpha_manual(values = c(0.4,0.8,0.8)) +
    ggtitle(title) +
    guides(alpha=F, size=F) +
    guides(fill = guide_legend(override.aes=list(size=3))) +
    theme_classic() +
    theme(
      legend.position = c(.1,.2),
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.3)),
      axis.text = element_text(color="black")
    )
  
  # Add labels to the top "nlabel "significant hits
  if (nlabel>0) {
    gg <- gg + ggrepel::geom_text_repel(data=head(dt,n=nlabel), 
              aes(x=x, y=y, label=feature_name), size=4)
  }  
  return(gg)
}


volcano_plot_scattermore <- function(diff_obj, mode = "mu", diff_feat_idx = NULL, evidence_thresh = 0.9,
                         change_thresh = log(1.5), size = 1.6, xlab = "Log odds ratio (A/B) change",
                         ylab = "Posterior evidence probability", title = NULL) {
  # Plot also
  if (!is.null(diff_feat_idx)) {
      diff_obj$true_diff <- FALSE
      diff_obj$true_diff[diff_feat_idx] <- TRUE
  }
  #diff_obj$true_diff <- factor(diff_obj$true_diff, levels = c("FALSE", "TRUE"))

  if (mode == "mu") {
    gg <- ggplot(data = diff_obj, aes(x = mean_LOR, y = mean_diff_tail_prob))
    if (!is.null(diff_feat_idx)) {
      gg <- gg +
        geom_scattermore(aes(color = mean_diff_test), pointsize = size) +
        scale_shape_manual(name = "True change", values = c(16, 3))
    } else {
      gg <- gg + 
        geom_scattermore(aes(color = mean_diff_test), pointsize = size)
    }
    
  } else if (mode == "gamma") {
    gg <- ggplot(data = diff_obj, aes(x = disp_LOR, y = disp_diff_tail_prob))
    if (!is.null(diff_feat_idx)) {
      gg <- gg +
        geom_scattermore(aes(color = disp_diff_test), pointsize = size) +
        scale_shape_manual(name = "True change", values = c(16, 3))
    } else {
      gg <- gg + 
        geom_scattermore(aes(color = disp_diff_test), pointsize = size)
    }
    
  } else if (mode == "epsilon") {
    gg <- ggplot(data = diff_obj, aes(x = res_disp_LOR, y = res_disp_diff_tail_prob))
    if (!is.null(diff_feat_idx)) {
      gg <- gg +
        geom_scattermore(aes(color = res_disp_diff_test), pointsize = size) +
        scale_shape_manual(name = "True change", values = c(16, 3))
    } else {
      gg <- gg + 
        geom_scattermore(aes(color = res_disp_diff_test), pointsize = size)
    }
  }
  
    # Add N= labels
    # A_hits <- dt[disp_diff_test==TRUE,.N]
    # B_hits <- dt[hvf==FALSE,.N]
    # all <- dt[,.N]
    # p <- p + annotate("text", x=0.50, y=ylim, size=8, label=sprintf("N=%d (%.2f%%)",sig_hits,100*(sig_hits/nosig_hits)))
  
  
  gg <- gg +
    geom_hline(yintercept = evidence_thresh, color = "orange", linetype = "dashed") +
    geom_vline(xintercept = c(-change_thresh, change_thresh), color = "dodgerblue4", linetype = "dashed") +
    scale_color_manual(name = "Hits", values = c("#595959", "lightpink3", "darkolivegreen3", "#595959", "#595959")) +
    xlab(xlab) + ylab(ylab) + 
    ggtitle(title) +
    theme_classic()
  
  return(gg)
}


MA_plot <- function(diff_obj, mode = c("gamma","epsilon"), title = NULL, 
                    xlab = "Mean methylation", ylab = "Log odds ratio (A/B) change", nlabel = 0) {
  
  mode <- match.arg(mode)
  
  dt <- diff_obj %>% as.data.table
  if (mode == "gamma") {
    dt %>% setnames(c("overall_mean","disp_LOR","disp_diff_test"),c("x","y","test"))
    dt[test=="ExcludedFromTesting",test:="NoDiff"]
    # x_text <- 3
  } else if (mode == "epsilon") {
    dt %>% setnames(c("overall_mean","res_disp_LOR","res_disp_diff_test"),c("x","y","test"))
    # x_text <- 1.5
  }
  dt %>% .[,abs_x:=abs(x)] %>% setorder(-y,-x) %>% .[,abs_x:=NULL]
  
  gg <- ggplot(dt, aes(x = x, y = y, fill = test)) +
    geom_jitter(aes(size = test, alpha = test), shape=21, stroke=0.15, height=0.01, width=0)
  
  # Add N= labels
  # foo <- unique(dt$test); foo <- foo[foo!="NoDiff"]
  # A_hits <- dt[test==foo[1],.N]
  # B_hits <- dt[test==foo[2],.N]
  # gg <- gg + annotate("text", x=-x_text, y=1.1, size=4.5, label=sprintf("N=%d (%.2f%%)",A_hits,100*(A_hits/dt[,.N])))
  # gg <- gg + annotate("text", x=x_text, y=1.1, size=4.5, label=sprintf("N=%d (%.2f%%)",B_hits,100*(B_hits/dt[,.N])))
  
  gg <- gg +
    geom_hline(yintercept = 0, color = "orange", linetype = "dashed") +
    # coord_cartesian(ylim=c(0,1.15)) +
    xlab(xlab) + ylab(ylab) + 
    scale_fill_manual(name = "Hits", values = c("#595959", "lightpink3", "darkolivegreen3","blue")) +
    scale_size_manual(values = c(0.6,1.75,1.75)) +
    scale_alpha_manual(values = c(0.4,0.8,0.8)) +
    ggtitle(title) +
    guides(alpha=F, size=F) +
    guides(fill = guide_legend(override.aes=list(size=3))) +
    theme_classic() +
    theme(
      # legend.position = c(.1,.2),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.3)),
      axis.text = element_text(color="black")
    )
  
  # Add labels to the top "nlabel "significant hits
  if (nlabel>0) {
    gg <- gg + ggrepel::geom_text_repel(data=head(dt,n=nlabel), 
                                        aes(x=x, y=y, label=feature_name), size=4)
  }  
  return(gg)
}

MA_plot_scattermore <- function(diff_obj, mode = "mu", diff_feat_idx = NULL, title = NULL) {
  
  # Store data in tibble
  # res_cor <- tibble(epsilon_diff = matrixStats::colMedians(fit_obj_A$posterior$epsilon -
  #                     fit_obj_B$posterior$epsilon ),
  #                   gamma_diff = matrixStats::colMedians(logitnorm::logit(fit_obj_A$posterior$gamma) -
  #                     logitnorm::logit(fit_obj_B$posterior$gamma) ),
  #                   mu_diff = matrixStats::colMedians(logitnorm::logit(fit_obj_A$posterior$mu) -
  #                     logitnorm::logit(fit_obj_B$posterior$mu) ),
  #                   mu_overall = (matrixStats::colMedians(fit_obj_A$posterior$mu) +
  #                     matrixStats::colMedians(fit_obj_B$posterior$mu) ) / 2)
  # if (!is.null(diff_feat_idx)) {
  #   res_cor$diff <- FALSE
  #   res_cor$diff[diff_feat_idx] <- TRUE
  # }
  # 
  
  if (mode == "mu") {
    gg <- ggplot(diff_obj, aes(x = overall_mean, y = mean_LOR)) +
      scattermore::geom_scattermore(aes(color = mean_diff_test), pointsize=1.6) +
      xlab("Mean methylation average")
    
  } else if (mode == "gamma") {
    gg <- ggplot(diff_obj, aes(x = overall_disp, y = disp_LOR)) +
      scattermore::geom_scattermore(aes(color = disp_diff_test), pointsize = 1.6) +
      xlab("Overdispersion average")
    
  } else if (mode == "epsilon") {
    gg <- ggplot(diff_obj, aes(x = overall_res_disp, y = res_disp_LOR)) +
      scattermore::geom_scattermore(aes(color = res_disp_diff_test), pointsize = 1.6) +
    xlab("Res overdispersion average")
  }

  gg <- gg +
    ylab("Log odds ratio (A/B) change") +
    ggtitle(title)
  
  if (!is.null(diff_feat_idx)) {
    gg <- gg +
      geom_point(aes(color = diff), size = 0.8) +
      scale_color_manual(values = c('#595959', 'red')) +
      geom_hline(yintercept = 0, color = "orange", linetype = "dashed")
  } else {
    gg <- gg +
      geom_hline(yintercept = 0, color = "orange", linetype = "dashed") +
      scale_color_manual(name = "Hits", values = c("#595959", "lightpink3", "darkolivegreen3", "#595959", "#595959"))
  }
  return(gg)
}