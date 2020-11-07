# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: MIT
# --------------------------------------------------- #
remove(list = ls())
library(tidyverse)
library(colorspace)
library(fields)

load("risks/analysis/results/risks_results_20200930.Rdata")
path <- "risks/analysis/figures/"

# ------------------------------------------
# Define ggplot themes to be used in all plots
theme1 <- function(){
  theme(strip.text = element_text(colour = 1, size = 15),
        panel.background = element_rect(fill = 'grey96'),
        legend.background = element_rect(fill = NA),
        legend.key.size   = unit(0.8, "cm"),
        legend.text.align = 0,
        legend.text  = element_text(colour = 1, size = 12),
        legend.title = element_text(colour = 1, size = 12),
        axis.text    = element_text(colour = 1, size = 12),
        axis.title   = element_text(colour = 1, size = 15),
        axis.title.x = element_text(
          margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.line.x  = element_line(linetype = 1),
        axis.ticks   = element_line(size = 1.2),
        panel.spacing = unit(0.5, "cm"),
        plot.margin   = unit(c(0.1, 1, 0.1, 0.1), "cm")
  )
}

theme2 <- function() {
  theme(strip.text = element_text(size = 12),
        legend.position = 'top',
        legend.key.width = unit(2.0, 'cm'),
        legend.key.size  = unit(0.5, "cm"),
        legend.background = element_rect(fill = NA),
        panel.spacing = unit(0.1, "cm"))
}

# ------------------------------------------
# Figure 1 & 2 - Linearity
plot1 <- function(X){
  xs   <- seq(0, 110, by = 10)
  xLab <- expression(paste('Life Expectancy at Birth, ', e[0]))
  yLab <- expression(paste('Age-Specific Death Rates, ', m[x]))
  xBreaks <- log(c(72, 78, 84))
  yBreaks <- log(c(0.0001, 0.001, 0.01, 0.1, 1))
  D    <- X %>% filter(Age %in% xs, mx > 0.00001) 
  
  r = NULL
  for (i in xs) {
    A <- X %>% filter(Age == i)
    log_mx <- log(A$mx)
    L <- log_mx > -100
    log_mx <- log_mx[L]
    log_e <- log(A$e0)[L]
    r <- c(r, signif(cor(log_mx, log_e), 2))
  }
  
  R = data.frame(Age = xs, r = r)
  
  ggplot(D, aes(x = log(e0), y = log(mx), color = Year)) +
    geom_point(size = 0.7) +
    stat_smooth(method = "lm", col = "red", size = 0.7) +
    facet_wrap(~Age, nrow = 2) +
    scale_color_gradient(name = "Year\n", low = "black", high = "cyan") +
    scale_x_continuous(name = xLab, breaks = xBreaks,
                       labels = exp(xBreaks)) +
    scale_y_continuous(name = yLab, breaks = yBreaks,
                       labels = paste0(exp(yBreaks)*100, "%")) +
    geom_text(data = R, 
              inherit.aes = FALSE, 
              size = 3.5, 
              color = "red",
              hjust = "inward",
              aes(x = log(71), y = log(.00005), label = paste0(r))) +
    theme1() + 
    theme2()
}

plot2 <- function(X){
  xs   <- seq(65, 110, by = 5)
  xLab <- expression(paste('Life Expectancy at Age 65, ', e[65]))
  yLab <- expression(paste('Age-Specific Death Rates, ', m[x]))
  xBreaks <- log(seq(15, 24, 3))
  yBreaks <- log(c(0.001, 0.01, 0.1, 1))
  D    <- X %>% filter(Age %in% xs, mx > 0.00001)
  
  r = NULL
  for (i in xs) {
    A <- X %>% filter(Age == i)
    log_mx <- log(A$mx)
    L <- log_mx > -100
    log_mx <- log_mx[L]
    log_e <- log(A$e0)[L]
    r <- c(r, signif(cor(log_mx, log_e), 2))
  }
  
  R = data.frame(Age = xs, r = r)
  
  ggplot(D, aes(x = log(e65), y = log(mx), color = Year)) +
    geom_point(size = 0.7) +
    stat_smooth(method = "lm", col = 1, size = 0.7) +
    facet_wrap(~Age, nrow = 2, scales = "fixed") +
    scale_color_gradient(name = "Year\n",low = "blue", high = "red") +
    scale_x_continuous(name = xLab, breaks = xBreaks,
                       labels = exp(xBreaks)) +
    scale_y_continuous(name = yLab, breaks = yBreaks,
                       labels = paste0(exp(yBreaks)*100, "%")) +
    geom_text(data = R, 
              inherit.aes = FALSE, 
              size = 3.5, 
              color = 1,
              hjust = "inward",
              aes(x = log(14), y = log(.004), label = paste0(r))) +
    theme1() + 
    theme2()
}


p1 = plot1(Plot1Data)
p1
ggsave(paste0(path, "Figure1-1.pdf"),
       plot= p1,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 5)

p2 = plot2(Plot1Data)
p2
ggsave(paste0(path, "Figure2-1.pdf"),
       plot= p2,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 5)

# ------------------------------------------
# Figure 3 - Estimates
plot3 <- function(X){
  ggplot(X, aes(x = x, y = value, color = Parameters, size = Parameters) ) +
    geom_line() +
    facet_grid(coeff ~ country, scales = "free", switch = "y",
               labeller = labeller(.rows = label_parsed)) +
    scale_colour_manual(name = '', values = c(4, 2)) +
    scale_size_manual(name = '', values = c(1, .6)) +
    scale_x_continuous(name = "Age, x", breaks = seq(0, 120, 30)) + 
    scale_y_continuous(name = "Estimated Parameters") + 
    theme1() +
    theme(legend.justification = c(1, 1),
          legend.position = c(1, 0.5),
          strip.text = element_text(size = 13))
}

p3 = plot3(parEstimate)
p3
ggsave(paste0(path, "Figure3-1.pdf"),
       plot= p3,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 5)

# ------------------------------------------
# Figure 4 - death rates
plot4 <- function(Y, M, K){
  M <- M %>% filter(year == Y)
  K <- K %>% mutate(ex = round(ex, 2), k = round(k, 2)) %>% filter(year == Y)
  xLab <- "Age, x"
  yLab <- expression(paste('Age-Specific Death Rates, ', m['x,2014']))
  yBreaks <- log(c(0.0001, 0.001, 0.01, 0.1, 1))
  
  ggplot(M, aes(x = x, y = log(mx), color = type)) +
    geom_line(aes(linetype=type), size = 1) +
    geom_point(aes(size = type)) +
    facet_wrap(~ country, ncol = 2, scales = "free") +
    scale_linetype_manual(name = '', values = c(1, NA)) +
    scale_colour_manual(name = '', values = c(3, 1)) +
    scale_size_manual(name = '', values = c(0.01, 1.5)) +
    scale_x_continuous(name = xLab, breaks = seq(0, 120, 30)) +
    scale_y_continuous(name = yLab, breaks = yBreaks,
                       labels = paste0(exp(yBreaks)*100, "%")) +
    geom_text(data = K, inherit.aes = FALSE, size = 4, hjust = "inward",
              aes(x = 0, y = log(.5), 
                  label = paste0("e0 =   ", ex, "\nk   =  ", k) )) +
    theme1() +
    theme(legend.justification = c(1, 0),
          legend.position = c(.99, 0.02),
          panel.spacing.x = unit(1.5, "cm"))
}

p4 = plot4(Y = 2018, M = mxEstimate, K = kEstimate)
p4
ggsave(paste0(path, "Figure4-1.pdf"),
       plot= p4,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 8)

# ------------------------------------------
# Figure 5 - Errors
plot5 <- function(E){
  breaks_err <- c(seq(0, 4.5, length.out = 10))
  xBreaks    <- c(1991, seq(1994, 2014, by = 4), 2018)
  myPalette  <- c('#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b', 1)
  legendN <- "Mean Relative\nError (%) \nin log(mx)"
  E2 <- E$E2 %>% mutate(err_cut = cut(MARLE, breaks = breaks_err))
  
  ggplot(E2, aes(x = year, y = country, fill = err_cut)) +
    geom_tile( colour = "white", size = 1.3) +
    scale_x_continuous(name = "", breaks = xBreaks) +
    scale_y_discrete(name = "") +
    scale_fill_manual(values = myPalette, name = legendN) +
    coord_equal(xlim = range(E$E2$year)) +
    theme1() +
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = 'top')
}
p5 = plot5(Errors)
ggsave(paste0(path, "Figure5-1.pdf"),
       plot= p5,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 2.8)

# ------------------------------------------
# Figure 6 - Relative errors
plot6 <- function(E) {
  breaks_err <- c(-25, -20, seq(-10, 10, length.out = 11), 20)
  xBreaks    <- c(1996, 2008, 2018)
  col2 <- c('#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026')
  col1 <- rev(c('#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#6e016b', 1))
  myPalette <- c(col1, col2)
  
  E1 <- E$E1 %>% mutate(err_cut = cut(pmax(RLE, quantile(RLE, 0.0005)), 
                                      breaks = breaks_err))
  ggplot(E1, aes(x = year, y = x, fill = err_cut)) +
    geom_tile() +
    facet_wrap(~ country, ncol = 4) +
    scale_x_continuous(name = 'Year, t', breaks = xBreaks, expand = c(0, 0)) +
    scale_y_continuous(name = 'Age, x', expand = c(0.01, 0.01)) +
    scale_fill_manual(name = "Relative Error (%)", values = myPalette,
                      guide = guide_legend(title.position = "top",
                                           label.position = "bottom",
                                           label.hjust = 0.5,
                                           nrow = 1)) +
    theme1() +
    theme(strip.text = element_text(size = 14),
          panel.background = element_rect(fill = NA, colour = NA),
          legend.background = element_rect(fill = NA),
          legend.key.width = unit(1.4, 'cm'),
          legend.key.size  = unit(0.4, "cm"),
          legend.position = "top")
}

p6 = plot6(Errors)
ggsave(paste0(path, "Figure6-1.pdf"),
       plot= p6,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 4.5)

# ------------------------------------------
# Figure 7 - Forecasts
plot7 <- function(M, K, ci) {
  T1 <- K %>% filter(quantile == ci[1]) %>% mutate(ex=round(ex,1), k=round(k,1))
  T2 <- K %>% filter(quantile == ci[2]) %>% mutate(ex=round(ex,1), k=round(k,1))
  T3 <- K %>% filter(quantile == ci[3]) %>% mutate(ex=round(ex,1), k=round(k,1))
  xLab <- "Age, x"
  yLab <- expression(paste("Predicted Death Rates, ", m['x,2040']))
  yBreaks <- log(c(0.0001, 0.001, 0.01, 0.1, 1))
  sz = 3.7
  
  ggplot(M, aes(x = x, y = log(mx), group = model)) +
    geom_ribbon(aes(ymin = log(mx_lw), ymax = log(mx_up), fill = model)) +
    geom_line(aes(color = model), size = 0.7) +
    facet_wrap(~ country, ncol = 2, scales = "free") +
    scale_colour_manual(name="", values = c('#e41a1c', 'blue')) +
    scale_fill_manual(name="", values = alpha(c('#d95f02', 'blue'), c(.7,.3))) +
    scale_x_continuous(name = xLab, breaks = seq(0, 120, 30)) +
    scale_y_continuous(name = yLab, breaks = yBreaks,
                       labels = paste0(exp(yBreaks)*100, "%")) +
    # Text on the plots
    geom_text(data = T1, inherit.aes = FALSE, size = sz, hjust = "inward",
              aes(x = 1, y = log(.6), 
                  label = paste0('Quantile    e0 =   ', ex, 
                                 '\n0.5%         k   = ', k))) +
    geom_text(data = T2, inherit.aes = FALSE, size = sz, hjust = "inward",
              aes(x = 1, y = log(.1), 
                  label = paste0('Median      e0 =   ', ex, 
                                 '\nprediction  k   = ', k))) +
    geom_text(data = T3, inherit.aes = FALSE, size = sz, hjust = "inward",
              aes(x = 1, y = log(.015), 
                  label = paste0('Quantile    e0 =   ', ex, 
                                 '\n99.5%       k   = ', k))) +
    theme1() + 
    theme(legend.justification = c(1, 0),
          legend.position = c(.99, 0.02),
          panel.spacing.x = unit(1.5, "cm"))
}

p7 = plot7(M = mxForecast, K = kForecast, ci = quant)
ggsave(paste0(path, "Figure7-1.pdf"),
       plot= p7,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 8)

# ------------------------------------------
# Figure 8 - Kannisto model
plot8 <- function(X){
  xLab <- "Age, x"
  yLab <- expression(paste('Old age mortality, ', m['x,2016']))
  
  ggplot(X, aes(x = Age, y = mx, group = Type, size = Type)) +
    geom_point() +
    geom_line(aes(color = Type, linetype = Type), size = 1) +
    facet_wrap(~ country, ncol = 2, scales = "free") +
    geom_vline(xintercept = c(80, 95), linetype = 3) +
    annotate("rect", xmin = 80, xmax = 95, ymin = 1, ymax = 0, alpha = .2) +
    annotate("text", x = 87.5, y = 0.7, label = "Fitted \n age-range" ) +
    scale_linetype_manual(name = "", values = c(2, 1)) +
    scale_color_manual(name = "", values = c(NA, 2)) +
    scale_size_manual(name = "", values = c(1.5, .01)) +
    scale_x_continuous(name = xLab) +
    scale_y_continuous(name = yLab) +
    coord_cartesian(ylim = c(0,1)) + 
    theme1() +
    theme(legend.justification = c(1, 0),
          legend.position = c(.99, 0.02),
          panel.spacing.x = unit(1.5, "cm"))
}

p8 = plot8(dataFig_8)
ggsave(paste0(path, "Figure8-1.pdf"),
       plot= p8,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 8)


# ------------------------------------------
# Figure 9 - LSE vs MLE
plot9 <- function(M1, M2){
  yr = "2018"
  tbl <- rbind(
    data.frame(age = x_LF, value = log(fitted(M1)[, yr]), 
               type = "log(m['x,2018'])", method = "OLS+SVD"),
    data.frame(
      age = x_LF, value = coef(M1)$bx, type = "beta[x]", method = "OLS+SVD"),
    data.frame(
      age = x_LF, value = coef(M1)$vx, type = "nu[x]", method = "OLS+SVD"),
    data.frame(
      age = x_LF, value = log(fitted(M2)[, yr]), type = "log(m['x,2018'])", 
               method = "MLE"),
    data.frame(
      age = x_LF, value = coef(M2)$bx, type = "beta[x]", method = "MLE"),
    data.frame(
      age = x_LF, value = coef(M2)$vx, type = "nu[x]", method = "MLE"))
  
  ggplot(tbl, aes(x = age, y = value, color = method)) +
    geom_line() +
    facet_wrap(~ type, ncol = 3, scales = "free", 
               labeller = labeller(.cols = label_parsed)) +
    scale_x_continuous(name = "Age, x", breaks = seq(0, 120, 30)) + 
    scale_y_continuous(name = "") +
    scale_colour_manual(name = '', values = c(3, 4)) +
    scale_size_manual(name = '', values = c(1, .6)) +
    theme1() +
    theme(legend.justification = c(1, 1),
          legend.position = c(0.15, 1.03),
          strip.text = element_text(size = 13),
          panel.spacing.x = unit(1, "cm"))
}

p9 = plot9(EW_LSE, EW_MLE)
ggsave(paste0(path, "Figure9-1.pdf"),
       plot= p9,
       device = "pdf",
       scale = 1,
       width = 10,
       height = 4)

# ------------------------------------------
# Figure 10 - vx rotation

plot10 <- function(x_LF, x_rot, vx_EW, vx_FR, vx_SW, vx_US, Cnt1){
  n  <- ncol(vx_EW)
  cl <- c(1, diverge_hcl(n - 1))
  ran <- range(vx_EW, vx_FR, vx_SW, vx_US)
  par(mfrow = c(2, 2), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  # Plot 1
  matplot(x_LF, vx_EW, t = "n", main = Cnt1[1], 
          ylab = expression(nu[x]), xlab = "", ylim = ran)
  grid()
  matlines(x_LF, vx_EW, col = cl, t = "l", lwd = 1.2, lty = 1)
  lines(x_LF, vx_EW[, 1], col = cl[1], lwd = 1.5)
  lines(x_LF, vx_EW[, n], col = cl[n], lwd = 2)
  legend('topright', col = cl[c(1, n)], lty = c(1, 1, NA, NA), 
         lwd = 2, bty = 'n', cex = 1.6,
         legend = c(expression(paste(nu[x], " Original")),
                    expression(paste(nu[x], " Ultimate"))))
  # Plot 2
  matplot(x_LF, vx_FR, t = "n", main = Cnt1[2], 
          ylab = "", xlab = "", ylim = ran)
  grid()
  matlines(x_LF, vx_FR, col = cl, t = "l", lwd = 1.2, lty = 1)
  lines(x_LF, vx_FR[, 1], col = cl[1], lwd = 1.5)
  lines(x_LF, vx_FR[, n], col = cl[n], lwd = 2)
  arrows(71, 0.0145, 71, 0.010, col = 1, lwd = 2, angle = 15)
  arrows(26, 0.0055, 26, 0.01, col = 1, lwd = 2, angle = 15)
  legend('topright', legend = expression(paste(e[0], " Level")),
         bty = 'n', cex = 1.6)
  image.plot(smallplot = c(.55, .91, .7, .75), horizontal = TRUE,
             legend.only = TRUE, zlim = range(x_rot),
             col = cl, cex.lab = 0.8, nlevel = n)
  # Plot 3
  matplot(x_LF, vx_SW, t = "n", main = Cnt1[3], 
          ylab = expression(nu[x]), xlab = "Age, x", ylim = ran)
  grid()
  matlines(x_LF, vx_SW, col = cl, t = "l", lwd = 1.2, lty = 1)
  lines(x_LF, vx_SW[, 1], col = cl[1], lwd = 1.5)
  lines(x_LF, vx_SW[, n], col = cl[n], lwd = 2)
  # Plot 4
  matplot(x_LF, vx_US, t = "n", main = Cnt1[4], 
          ylab = "", xlab = "Age, x", ylim = ran)
  grid()
  matlines(x_LF, vx_US, col = cl, t = "l", lwd = 1.2, lty = 1)
  lines(x_LF, vx_US[, 1], col = cl[1], lwd = 1.5)
  lines(x_LF, vx_US[, n], col = cl[n], lwd = 2)
}  

p10 = plot10(x_LF, x_rot, vx_EW, vx_FR, vx_SW, vx_US, Cnt1)

pdf(paste0(path, "Figure10-1.pdf"),
    width = 10,
    height = 10)
plot10(x_LF, x_rot, vx_EW, vx_FR, vx_SW, vx_US, Cnt1)
dev.off()


