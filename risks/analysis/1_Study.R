# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: MIT
# --------------------------------------------------- #
remove(list = ls())
set.seed(1234)
library(MortalityEstimate)  #v0.9.10
library(MortalityLaws)      #v1.8.5
library(StMoMo)             #v0.4.1
library(tidyverse)          #v1.3.0

# ------------------------------------------
# Analysis set-up
# DATA 
sex    <- "Female"
model1 <- 'Linear-Link'
model2 <- 'Lee-Carter'
model3 <- 'kannisto'
# Scenario 1 - Backtesting
years_O <- 1965:1990 # Years LL fit
years_E <- 1991:2018 # Years LL estimate
years_T <- c(years_O, years_E) # All years
# Scenario 2 - Forecasting
years_O2 <- 1980:2018 # Years
years_FC <- 2019:2040 # Years to forecast
nsim     <- 1000
quant    <- c(0.5, 50, 99.5)
# vx rotation
e0_start     <- 80
e0_threshold <- 75
e0_ultimate  <- 102
# ages
x_LF   <- 0:120     # Ages LL fit
x_KF   <- 80:95     # Ages Kannisto fit
x_KE   <- 85:120    # Ages Kannisto predict
x_KE2  <- 96:120    # Ages Kannisto predict
x_E    <- 0:100     # Ages for computing errors
x_LC   <- 0:95      # Ages for Lee-Carter fit
x_rot  <- seq(e0_start, e0_ultimate) # Levels of e0 to use in rotation plot
# Country codes
Cnt1 <- c("ENGLAND & WALES", "FRANCE", "SWEDEN", "USA")
Cnt2 <- c("GBRTENW", "FRATNP", "SWE", "USA") #
Cnt3 <- c("EW", "FR", "SW", "US") #
N <- length(Cnt1)
# Data for 4 countries
mx_EW <- HMD4mx$GBRTENW
mx_FR <- HMD4mx$FRATNP
mx_SW <- HMD4mx$SWE
mx_US <- HMD4mx$USA
  
# ------------------------------------------
# BACKTESTING  
  
# Kannisto model (see K1.. K3): extend mortality curve up to 120
for (i in 1:N) {
  Cn  <- Cnt3[i]
  mxi <- get(paste0("mx_", Cn))
  mdl <- MortalityLaw(x = x_KF, 
                      mx = mxi[paste(x_KF), ], 
                      law = "kannisto", 
                      scale.x = TRUE, 
                      show = FALSE) # Fit Kannisto
  prd <- rbind(mxi[paste(0:84), ], predict(mdl, x = x_KE)) #predict
  # assign(paste0("K_", Cn), mdl)
  assign(paste0("Mx_", Cn), prd)  
  remove(i, Cn, mxi, mdl, prd)
}

# Fit Linear-Link (see F1...F6)
# Scenario 1 - smooth
for (j in 1:N) {
  Cn  <- Cnt3[j]
  mxi <- get(paste0("Mx_", Cn))[, paste0(years_O)]
  L1  <- LinearLink(x = x_LF, 
                    mx = mxi, 
                    y = years_O, 
                    country = Cnt1[j], 
                    theta = 0)
  L2  <- LinearLink(x = x_LF, 
                    mx = mxi, 
                    y = years_O, 
                    country = Cnt1[j], 
                    theta = 0, 
                    use.smooth = FALSE)
  assign(paste0("F_", Cn, 1), L1)
  assign(paste0("F_", Cn, 2), L2)
  remove(j, Cn, mxi, L1, L2)
}

# Compute observed life tables (LT1)
R1 <- expand.grid(country = Cnt3, year = years_T, stringsAsFactors = FALSE)
LT_O <- NULL 

for (j in 1:nrow(R1)) {
  # j = 201
  # print(j)
  Cn  <- R1[j, 1]
  yj  <- R1[j, 2]
  mxi <- get(paste0("mx_", Cn))[paste(x_E), paste0(yj)]
  LT  <- LifeTable(x = x_E, mx = mxi)$lt[, c(1:3, 10)]
  LT  <- data.frame(country = Cn, year = yj, LT)
  LT_O <- rbind(LT_O, LT)
  remove(j, Cn, yj, mxi, LT)
}
LT_O$type <- 'Observed'

# Predict life tables between 1991 and 2018
kEstimate <- LT_O %>% 
  filter(x == 0, year %in% years_E) %>% 
  mutate(k = NA)
LT_E <- NULL 

for (j in 1:nrow(kEstimate)) {
  Cn  <- kEstimate[j, 1]
  yj  <- kEstimate[j, 2]
  exj <- kEstimate[j, "ex"]
  Fx  <- get(paste0("F_", Cn, 1))
  Px  <- LinearLinkLT(object = Fx, 
                      ex = exj, 
                      use.vx.rotation = TRUE, 
                      e0_threshold = e0_threshold)
  LT   <- data.frame(country = Cn, year = yj, Px$lt[, c(1:3, 10)])
  LT_E <- rbind(LT_E, LT)
  kEstimate$k[j] <- Px$k
  remove(j, Cn, yj, exj, Fx, Px, LT)
}

LT_E$type <- "Estimated"
levels(LT_O$country) <- 
  levels(LT_E$country) <- 
  levels(kEstimate$country) <- Cnt1

# Mortality curves: observed (LT_O), estimated (LT_E)
mxEstimate <- rbind(LT_O, LT_E) %>% filter(year %in% years_E)
mxEstimate$type <- as.factor(mxEstimate$type)

remove(LT_E, LT_O, mx_EW, mx_FR, mx_SW, mx_US, R1)

#-----------------------------------------
# Create a table with coefficients to be used in ggplots
parEstimate = NULL
for (k in 1:(N*2)) {
  env.objects = ls()
  mdls <- env.objects[grep("F_", env.objects)]
  mdl  <- get(mdls[k])
  cf   <- coef(mdl)
  st   <- ifelse(mdl$input$use.smooth, "Smooth", "Raw")
  out  <- with(mdl$input, 
               data.frame(
                 bx = cf$bx, 
                 vx = cf$vx, 
                 country = country, 
                 x = x, 
                 Parameters = st))
  parEstimate <- rbind(parEstimate, out)
  remove(k, env.objects, mdls, mdl, cf, st, out)
}
parEstimate <- parEstimate %>% 
  gather(., key = coef, value = value, -(country:Parameters))
parEstimate$coeff <- factor(parEstimate$coef, labels = c("beta[x]", "nu[x]"))

# ----------------------------------------------
# vx rotation
vx <- matrix(NA, nrow = length(x_LF), ncol = length(x_rot))
vx_EW = vx_FR = vx_SW = vx_US <- vx
for (i in 1:length(x_rot)) {
  vx_EW[, i] <- LinearLinkLT(object = F_EW1, ex = x_rot[i], 
                             use.vx.rotation = TRUE, 
                             e0_threshold = e0_threshold)$vx
  vx_FR[, i] <- LinearLinkLT(object = F_FR1, ex = x_rot[i], 
                             use.vx.rotation = TRUE, 
                             e0_threshold = e0_threshold)$vx
  vx_SW[, i] <- LinearLinkLT(object = F_SW1, ex = x_rot[i], 
                             use.vx.rotation = TRUE, 
                             e0_threshold = e0_threshold)$vx
  vx_US[, i] <- LinearLinkLT(object = F_US1, ex = x_rot[i], 
                             use.vx.rotation = TRUE, 
                             e0_threshold = e0_threshold)$vx
}
remove(i, vx)

# ------------------------------------------
# Compute_Errors  
compute_errors <- function(M, vsn = 0.00001) {
  E1 <- M %>% 
    mutate(mx_ = pmax(mx, 0.00001)) %>% 
    filter(x %in% x_E) %>% 
    select(-c(x.int, ex, mx)) %>% 
    spread(key = type, value = mx_) %>% 
    mutate(RLE = (log(Observed) - log(Estimated))/log(Observed) * 100, 
           ARLE = abs(RLE))
  
  E2 <- E1 %>% 
    group_by(country, year) %>% 
    summarise(MARLE = mean(ARLE, na.rm = TRUE)) %>% 
    data.frame()
  return(list(E1 = E1, E2 = E2))
}

Errors <- compute_errors(mxEstimate)
rel_errors <- round(range(Errors$E2$MARLE), 1)
  
# ------------------------------------------
# Lee-Carter Forecast  
  
# Lee-Carter model (1992)
constLC <- function(ax, bx, kt, b0x, gc, wxt, ages){
  c1 <- mean(kt[1, ], na.rm = TRUE)
  c2 <- sum(bx[, 1], na.rm = TRUE)
  list(ax = ax + c1 * bx, bx = bx / c2, kt = c2 * (kt - c1))
}

fun_ex <- function(Z) LifeTable(x = x_LF, mx = Z)$lt$ex

LCforecast_ex = LCforecast_mx <- data.frame()
# k = 1
# i = 1
for (k in 1:N) {
  # Data
  Cn  <- Cnt3[k]
  mxt <- get(paste0("Mx_", Cn))
  Dxt <- mxt[paste(x_LC), paste(years_O2)] * 10^5
  Ext <- Dxt * 0 + 10^5
  yrt <- as.numeric(colnames(Dxt))
  #--------------------------------------------------------------------------
  # Defining the models
  wxt <- genWeightMat(ages = x_LC, years = yrt, clip = 3) # weighting matrix
  h   <- max(years_FC) - max(yrt)
  LC  <- StMoMo(link = "logit", staticAgeFun = TRUE, periodAgeFun = "NP",
                constFun = constLC)
  LCfit <- fit(LC, Dxt = Dxt, Ext = Ext, ages = x_LC, years = yrt,
               ages.fit = x_LC, wxt = wxt)
  LCfor <- forecast(LCfit, h = h, jumpchoice = "actual")
  LCsim <- simulate(LCfit, nsim = nsim, h = h, jumpchoice = "actual")
  
  LC_mx = LC_ex <- data.frame()
  for (i in 1:length(quant)) {
    LC_q  <- apply(LCsim$rates, c(1, 2), quantile, probs = quant[i]/100)
    Kan_q <- MortalityLaw(x = x_KF, 
                          mx = LC_q[paste(x_KF), ], 
                          law = model3, 
                          scale.x = TRUE, 
                          show = FALSE)
    J1 <- data.frame(rbind(LC_q, predict(Kan_q, x_KE2)))
    J2 <- data.frame(apply(J1, 2, fun_ex))
    colnames(J1) <- colnames(J2) <- years_FC
    
    J1 <- J1 %>% 
      mutate(x = x_LF) %>%
      reshape(., direction = 'long', 
              idvar = 'x', 
              varying = paste(years_FC),
              v.names = 'mx', 
              timevar = 'year', 
              times = years_FC) %>%
      mutate(quantile = quant[i], 
             country = Cn, 
             model = model2)
    
    J2 <- J2 %>% mutate(x = x_LF) %>%
      reshape(., direction = 'long', 
              idvar = 'x', 
              varying = paste(years_FC),
              v.names = 'ex', 
              timevar = 'year', 
              times = years_FC) %>%
      mutate(quantile = quant[i], 
             country = Cn, 
             model = model2)
    
    LC_mx <- rbind(LC_mx, J1)
    LC_ex <- rbind(LC_ex, J2)
    remove(i, J1, J2)
  }
  
  LCforecast_mx <- rbind(LCforecast_mx, LC_mx)
  LCforecast_ex <- rbind(LCforecast_ex, LC_ex)
  remove(k, Cn, mxt, Dxt, Ext, wxt, h, yrt, LC, LC_mx, LC_ex, LC_q, Kan_q,
         LCfit, LCfor, LCsim)
}

LCforecast_ex <- LCforecast_ex %>% filter(year == max(years_FC), x == 0)
LCforecast_mx <- LCforecast_mx %>% filter(year == max(years_FC))
LCforecast_mx$country <- factor(LCforecast_mx$country, labels = Cnt1)
  
# ------------------------------------------
# Linear-Link Forecast
  
# Fit Linear-Link - smooth
for (j in 1:N) {
  Cn  <- Cnt3[j]
  mxi <- get(paste0("Mx_", Cn))[, paste0(years_O2)]
  mdl <- LinearLink(x = x_LF, 
                    mx = mxi, 
                    y = years_O2, 
                    country = Cnt1[j], 
                    theta = 0, 
                    use.smooth = FALSE)
  assign(paste0("F_", Cn, 3), mdl)
  remove(j, Cn, mxi, mdl)
}

LLforecast_mx <- kForecast <- NULL
for (i in 1:nrow(LCforecast_ex)) {
  Cn <- LCforecast_ex[i, 'country']
  Yn <- LCforecast_ex[i, 'year']
  Qn <- LCforecast_ex[i, 'quantile']
  en <- LCforecast_ex[i, 'ex']
  Fx <- get(paste0('F_', Cn, 3))
  Px <- LinearLinkLT(Fx, en, use.vx.rotation = TRUE)
  LLforecast_mx <- rbind(LLforecast_mx, 
                         data.frame(
                           country = Cn, 
                           year = Yn, 
                           quantile = Qn, 
                           Px$lt[, c('x', 'mx')], 
                           model = model1))
  kForecast <- rbind(kForecast, 
                     data.frame(
                       country = Cn, 
                       year = Yn, 
                       quantile = Qn, 
                       k = Px$k, 
                       ex = Px$lt$ex[1]))
}

LLforecast_mx$country <- factor(LLforecast_mx$country, labels = Cnt1)
kForecast$country <- factor(kForecast$country, labels = Cnt1)

mxFC <- rbind(LCforecast_mx, LLforecast_mx)
mxForecast       <- mxFC[mxFC$quantile == quant[2], ]
mxForecast$mx_lw <- mxFC[mxFC$quantile == quant[1], "mx"]
mxForecast$mx_up <- mxFC[mxFC$quantile == quant[3], "mx"]

remove(i, Cn, Yn, Qn, en, Fx, Px, mxFC, 
       LCforecast_ex, LCforecast_mx, LLforecast_mx)
remove(F_EW1, F_EW2, F_EW3, F_FR1, F_FR2, F_FR3, 
       F_SW1, F_SW2, F_SW3, F_US1, F_US2, F_US3)

# ------------------------------------------ 
# Least square vs Maximum Likelihood estimation  

# Comparison between LSE and MLE estimation methods
EW_LSE <- LinearLink(x = x_LF, 
                     mx = Mx_EW[, paste(years_O2)], 
                     years_O2, 
                     theta = 0, 
                     method = "LSE", 
                     use.smooth = FALSE)

EW_MLE <- LinearLink(x = x_LF, 
                     mx = Mx_EW[, paste(years_O2)], 
                     years_O2, 
                     theta = 0, 
                     method = "MLE", 
                     use.smooth = FALSE)
  
# ------------------------------------------
# Get ready the data for plot 1 and 2  

# HMD Dx and Ex to compute raw mx
load("risks/analysis/data/HMD_Dx_20200930.Rdata") 
load("risks/analysis/data/HMD_Ex_20200930.Rdata") 
HMD_mx <- cbind(HMD_Dx$data[, 1:3], HMD_Dx$data[, 4:6]/HMD_Ex$data[, 4:6])
Cnt4   <- HMD_Dx$input$countries # all countries here
EDxt   <- HMD_mx[HMD_mx$Year %in% years_O2, 1:4] %>% 
  rename(mx = Female) %>% 
  mutate(ex = NA, e0 = NA, e65 = NA)
EDxt[is.na(EDxt$mx), ]$mx <- 1

Plot1Data <- NULL
for (i in 1:length(Cnt4)) {
  DE <- EDxt[EDxt$country == Cnt4[i], ]
  yr <- unique(DE$Year) 
  
  for (j in 1:length(yr)) {
    tab     <- DE[DE$Year == yr[j], ]
    tab$ex  <- with(tab, LifeTable(x = Age, mx = mx)$lt$ex)
    tab     <- tab[tab$Age %in% seq(0, 110, by = 5), ]
    tab$e0  <- tab[1, "ex"]
    tab$e65 <- tab[14, "ex"]
    Plot1Data <- rbind(Plot1Data, tab)
  }
}
dataFig_1_2 <- Plot1Data[complete.cases(Plot1Data), ]
remove(i, j, DE, yr, tab, EDxt, HMD_Dx, HMD_Ex)

# ------------------------------------------
# Generate results for plot 7  
yr = 2016
K_E <- rbind(
  data.frame(Age = x_LF, country = Cnt1[1], mx = Mx_EW[, paste(yr)]),
  data.frame(Age = x_LF, country = Cnt1[2], mx = Mx_FR[, paste(yr)]),
  data.frame(Age = x_LF, country = Cnt1[3], mx = Mx_SW[, paste(yr)]),
  data.frame(Age = x_LF, country = Cnt1[4], mx = Mx_US[, paste(yr)])) %>%
  filter(Age %in% x_KE) %>% mutate(Type = "Predicted")

K_O <- HMD_mx[, 1:4] %>% 
  filter(Year == yr, Age >= 75, country %in% Cnt2) %>% 
  mutate(Type = "Observed") %>% 
  select(-Year) %>% 
  rename(mx = Female)

K_O$country <- factor(K_O$country, levels = Cnt2)
K_O$country <- factor(K_O$country, labels = Cnt1)
dataFig_8 <- rbind(K_O, K_E)

remove(HMD_mx, K_O, K_E, Mx_EW, Mx_SW, Mx_FR, Mx_US)
  
# ------------------------------------------
# SAVE the results  
save(list = ls(), file = "risks/analysis/results/risks_results_20200930.Rdata")
  
  