# R-BIBLIOTEKER.
library(quantmod)
library(moments)
library(VAR.etp)
library(vars)
library(MASS)
library(tseries)
library(car)
library(zoo)
library(xts)

# INDLÆS DATA. BENYT: "Samling.csv".
input <- zoo::read.csv.zoo('C:/Users/AKF/Documents/Matematik-Okonomi/5. Ar/2. Semester/Master-Thesis/Input/Samling.csv', format = "%Y%m%d")

# MANIPULATION AF DATA.

# AKTIVKLASSER
INF    <- log(1+input$cpiret)
TB     <- log(1+input$t90ret)
NET_TB <- TB - INF
AKT    <- log(1 + input$vwretd) - TB
S_OBL  <- log(1 + input$b10ret) - TB
V_OBL  <- log(1 + input$corpr) - TB

# FINANSIELLE INDIKATORER
DP <- input$dp
EP <- log(as.numeric(input$vwindx)) - log(input$earn)

BM      <- input$bm
AKT_VAR <- log(1 + input$svar)
SMB     <- input$smb
HML     <- input$hml
EP_J    <- EP + (var(EP, na.rm = TRUE)/2)
DP_J    <- DP + (var(DP, na.rm = TRUE)/2)


# TERM STRUCTURE
T_SPREAD <- input$t10y/100 - input$t90y/100
Y_SPREAD <- input$t5y/100 - input$t90y/100
C_SPREAD <- input$baa/100 - input$t90y/100
D_SPREAD <- input$baa/100 - input$aaa/100

# MAKROØKONOMISKE 
FR <- input$fedfund/100

# JUSTERING
TB_J     <- TB + (var(TB, na.rm = TRUE)/2)
NET_TB_J <- NET_TB + (var(NET_TB, na.rm = TRUE)/2)
AKT_J    <- AKT + (var(AKT, na.rm = TRUE)/2)
S_OBL_J  <- S_OBL + (var(S_OBL, na.rm = TRUE)/2)
V_OBL_J  <- V_OBL + (var(V_OBL, na.rm = TRUE)/2)

# BESKRIVENDE STATISTIK

## AKTIVKLASSER
DATA   <- data.frame(NET_TB, AKT, S_OBL, V_OBL)
DATA_J <- data.frame(NET_TB_J, AKT_J, S_OBL_J, V_OBL_J)

ADJ_DATA_MEAN <- colMeans(DATA_J)
DATA_STD      <- apply(DATA, 2, sd)
DATA_QUANTILE <- apply(DATA, 2, quantile)
DATA_SKEWNESS <- apply(DATA, 2, skewness)
DATA_KURTOSIS <- apply(DATA, 2, kurtosis)

DESCRIPTIVE <- data.frame(adj_mean=ADJ_DATA_MEAN,
                          sd=DATA_STD,
                          sr=ADJ_DATA_MEAN/DATA_STD,
                          skew=DATA_SKEWNESS,
                          kurt=DATA_KURTOSIS,
                          quantile=t(DATA_QUANTILE),
                          ac=unlist(lapply(1:ncol(DATA), function(x) acf(DATA[,x], plot = FALSE, lag.max = 1)$acf))[c(rep(FALSE,1),TRUE)])
DESCRIPTIVE$sr[1] <- NA

JB <- data.frame(teststat = unlist(lapply(1:ncol(DATA), function(x) jarque.bera.test(DATA[,x])$statistic)),
                 p = unlist(lapply(1:ncol(DATA), function(x) jarque.bera.test(DATA[,x])$p.value)))

## TILSTANDSVARIABLE
DATA_T   <- xts::xts(merge(DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(TB))
DATA_T_J <- xts::xts(merge(DP_J, EP_J, TB_J), index(TB_J))

ADJ_DATA_T_MEAN <- colMeans(DATA_T_J)
DATA_T_MEAN     <- colMeans(DATA_T)
DATA_T_STD      <- apply(DATA_T, 2, sd)
DATA_T_QUANTILE <- apply(DATA_T, 2, quantile)
DATA_T_SKEWNESS <- apply(DATA_T, 2, skewness)
DATA_T_KURTOSIS <- apply(DATA_T, 2, kurtosis)

DESCRIPTIVE_T <- data.frame(mean=DATA_T_MEAN,
                            sd=DATA_T_STD,
                            skew=DATA_T_SKEWNESS,
                            kurt=DATA_T_KURTOSIS,
                            quantile=t(DATA_T_QUANTILE),
                            ac=unlist(lapply(1:ncol(DATA_T), function(x) acf(DATA_T[,x], plot = FALSE, lag.max = 1)$acf))[c(rep(FALSE,1),TRUE)])


DESCRIPTIVE_T$mean[1:2] <- ADJ_DATA_T_MEAN[1:2]
DESCRIPTIVE_T$mean[7]   <- ADJ_DATA_T_MEAN[3]

JB_T <- data.frame(teststat = unlist(lapply(1:ncol(DATA_T), function(x) jarque.bera.test(DATA_T[,x])$statistic)),
                   p = unlist(lapply(1:ncol(DATA_T), function(x) jarque.bera.test(DATA_T[,x])$p.value)))

# STATIONARITET
data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

STAT_TABLE <- data.frame(teststat_0=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=0)$statistic)),
                         p_0=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=0)$p.value)),
                         teststat_1=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=1)$statistic)),
                         p_1=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=1)$p.value)),
                         teststat_2=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=2)$statistic)),
                         p_2=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=2)$p.value)),
                         teststat_3=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=3)$statistic)),
                         p_3=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=3)$p.value)),
                         teststat_4=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=4)$statistic)),
                         p_4=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=4)$p.value)),
                         teststat_5=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=5)$statistic)),
                         p_5=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=5)$p.value)))

data_diff <- xts::xts(merge(diff(DP), diff(BM), diff(TB), diff(FR)), index(diff(DP)))

STAT_TABLE_DIFF <- data.frame(teststat_0=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=0)$statistic)),
                              p_0=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=0)$p.value)),
                              teststat_1=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=1)$statistic)),
                              p_1=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=1)$p.value)),
                              teststat_2=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=2)$statistic)),
                              p_2=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=2)$p.value)),
                              teststat_3=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=3)$statistic)),
                              p_3=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=3)$p.value)),
                              teststat_4=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=4)$statistic)),
                              p_4=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=4)$p.value)),
                              teststat_5=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=5)$statistic)),
                              p_5=unlist(lapply(1:ncol(data_diff), function(x) adf.test(data_diff[,x], k=5)$p.value)))

STAT_TABLE <- data.frame(teststat_0=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=0)$statistic)),
                         p_0=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=0)$p.value)),
                         teststat_1=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=1)$statistic)),
                         p_1=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=1)$p.value)),
                         teststat_2=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=2)$statistic)),
                         p_2=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=2)$p.value)),
                         teststat_3=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=3)$statistic)),
                         p_3=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=3)$p.value)),
                         teststat_4=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=4)$statistic)),
                         p_4=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=4)$p.value)),
                         teststat_5=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=5)$statistic)),
                         p_5=unlist(lapply(1:ncol(data), function(x) adf.test(data[,x], k=5)$p.value)))

STAT_TABLE_KPSS <- data.frame(teststat_l=unlist(lapply(1:ncol(data), function(x) kpss.test(data[,x], null = "L")$statistic)),
                              p_l=unlist(lapply(1:ncol(data), function(x) kpss.test(data[,x], null = "L")$p.value)),
                              teststat_t=unlist(lapply(1:ncol(data), function(x) kpss.test(data[,x], null = "T")$statistic)),
                              p_t=unlist(lapply(1:ncol(data), function(x) kpss.test(data[,x], null = "T")$p.value)))

STAT_TABLE_PP <- data.frame(teststat_a=unlist(lapply(1:ncol(data), function(x) pp.test(data[,x], type = "Z(alpha)")$statistic)),
                            p_a=unlist(lapply(1:ncol(data), function(x) pp.test(data[,x], type = "Z(alpha)")$p.value)),
                            teststat_at=unlist(lapply(1:ncol(data), function(x) pp.test(data[,x], type = "Z(t_alpha)")$statistic)),
                            p_at=unlist(lapply(1:ncol(data), function(x) pp.test(data[,x], type = "Z(t_alpha)")$p.value)))

# UNIVARIATE
FIT_AKT <- lapply(5:ncol(data), function(x) lm(as.numeric(xts::lag.xts(data$AKT,-1, na.pad = FALSE)) ~ as.numeric(head(data[,x], -1)), data=data))
FIT_AKT_SUMMARY <- lapply(FIT_AKT, summary)

FIT_TABLE <- data.frame(coef=sapply(FIT_AKT, coef)[2,],
                        cf=data.frame(unlist(lapply(1:length(FIT_AKT), function(x) coefci(FIT_AKT[[x]], vcov. = NeweyWest(FIT_AKT[[x]]), parm = 2)))[c(rep(TRUE,1),FALSE)], unlist(lapply(1:length(FIT_AKT), function(x) coefci(FIT_AKT[[x]], vcov. = NeweyWest(FIT_AKT[[x]]), parm = 2)))[c(rep(FALSE,1),TRUE)]),
                        sd=unlist(lapply(1:length(FIT_AKT), function(x) sqrt(diag(NeweyWest(FIT_AKT[[x]])))[2])),
                        t=unlist(lapply(1:length(FIT_AKT), function(x) coeftest(FIT_AKT[[x]], vcov. = NeweyWest(FIT_AKT[[x]]))[6])),
                        p=unlist(lapply(1:length(FIT_AKT), function(x) coeftest(FIT_AKT[[x]], vcov. = NeweyWest(FIT_AKT[[x]]))[8])),
                        rsq=unlist(lapply(1:length(FIT_AKT), function(x) FIT_AKT_SUMMARY[[x]]$r.squared)))

FIT_S_OBL <- lapply(5:ncol(data), function(x) lm(as.numeric(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE)) ~ as.numeric(head(data[,x], -1)), data=data))
FIT_S_OBL_SUMMARY <- lapply(FIT_S_OBL, summary)

FIT_S_TABLE <- data.frame(coef=sapply(FIT_S_OBL, coef)[2,],
                          cf=data.frame(unlist(lapply(1:length(FIT_S_OBL), function(x) coefci(FIT_S_OBL[[x]], vcov. = NeweyWest(FIT_S_OBL[[x]]), parm = 2)))[c(rep(TRUE,1),FALSE)], unlist(lapply(1:length(FIT_S_OBL), function(x) coefci(FIT_S_OBL[[x]], vcov. = NeweyWest(FIT_S_OBL[[x]]), parm = 2)))[c(rep(FALSE,1),TRUE)]),
                          sd=unlist(lapply(1:length(FIT_S_OBL), function(x) sqrt(diag(NeweyWest(FIT_S_OBL[[x]])))[2])),
                          t=unlist(lapply(1:length(FIT_S_OBL), function(x) coeftest(FIT_S_OBL[[x]], vcov. = NeweyWest(FIT_S_OBL[[x]]))[6])),
                          p=unlist(lapply(1:length(FIT_S_OBL), function(x) coeftest(FIT_S_OBL[[x]], vcov. = NeweyWest(FIT_S_OBL[[x]]))[8])),
                          rsq=unlist(lapply(1:length(FIT_S_OBL), function(x) FIT_S_OBL_SUMMARY[[x]]$r.squared)))

FIT_V_OBL <- lapply(5:ncol(data), function(x) lm(as.numeric(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE)) ~ as.numeric(head(data[,x], -1)), data=data))
FIT_V_OBL_SUMMARY <- lapply(FIT_V_OBL, summary)

FIT_V_TABLE <- data.frame(coef=sapply(FIT_V_OBL, coef)[2,],
                          cf=data.frame(unlist(lapply(1:length(FIT_V_OBL), function(x) coefci(FIT_V_OBL[[x]], vcov. = NeweyWest(FIT_V_OBL[[x]]), parm = 2)))[c(rep(TRUE,1),FALSE)], unlist(lapply(1:length(FIT_V_OBL), function(x) coefci(FIT_V_OBL[[x]], vcov. = NeweyWest(FIT_V_OBL[[x]]), parm = 2)))[c(rep(FALSE,1),TRUE)]),
                          sd=unlist(lapply(1:length(FIT_V_OBL), function(x) sqrt(diag(NeweyWest(FIT_V_OBL[[x]])))[2])),
                          t=unlist(lapply(1:length(FIT_V_OBL), function(x) coeftest(FIT_V_OBL[[x]], vcov. = NeweyWest(FIT_V_OBL[[x]]))[6])),
                          p=unlist(lapply(1:length(FIT_V_OBL), function(x) coeftest(FIT_V_OBL[[x]], vcov. = NeweyWest(FIT_V_OBL[[x]]))[8])),
                          rsq=unlist(lapply(1:length(FIT_V_OBL), function(x) FIT_V_OBL_SUMMARY[[x]]$r.squared)))


# MULTIPLE REGRESSION 
data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

data_mult_akt <- as.xts(merge(xts::lag.xts(data$AKT, -1, na.pad = FALSE), head(data, -1)))
data_mult_akt <- subset(data_mult_akt, select = -AKT.1)

data_mult_s <- as.xts(merge(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE), head(data, -1)))
data_mult_s <- subset(data_mult_s, select = -S_OBL.1)

data_mult_v <- as.xts(merge(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE), head(data, -1)))
data_mult_v <- subset(data_mult_v, select = -V_OBL.1)

# KITCHEN SINK
fit_AKT <- lm(AKT ~ DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_mult_akt)

step_AKT_model <- step(fit_AKT, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=AKT ~ TB, upper=fit_AKT))

summary(step_AKT_model)

FIT_M_TABLE <- data.frame(coef=coef(fit_AKT),
                          cf=data.frame(coefci(fit_AKT, vcov. = NeweyWest(fit_AKT))[,1],
                                        coefci(fit_AKT, vcov. = NeweyWest(fit_AKT))[,2]),
                          sd=sqrt(diag(NeweyWest(fit_AKT))),
                          t=coeftest(fit_AKT, vcov. = NeweyWest(fit_AKT))[,3],
                          p=coeftest(fit_AKT, vcov. = NeweyWest(fit_AKT))[,4])

FIT_M_STEP_TABLE <- data.frame(coef=coef(step_AKT_model),
                               cf=data.frame(coefci(step_AKT_model, vcov. = NeweyWest(step_AKT_model))[,1],
                                             coefci(step_AKT_model, vcov. = NeweyWest(step_AKT_model))[,2]),
                               sd=sqrt(diag(NeweyWest(step_AKT_model))),
                               t=coeftest(step_AKT_model, vcov. = NeweyWest(step_AKT_model))[,3],
                               p=coeftest(step_AKT_model, vcov. = NeweyWest(step_AKT_model))[,4])

fit_S_OBL <- lm(S_OBL ~ DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_mult_s)

step_S_OBL_model <- step(fit_S_OBL, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=S_OBL ~ TB, upper=fit_S_OBL))

summary(step_S_OBL_model)

FIT_M_S_TABLE <- data.frame(coef=coef(fit_S_OBL),
                            cf=data.frame(coefci(fit_S_OBL, vcov. = NeweyWest(fit_S_OBL))[,1],
                                          coefci(fit_S_OBL, vcov. = NeweyWest(fit_S_OBL))[,2]),
                            sd=sqrt(diag(NeweyWest(fit_S_OBL))),
                            t=coeftest(fit_S_OBL, vcov. = NeweyWest(fit_S_OBL))[,3],
                            p=coeftest(fit_S_OBL, vcov. = NeweyWest(fit_S_OBL))[,4])

FIT_M_S_STEP_TABLE <- data.frame(coef=coef(step_S_OBL_model),
                                 cf=data.frame(coefci(step_S_OBL_model, vcov. = NeweyWest(step_S_OBL_model))[,1],
                                               coefci(step_S_OBL_model, vcov. = NeweyWest(step_S_OBL_model))[,2]),
                                 sd=sqrt(diag(NeweyWest(step_S_OBL_model))),
                                 t=coeftest(step_S_OBL_model, vcov. = NeweyWest(step_S_OBL_model))[,3],
                                 p=coeftest(step_S_OBL_model, vcov. = NeweyWest(step_S_OBL_model))[,4])


fit_V_OBL <- lm(V_OBL ~ DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_mult_v)

step_V_OBL_model <- step(fit_V_OBL, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=V_OBL ~ TB, upper=fit_V_OBL))

summary(step_V_OBL_model)

FIT_M_V_TABLE <- data.frame(coef=coef(fit_V_OBL),
                            cf=data.frame(coefci(fit_V_OBL, vcov. = NeweyWest(fit_V_OBL))[,1],
                                          coefci(fit_V_OBL, vcov. = NeweyWest(fit_V_OBL))[,2]),
                            sd=sqrt(diag(NeweyWest(fit_V_OBL))),
                            t=coeftest(fit_V_OBL, vcov. = NeweyWest(fit_V_OBL))[,3],
                            p=coeftest(fit_V_OBL, vcov. = NeweyWest(fit_V_OBL))[,4])



FIT_M_V_STEP_TABLE <- data.frame(coef=coef(step_V_OBL_model),
                                 cf=data.frame(coefci(step_V_OBL_model, vcov. = NeweyWest(step_V_OBL_model))[,1],
                                               coefci(step_V_OBL_model, vcov. = NeweyWest(step_V_OBL_model))[,2]),
                                 sd=sqrt(diag(NeweyWest(step_V_OBL_model))),
                                 t=coeftest(step_V_OBL_model, vcov. = NeweyWest(step_V_OBL_model))[,3],
                                 p=coeftest(step_V_OBL_model, vcov. = NeweyWest(step_V_OBL_model))[,4])

# MODELSELEKTIONS
data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

## RISIKOFRI
data_var_rf <- as.xts(merge(xts::lag.xts(data$NET_TB, -1, na.pad = FALSE), head(data, -1)))

fit_var_rf  <- lm(NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_rf)

step_var_rf_model <- step(fit_var_rf, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + TB, upper=fit_var_rf))

summary(step_var_rf_model)

prod_var_rf <- lm(NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_rf)

car::vif(prod_var_rf)

## AKTIER
data_var_akt <- as.xts(merge(xts::lag.xts(data$AKT, -1, na.pad = FALSE), head(data, -1)))

fit_var_akt  <- lm(AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_akt)

step_var_akt_model <- step(fit_var_akt, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + TB, upper=fit_var_akt))

summary(step_var_akt_model)

prod_var_akt <- lm(AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_akt)

car::vif(prod_var_akt)

## STATS
data_var_s <- as.xts(merge(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE), head(data, -1)))

fit_var_s  <- lm(S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_s)

step_var_s_model <- step(fit_var_s, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + TB, upper=fit_var_s))

summary(step_var_s_model)

prod_var_s <- lm(S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_s)

car::vif(prod_var_s)

## VIRKSOMHEDER
data_var_v <- as.xts(merge(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE), head(data, -1)))

fit_var_v  <- lm(V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_v)

step_var_v_model <- step(fit_var_v, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + TB, upper=fit_var_v))

summary(step_var_v_model)

prod_var_v <- lm(V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + BM + SMB + TB + Y_SPREAD, data=data_var_v)

car::vif(prod_var_v)

data_v <- as.ts(data)

# VECTOR AUTOREGRESSIVE MODEL
data_sub <- subset(data_v, select = -c(DP, EP, AKT_VAR, HML, T_SPREAD, C_SPREAD, D_SPREAD, FR))

var <- VAR(data_sub, p = 1)

summary(var)

coef <- Bcoef(var)
cov  <- summary(var)$covres
cor  <- summary(var)$corres

roots <- roots(var)

var.p <- VAR.Pope(data_sub, p = 1)
coef <- var.p$coef
cov <- var.p$sigu

# PORTEFØLJE ALLOKERING.
coef

Phi0_temp <- coef[ , ncol(coef)]
Phi1_temp <- coef[ , 1:(ncol(coef)-1)]

z  <- data_sub[ , 1:(ncol(data_sub))]
z  <- as.matrix(z)
z2 <- t(z)

port_calc <- function(gamma = 5, K = 100, n= 4){
  m <- nrow(Phi1_temp)
  
  H <- diag(x = 1, m, m)
  
  H1 <- H[1,]
  H1 <- t(as.matrix(H1))
  
  Hx <- H[2:n, ]
  Hx <- as.matrix(Hx)
  
  Phi0 <- Phi0_temp
  Phi0 <- as.matrix(Phi0)
  
  Phi1 <- Phi1_temp
  Phi1 <- as.matrix(Phi1)
  
  Phi0X <- Hx %*% Phi0
  Phi1X <- Hx %*% Phi1
  Phi01 <- H1 %*% Phi0
  Phi11 <- H1 %*% Phi1
  
  Sigma   <- cov
  SigmaXX <- Hx %*% Sigma %*% t(Hx)
  SigmaX2 <- diag(SigmaXX)
  SigmaX  <- Hx %*% Sigma
  SigmaX1 <- Hx %*% Sigma %*% t(H1)
  Sigma11 <- H1 %*% Sigma %*% t(H1)
  
  A0     <- matrix(0, nrow = n-1, ncol = 1)
  A1     <- matrix(0, nrow = n-1, ncol = m)
  B1     <- matrix(0, nrow = 1, ncol = m)
  B2     <- matrix(0, nrow = m, ncol = m)
  Lambda <- matrix(0, nrow = m, ncol = m)
  Gamma  <- matrix(0, nrow = m, ncol = m)
  Xi     <- matrix(0, nrow = m, ncol = m)
  
  for (i in 1:(K-1)) {
    if (i==1) {
      # A0(1).
      A0 <- (1/gamma) * solve(SigmaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1)
      
      # A1(1).
      A1 <- (1/gamma) * solve(SigmaXX) %*% Phi1X
      
      # B1(1).
      B1 <- Phi11 + t(A0) %*% (Phi1X - gamma * SigmaXX %*% A1) + t(Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1) %*% A1
      
      # B2(1).
      B2 <- t(A1) %*% (Phi1X - gamma/2 * SigmaXX %*% A1)
    } else {
      
      # Gem sidste periode A0 og A1.
      A0a <- A0
      A1a <- A1
      
      # Omega.
      Omega   <- solve(solve(Sigma)-2 * (1 - gamma) * B2)
      OmegaXX <- Hx %*% Omega %*% t(Hx)
      OmegaX1 <- Hx %*% Omega %*% t(H1)
      OmegaX  <- Hx %*% Omega
      
      # A0.
      A0 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * (OmegaX1 + OmegaX %*% (t(B1) + 2 * B2 %*% Phi0)))
      
      # A1.
      A1 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi1X + 2 * (1 - gamma) * OmegaX %*% B2 %*% Phi1)
      
      # Lambda.
      Lambda <- B2 %*% Omega %*% t(B2)
      
      # Gamma.
      Gamma <- 2 * B2 %*% Omega
      
      # Xi.
      Xi  <- t(Gamma)
      Xix <- Hx %*% Xi
      Xi1 <- H1 %*% Xi
      
      # Beregn B1 og B2 for næste periode.
      
      # B1.
      B1 <- Phi11 + t(A0a) %*% (Phi1X - (SigmaXX - (1 - gamma) * OmegaXX) %*% A1a) + t(Phi0X + 1/2 * SigmaX2 + (1 - gamma) * OmegaX1) %*% A1a + (B1 + 2 * t(Phi0) %*% B2) %*% Phi1 + (1 - gamma) * (B1 %*% t(OmegaX) + t(Phi0) %*% t(Xix)) %*% A1a + (1 - gamma) * (4 * t(Phi0) %*% Lambda + Xi1 + t(A0a) %*% Xix + B1 %*% Xi) %*% Phi1
      
      # B2.
      B2 <- t(A1a) %*% (Phi1X - gamma/2 * SigmaXX %*% A1a) + t(Phi1) %*% (B2 + 2 * (1 - gamma) * Lambda) %*% Phi1 + (1 - gamma) * t(Phi1) %*% t(Xix) %*% A1a
    }
  }
  
  if (K > 1) {
    # A0.
    A0 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * (OmegaX1 + OmegaX %*% (t(B1) + 2 * B2 %*% Phi0)))
    
    # A1.
    A1 <- solve(SigmaXX - (1 - gamma) * OmegaXX) %*% (Phi1X + 2 * (1 - gamma) * OmegaX %*% B2 %*% Phi1)
    
  } else {
    # A0.
    A0 <- 1/gamma * solve(SigmaXX) %*% (Phi0X + 1/2 * SigmaX2 + (1 - gamma) * SigmaX1)
    
    # A1.
    A1 <- 1/gamma * solve(SigmaXX) %*% Phi1X
  }
  
  alphafixedT <- matrix(0, nrow = n-1, ncol = ncol(z2))
  
  for (t in (2:ncol(z2))) {
    alphafixedT[, t] <- A0 + A1 %*% z2[, t-1]
  }
  
  MyopicD <- matrix(0, nrow = n-1, ncol = ncol(z2)) 
  
  for (K in (2:ncol(z2))) {
    MyopicD[, K] <- (1/gamma) * solve(SigmaXX) %*% (Hx %*% (Phi0 + Phi1 %*% z2[, K-1]) + 1/2 * SigmaX2 + (1-gamma) * SigmaX1)
  }
  
  HedgingD  <- alphafixedT - MyopicD
  
  TangencyP <- matrix(0, nrow = n-1, ncol = ncol(z2))
  
  for (t in (2:ncol(z2))){
    TangencyP[, t] <- solve(SigmaXX) %*% Hx %*% (Phi0 + Phi1 %*% z2[, t-1]) + 1/2 * SigmaX2
  }
  
  GMV             <- -solve(SigmaXX)%*%SigmaX1
  
  alphafixedT     <- rbind(1-colSums(alphafixedT), alphafixedT)
  
  res_alphafixedT <- matrix(0, nrow = n, ncol = ncol(z2))
  
  for (t in (2:ncol(z2))) {
    res_alphafixedT[, t] <- c(max(0, alphafixedT[1, t])/sum(max(0, alphafixedT[1,t]),
                                                            max(0, alphafixedT[2,t]),
                                                            max(0, alphafixedT[3,t]),
                                                            max(0, alphafixedT[4,t])),
                              
                              max(0, alphafixedT[2, t])/sum(max(0, alphafixedT[1,t]),
                                                            max(0, alphafixedT[2,t]),
                                                            max(0, alphafixedT[3,t]),
                                                            max(0, alphafixedT[4,t])),
                              
                              max(0, alphafixedT[3, t])/sum(max(0, alphafixedT[1,t]),
                                                            max(0, alphafixedT[2,t]),
                                                            max(0, alphafixedT[3,t]),
                                                            max(0, alphafixedT[4,t])),
                              
                              max(0, alphafixedT[4, t])/sum(max(0, alphafixedT[1,t]),
                                                            max(0, alphafixedT[2,t]),
                                                            max(0, alphafixedT[3,t]),
                                                            max(0, alphafixedT[4,t]))
    )
  }
  
  return(list(alphafixedT=alphafixedT,
              MyopicD=MyopicD,
              HedgingD=HedgingD,
              TangencyP=TangencyP,
              GMV=GMV,
              res_alphafixedT=res_alphafixedT)
  )
}

dynport <- port_calc(gamma = 5, K = 100, n = 4)

KOMP <- data.frame(Tan = c(1-sum(rowMeans(dynport$TangencyP)),
                         rowMeans(dynport$TangencyP)),
                   GMV = c(1-sum(dynport$GMV),
                         dynport$GMV),
                   Myopic = c(1-sum(rowMeans(dynport$MyopicD)),
                            rowMeans(dynport$MyopicD)),
                   IHD = c(1-sum(rowMeans(dynport$HedgingD)), rowMeans(dynport$HedgingD))
)

BIG_TABLE   <- matrix(0, nrow = 16, ncol = 8)
BIG_R_TABLE <- matrix(0, nrow = 16, ncol = 8)

gammas   <- c(2, 5, 10 , 20)
horizons <-  c(1, 4, 8, 20, 40, 60, 80, 100)
k        <- 1

for (i in c(1, 5, 9, 13)) {
  for (j in (1:length(horizons))) {
    BIG_TABLE[i:(i+3), j]   <- rowMeans(port_calc(gamma = gammas[k], K = horizons[j])$alphafixedT)
    BIG_R_TABLE[i:(i+3), j] <- rowMeans(port_calc(gamma = gammas[k], K = horizons[j])$res_alphafixedT)
  }
  
  k <- k+1
}

HOR_ANA_VAL <- c(1, 3:100)

HOR_ANA   <- matrix(0, nrow = 4, ncol = length(HOR_ANA_VAL))
HOR_ANA_R <- matrix(0, nrow = 4, ncol = length(HOR_ANA_VAL))

for (i in (1:length(HOR_ANA_VAL))) {
  a <- port_calc(gamma = 5, K = HOR_ANA_VAL[i])
  
  HOR_ANA[, i]   <- rowMeans(a$alphafixedT)
  HOR_ANA_R[, i] <- rowMeans(a$res_alphafixedT)
}

GAM_ANA_VAL <- seq(1, 50, 0.1)

GAM_ANA     <- matrix(0, nrow = 4, ncol = length(GAM_ANA_VAL))
GAM_ANA_R   <- matrix(0, nrow = 4, ncol = length(GAM_ANA_VAL))
GAM_ANA_MYO <- matrix(0, nrow = 4, ncol = length(GAM_ANA_VAL))

for (i in (1:length(GAM_ANA_VAL))) {
  b <- port_calc(gamma = GAM_ANA_VAL[i], K = 100)
  
  GAM_ANA[, i]     <- rowMeans(b$alphafixedT)
  GAM_ANA_R[, i]   <- rowMeans(b$res_alphafixedT)
  
  GAM_ANA_MYO[, i] <- rowMeans(rbind(1-colSums(b$MyopicD), b$MyopicD))
}
