library(quantmod)
library(moments)
library(VAR.etp)
library(vars)
library(MASS)
library(tseries)

# LOAD DATA.
input <- zoo::read.csv.zoo('C:/Users/AKF/Documents/Matematik-Okonomi/5. Ar/2. Semester/Master-Thesis/Input/Samling.csv', format = "%Y%m%d")
#input <- as.xts(input)

# MANIPULATING DATA.

# AKTIVKLASSER
INF <- log(1+input$cpiret)
TB  <- log(1+input$t90ret)
NET_TB <- TB - INF
AKT    <- log(1 + input$vwretd) - TB
S_OBL  <- log(1 + input$b10ret) - TB
V_OBL  <- log(1 + input$corpr) - TB

# FINANSIELLE INDIKATORER
DP <- input$dp
EP <- log(as.numeric(input$vwindx)) - log(input$earn)

BM <- input$bm
AKT_VAR <- log(1 + input$svar)
SMB <- input$smb
HML <- input$hml
EP_J <- EP + (var(EP, na.rm = TRUE)/2)
DP_J <- DP + (var(DP, na.rm = TRUE)/2)


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
DATA <- data.frame(NET_TB, AKT, S_OBL, V_OBL)
DATA_J <- data.frame(NET_TB_J, AKT_J, S_OBL_J, V_OBL_J)

ADJ_DATA_MEAN <- colMeans(DATA_J)
DATA_STD <- apply(DATA, 2, sd)
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
DATA_T <- xts::xts(merge(DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(TB))

DATA_T_J <- xts::xts(merge(DP_J, EP_J, TB_J), index(TB_J))

ADJ_DATA_T_MEAN <- colMeans(DATA_T_J)
DATA_T_MEAN <- colMeans(DATA_T)
DATA_T_STD <- apply(DATA_T, 2, sd)
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
DESCRIPTIVE_T$mean[7] <- ADJ_DATA_T_MEAN[3]

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







# VAR

data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

## RISIKOFRI
data_var_rf <- as.xts(merge(xts::lag.xts(data$NET_TB, -1, na.pad = FALSE), head(data, -1)))

fit_var_rf  <- lm(NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_rf)

step_var_rf_model <- step(fit_var_rf, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + TB, upper=fit_var_rf))

summary(step_var_rf_model)

prod_var_rf <- lm(NET_TB ~ NET_TB.1 + AKT + S_OBL + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_rf)

car::vif(prod_var_rf)

step_table_rf <- data.frame(dp=c(-2516.1, -2519.7, -2523.5, -2525.3, NA),
                            pe=c(-2520.2, NA, NA, NA, NA),
                            bm=c(-2497.1, -2499.3, -2502.5, -2505.2, -2499.1),
                            avar=c(-2515.0, -2518.0, -2521.5, -2520.9,  -2522.5),
                            hml=c(-2517.7, -2521.5, -2525.0, NA, NA),
                            smb=c(-2513.6, -2517.3, -2521.0, -2522.7, -2521.7),
                            ts=c(-2519.7, -2523.5, NA, NA, NA),
                            ys=c(-2514.8, -2518.7, -2514.5, -2516.3, -2510.2),
                            cs=c(-2517.9, -2521.7, -2511.8, -2514.2, -2507.9),
                            ds=c(-2503.7, -2507.5, -2495.9, -2498.9, -2497.0),
                            fr=c(-2516.4, -2520.2, -2523.7, -2525.1, -2524.7))
step_table_rf <- cbind(min=apply(step_table_rf, 1, FUN=min, na.rm=TRUE), step_table_rf)


## AKTIER
data_var_akt <- as.xts(merge(xts::lag.xts(data$AKT, -1, na.pad = FALSE), head(data, -1)))

fit_var_akt  <- lm(AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_akt)

step_var_akt_model <- step(fit_var_akt, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + TB, upper=fit_var_akt))

summary(step_var_akt_model)

prod_var_akt <- lm(AKT ~ NET_TB + AKT.1 + S_OBL + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_akt)

car::vif(prod_var_akt)

step_table_akt <- data.frame(dp=c(-1264.3, NA, NA, NA, NA, NA, NA, NA, NA),
                            pe=c(-1264.0, -1267.7, -1271.3, NA, NA, NA, NA, NA, NA),
                            bm=c(-1264.1, -1267.3, -1269.8, -1269.4, -1271.3, -1271.7, -1272.7,-1271.8,  -1275.1 ),
                            avar=c(-1262.8, -1266.5, -1269.8,-1273.5, NA, NA, NA, NA, NA ),
                            hml=c(-1262.3, -1266.0, -1269.6,-1272.9, -1274.0, NA, NA, NA, NA ),
                            smb=c(-1259.4,-1263.0, -1266.5, -1269.7,  -1272.5, -1272.8,-1272.2, -1273.5, -1276.4 ),
                            ts=c(-1260.2,-1263.1, -1266.0, -1268.7, -1271.7, -1273.2, -1276.3, NA, NA ),
                            ys=c(-1261.3, -1264.9, -1268.7, -1271.6, -1273.9, -1274.9, NA, NA, NA),
                            cs=c(-1261.8, -1264.5, -1264.2, -1266.7, -1270.5, -1272.3,-1275.9, -1279.8, NA ),
                            ds=c( -1264.2,  -1267.9, NA, NA, NA, NA, NA, NA, NA),
                            fr=c(-1258.8, -1262.4, -1266.1,-1269.5, -1272.2, -1272.8, -1274.0, -1272.1, -1275.5))
step_table_akt <- cbind(min=apply(step_table_akt, 1, FUN=min, na.rm=TRUE), step_table_akt)

## STATS
data_var_s <- as.xts(merge(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE), head(data, -1)))

fit_var_s  <- lm(S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_s)

step_var_s_model <- step(fit_var_s, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + TB, upper=fit_var_s))

summary(step_var_s_model)

prod_var_s <- lm(S_OBL ~ NET_TB + AKT + S_OBL.1 + V_OBL + BM + SMB + TB + Y_SPREAD, data=data_var_s)

car::vif(prod_var_s)

step_table_s <- data.frame(dp=c(-1710.4,-1714.2,  -1717.8, -1721.7, -1721.9, -1725.8, NA ),
                             pe=c(-1711.0, -1714.8, -1718.6,-1722.2, NA, NA, NA  ),
                             bm=c(-1711.9,-1715.7, NA, NA, NA, NA, NA ),
                             avar=c(-1705.2,-1708.9, -1712.6, -1716.3, -1717.0, -1717.9, -1720.2 ),
                             hml=c( -1712.1, NA, NA, NA, NA, NA, NA),
                             smb=c(-1708.0, -1711.7, -1714.9,-1718.7, -1721.7, -1723.9,  -1725.1 ),
                             ts=c(-1711.9, -1715.7, -1719.2,  NA, NA, NA, NA),
                             ys=c(-1707.8,-1711.6,  -1715.1,-1694.6,-1695.9,  -1691.6,-1695.2  ),
                             cs=c(-1702.1, -1705.9, -1709.3,-1696.9, -1699.6, -1681.0, -1684.1 ),
                             ds=c(-1711.0,-1714.8, -1718.1, -1721.9, -1724.8, NA, NA ),
                             fr=c(-1687.6,-1691.4,  -1694.2, -1697.9, -1701.3,-1702.7, -1704.9  ))
step_table_s <- cbind(min=apply(step_table_s, 1, FUN=min, na.rm=TRUE), step_table_s)

## VIRKSOMHEDER
data_var_v <- as.xts(merge(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE), head(data, -1)))

fit_var_v  <- lm(V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + DP + EP + BM + AKT_VAR + HML + SMB + TB + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data=data_var_v)

step_var_v_model <- step(fit_var_v, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + TB, upper=fit_var_v))

summary(step_var_v_model)

prod_var_v <- lm(V_OBL ~ NET_TB + AKT + S_OBL + V_OBL.1 + BM + SMB + TB + Y_SPREAD, data=data_var_v)

car::vif(prod_var_v)

step_table_v <- data.frame(dp=c(-1607.8, -1611.5, -1615.4, -1618.2,-1620.0, -1623.1,NA ),
                             pe=c(-1608.6, -1612.4, NA, NA, NA, NA, NA),
                             bm=c(-1608.6, -1612.4,  -1616.2, NA, NA ,NA , NA),
                             avar=c(-1603.4, -1607.3, -1610.8,  -1614.5, -1618.4, -1619.8, -1621.2),
                             hml=c( -1607.6, -1611.5, -1615.3, -1619.1,-1622.8, NA, NA ),
                             smb=c(-1603.0, -1606.7,-1610.6,  -1614.2, -1618.0, -1621.1, -1621.4),
                             ts=c(-1608.6, NA, NA, NA, NA, NA, NA),
                             ys=c(-1606.7,-1604.2, -1607.8, -1611.6, -1607.3, -1610.6, -1614.2 ),
                             cs=c(-1605.8, -1602.3, -1606.1,-1609.8,  -1590.1, -1593.0,  -1596.7 ),
                             ds=c(-1608.4,-1612.2, -1616.0, -1619.8, NA, NA, NA ),
                             fr=c(-1584.4,-1588.2,  -1592.0, -1595.4,-1599.2,   -1602.2,-1603.8  ))
step_table_v <- cbind(min=apply(step_table_v, 1, FUN=min, na.rm=TRUE), step_table_v)





data_v <- as.ts(data)
# 
data_sub <- subset(data_v, select = -c(DP, EP, AKT_VAR, HML, T_SPREAD, C_SPREAD, D_SPREAD, FR))
# 
var <- VAR(data_sub, p = 1)
var.p <- VAR.Pope(data_sub, p = 1)

coef <- Bcoef(var)

coef1 <- coef(var)
cov <- summary(var)$covres
cor <- summary(var)$corres

test <- summary(var)

VAR_TABLE <- data.frame(coef,
                        rsq=unlist(lapply(1:8, function(x) test$varresult[[x]]$adj.r.squared)))

VAR_TABLE <- format(round(VAR_TABLE, digits = 3), nsmall = 3)

for (j in c(1:8)) {
  for (i in c(1:9)){
    if (abs(test$varresult[[j]]$coefficients[,3][i])>=qnorm(0.975)) {
      VAR_TABLE[j,i] <- paste0("\\textbf{", VAR_TABLE[j,i],"}")
    }
    
  }
}

j <- 0

for (i in c(1:8)) {

  VAR_TABLE <- rbind(VAR_TABLE[(1):(1+j),],
                     format(round(var.p$coef[i,],3), nsmall = 3),
                     VAR_TABLE[-(1:(1+j)),])

  VAR_TABLE <- rbind(VAR_TABLE[1:(2+j),],
                     format(round(test$varresult[[i]]$coefficients[,3],3), nsmall = 3),
                     VAR_TABLE[-(1:(2+j)),])

  j <- j + 3
}

j <- 1

for (i in seq(2,23,3)) {
  
  VAR_TABLE[i, 10] <- format(round(test$varresult[[j]]$fstatistic[1], 3), nsmall = 3)
  
  j <- j + 1
}

for (i in seq(1,24,3)){

VAR_TABLE[i+1,1:9] <- gsub("\\s", "", paste0("(", VAR_TABLE[i+1,], ")"))

VAR_TABLE[i+2,1:9] <- gsub("\\s", "", paste0("[", VAR_TABLE[i+2,], "]"))

VAR_TABLE[i+1,10] <- gsub("\\s", "", paste0("$\\langle", VAR_TABLE[i+1,10], "\\rangle$"))

}

VAR_TABLE <- cbind(NAME=c("$r_t^{\\text{rf}}$","","", "$rx_t^{\\text{a}}$","","", "$rx_t^{\\text{s}}$","","", "$rx_t^{\\text{v}}$","","", '$x_t^{\\text{bm}}$',"","", '$x_t^{\\text{smb}}$',"","", '$x_t^{\\text{b}}$',"","", '$x_t^{\\text{ys}}$',"",""), VAR_TABLE)

VAR_TABLE <- VAR_TABLE[,c(1,10,2,3,4,5,6,7,8,9,11)]

for (i in seq(1,24,3)){
  
  VAR_TABLE[i+2,11] <- NA
}

coef <- var.p$coef
cov <- var.p$sigu

recessions.df = read.table(textConnection(
  "Peak, Trough
1857-06-01, 1858-12-01
1860-10-01, 1861-06-01
1865-04-01, 1867-12-01
1869-06-01, 1870-12-01
1873-10-01, 1879-03-01
1882-03-01, 1885-05-01
1887-03-01, 1888-04-01
1890-07-01, 1891-05-01
1893-01-01, 1894-06-01
1895-12-01, 1897-06-01
1899-06-01, 1900-12-01
1902-09-01, 1904-08-01
1907-05-01, 1908-06-01
1910-01-01, 1912-01-01
1913-01-01, 1914-12-01
1918-08-01, 1919-03-01
1920-01-01, 1921-07-01
1923-05-01, 1924-07-01
1926-10-01, 1927-11-01
1929-08-01, 1933-03-01
1937-05-01, 1938-06-01
1945-02-01, 1945-10-01
1948-11-01, 1949-10-01
1953-07-01, 1954-05-01
1957-08-01, 1958-04-01
1960-04-01, 1961-02-01
1969-12-01, 1970-11-01
1973-11-01, 1975-03-01
1980-01-01, 1980-07-01
1981-07-01, 1982-11-01
1990-07-01, 1991-03-01
2001-03-01, 2001-11-01
2007-12-01, 2009-06-01"), sep=',',
  colClasses=c('Date', 'Date'), header=TRUE)

recessions.trim = subset(recessions.df, Peak >= min(index(data)) )

# a: AKT.1 BM AKT_VAR SMB TB Y_SPREAD
# s: S_OBL.1 AKT_VAR SMB TB Y_SPREAD
# v: V_OBL.1 AKT_VAR SMB TB Y_SPREAD

# data_test <- subset(data_v, select = c(BM, AKT_VAR, SMB, TB, Y_SPREAD))
# 
# test <- lm(AKT ~ AKT.1 + BM + AKT_VAR + SMB + TB + Y_SPREAD, data=data_var_akt)
# extractAIC(step_var_akt_model)
# extractAIC(test)
