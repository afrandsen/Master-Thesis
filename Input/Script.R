library(quantmod)
library(moments)

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
EP <- log(input$earn) - log(as.numeric(input$vwindx))

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
                          quantile=t(DATA_QUANTILE))
DESCRIPTIVE$sr[1] <- NA

## TILSTANDSVARIABLE
DATA_T <- xts::xts(merge(TB, EP, DP, BM, AKT_VAR, SMB, HML, FR, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD), index(TB))

DATA_T_J <- xts::xts(merge(TB_J, EP_J, DP_J), index(TB_J))

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
                          quantile=t(DATA_T_QUANTILE))


DESCRIPTIVE_T$mean[1:3] <- ADJ_DATA_T_MEAN

data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, TB, EP, DP, BM, AKT_VAR, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR, HML, SMB), index(NET_TB))

# UNIVARIATE
FIT_AKT <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$AKT,-1) ~ data[,x], data=data))
FIT_AKT_SUMMARY <- lapply(FIT_AKT, summary)

FIT_TABLE <- data.frame(coef=sapply(FIT_AKT, coef)[2,],
                        sd=sapply(FIT_AKT_SUMMARY, coef)[4,],
                        t=sapply(FIT_AKT_SUMMARY, coef)[6,],
                        p=sapply(FIT_AKT_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_AKT), function(x) FIT_AKT_SUMMARY[[x]]$r.squared)))

FIT_S_OBL <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$S_OBL, -1) ~ data[,x], data=data))
FIT_S_OBL_SUMMARY <- lapply(FIT_S_OBL, summary)

FIT_S_TABLE <- data.frame(coef=sapply(FIT_S_OBL, coef)[2,],
                        sd=sapply(FIT_S_OBL_SUMMARY, coef)[4,],
                        t=sapply(FIT_S_OBL_SUMMARY, coef)[6,],
                        p=sapply(FIT_S_OBL_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_S_OBL), function(x) FIT_S_OBL_SUMMARY[[x]]$r.squared)))

FIT_V_OBL <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$V_OBL, -1) ~ data[,x], data=data))
FIT_S_OBL_SUMMARY <- lapply(FIT_V_OBL, summary)

FIT_V_TABLE <- data.frame(coef=sapply(FIT_V_OBL, coef)[2,],
                        sd=sapply(FIT_S_OBL_SUMMARY, coef)[4,],
                        t=sapply(FIT_S_OBL_SUMMARY, coef)[6,],
                        p=sapply(FIT_S_OBL_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_V_OBL), function(x) FIT_S_OBL_SUMMARY[[x]]$r.squared)))

# MULTIPLE REGRESSION 






# # KITCHEN SINK
# fit_AKT <- lm(lag(AKT, -1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)
# fit_S_OBL <- lm(lag(S_OBL, -1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)
# fit_V_OBL <- lm(lag(V_OBL, -1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)
# 
# m_fit <- lm(cbind(AKT, S_OBL, V_OBL) ~ lag(DP, 1) + lag(BM, 1) + lag(AKT_VAR, 1) + lag(T_SPREAD, 1) + lag(Y_SPREAD, 1) + lag(C_SPREAD, 1) + lag(D_SPREAD, 1) + lag(FR, 1), data = data)
# 
# 
# library(vars)
# 
# var <- VAR(data, p = 1, ic = 'AIC')
