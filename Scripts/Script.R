library(quantmod)
library(moments)
library(VAR.etp)
library(vars)
library(MASS)

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
                          quantile=t(DATA_QUANTILE))
DESCRIPTIVE$sr[1] <- NA

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
                          quantile=t(DATA_T_QUANTILE))


DESCRIPTIVE_T$mean[1:2] <- ADJ_DATA_T_MEAN[1:2]
DESCRIPTIVE_T$mean[7] <- ADJ_DATA_T_MEAN[3]

data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, EP, BM, AKT_VAR, HML, SMB, TB, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

# UNIVARIATE
FIT_AKT <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$AKT,-1, na.pad = FALSE) ~ head(data[,x], -1), data=data))
FIT_AKT_SUMMARY <- lapply(FIT_AKT, summary)

FIT_TABLE <- data.frame(coef=sapply(FIT_AKT, coef)[2,],
                        sd=sapply(FIT_AKT_SUMMARY, coef)[4,],
                        t=sapply(FIT_AKT_SUMMARY, coef)[6,],
                        p=sapply(FIT_AKT_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_AKT), function(x) FIT_AKT_SUMMARY[[x]]$r.squared)))

FIT_S_OBL <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE) ~ head(data[,x], -1), data=data))
FIT_S_OBL_SUMMARY <- lapply(FIT_S_OBL, summary)

FIT_S_TABLE <- data.frame(coef=sapply(FIT_S_OBL, coef)[2,],
                        sd=sapply(FIT_S_OBL_SUMMARY, coef)[4,],
                        t=sapply(FIT_S_OBL_SUMMARY, coef)[6,],
                        p=sapply(FIT_S_OBL_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_S_OBL), function(x) FIT_S_OBL_SUMMARY[[x]]$r.squared)))

FIT_V_OBL <- lapply(5:ncol(data), function(x) lm(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE) ~ head(data[,x], -1), data=data))
FIT_S_OBL_SUMMARY <- lapply(FIT_V_OBL, summary)

FIT_V_TABLE <- data.frame(coef=sapply(FIT_V_OBL, coef)[2,],
                        sd=sapply(FIT_S_OBL_SUMMARY, coef)[4,],
                        t=sapply(FIT_S_OBL_SUMMARY, coef)[6,],
                        p=sapply(FIT_S_OBL_SUMMARY, coef)[8,],
                        rsq=unlist(lapply(1:length(FIT_V_OBL), function(x) FIT_S_OBL_SUMMARY[[x]]$r.squared)))

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

# MULTIPLE REGRESSION 

data_mult_akt <- as.xts(merge(xts::lag.xts(data$AKT, -1, na.pad = FALSE), head(data, -1)))
data_mult_akt <- subset(data_mult_akt, select = -AKT.1)

data_mult_s <- as.xts(merge(xts::lag.xts(data$S_OBL, -1, na.pad = FALSE), head(data, -1)))
data_mult_s <- subset(data_mult_s, select = -S_OBL.1)

data_mult_v <- as.xts(merge(xts::lag.xts(data$V_OBL, -1, na.pad = FALSE), head(data, -1)))
data_mult_v <- subset(data_mult_v, select = -V_OBL.1)

# KITCHEN SINK
fit_AKT <- lm(AKT ~ ., data=data_mult_akt)

test <- lm(AKT~NET_TB + S_OBL + V_OBL + TB + BM + T_SPREAD + FR, data = data_mult_akt)

step_AKT_model <- step(fit_AKT, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=AKT ~ NET_TB + S_OBL + V_OBL + TB, upper=fit_AKT))

summary(step_AKT_model)

fit_S_OBL <- lm(S_OBL ~ ., data=data_mult_s)

step_S_OBL_model <- step(fit_S_OBL, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=S_OBL ~ NET_TB + AKT + V_OBL + TB, upper=fit_S_OBL))

summary(step_S_OBL_model)

fit_V_OBL <- lm(V_OBL ~ ., data=data_mult_v)

step_V_OBL_model <- step(fit_V_OBL, k = qchisq(0.05,df = 1, lower.tail = F), scope = list(lower=V_OBL ~ NET_TB + AKT + S_OBL + TB, upper=fit_V_OBL))

summary(step_V_OBL_model)


# m_fit <- lm(cbind(AKT, S_OBL, V_OBL) ~ lag(DP, 1) + lag(BM, 1) + lag(AKT_VAR, 1) + lag(T_SPREAD, 1) + lag(Y_SPREAD, 1) + lag(C_SPREAD, 1) + lag(D_SPREAD, 1) + lag(FR, 1), data = data)
# 

# VAR
data_v <- as.ts(data)
# 
data_sub <- subset(data_v, select = -c(DP, EP, AKT_VAR, HML, Y_SPREAD, C_SPREAD, D_SPREAD, FR))
# 
var <- VAR(data_sub, p = 1)

coef <- Bcoef(var)
coef1 <- coef(var)
cov <- summary(var)$covres
cor <- summary(var)$corres

# fileName_var <- 'varAM.xlsx'
# 
# excel_var <- createWorkbook(fileName_var)
# 
# addWorksheet(excel_var,'Coef')
# addWorksheet(excel_var,'CoefStdTp')
# addWorksheet(excel_var,'Cov')
# addWorksheet(excel_var,'Corr')
# 
# 
# writeData(excel_var, sheet = 1, coef)
# writeData(excel_var, sheet = 2, coef1)
# writeData(excel_var, sheet = 3, cov)
# writeData(excel_var, sheet = 4, cor)
# 
# saveWorkbook(excel_var, fileName_var, overwrite = T)

# VAR.select(data_sub,p=1)


# var <- VAR.etp::VAR.est(data, p=1)
# 
# pope_var <- VAR.Pope(data, p = 1)
