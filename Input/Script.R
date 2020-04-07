library(quantmod)

# LOAD DATA.
input <- zoo::read.csv.zoo('Input/Samling.csv', format = "%Y%m%d")
input <- as.xts(input)

# MANIPULATING DATA.

# AKTIVKLASSER
INF <- log(1+input$cpiret)

TB  <- log(1+input$t90ret)
NET_TB <- TB - INF
AKT    <- log(1 + input$vwretd) - TB
S_OBL  <- log(1 + input$b10ret) - TB
V_OBL  <- log(1 + input$corpr) - TB

# FINANSIELLE INDIKATORER
N <- length(input$b10ret)

div <- (1 + as.numeric(input$vwretd[2:N])) * as.numeric(input$vwindx[1:(N-1)]) - as.numeric(input$vwindx[2:N])

logDiv <- c()

for (i in c(4:length(div))) {
  logDiv[i - 3] <- log(sum(div[(i-3):i]))
}

divPrice <- logDiv - log(as.numeric(input$vwindx[5:N]))

DP <- xts::xts(divPrice, index(input$vwretd[5:N]))
BM <- input$bm
AKT_VAR <- log(1 + input$svar)

# TERM STRUCTURE
T_SPREAD <- input$t10y/100 - input$t90y/100
Y_SPREAD <- input$t5y/100 - input$t90y/100
C_SPREAD <- input$baa/100 - input$t90y/100
D_SPREAD <- input$baa/100 - input$aaa/100

# MAKROØKONOMISKE
FR <- input$fedfund/100

# JUSTERING
TB_J <- TB + (var(TB, na.rm = TRUE)/2)
NET_TB_J <- NET_TB + (var(NET_TB, na.rm = TRUE)/2)
AKT_J <- AKT + (var(AKT, na.rm = TRUE)/2)
S_OBL_J <- S_OBL + (var(S_OBL, na.rm = TRUE)/2)
V_OBL_J <- V_OBL + (var(V_OBL, na.rm = TRUE)/2)

data <- xts::xts(merge(NET_TB, AKT, S_OBL, V_OBL, DP, BM, AKT_VAR, T_SPREAD, Y_SPREAD, C_SPREAD, D_SPREAD, FR), index(NET_TB))

data <- na.omit(data)

# KITCHEN SINK
fit_AKT <- lm(lag(AKT, 1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)
fit_S_OBL <- lm(lag(S_OBL, 1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)
fit_V_OBL <- lm(lag(V_OBL, 1) ~ DP + BM + AKT_VAR + T_SPREAD + Y_SPREAD + C_SPREAD + D_SPREAD + FR, data = data)

m_fit <- lm(cbind(AKT, S_OBL, V_OBL) ~ lag(DP, 1) + lag(BM, 1) + lag(AKT_VAR, 1) + lag(T_SPREAD, 1) + lag(Y_SPREAD, 1) + lag(C_SPREAD, 1) + lag(D_SPREAD, 1) + lag(FR, 1), data = data)


library(vars)

var <- VAR(data, p = 1)




# LOAD DATA

rf <- read.csv('Input/Riskfree.csv')
stock <- read.csv('Input/Stock.csv')
yield <- read.csv('Input/Yield.csv')

rf

(mean(log(1+rf$t90ret), na.rm = TRUE)+(var(log(1+rf$t90ret), na.rm=TRUE)/2))*12*100

# DATA MANIPULATION

## RISKFREE
date <- zoo::as.Date(as.character(rf$caldt), "%Y%m%d")
tbill90 <- as.numeric(rf$t90ret)

rf <- xts::xts(order.by = date, tbill90)

log_rf <- log(1+rf)

## STOCK
stockret <- as.numeric(stock$vwretd)

stock <- xts::xts(order.by = date, stockret)

log_stock <- log(1+stock)

plot(log(1+stock$vwretd), type = 'l')
