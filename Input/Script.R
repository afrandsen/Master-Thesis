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
