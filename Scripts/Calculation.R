library(quantmod)

input <- zoo::read.csv.zoo('C:/Users/AKF/Documents/Matematik-Okonomi/5. Ar/2. Semester/Master-Thesis/Input/test.csv', format = "%Y%m%d")
#input <- as.xts(input)

N <- length(input$vwretd)
div <- (1 + as.numeric(input$vwretd[2:N])) * as.numeric(input$vwindx[1:(N-1)]) - as.numeric(input$vwindx[2:N])
logDiv <- c()
for (i in c(4:length(div))) {
  logDiv[i - 3] <- log(sum(div[(i-3):i]))
}
divPrice <- logDiv - log(as.numeric(input$vwindx[5:N]))

DP <- xts::xts(divPrice, index(input$vwretd[5:N]))
DP <- as.zoo(DP)

DP

write.csv(DP, 'C:/Users/AKF/Documents/Matematik-Okonomi/5. Ar/2. Semester/Master-Thesis/Input/ny.csv')
