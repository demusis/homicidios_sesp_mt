library(dplyr)
library(factoextra)
library(forecast)
library(imputeTS)
library(lubridate)
library(xts)
library(readr)
library(readxl)
library(writexl)
library(zoo)

setwd("~/INMET")

INMET_Cuiaba <- read.csv2("~/INMET/INMET_Cuiaba.CSV")

# Função para substituir -9999 por NA
substituir_valores <- function(x) {
  x[x == -9999] <- NA
  return(x)
}

# Aplicando a função a todas as colunas e linhas do dataframe
INMET_Cuiaba <- INMET_Cuiaba %>%
  mutate(across(everything(), ~substituir_valores(.)))

# Combinar colunas de data e hora em uma única coluna POSIXct
INMET_Cuiaba$datahora <- as.POSIXct(paste(INMET_Cuiaba$DATA..YYYY.MM.DD., 
                                          INMET_Cuiaba$HORA..UTC.), 
                                    format = "%d/%m/%Y %H:%M")

# Seleciona as variáveis a serem analisadas
INMET_Cuiaba <- INMET_Cuiaba[, c(3, 8, 16, 20)]

# Cria objeto xts com preenchimento de falhas
inicio <- first(INMET_Cuiaba$datahora)
fim <- last(INMET_Cuiaba$datahora)
  
# Gerando sequência de horas
horas <- seq(from = inicio, to = fim, by = "hour")
  
# Criando um objeto do tipo ts (xts)
INMET_xts <- xts(INMET_Cuiaba[, -ncol(INMET_Cuiaba)], 
                 order.by = horas)
  
# Preencher valores NA usando filtro de Kalman
INMET_xts <- na_kalman(INMET_xts,
                       model = "auto.arima", 
                       smooth = TRUE)

# Agregando as colunas em médias mensais
INMET_mensal_xts <- apply.monthly(INMET_xts, FUN = mean)

# Carregar dados de homicídios em Cuiabá + VG 
homicidios <- read_excel("Homicidios dolosos (CBA e VG) 2012-2022.xlsx")
colnames(homicidios) <- c("mes", "ano", 'homicidios')

meses <- c("JANEIRO", "FEVEREIRO", "MARÇO", "ABRIL", "MAIO", "JUNHO", "JULHO", "AGOSTO", "SETEMBRO", "OUTUBRO", "NOVEMBRO", "DEZEMBRO")
homicidios$mes <- match(homicidios$mes, meses)

# Combinar as colunas de ano e mês para criar uma coluna de datas
homicidios$data <- as.Date(paste(homicidios$ano, 
                                 homicidios$mes, 
                                 1, 
                                 sep = "-"))

# Converter o dataframe em um objeto xts
homicidios_mensal_xts <- as.xts(homicidios$homicidios, 
                                order.by = as.Date(homicidios$data))


# ---------------------------------------------

resultados <- data.frame(covariaveis = character(), 
                         aic = numeric(), 
                         aicc = numeric(), 
                         bic = numeric(), 
                         mase = numeric(), 
                         stringsAsFactors = FALSE)

# Dividir os dados em conjuntos de treinamento (70%) e teste (30%)
nh <- nrow(homicidios_mensal_xts)
n_treino <- floor(0.7 * nh)
homicidios_mensal_xts_treino <- homicidios_mensal_xts[1:n_treino]
homicidios_mensal_xts_teste <- homicidios_mensal_xts[(n_treino + 1):nh]

ajuste <- function(aux = NA) {
  # Dividir os dados em conjuntos de treinamento (70%) e teste (30%)
  INMET_treino <- as.matrix(INMET_mensal_xts[, aux])[1:n_treino,]
  INMET_teste <- as.matrix(INMET_mensal_xts[, aux])[(n_treino=1):nh,]
  
  # Testar os modelos
  if (any(is.na(aux))) {
    res <- auto.arima(
      homicidios_mensal_xts,
      lambda = "auto",
      stepwise = FALSE,
      approximation = FALSE,
      seasonal = TRUE,
      allowdrift = TRUE,
      trace = FALSE
    )  
    res2 <- auto.arima(
      homicidios_mensal_xts_treino,
      lambda = "auto",
      stepwise = FALSE,
      approximation = FALSE,
      seasonal = TRUE,
      allowdrift = TRUE,
      trace = FALSE
    )
    
  } else {
    res <- auto.arima(
      homicidios_mensal_xts,
      xreg = as.matrix(INMET_mensal_xts[, aux]),
      lambda = "auto",
      stepwise = FALSE,
      approximation = FALSE,
      seasonal = TRUE,
      allowdrift = TRUE,
      trace = FALSE
    )
    res2 <- auto.arima(
      homicidios_mensal_xts_treino,
      xreg = INMET_treino,
      lambda = "auto",
      stepwise = FALSE,
      approximation = FALSE,
      seasonal = TRUE,
      allowdrift = TRUE,
      trace = FALSE
    )
  }
  
  df2 <- data.frame(matrix(nrow = 1,
                           ncol = ncol(resultados)))
  names(df2) <- names(resultados)
  df2[1, 1] <- paste(aux,
                     collapse = ' ')
  df2[1, 2] <- res$aic
  df2[1, 3] <- res$aicc
  df2[1, 4] <- res$bic
  
  # Previsão do conjunto de teste
  if (any(is.na(aux))) {
    prev_teste <- forecast(res2, 
                           h = nrow(homicidios_mensal_xts_teste))
  } else {
    prev_teste <- forecast(res2, 
                           xreg = INMET_teste, 
                           h = nrow(homicidios_mensal_xts_teste))
  }
  
  # Calcular o MASE
  df2[1, 5] <- accuracy(prev_teste, homicidios_mensal_xts_teste)[2,6]
  
  return(df2)
}

# Obter todas as combinações possíveis das colunas
num_colunas <- ncol(INMET_mensal_xts)

resultados <- ajuste(aux=NA)
# Loop para iterar todas as combinações possíveis de um a n
for (i in 1:num_colunas) {
  combinacoes <- combn(1:num_colunas, i)
  for (j in 1:ncol(combinacoes)) {
    resultados <- rbind(resultados, 
                        ajuste(combinacoes[, j]))
  }
}

resultados <- na.omit(resultados)
write_xlsx(resultados, 'resultados.xlsx')

homicidios_final <- auto.arima(
  homicidios_mensal_xts,
  xreg = as.matrix(INMET_mensal_xts[, 2]),
  lambda = "auto",
  stepwise = FALSE,
  approximation = FALSE,
  seasonal = TRUE,
  allowdrift = TRUE,
  trace = FALSE
)

regressor_final <- auto.arima(
  INMET_mensal_xts[,2],
  lambda = "auto",
  stepwise = FALSE,
  approximation = FALSE,
  seasonal = TRUE,
  allowdrift = TRUE,
  trace = FALSE
)

regressor_final_est <- as.matrix(forecast(regressor_final, h=3)$mean[1:3])

colnames(regressor_final_est) <- colnames(as.matrix(INMET_mensal_xts[, 2]))
rownames(regressor_final_est) <- c("2023-03-31 23:00:00", 
                                   "2023-04-30 23:00:00",
                                   "2023-05-31 23:00:00")


# Predições
predicoes <- forecast(homicidios_final, xreg = regressor_final_est, h=1)

autoplot(predicoes)

# Plote o gráfico com a série original e estimativas

homicidios_estimados <- xts(homicidios_final$fitted, 
                            order.by = homicidios$data)

homicidios_est_obs <- cbind(homicidios_mensal_xts, homicidios_estimados)
plot(homicidios_mensal_xts, 
     main = "Série Original e Estimativas", 
     xlab = "Data", 
     ylab = "Homicídios", 
     col = "blue", 
     lwd = 2, 
     type = "l", 
     xaxt = "n")
lines(homicidios_estimados, 
      col = "red", 
      lwd = 2, 
      lty = 2)

