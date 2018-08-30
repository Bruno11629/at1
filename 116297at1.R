#O codigo presume que a pasta at1, possui o mesmo conteudo 
#e mesma ordenação daquela disponivel na descricao da atividade

#lendo os arquivos
library(tidyverse)
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read_csv)
dados= bind_rows(myfiles[-7])
centros= myfiles[[7]]

#limpando os dados 
dados= dados[-c(4:19)]
dados= dados[-3]
dados$data= as.Date(dados$data)
