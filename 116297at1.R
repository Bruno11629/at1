#O codigo presume que a pasta at1, possui o mesmo conteudo 
#e mesma ordenacao daquela disponivel na descricao da atividade

#setup
install.packages("tidyverse")
install.packages("geosphere")

library(tidyverse)
library(geosphere)

#lendo os arquivos
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read_csv)
dados= bind_rows(myfiles[-7]) 
centros= myfiles[[7]]

#tratando os dados
dados= dados[-c(4:19)]
dados= dados[-3]
dados$data= as.Date(dados$data,format = "%d/%m/%Y")
dados= dados %>% filter(data >= "2017-09-03"& data < "2018-08-29")

#construindo tabela

conta_chuva= function(dado){
faltantes=dado %>% filter(precipitacao== "////") %>% nrow()
if(faltantes== nrow(dado)){
  return(NULL)
}
  else{
    prop_falt= faltantes/nrow(dado)
    dadovalue= dado  %>% filter(precipitacao!="////")
    dadovalue$precipitacao = as.numeric(dadovalue$precipitacao)
    soma= sum(dadovalue$precipitacao)/(1-prop_falt)
    return(soma)
}
}


corta_dados= function(dados){
idx= 1
i=1
dados= dados%>% arrange(data)
data= dados$data[1]
chuva= NULL
while (idx <= nrow(dados)){
  while ((dados$data[i]-dados$data[idx])<10 & i<=nrow(dados)) {
    i=i+1
    }
  dado10= dados %>% slice(idx:(i-1))
  guarda.chuva= conta_chuva(dado10)
  if(!is.null(guarda.chuva)){
    data= c(data,dados$data[idx])
    chuva= c(chuva,guarda.chuva)
  }
  idx=i
}
saida= as.data.frame(cbind(as.character(data[-1]),chuva))
names(saida)= c("data","precipitacao")
return(saida)
}

#Correcao nos dados para comportar datas nao conformes
dados= dados %>% filter(!(codigo_estacao == "A728" & data== "2017-12-21") )
dados= dados %>% filter(!(codigo_estacao == "A744" & data < "2017-12-22") )
dados= dados %>% filter(!(codigo_estacao == "A771" & data < "2018-03-22") )

#Aplicando funcoes
A530= dados %>% filter(codigo_estacao== "A530") %>% corta_dados()
A701= dados %>% filter(codigo_estacao== "A701") %>% corta_dados()
A706= dados %>% filter(codigo_estacao== "A706") %>% corta_dados()
A711= dados %>% filter(codigo_estacao== "A711") %>% corta_dados()
A713= dados %>% filter(codigo_estacao== "A713") %>% corta_dados()
A715= dados %>% filter(codigo_estacao== "A715") %>% corta_dados()
A726= dados %>% filter(codigo_estacao== "A726") %>% corta_dados()
A728= dados %>% filter(codigo_estacao== "A728") %>% corta_dados()
A738= dados %>% filter(codigo_estacao== "A738") %>% corta_dados()
A739= dados %>% filter(codigo_estacao== "A739") %>% corta_dados()
A744= dados %>% filter(codigo_estacao== "A744") %>% corta_dados()
A755= dados %>% filter(codigo_estacao== "A755") %>% corta_dados()
A765= dados %>% filter(codigo_estacao== "A765") %>% corta_dados()
A771= dados %>% filter(codigo_estacao== "A771") %>% corta_dados()

Y_matriz= full_join(A530,A701,by= "data")
Y_matriz= full_join(Y_matriz,A706,by="data")
Y_matriz= full_join(Y_matriz,A711,by="data")
Y_matriz= full_join(Y_matriz,A713,by="data")
Y_matriz= full_join(Y_matriz,A715,by="data")
Y_matriz= full_join(Y_matriz,A726,by="data")
Y_matriz= full_join(Y_matriz,A728,by="data")
Y_matriz= full_join(Y_matriz,A738,by="data")
Y_matriz= full_join(Y_matriz,A739,by="data")
Y_matriz= full_join(Y_matriz,A744,by="data")
Y_matriz= full_join(Y_matriz,A755,by="data")
Y_matriz= full_join(Y_matriz,A765,by="data")
Y_matriz= full_join(Y_matriz,A771,by="data")

names(Y_matriz)= c("Data","A530","A701","A706","A711","A713","A715","A726","A728","A738","A739",
                   "A744","A755","A765","A771")

Y_vec= gather(Y_matriz,key = "codigo_estacao", value = "Precipitacao",-Data)
Y_vec$Precipitacao= as.numeric(Y_vec$Precipitacao)
centros= centros %>% arrange(Estação)
centros$codigo_estacao= c("A755", "A765","A744","A530","A706","A738","A739","A726","A711",
                          "A715","A771","A701","A713","A728")
centros= centros %>% arrange(codigo_estacao)




f.theta= function(theta,Y_vec,centros){
  y=Y_vec$Precipitacao[!is.na(Y_vec$Precipitacao)]
  mu=theta[1]
  sigma2=theta[2]
  phi.t=theta[3]
  phi.s=theta[4]
  sigma.T= matrix(nrow=36,ncol = 36)
  for(i in 1:36){
    for(j in 1:36){
     sigma.T[i,j]= sigma2*exp(-abs(i-j)/phi.t) 
    }
  }
  sigma.S= matrix(nrow = 14, ncol = 14)
  for(i in 1:14){
    for(j in 1:14){
      sigma.S[i,j]= exp(-distHaversine(c(centros$Longitude[i],centros$Latitude[i]),
                               c(centros$Longitude[j],centros$Latitude[j]))/phi.s )
    }
  }
  sigma.T.inv = solve(sigma.T)
  sigma.S.inv= solve(sigma.S)
  kron= sigma.S.inv %x% sigma.T.inv
  kron= kron[!is.na(Y_vec$Precipitacao),!is.na(Y_vec$Precipitacao)]
  L= -(-504*log(2*pi) - 18*log(det(sigma.T)) - 7*log(det(sigma.S)) -
    0.5*t((y- mu))%*%kron%*%(y- mu))
  return(L)
}

mu=mean(Y_vec$Precipitacao[!is.na(Y_vec$Precipitacao)])
phi.t=0.1
phi.s=1000
sigma2=1000

theta_hat=optim(c(mu,sigma2,phi.t,phi.s), f.theta,Y_vec=Y_vec,centros= centros)


