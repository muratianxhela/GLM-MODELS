# GLM-MODELS

Utilizzando il software statistico R, vengono forniti gli strumenti per l'analisi dei dati tramite i diversi modelli parametrici e semiparametrici, usando il dataset "Datasets for Use with Salvan, Sartori and Pace (2020)", in particolare in questo progetto si è considerato il dataset clotting and beetels. 

**Principali argomenti di glm**:

 - formula: come per lm;
 - family: che può essere binomial, gaussian, Gamma, Poisson,...
 -  link: come argomento di family (default è il legame canonico)

*Principali funzioni applicabili all'oggetto glm*:

  - summary: per ottenere il risultato di un oggetto glm
  -  confint: per ottenere gli IC per i coefficienti di regressione
  - anova: per ocnfrontare modelli annidati tramite test di devianza
  - plot: per ottenere l'analisi grafica dei residui
  -  fitted: per ottenere i valori stimati da eta_i.hat
  - predict: per ottenere i valori del predittore lineare
  - residuals: per ottenere i resudui del modello
  - rstandard: per ottenere i residui standardizzati
  - qresiduals: (libreria 'statmod') per ottenere i residui-quantile

```
setwd("C:\\Users\\39389\\CartellaDati")

clotting <- read.table("clotting.dat") 
attach(clotting) 
```
ANALISI GRAFICA
```
plot(log(u),tempo,type="n")
points(log(u[lotto == "uno"]),tempo[lotto == "uno"])
points(log(u[lotto != "uno"]),tempo[lotto != "uno"],pch=19,col=2)
```
ADATTAMENTO DEL MODELLO GLM GAMMA-funzione di legame: legame canonico.
```
clotting.glm <- glm(tempo~log(u)+lotto,family=Gamma("inverse"))
summary(clotting.glm)
```
Entrambe le variabili esplicative risultano significative. Adesso possiamo introdure dell'effetto di interazione tra le due esplicative.
```
clotting.glm2 <- update(clotting.glm,.~.+log(u):lotto) 
summary(clotting.glm2)
```
L'effetto di interazione è fortemente significativo, il secondo modello risulta preferibile.

Conferma ulteriore: CONFRONTO TRA LE DEVIANZE ('anova')
```
anova(clotting.glm,clotting.glm2,test="F")
```
Il modello sembra fornire un buon adattamento ai dati osservati. 

Grafico
```
par(mfrow=c(1,2),pty="s") 
plot(log(u),tempo,type="n",ylim=range(fitted(clotting.glm2),tempo)) 
points(log(u[lotto == "uno"]),tempo[lotto == "uno"])
points(log(u[lotto != "uno"]),tempo[lotto != "uno"],pch=19,col=2)
points(log(u),fitted(clotting.glm2),col=3,pch=4) 
plot(tempo,fitted(clotting.glm2),xlab="valori osservati", ylab="valori stimati")
abline(0,1,col=2)
```
Provi a stimare un modello gamma ocn funzione di legame logaritmica e si confronti. L'adattamento con quello ottenuto attraverso il legame canonico.
```
  clotting.glmE1 <- glm(tempo~log(u)+lotto,family=Gamma("log"))
  summary(clotting.glmE1)
```
 Entrambe le variabili esplicative risultano significative; p-value altamente. Introduzione dell'effetto di interazione tra le due esplicative.
 ```
  clotting.glmE2 <- update(clotting.glmE1,.~.+log(u):lotto) 
  summary(clotting.glmE2)
  ```
L'effetto di interazione NON è fortemente significativo, il primo modello risulta preferibile.

Conferma ulteriore: CONFRONTO TRA LE DEVIANZE ('anova')
  ```
  anova(clotting.glmE1,clotting.glmE2,test="F")
  ```
 Il modello sembra fornire un buon adattamento ai dati osservati.

Grafico 
```
  par(mfrow=c(1,2),pty="s") 
  plot(log(u),tempo,type="n",ylim=range(fitted(clotting.glmE1),tempo)) 
  points(log(u[lotto == "uno"]),tempo[lotto == "uno"])
  points(log(u[lotto != "uno"]),tempo[lotto != "uno"],pch=19,col=2)
  points(log(u),fitted(clotting.glmE1),col=3,pch=4) 
  plot(tempo,fitted(clotting.glmE1),xlab="valori osservati", ylab="valori stimati")
  abline(0,1,col=2)
  
detach(clotting)  
```
**SECOND DATASET- dataset="chimps.dat"**
```
library(SMPracticals)
data(chimps)
attach(chimps)
head(chimps)  
```
STIMA DEL MODELLO
Modello Gamma con predittore additivo (senza termine di interazione), con legame canonico (inversa).
```
chimps.glm<-glm(y~chimp+word,family=Gamma)
summary(chimps.glm)
```
N.B. viene calcolato anche il parametro di dispersioneche tra l'altro può essere calcolato manualmente anche così:
```
muhat<-fitted(chimps.glm)
phihat<-sum((y-muhat)^2/muhat^2)/chimps.glm$df.residual
phihat
```
Tre possibili modi di confroto della significaatività di un paramentro:

-  TEST TRA DEVIANZE
```
D0<-chimps.glm$null.deviance     #Devianza nulla del glm
DC<-chimps.glm$deviance          #Devianza residua del glm
df<-chimps.glm$df.null-chimps.glm$df.residual #gradi di libertà (differenza)
t<-(D0-DC)/phihat      #test statistico sulla devianza
alphaoss<- 1-pchisq(t,df)   #livello di significatività osservato
alphaoss
```
Il modello stimato spiega i dati significativamente meglio del modello nullo (di ccs)

-  TEST RAPPORTO DI VEROSIMIGLIANZA
Assenza dell'effetto parola, avendo tenuto conto solo dell'effetto scimpanzè)
```
anova(chimps.glm, test="Chisq")
```
Viene rifiutata l'ipotesi nulla (H0: tutti i coeff =0)

- CONFRONTO TRA DUE MODELLI (base e complesso)
```
chimps.glm0<-glm(y~chimp, family=Gamma)
anova(chimps.glm0,chimps.glm, test="Chisq")
```
- NOTA_1: l'ultimo risultato è preferibile perche, il risultato 'anova' applicato a un singolo glm. Dipende dell'ordine secondo cui le variabili appaiono nella formula del predittore lineare.
- NOTA_2: nel modello gamma, dove l'ignoto parametro di dispersione deve essere stimato, la distribuzione nulla del test del log rapporto di verosimiglianza può essere approssimato da una distribuzione F.
```
anova(chimps.glm, test="F")
```
Con F i livelli di significatività osservati tendono ad essere più grandi di quelli del X^2. In questo caso la soluzione è più CONSERVATIVA. Lo stesso discorso vale se si utilitta la t (default) rispetto alla normale nei test di Wald sui singoli coefficienti.

Funzione di legame logaritmo:
```
chimps.glm1<-glm(y~chimp+word,family=Gamma(link=log))
summary(chimps.glm1)
```
Interpretazione Intercetta: il logaritmo della media del tempo di apprendimento della parola da parte del primo scimpanzè. Più coeff. significativamente diversi da zero, AIC più basso => modello leggermente preferibile.
```
par(mfrow=c(2,2))
plot(chimps.glm,1:4) #legame canonico
plot(chimps.glm1,1:4) #legame logartmico

detach(chimps)
par(mfrow=c(1,1))
```
**THIRD DATASET-dataset= "cement.dat"**
```
cement<-read.table("cement.dat",header=TRUE)
attach(cement)
```
ANALISI PRELIMINARI
```
str(cement)
plot(tempo,resistenza)
plot(1/tempo,resistenza)
par(mfrow=c(1,2))
plot(tempo,resistenza)
plot(1/tempo,resistenza)
```
Come si pò vedere dai resultati il reciproco del tempo di indurimento sembra più lineare

ADATTAMENTO DEL MODELLO:

Modello Gamma con legame canonico sul reciproco del tempo
```
cement.glm<-glm(resistenza~I(1/tempo),family = Gamma("inverse"))
summary(cement.glm)
```
ANALISI DEI RESIDUI
```
par(mfrow=c(2,2))
plot(cement.glm,1:4)
```
Mostra una possibile esistenza di una relazione quadratica.

Modello con l'aggiunta di un termine quadratico della variabile esplicativa
```
cement.glm1<-update(cement.glm,.~.+I(1/tempo^2),family=Gamma("inverse"))
summary(cement.glm1)
```
ANALISI DEI RESIDUI
```
plot(cement.glm1,1:4)
```
CONFRONTO TRA MODELLI
```
anova(cement.glm,cement.glm1)
```
Il secondo modello, già dall'analisi dei residui sembra migliorare nell'adattamento anche con anova.

Modello: (alternativa) con diversa funzione di legame (logaritmica).
```
cement.glm2<-glm(resistenza~I(1/tempo),family= Gamma("log"))
summary(cement.glm2)
```
ANALISI DEI RESIDUI
```
plot(cement.glm2,1:4)
```
CONFRONTO TRA MODELLI:in termini di vicinaza tra i valori predetti ai valori osservati:
```
par(mfrow=c(1,1))
plot(tempo,resistenza)
lines(fitted(cement.glm)~tempo, col=2)
lines(fitted(cement.glm1)~tempo, col=3)
```
Visto che la curva verde sembra più vicina ai dati e si riferisce ad un modello più parsimonioso, il modello Gamma con funzione di legame logaritmica risulta preferibile.

```
par(mfrow=c(1,1))
plot(sqrt(1/tempo),1/resistenza)
cement.lm<-lm(I(1/resistenza)~I(sqrt(1/tempo)))
summary(cement.lm)
tempolm<-sqrt(1/tempo)
resistenzalm<-1/resistenza
cement.lm1<-lm(resistenzalm~tempolm)
summary(cement.lm1)
par(mfrow=c(2,2))
plot(cement.lm1,1:4)

par(mfrow=c(1,1))
plot(tempo,resistenza)
lines(fitted(cement.glm)~tempo, col=2)
lines(fitted(cement.glm1)~tempo, col=3)
lines(fitted(cement.lm)~tempo,col="blue")
```
LAST DATASET-Beetles
```
Beetles<-read.table("beetles.dat", header=T)
colnames(Beetles)<-c("num", "uccisi", "logdose")
Beetles
attach(Beetles)
```
y<-uccisi/num
y
```
- Var risp= prop successi
- Variaz dose ha effetto e descrivere qdt effetto
- Modello saturo: ( distr binomiali, vet y= vet smv nel mod saturo)

Analisi grafiche
```
plot(y~logdose)
```
- Andamento  monotono cresc (c'è effetto dose);
- Andamento curva nn lineare-> usare funz legame diversa da identità, anche perchè funz legame deve essere coerente.

Analisi esplorativa
Sepre obiettivo se e che tipo di rel c'è tra risp (prop uccisi) ed esplic e lo vediamo su diverse scale x scegliere fnz legame
Logit: anzichè prendere direttamente y prendo trasformazione
```
logit<-log(y/(1-y))
logit  #ultimo è infinito, nn dà errore

plot(logdose,logit)
```
Il num di pti nn è 8 ma 7 xk ultimo pto è = oo tasformazione logistica empirica (correzione asintoticamente nn rilevante)
```
logitc<-log((y+0.5/num)/(1-y+0.5/num))
logitc

plot(logdose,logitc)
```
`Logit circa 5= prob quote circa 5-> exp(5)` cioè numero molto grande della proporz di successi su insuccessi. Vediamo quanto lineare è la relazione (analisi grafica)
```
abline(lm(logitc~logdose), col=2)
```
Discreto adattamento
```
summary(lm(logitc~logdose))$r.squared # R-squared:  0.96
```
Probit:
```
qnorm(y)
```
- Val di x per cui qnorm=1 infinito -> serve correzione
- Aggiusti un po' le y in modo che il qnorm ti venga un numero finito
- qnorm(0)=-inf

Probit con correzione empirica
```
yc<-y+0.1/num #correggiamo tt in un verso poi aggiustiamo
yc
yc[y>0.5]<-y[y>0.5]-0.1/num[y>0.5]
yc
y

probitc<-qnorm(yc)
probitc #tt valori strett finiti

plot(logdose, logitc)
abline(lm(logitc~logdose), col=1)
points(logdose,probitc, pch=2, col=2)
abline(lm(probitc~logdose), col=2)
```
Non si capisce ad occhio la migliore
```
summary(lm(probitc~logdose))$r.squared
```
- Adattamento stesso livello.
Faccio stessa cosa con correzioni 0 e 1 con le loglog complementare e cauchy
```
loglogc<-log(-log(1-yc))
loglogc
plot(logdose, loglogc)
cauchy<-tan(pi*(yc-0.5))

plot(logdose, logitc)
abline(lm(logitc~logdose), col=1)
points(logdose,probitc, pch=2, col=2)
abline(lm(probitc~logdose), col=2)
points(logdose,loglogc, pch=3, col=3)
abline(lm(loglogc~logdose), col=3)
points(logdose,cauchy, pch=4, col=4)
abline(lm(cauchy~logdose), col=4)

summary(lm(loglogc~logdose))$r.squared #0.99 bene
summary(lm(cauchy~logdose))$r.squared # 0.3637069
```
Come si può vedere dai resultati in questo caso il modelllo loglogcè il migliore per i nostri dati 

Modelli lin generalizzzati


**Regressione logistica**
Abbiamo 2 possibilità con dati raggruppati
Beetles.glm<-glm(y~logdose, family=binomial, weights=num)
summary(Beetles.glm)
#equivalente:
Beetles.glm1<-glm(cbind(uccisi, num-uccisi)~logdose, family=binomial)
summary(Beetles.glm1)
#esempio 3.1 null deviance e residual deviance
#dati raggrupp: posso usare dev res cm ind bontà di adattam, è trv modello corrente vs modllo saturo
#gdl= 8-2=6 gdl

#provo a calc con de res il pval del test del mod corr risp a od saturo

1-pchisq(11.232, 6)
#pval viene  0.08146544 qundi nn c'è evidenza forte contro il modello corrente verso il modello saturo

#ricalcolo a mano la matr di cov dello stimatore betahat e di ottenere per via diretta std err beta
#vedi es a pag 169
#XtWx
#W= v(mui)*ai(phi)

#w_i=m_i*pi_i*(1-pi_i) con pi_i stimati (mui cappello)ottenuti da fitted(Beetles.glm)


hatpi<-fitted(Beetles.glm)
W<-diag(num*hatpi*(1-hatpi))
X=model.matrix(Beetles.glm) #dà matr X, 8x2
X

varbeta<-solve(t(X)%*%W%*%X)
stdbeta<-sqrt(diag(inversa))
#altrimenti 
vcov(Beetles.glm) #viene l'inversa, cioè la matr di var


#predittore lineare e IC

hatbeta<-coef(Beetles.glm)
hateta<-t(c(1, 1.8))%*%hatbeta  #etacappello
hateta
#varianza=var trasf lin

varhateta<-t(c(1, 1.8))%*% varbeta %*% c(1, 1.8)
sqrt(varhateta)


#IC per eta con c=1.8

0.961318+c(-1, 1) * qnorm(0.975) * 0.1450564
predict(Beetles.glm, newdata=data.frame(logdose=1.8), se=T) 
#ic per pi(1.8). 2 metodi: da ic per eta a ic per pi (e^eta/1-e^eta)
#provo con predict
#metodo 1:trasformare con l'inversa della logit l'ic per eta


#2: usare il metodo delta:usare dirett predict:
predict(Beetles.glm, newdata=data.frame(logdose=1.8), se=T, type="response") #ic fatto su scala della risposta
#da solo:adatto modello probit e loglogc
Beetles.glm2<-glm(y~logdose, family=binomial(link=probit), weights=num) #se uso y devo mettere pesi
#rel tra coeff modello logit = 1.6 * probit

summary(Beetles.glm2)$coef
19.73*1.6 #circa 31

plot(logdose,y)
lines(logdose, hatpi, lty=2)  #associare ai pti osservati sulla scala originaria (logsdose) i valori predetti (linea=modello stimato)
#posso sovrapporre pti stimati con modello probit

hatpi2<-fitted(Beetles.glm2)
lines(logdose, hatpi2, lty=3, col=2)
#a livello grafico diff in term di aic poca


#valutare capacità predittiva con curva ROC

library(pROC)

#nella disp fatta con loglogc, noi probit
#carichiamo dati non raggruppati

Beetles10<-read.table("beetles10.dat", header=T)
colnames(Beetles10)<-c("logdose10", "ucciso10")
attach(Beetles10)
Beetles10.glm2<-glm(ucciso10~logdose10, family=binomial(probit))
pihat210<-fitted(Beetles10.glm2)
#lunghezza 481 ma ne ottengo di 8 tipi

#controllo che summary dia stesse stime di dati nìraggruppati

plot(roc(ucciso10, pihat210), print.auc=TRUE)
#stampa anche valore area sotto la curva: 0.91 va da 0 e 1, piu grande è migliore capacità predittiva, è tasso veri positivi/falsi positivi


