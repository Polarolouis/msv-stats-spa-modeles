
#==================================================================================
# QuadratsFonctions : fonctions qui calcule des stqtistiques à partir de quadrats
#                     rquadrat
#                     npquadrat
#                     npresquadrat
#                     plotquadrat
#                     testquadrat
#                     np2quadrat
#==================================================================================

rquadrat = function(n,Dlimx,Dlimy,qlimx,qlimy){
# tirage de n quadrats de taille qlimx*qlimy dans la fenetre Dlimx*Dlimy
#-----------------------------------------------------------------------
           x=runif(n,Dlimx[1],Dlimx[2]-qlimx)
           y=runif(n,Dlimy[1],Dlimy[2]-qlimy)
           quad = list()
           for (ifen in 1:n){
               quad[[ifen]] = owin(c(x[ifen],x[ifen]+qlimx),c(y[ifen],y[ifen]+qlimy))
           }
           quad
}

npquadrat = function(quad,X,Y){
# comptage des couples(X,Y) dans les quadrats quad
#-------------------------------------------------
            n = length(quad)
            npq = vector()
            for (ifen in 1:n){
                npq[ifen] = sum(inside.owin(X,Y,quad[[ifen]]))
            }
            npq
}

npresquadrat = function(quad,X,Y){
# indicateur de presence  des couples(X,Y) dans les quadrats quad
#----------------------------------------------------------------
            n = length(quad)
            npresq = vector()
            for (ifen in 1:n){
                npresq[ifen] = 1-prod(!inside.owin(X,Y,quad[[ifen]]))
            }
            npresq
}

plotquadrat = function(quad){
# dessin des quadrats quad sur un dessin existant
#-----------------------------------------------
              n = length(quad)
              for (ifen in 1:n){
                plot(quad[[ifen]],add=TRUE)
              }
}

testquadrat = function(PP,n=30,qlimx=1,qlimy=1){
# test de repartition aleatoire du processus PP, al'aide de n quadrats de taille qlimx*qlimy
#-------------------------------------------------------------------------------------------
              aireq = qlimx*qlimy
              aireD = diff(PP$window$xrange)*diff(PP$window$yrange)
              lambda = length(PP$x)/aireD
              quad = rquadrat(n,PP$window$xrange,PP$window$yrange,qlimx,qlimy)
              npq = npquadrat(quad,PP$x,PP$y)
              obstmp = table(npq)
              nb = as.integer(names(obstmp))
              eff = as.integer(obstmp)
              obs = list()
              for (iclas in 1:(max(nb)+1)){
                  classe = iclas-1
                  effectif = 0
                  if (sum(names(obstmp)==(iclas-1))==1) effectif = eff[names(obstmp)==(iclas-1)]
                  obs[[iclas]] = list(classe=classe,effectif=effectif)
              }
              effectif=vector()
              for (iclas in 1:length(obs)) effectif[iclas] = obs[[iclas]]$effectif
              
              while (any(effectif <5)){
                    iclas = min(which(effectif<5))
                    tmp=list()
                    if (iclas==1) {
                       classe = c(obs[[1]]$classe,obs[[2]]$classe)
                       eff = obs[[1]]$effectif + obs[[2]]$effectif
                       tmp[[1]] = list(classe=classe,effectif=eff)
                       for (iclas in 2:(length(obs)-1)){
                           tmp[[iclas]] = obs[[iclas+1]]  }          
                    } else if (iclas==length(obs)) {
                       classe = c(obs[[iclas-1]]$classe,obs[[iclas]]$classe)
                       eff = obs[[iclas-1]]$effectif + obs[[iclas]]$effectif
                       tmp[[iclas-1]] = list(classe=classe,effectif=eff)
                       for (iclas in 1:(length(obs)-2)){
                           tmp[[iclas]] = obs[[iclas]]  }  
                    } else {
                          ivois = iclas + 2*which.min(c(effectif[iclas-1],effectif[iclas+1]))-3
                          classe = c(obs[[ivois]]$classe,obs[[iclas]]$classe)
                          eff = obs[[ivois]]$effectif + obs[[iclas]]$effectif
                          obs[[ivois]] = list(classe=classe,effectif=eff)
                          for (icl in 1:(iclas-1)){
                              tmp[[icl]] = obs[[icl]]  }  
                          for (icl in (iclas):(length(obs)-1)){
                              tmp[[icl]] = obs[[icl+1]]  } 
                    }
                    effectif=vector()
                    obs = tmp
                    for (iclas in 1:length(obs)) effectif[iclas] = obs[[iclas]]$effectif
              }
              prob = vector()
              for (iclas in 1:(length(obs)-1)){
                  classe = obs[[iclas]]$classe
                  prob[iclas] = sum(dpois(classe,lambda*aireq))
              }
              prob[length(obs)] = 1-sum(prob)
              theo = n*prob
              DChi2 = sum((effectif-theo)^2/theo)
              dl = length(obs)-2
              pvaleur = 1-pchisq(DChi2,dl)
              tt = list(DChi2=DChi2,dl=dl,pvaleur=pvaleur,theoriques = theo,observes = effectif)
              return(tt)
}

np2quadrat = function(quad,PP1,PP2){
# comptage des processus PP1 et PP2 dans les quadrats quad
#-------------------------------------------------
            n = length(quad)
            np2q = matrix(0,n,4)
            for (ifen in 1:n){
                np2q[ifen,1] = (1-prod(!inside.owin(PP1$x,PP1$y,quad[[ifen]])))*(1-prod(!inside.owin(PP2$x,PP2$y,quad[[ifen]])))
                np2q[ifen,2] = (1-prod(!inside.owin(PP1$x,PP1$y,quad[[ifen]])))*(prod(!inside.owin(PP2$x,PP2$y,quad[[ifen]])))
                np2q[ifen,3] = (prod(!inside.owin(PP1$x,PP1$y,quad[[ifen]])))*(1-prod(!inside.owin(PP2$x,PP2$y,quad[[ifen]])))
                np2q[ifen,4] = (prod(!inside.owin(PP1$x,PP1$y,quad[[ifen]])))*(prod(!inside.owin(PP2$x,PP2$y,quad[[ifen]])))
            }
            tab2q = matrix(apply(np2q,2,sum),2,2)
            colnames(tab2q)=c('pres 1' ,'abs 1')
            rownames(tab2q)=c('pres 2' ,'abs 2')
            tab2q
}

Pielou = function(PPm){
# test du Khi2 sur le tableau des distances et indice de Pielou pour un processus a 2 marques
#------------------------------------------------------------
  coords = cbind(PPm$x,PPm$y)
  Dist = as.matrix(dist(coords))
  diag(Dist)=max(Dist)
  marquePV = PPm$marks[apply(Dist,1,which.min)]
  tabPV = table(PPm$marks,marquePV)
  test = chisq.test(tabPV)
  indice = 1 - sum(tabPV)*(tabPV[1,2]+tabPV[2,1])/(sum(PPm$marks==1)*sum(marquePV==2)+sum(PPm$marks==2)*sum(marquePV==1))
  retour = list(test=test,indice =indice)
  return(retour)
}
