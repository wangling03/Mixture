GenBeta <- function (P, w, betatype, vartype) {
  
  if (tolower(betatype) !=tolower("Uninormal") && tolower(betatype) !=tolower("Binormal"))
   {message("Input the wrong betatype name")
    return()}
  
  if (tolower(vartype) !=tolower("v1") && tolower(vartype) !=tolower("v2")
      && tolower(vartype) !=tolower("v3"))
  {message("Input the wrong vartype name")
    return()}
  
  cpct <- 0.4
  beta=matrix(0,1,P)
  
  bimodalnormal <- function (n,cpct, mu1, mu2, sig1, sig2) {
    y0 <- rnorm(n,mean=mu1, sd = sig1)
    y1 <- rnorm(n,mean=mu2, sd = sig2)
    
    flag <- rbinom(n,size=1,prob=cpct)
    y <- y1*(1 - flag) + y0*flag 
  }

switch(vartype, 
       "v1"={ unimean=5
              univar=1
              bimean1=3
              bivar1=1
              bimean2=8
              bivar2=1
             },
       "v2"={ unimean=10
              univar=2
              bimean1=6
              bivar1=2
              bimean2=16
              bivar2=2
             },
       "v3"={ unimean=10
              univar=3
              bimean1=6
              bivar1=3
              bimean2=16
              bivar2=3
              }
)

  switch(betatype, 
         "Uninormal"={J<-rbern(P, w) 
                      for (j in 1:P){
                        if (J[j]==1) {beta[j]<-0}
                           else {beta[j]<-rnorm(1,unimean,abs(sqrt(univar)))} 
                       }
         },
         "Binormal"={ J<-rbern(P, w) 
                      for (j in 1:P){
                        if (J[j]==1) {beta[j]<-0}
                         else {beta[j]<-bimodalnormal(1,cpct,bimean1,bimean2,abs(sqrt(bivar1)),abs(sqrt(bivar2)))}
                    }
         }
  )
  return(beta)  
}