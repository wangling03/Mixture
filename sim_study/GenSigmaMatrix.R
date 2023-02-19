GenSigmaMatrix <- function (N,P, type) {

if (tolower(type) !=tolower("equalcorr4") && tolower(type) !=tolower("equalcorr8")
    && tolower(type) !=tolower("FanandSong") && tolower(type) !=tolower("IndepCauchy")
    && tolower(type) !=tolower("Nonlinear") && tolower(type) !=tolower("Threefactor")
    && tolower(type) !=tolower("Twofactor"))
        
    {message("Input the wrong type name")
     return()}
        
err=matrix(0,N,P)

switch(type, 
       "equalcorr4"={ rho=0.4
                     Sigma <- matrix(rep(rho,P*P), nrow=P)
                     diag(Sigma)<- 1
                    },
       "equalcorr8"={ rho=0.8
                      Sigma <- matrix(rep(rho,P*P), nrow=P)
                      diag(Sigma)<- 1
                    },
       "FanandSong"={ numb=100
                      for (i in 1:(P-numb))
                        {err[,i]=rnorm(N,0,1)}
         
                         e_temp1=err[,1:10]
                         e_temp1[,seq(2, 10, 2)]=-e_temp1[,seq(2, 10, 2)]
                         e_temp2=rowSums(e_temp1)/5
         
                       for (i in (P-numb+1):P)
                         {e_temp3=sqrt(1-(10/25))*rnorm(N,0,1)
                          err[,i]=e_temp2+e_temp3}
         
                      Sigma <- cor(err)
                      },
       "IndepCauchy"= { for (i in 1:P)
                         {err[,i]=rcauchy(N,0,1)}
                       Sigma <- cor(err)
                       },
       "Nonlinear"= { W1=rnorm(1,0,1)
                      W2=rnorm(1,0,1)
                      for (i in 1:P)
                        {ro1=runif(N,-1,1)
                         ro2=runif(N,-1,1)
                         err[,i]=sin(ro1*W1)+sign(ro2)*exp(abs(ro2)*W2)+rnorm(N,0,1)}
                      Sigma <- cor(err)
                     },
       "Threefactor"= {W1=rnorm(1,-2,1)
                       W2=rnorm(1,1,1)
                       W3=rnorm(1,4,1)
       
                       for (i in 1:P)
                        {err[,i]=runif(N,-1,1)*W1+runif(N,-1,1)*W2+runif(N,-1,1)*W3+rnorm(N,0,1)}
                       Sigma <- cor(err)
                     },
       "Twofactor"= {W1=rnorm(1,0,1)
                     W2=rnorm(1,0,1)
                     for (i in 1:P)
                       {err[,i]=runif(N,-1,1)*W1+runif(N,-1,1)*W2+rnorm(N,0,1)}
                     Sigma <- cor(err)
                     }
)

return(Sigma)
}
