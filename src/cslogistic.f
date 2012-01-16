c=======================================================================                  
c=======================================================================                  
c     BAYESIAN CONDITIONALLY SPECIFIED LOGISTIC REGRESSION
c=======================================================================                  
c     cslogistic.f is a fortram code to fit a conditionally
c     specified logistic regression model for multivariate clustered
c     binary data in a Frequentist and in a Bayesian way. 
c     For the Bayesian analysis a multivariate normal prior is used.
c     
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
c      Fax  : +56-2-3547729  Email: atjara@uc.cl
c
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/
c      Fax  : +56-2-3547729  Email: mjgarcia@uc.cl

c     This software is distributed under the terms of the GNU GENERAL
c     PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
c     file for more information.
c
c     Copyright (C) 2005-2012 Alejandro Jara and Maria Jose Garcia-Zattera
c=======================================================================                  
c=======================================================================                  

c=======================================================================                  
      subroutine mha(mcmc,burnin,thin,nstore,seed1,seed2,seed3,propv,
     &               p,dima,dimen,n,nrec,y,x,perm,
     &               beta,betac,beta_pm,beta_pv,storemat,arate,
     &               workm1,workm2,workm3,workv1,workv2,workv3,
     &               iflag)
c=======================================================================
c     Perform MCMC algorithm for a conditional specified logistic
c     regression model. Equal effect of intercept and covariates.
c     A.J.V., 2005
      implicit none
      integer i,j,total_iter,dimen,n2,nx,mx
      integer mcmc,burnin,thin,nstore
      integer accepts,counter,seed1,seed2,seed3
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(dimen)
      real*8 x(nrec,p),propv(dimen,dimen)
      real*8 beta(dimen),betac(dimen)
      real*8 beta_pm(dimen),beta_pv(dimen,dimen)
      real*8 logpost_cur,logpost_can,ratio
      real*8 storemat(nstore,dimen),arate
      real*8 workm1(dimen*(dimen+1)/2),workm2(dimen,dimen)
      real*8 workm3(dimen,dimen)
      real*8 workv1(dimen),workv2(dimen),workv3(dimen)
      real*8 u
      real runif

c++++ initialize parameters

      accepts=0
      counter=0
      total_iter=mcmc+burnin

c++++ set random number generator

      call setrand(seed1,seed2,seed3)
       
c++++ calculate the permutations of vector response

      n2=2**n 
      do i=1,n2
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      call permuta(n,perm)     
 
c++++ evaluate log-posterior for current values of parameters

      call logposta(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &              logpost_cur,workv2,workm2,workm3,workv3,iflag)

c++++ start the Metropolis algorithm

      do i=1,total_iter

c+++++++ check if the user has requested an interrupt

         call rchkusr()

c+++++++ sample candidates
         
         call rmvnorm(dimen,beta,propv,workm1,workv1,betac)
      
c+++++++ evaluate log-posterior for candidate values of parameters

         call logposta(p,dima,betac,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                 logpost_can,workv2,workm2,workm3,workv3,iflag)

c+++++++ aceptation step

         ratio=exp(logpost_can-logpost_cur)
         u=runif()
         
         if(u.lt.ratio)then
            do j=1,dimen
               beta(j)=betac(j)
            end do
            
            logpost_cur=logpost_can
            
            accepts=accepts+1
         end if

c+++++++ save the samples
         
         if(i.gt.burnin) then
            nx=i/thin
	    mx=nx*thin
            if(mx.eq.i)then
               counter=counter+1 
               do j=1,dimen
                  storemat(counter,j)=beta(j)
               end do
            end if   
         end if
         
      end do
      
      arate=float(accepts)/float(total_iter)
      
      return
      end


c=======================================================================                  
      subroutine logposta(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                    logpost,vv,a,sigma2,res,iflag)
c=======================================================================
c     evaluate the log-posterior of a conditional specified logistic
c     regression model. Equal effect.
c     A.J.V., 2005
      implicit none
      integer i,j,k,l,dimen,n2,counter 
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(p+dima)
      real*8 beta(p+dima)
      real*8 beta_pm(p+dima),beta_pv(p+dima,p+dima)
      real*8 logprior
      real*8 x(nrec,p),loglike
      real*8 aa,bb,cc,xb,bi,ci,cci,logpost
      real*8 vv(p+dima),a(p+dima,p+dima)
      real*8 sigma2(p+dima,p+dima),res(p+dima)
      
      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n
      dimen=p+dima
      
c++++ Likelihood
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               do l=1,p
                  xb=xb+beta(l)*x(i,l)               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            counter=0
            do k=1,n
               do l=k+1,n
                  counter=counter+1
                  bi=bi+beta(p+counter)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
         aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            do k=1,p
               xb=xb+beta(k)*x(i,k)               
            end do
            bb=bb+xb*y(i,j)
         end do   

         counter=0
         do j=1,n
            do k=j+1,n
               counter=counter+1
               cc=cc+beta(p+counter)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=-aa+bb+cc

c++++ Prior
      call dmvn(dimen,beta,beta_pm,beta_pv,logprior,vv,a,sigma2,res,
     &          iflag)   

c++++ Posterior
      logpost=loglike+logprior      
      
      return
      end


c=======================================================================                  
      subroutine mhb(mcmc,burnin,thin,nstore,seed1,seed2,seed3,propv,
     &               p,dima,dimen,n,nrec,y,x,perm,
     &               beta,betac,beta_pm,beta_pv,storemat,arate,
     &               workm1,workm2,workm3,workv1,workv2,workv3,
     &               iflag)
c=======================================================================
c     Perform MCMC algorithm for a conditional specified logistic
c     regression model. Different effect
c     A.J.V., 2005
      implicit none
      integer i,j,total_iter,dimen,n2,nx,mx
      integer mcmc,burnin,thin,nstore
      integer accepts,counter,seed1,seed2,seed3
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(dimen)
      real*8 x(nrec,p),propv(dimen,dimen)
      real*8 beta(dimen),betac(dimen)
      real*8 beta_pm(dimen),beta_pv(dimen,dimen)
      real*8 logpost_cur,logpost_can,ratio
      real*8 storemat(nstore,dimen),arate
      real*8 workm1(dimen*(dimen+1)/2),workm2(dimen,dimen)
      real*8 workm3(dimen,dimen)
      real*8 workv1(dimen),workv2(dimen),workv3(dimen)
      real*8 u
      real runif

c++++ initialize parameters

      accepts=0
      counter=0
      total_iter=mcmc+burnin

c++++ set random number generator

      call setrand(seed1,seed2,seed3)
       
c++++ calculate the permutations of vector response

      n2=2**n 
      do i=1,n2
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      call permuta(n,perm)     
 
c++++ evaluate log-posterior for current values of parameters

      call logpostb(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &              logpost_cur,workv2,workm2,workm3,workv3,iflag)

c++++ start the Metropolis algorithm

      do i=1,total_iter

c+++++++ check if the user has requested an interrupt

         call rchkusr()

c+++++++ sample candidates
         
         call rmvnorm(dimen,beta,propv,workm1,workv1,betac)
      
c+++++++ evaluate log-posterior for candidate values of parameters

         call logpostb(p,dima,betac,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                 logpost_can,workv2,workm2,workm3,workv3,iflag)

c+++++++ aceptation step

         ratio=exp(logpost_can-logpost_cur)
         u=runif()
         
         if(u.lt.ratio)then
            do j=1,dimen
               beta(j)=betac(j)
            end do
            
            logpost_cur=logpost_can
            
            accepts=accepts+1
         end if

c+++++++ save the samples
         
         if(i.gt.burnin) then
            nx=i/thin
	    mx=nx*thin
            if(mx.eq.i)then
               counter=counter+1 
               do j=1,dimen
                  storemat(counter,j)=beta(j)
               end do
            end if   
         end if
         
      end do
      
      arate=float(accepts)/float(total_iter)
      
      return
      end


c=======================================================================                  
      subroutine logpostb(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                    logpost,vv,a,sigma2,res,iflag)
c=======================================================================
c     evaluate the log-posterior of a conditional specified logistic
c     regression model. Different effect.
c     A.J.V., 2005
      implicit none
      integer i,j,k,l,dimen,n2,counter 
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(n*p+dima)
      real*8 beta(n*p+dima)
      real*8 beta_pm(n*p+dima),beta_pv(n*p+dima,n*p+dima)
      real*8 logprior
      real*8 x(nrec,p),loglike
      real*8 aa,bb,cc,xb,bi,ci,cci,logpost
      real*8 vv(n*p+dima),a(n*p+dima,n*p+dima)
      real*8 sigma2(n*p+dima,n*p+dima),res(n*p+dima)
      
      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n
      dimen=n*p+dima
      
c++++ Likelihood
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               do l=1,p
                  xb=xb+beta((k-1)*p+l)*x(i,l)               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            counter=0
            do k=1,n
               do l=k+1,n
                  counter=counter+1
                  bi=bi+beta(n*p+counter)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
         aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            do k=1,p
               xb=xb+beta((j-1)*p+k)*x(i,k)               
            end do
            bb=bb+xb*y(i,j)
         end do   

         counter=0
         do j=1,n
            do k=j+1,n
               counter=counter+1
               cc=cc+beta(n*p+counter)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=-aa+bb+cc

c++++ Prior
      call dmvn(dimen,beta,beta_pm,beta_pv,logprior,vv,a,sigma2,res,
     &          iflag)   

c++++ Posterior
      logpost=loglike+logprior      
      
      return
      end


c=======================================================================                  
      subroutine mhc(mcmc,burnin,thin,nstore,seed1,seed2,seed3,propv,
     &               p,dima,dimen,n,nrec,y,x,perm,
     &               beta,betac,beta_pm,beta_pv,storemat,arate,
     &               workm1,workm2,workm3,workv1,workv2,workv3,
     &               iflag)
c=======================================================================
c     Perform MCMC algorithm for a conditional specified logistic
c     regression model. Different intercept and common effect for the
c     covariates.
c     A.J.V., 2005
      implicit none
      integer i,j,total_iter,dimen,n2,nx,mx
      integer mcmc,burnin,thin,nstore
      integer accepts,counter,seed1,seed2,seed3
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(dimen)
      real*8 x(nrec,p),propv(dimen,dimen)
      real*8 beta(dimen),betac(dimen)
      real*8 beta_pm(dimen),beta_pv(dimen,dimen)
      real*8 logpost_cur,logpost_can,ratio
      real*8 storemat(nstore,dimen),arate
      real*8 workm1(dimen*(dimen+1)/2),workm2(dimen,dimen)
      real*8 workm3(dimen,dimen)
      real*8 workv1(dimen),workv2(dimen),workv3(dimen)
      real*8 u
      real runif

c++++ initialize parameters

      accepts=0
      counter=0
      total_iter=mcmc+burnin

c++++ set random number generator

      call setrand(seed1,seed2,seed3)
       
c++++ calculate the permutations of vector response

      n2=2**n 
      do i=1,n2
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      call permuta(n,perm)     
 
c++++ evaluate log-posterior for current values of parameters

      call logpostc(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &              logpost_cur,workv2,workm2,workm3,workv3,iflag)

c++++ start the Metropolis algorithm

      do i=1,total_iter

c+++++++ check if the user has requested an interrupt

         call rchkusr()

c+++++++ sample candidates
         
         call rmvnorm(dimen,beta,propv,workm1,workv1,betac)
      
c+++++++ evaluate log-posterior for candidate values of parameters

         call logpostc(p,dima,betac,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                 logpost_can,workv2,workm2,workm3,workv3,iflag)

c+++++++ aceptation step

         ratio=exp(logpost_can-logpost_cur)
         u=runif()
         
         if(u.lt.ratio)then
            do j=1,dimen
               beta(j)=betac(j)
            end do
            
            logpost_cur=logpost_can
            
            accepts=accepts+1
         end if

c+++++++ save the samples
         
         if(i.gt.burnin) then
            nx=i/thin
	    mx=nx*thin
            if(mx.eq.i)then
               counter=counter+1 
               do j=1,dimen
                  storemat(counter,j)=beta(j)
               end do
            end if   
         end if
         
      end do
      
      arate=float(accepts)/float(total_iter)
      
      return
      end


c=======================================================================                  
      subroutine logpostc(p,dima,beta,beta_pm,beta_pv,nrec,n,y,x,perm,
     &                    logpost,vv,a,sigma2,res,iflag)
c=======================================================================
c     evaluate the log-posterior of a conditional specified logistic
c     regression model. Different intercept and common covariates.
c     A.J.V., 2005
      implicit none
      integer i,j,k,l,dimen,n2,counter 
      integer p,dima,nrec,n,y(nrec,n),perm(2**n,n),iflag(n+(p-1)+dima)
      real*8 beta(n+(p-1)+dima)
      real*8 beta_pm(n+(p-1)+dima),beta_pv(n+(p-1)+dima,n+(p-1)+dima)
      real*8 logprior
      real*8 x(nrec,p),loglike
      real*8 aa,bb,cc,xb,bi,ci,cci,logpost
      real*8 vv(n+(p-1)+dima),a(n+(p-1)+dima,n+(p-1)+dima)
      real*8 sigma2(n+(p-1)+dima,n+(p-1)+dima),res(n+(p-1)+dima)
      
      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n
      dimen=n+(p-1)+dima
      
c++++ Likelihood
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               xb=xb+beta(k)
               do l=1,(p-1)
                  xb=xb+beta(n+l)*x(i,(l+1))               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            counter=0
            do k=1,n
               do l=k+1,n
                  counter=counter+1
                  bi=bi+beta(n+(p-1)+counter)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
         aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            xb=xb+beta(j) 
            do k=1,(p-1)
               xb=xb+beta(n+k)*x(i,(k+1))               
            end do
            bb=bb+xb*y(i,j)
         end do   

         counter=0
         do j=1,n
            do k=j+1,n
               counter=counter+1
               cc=cc+beta(n+(p-1)+counter)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=-aa+bb+cc

c++++ Prior
      call dmvn(dimen,beta,beta_pm,beta_pv,logprior,vv,a,sigma2,res,
     &          iflag)   

c++++ Posterior
      logpost=loglike+logprior      
      
      return
      end


c=======================================================================
      subroutine cloga(n,nrec,p,dima,y,x,beta,perm,loglike)
c=======================================================================      
c     evaluate the negative log-likelihood of a conditional specified 
c     logistic regression model. Equal effect.
c     A.J.V., 2005
      implicit none 
      integer i,j,k,l,count
      integer n,n2,nrec,p,dima,perm(2**n,n),y(nrec,n)
      real*8 x(nrec,p),beta(p+dima)
      real*8 aa,bb,cc,xb,bi,ci,cci,loglike

      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n 
      
      do i=1,2**n
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      
      call permuta(n,perm)      
      
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               do l=1,p
                  xb=xb+beta(l)*x(i,l)               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            count=0
            do k=1,n
               do l=k+1,n
                  count=count+1
                  bi=bi+beta(p+count)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
        aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            do k=1,p
               xb=xb+beta(k)*x(i,k)               
            end do
            bb=bb+xb*y(i,j)
         end do   

         count=0
         do j=1,n
            do k=j+1,n
               count=count+1
               cc=cc+beta(p+count)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=aa-(bb+cc)
      
      return
      end

c=======================================================================
      subroutine clogb(n,nrec,p,dima,y,x,beta,perm,loglike)
c=======================================================================      
c     evaluate the log-likelihood of a conditional specified logistic
c     regression model. Different effect.
c     A.J.V., 2005
      implicit none 
      integer i,j,k,l,count
      integer n,n2,nrec,p,dima,perm(2**n,n),y(nrec,n)
      real*8 x(nrec,p),beta(n*p+dima)
      real*8 aa,bb,cc,xb,bi,ci,cci,loglike
      
      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n 
      
      do i=1,2**n
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      
      call permuta(n,perm)      
      
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               do l=1,p
                  xb=xb+beta((k-1)*p+l)*x(i,l)               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            count=0
            do k=1,n
               do l=k+1,n
                  count=count+1
                  bi=bi+beta(n*p+count)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
        aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            do k=1,p
               xb=xb+beta((j-1)*p+k)*x(i,k)               
            end do
            bb=bb+xb*y(i,j)
         end do   

         count=0
         do j=1,n
            do k=j+1,n
               count=count+1
               cc=cc+beta(n*p+count)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=aa-(bb+cc)
      
      return
      end
   

c=======================================================================
      subroutine clogc(n,nrec,p,dima,y,x,beta,perm,loglike)
c=======================================================================      
c     evaluate the log-likelihood of a conditional specified logistic
c     regression model. Different Intercept and common covariates.
c     A.J.V., 2005
      implicit none 
      integer i,j,k,l,count
      integer n,n2,nrec,p,dima,perm(2**n,n),y(nrec,n)
      real*8 x(nrec,p),beta(n+(p-1)+dima)
      real*8 aa,bb,cc,xb,bi,ci,cci,loglike
      
      aa=0.d0
      bb=0.d0
      cc=0.d0
      bi=0.d0
      ci=0.d0
      cci=0.d0
      
      n2=2**n 
      
      do i=1,2**n
         do j=1,n
            perm(i,j)=0 
         end do
      end do
      
      call permuta(n,perm)      
      
      do i=1,nrec
         ci=0.d0
         do j=1,n2
            bi=0.d0
            do k=1,n
               xb=0.d0   
               xb=xb+beta(k)
               do l=1,(p-1)
                  xb=xb+beta(n+l)*x(i,(l+1))               
               end do
               bi=bi+xb*perm(j,k)
            end do   

            count=0
            do k=1,n
               do l=k+1,n
                  count=count+1
                  bi=bi+beta(n+(p-1)+count)*perm(j,k)*perm(j,l)
               end do
            end do

            ci=ci+dexp(bi)
         end do
        aa=aa+dlog(ci)

         do j=1,n
            xb=0.d0
            xb=xb+beta(j) 
            do k=1,(p-1)
               xb=xb+beta(n+k)*x(i,(k+1))               
            end do
            bb=bb+xb*y(i,j)
         end do   

         count=0
         do j=1,n
            do k=j+1,n
               count=count+1
               cc=cc+beta(n+(p-1)+count)*y(i,j)*y(i,k)
            end do
         end do         
      end do
      
      loglike=aa-(bb+cc)

      return
      end
   
c=======================================================================            
      subroutine permuta(n,a)
c=======================================================================
c     generate permutations of response variable
c     A.J.V., 2005      
      implicit none
      integer i,j,k,n,a(2**n,n),c1,numb,count,start
      
      numb=2**n
      do i=1,2**n
         do j=1,n
            a(i,j)=0
         end do
      end do
      
      do i=1,n
         c1=2**(n-i)
         count=numb/(2.d0*c1)
         start=c1
         do j=1,count
            do k=start+1,start+c1 
               a(k,i)=1
            end do
            start=start+2*c1
         end do
      end do
      return
      end
      

c=======================================================================                  
c=======================================================================                  
c     HELPPERS
c=======================================================================                  
c=======================================================================                  

c=======================================================================            
      subroutine setrand(seed1,seed2,seed3)
c=======================================================================            
c     This routine initialize the random number generator.
c     A.J.V., 2005
      integer seed1,seed2,seed3
      common/rrunif/iix,iiy,iiz
      iix=seed1
      iiy=seed2
      iiz=seed3
      return
      end

c=======================================================================                  
      real function runif()
c=======================================================================                  
c     This routine generates a uniform random variable using the
c     algorithm of Wichman and Hill (1982) 
c     A.J.V., 2005
      integer ix,iy,iz
      common/rrunif/ix,iy,iz
      ix=171*mod(ix,177)-2*(ix/177)
      iy=172*mod(iy,176)-35*(iy/176)
      iz=170*mod(iz,178)-63*(iz/178)
      if (ix.lt.0) ix=ix+30269
      if (iy.lt.0) iy=iy+30307
      if (iz.lt.0) iz=iz+30323
      runif=amod(float(ix)/30269.0+float(iy)/30307.0+
     +		    float(iz)/30323.0,1.0)
      return
      end
            
      
c=======================================================================            
      real function rnorm(mu,sd)
c=======================================================================            
c     This function generates a N(mu,sd^2) random values.
c     A.J.V., 2005
      real*8 mu,sd
      real u,runif
      u=runif()
      rnorm = mu+sd*ppnda(u)
      return
      end
      
      
c=======================================================================            
      real function ppnda(p)
c=======================================================================            
c     This function calculates the inverse normal distribution function
c     usin the algorithm of Beasley and Springer (1977)
c     A.J.V., 2005
      real zero,split,half,one
      real a0,a1,a2,a3,b1,b2,b3,b4,c0,c1,c2,c3,d1,d2
      real p,q,r
      zero=0.0e0
      half=0.5e0
      one=1.0e0
      split=0.42e0
      a0=2.50662823884e0
      a1=-18.61500062529e0
      a2=41.39119773534e0
      a3=-25.44106049637e0
      b1=-8.47351093090e0
      b2=23.08336743743e0
      b3=-21.06224101826e0
      b4=3.13082909833e0
      c0=-2.78718931138e0
      c1=-2.29796479134e0
      c2=4.85014127135e0
      c3=2.32121276858e0
      d1=3.54388924762e0
      d2=1.63706781897e0
      ifault=0
      q=p-half
      if (abs(q).gt.split) go to 101
      r=q*q
      ppnda=q*(((a3*r+a2)*r+a1)*r+a0)/
     +	  ((((b4*r+b3)*r+b2)*r+b1)*r+one)
      return
101   r=p
	  if (q.gt.zero) r=one-p
	  if (r.le.zero) go to 102
	  r=sqrt(-alog(r))
	  ppnda=(((c3*r+c2)*r+c1)*r+c0)/
     +	  ((d2*r+d1)*r+one)
	  if (q.lt.zero) ppnda=-ppnda
      return
102	  ifault=1
	  ppnda=99.9
      return
      end


c=======================================================================            
      subroutine rmvnorm(n,mean,sigma,work1,work2,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = SIGMA
c     A.J.V., 2005
      integer n
      real*8 mean(n),sigma(n,n),work1(n*(n+1)/2),work2(n),y(n)

      call cholesky(n,sigma,work1)
      call mvnchol(n,work1,mean,work2,y)
      return
      end      
      
      
c=======================================================================      
      subroutine mvnchol(n,l,mean,work,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = LL', i.e., L is the Cholesky
c     decomposition of the desired variance-covariance structure.
c     WORK is a double precision work vector of at least N elements.
c     The subroutine calls NORMAL for a vector of N iid normal(0,1)
c     deviates (stored in work).  The new variables are calculated as 
c     MEAN + L*WORK.
c     A.J.V., 2005
      implicit none
      integer n,i,j,jj
      real*8 l(n*(n+1)/2),mean(n),work(n),y(n)
      
      call normalvec(n,work)
      
      do i=1,n
         y(i)=mean(i)
      end do
      
      jj = 0
      
      do i = 1,n
         do j = i,n
            jj = jj + 1
            y(j) = y(j) + l(jj)*work(i)
         end do
      end do
      return
      end
 
      
c=======================================================================      
      subroutine normalvec(n,work)
c=======================================================================
c     generates a vector of normal variables
c     A.J.V., 2005
      implicit none
      integer i,n
      real*8 work(n)
      real rnorm
      
      do i=1,n
         work(i)=rnorm(0.d0,1.d0)
      end do
      return
      end

      
c=======================================================================      
      subroutine cholesky(n,a,l)
c=======================================================================      
c     Subroutine to do a Double precision Half stored CHOLesky
c     decomposition.  Calculate Cholesky decomposition of symmetric,
c     positive definite matrix A which is LOWER HALF-STORED.  The matrix
c     l is the output.
c     A.J.V., 2005
      implicit none
      integer n,i,ii,j,jj,k,kk
      real*8 a(n,n),l(n*(n+1)/2)
      real*8 aii,scal,rtemp

      jj=0
      do i=1,n
         do j=i,n
            jj=jj+1
            l(jj)=a(i,j)
         end do
      end do   
      
      ii = 1
      do i = 1,n-1
         aii = l(ii)
         if (aii.le.0.d0) then
            call intpr('matrix is not positive definite',-1,0,0)
            call intpr('in chol subroutine',-1,0,0)            
            return
         end if
         aii = sqrt(aii)
         l(ii) = aii
         scal = 1.d0/aii
         do j = ii+1,ii+n-i
            l(j) = scal*l(j)
         end do

         jj = 1
         do j = ii+1,ii+n-i
            if (l(j).ne.0.d0) then
               rtemp = -1.d0 * l(j)
               kk = ii + jj + n - i
               do k = j,ii+n-i
                  l(kk) = l(kk) + l(k)*rtemp
                  kk = kk + 1
               end do
            end if
            jj = jj + n - i - j + ii + 1
         end do
         ii = ii + n - i + 1
      end do
      aii = l(ii)
      if (aii.le.0.d0) then
          call intpr('matrix is not positive definite',-1,0,0)
          call intpr('in chol subroutine',-1,0,0)            
          return 
      end if
      aii = sqrt(aii)
      l(ii) = aii
      return
      end


c=======================================================================
      subroutine dmvn(n,x,mu,sigma,eval,vv,a,sigma2,res,iflag)        
c=======================================================================
c     return the log of a multivariate normal density
c     A.J.V., 2005
      implicit none
      integer n,i,j
      real*8 mu(n),sigma(n,n),x(n),vv(n)
      real*8 a(n,n),sigma2(n,n),res(n),det,sse,eval
      integer iflag(n)
      real*8 work1,work2,work3,tpi
      
      work1=0.d0
      work2=0.d0
      work3=0.d0
      det=0.d0
      sse=0.d0
	  
      det=0.d0
      
      tpi=6.283185307179586476925286766559d0
       
      work1=-(float(n)*log(tpi))
    
      do i=1,n
         do j=1,n
            a(i,j)=sigma(i,j)
         end do
      end do

      call invdet(a,n,sigma2,det,iflag,vv)

      work2=det
   
      do i=1,n
         res(i)=x(i)-mu(i)
      end do
   
      do i=1,n
         do j=1,n
            sse=sse+res(i)*sigma2(i,j)*res(j)          
         end do
      end do

      work3=sse
     
      eval=(work1-work2-work3)/2.d0
      
      return
      end


c=======================================================================
      subroutine invdet(a,n,ainv,detlog,indx,vv)
c=======================================================================
c     this routine returns the inverse of A in AINV and the log of the
c     determinate of A in DETLOG. N.B: log of determinate only exists
c     (as a real number) if A is a postive definate matrix.  This routine
c     assume that the matrix is a covariance matrix, and therefore it is
c     the dabs of a(j,j) is used.
c     This is from numerical recipies, pages 38 and 39. (The fortran
c     version).
c     This program calls the lu decomposition routines LUDCMP and LUBKSB,
c     see numerical recipies pages 35-37.
c     A.J.V., 2005
      implicit none   
      integer i,j,n,indx(n)
      real*8 a(n,n),ainv(n,n),vv(n),detlog
c
c  set up identity matrix
c
      do 6 i=1,n
	  do 5 j=1,n
	    ainv(i,j) = 0.0d0
 5	    continue
	  ainv(i,i) = 1.0d0
 6	  continue
c
c  decompose the matrix
c
      call dludcmp(a,n,indx,detlog,vv)
c
c  calculate the determinate
c
      detlog = 0.d0
      do 11 j=1,n
	    detlog = detlog + dlog(dabs(a(j,j)))
 11	    continue
c
c  find inverse by columns
c
      do 13 j=1,n
	    call dlubksb(a,n,indx,ainv(1,j))
 13	    continue
      return
      end


c=======================================================================
      subroutine dludcmp(a,n,indx,d,vv)
c=======================================================================
c     This is copied from Press, et. al. {\em Numerical Recipes:
c     Fortran version} pg 35-36
c
c     Given an N x N matrix A, with physical dimension NP, this routine
c     replaces it by the LU decomposition of a rowwise permutation of
c     itself. A and N are input. A is output, arranged as in equation
c     (2.3.14, in Press et. al., pg 34.);	INDX is an output vector which
c     records the row permutation effected by the partial pivoting; D is
c     output as +/-1.d0 depending on wheter the number of row interchanges
c     was even or odd, respectively.  This routine is used in combination
c     with DLUBDSB to solve linear equatoin or invert a matrix.
c     A.J.V., 2005
      implicit none 
      integer i,j,k,n
      real*8 tiny
      parameter (tiny=1.0d-22)
      real*8 a(n,n), vv(n),d,aamax,sum,dum
      integer indx(n),iimax
      
      iimax=1
      d=1.d0
      do 12 i=1,n
	    aamax=0.d0
	    do 11 j=1,n
	       if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11	       continue
	       if (aamax.eq.0.d0) then
                  call intpr('ERROR: inverting a zero matrix',-1,0,0)
                  call intpr('in dludcmp subroutine',-1,0,0)            
                  return 
               end if
	    vv(i)=1.d0/aamax
 12	    continue
	  do 19 j=1,n
	     do 14 i=1,j-1
	        sum=a(i,j)
	        do 13 k=1,i-1
		       sum=sum-a(i,k)*a(k,j)
 13		       continue
	         a(i,j)=sum
 14	         continue
	     aamax = 0.d0
	     do 16 i=j,n
	        sum=a(i,j)
	        do 15 k=1,j-1
		       sum = sum - a(i,k)*a(k,j)
 15		       continue
	        a(i,j)=sum
	        dum=vv(i)*dabs(sum)
	        if(dum.ge.aamax)then
		       iimax=i
		       aamax=dum
		       endif
 16	        continue
	    if(j.ne.iimax)then
	       do 17 k=1,n
		      dum=a(iimax,k)
		      a(iimax,k)=a(j,k)
		      a(j,k)=dum
 17		      continue
	       d=(-1.d0)*d
	       vv(iimax)=vv(j)
	       endif
	    indx(j)=iimax
	    if(a(j,j).eq.0.d0)a(j,j)=tiny
	    if(j.ne.n)then
	       dum=1.d0/a(j,j)
	       do 18 i=j+1,n
		      a(i,j)=a(i,j)*dum
 18		      continue
	       endif
 19	    continue
	  return
	  end


c=======================================================================
       subroutine dlubksb(a,n,indx,b)
c=======================================================================
c      This is copied from Press, et. al. {\em Numerical Recipes:
c      Fortran version} pg 35-36
c
c      Solves the set of N linear equations A*X=B.	
c      Here A is input, not as
c      the matrix A but rather as its LU decomposition, determined by the
c      routine LUDCMP. INDX is input as the permutation vector returned by
c      LUDCMP. B is input as the right-hand side vector B, and returns with
c      the solution vector X. A,N,NP and INDX are not modified by this
c      routine and can be left in place for successive calls with different
c      right hand sides B. This routine tades into account the possibility
c      that B will begin with many zero elements, so it is efficient for use
c      in matrix inversion.
c      A.J.V., 2005
       implicit none
       integer n,i,j,ii,ll
       real*8 a(n,n), b(n), sum
       integer indx(n)

       ii=0
        do 12 i=1,n
	    ll = indx(i)
	    sum=b(ll)
	    b(ll)=b(i)
	    if(ii.ne.0) then
	       do 11 j=ii,i-1
		  sum=sum-a(i,j)*b(j)
 11		  continue
	    else if (sum.ne.0.d0)then
	       ii=i
	    endif
	    b(i)=sum
 12	    continue
	  do 14 i=n,1,-1
	    sum=b(i)
	    do 13 j=i+1,n
	       sum=sum-a(i,j)*b(j)
 13	       continue
	    b(i)=sum/a(i,i)
 14	    continue
	  return
	  end


c=======================================================================      
       subroutine hpd(n,alpha,x,alow,aupp)
c=======================================================================       
c      computing 100(1-alpha)% HPD and credible intervals for x 
c      use Chen-Shao HPD Estimation Algorithm 
c      (see page 219 of Monte Carlo Methods in Bayesian Computation, 
c      Springer-Verlag, 2000) 
c      ming-hui chen
c      july 23, 2001 at wpi
c      input:
c            alpha: confidence level,  0 < alpha < 1
c            n = mcmc sample size
c            x(n): a univariate vector of mcmc sample 
c      output: 
c            (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c            (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible interval
c
       implicit real*8 (a-h,o-z)
       real*8 x(n)
       real*8 aupp(2),alow(2)

       whpd=0.d0
       aupp1=0.d0
       alow1=0.d0

       q1=(alpha/2.0d0)*float(n)
       q2=(1.0d0-alpha/2.0d0)*float(n)
       nq1=nint(q1)
       nq2=nint(q2)
       nq=nq2-nq1
       do 100 i=1,n-1
          do 110 j=i+1,n
             if (x(i) .gt. x(j)) then
                temp=x(i)
                x(i)=x(j)
                x(j)=temp
             end if
 110   continue
 100   continue
       do 120 j=1,n-nq
              pdiff1=x(j)
              pdiff2=x(j+nq)
              wb=pdiff2-pdiff1
              if (j .eq. 1) then
                 whpd=wb
                 aupp1=pdiff2
                 alow1=pdiff1
               else
                 if (whpd .gt. wb) then
                    whpd=wb
                    aupp1=pdiff2
                    alow1=pdiff1
                 end if
              end if
 120   continue
       alow(1)=alow1       
       aupp(1)=aupp1       
       alow(2)=x(nq1)
       aupp(2)=x(nq2)

       return
       end
      
           
