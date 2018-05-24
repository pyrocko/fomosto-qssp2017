      subroutine swavelet(tau,df,nf,wvf)
      implicit none
      integer*4 nf
      real*8 tau,df
      complex*16 wvf(nf)
c
      integer*4 l,n
      real*8 f,omi,x,dt0
      complex*16 alfa,beta,gamma,eta
c
      real*8 pi,pi2,eps
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-04/
c
c     for wavelet: normalized square half-sinus
c
      do l=1,nf
        f=df*dble(l-1)
        x=f*tau
        if(x.eq.0.d0)then
          wvf(l)=(1.d0,0.d0)
        else if(x.ge.1.d0-eps.and.x.le.1.d0+eps)then
          wvf(l)=dcmplx(-1.d0/x/(1+x),0.d0)
        else
          wvf(l)=dcmplx(0.d0,1.d0/(pi2*x*(1.d0+x)*(1.d0-x)))
     &             *(cdexp(dcmplx(0.d0,-pi2*x))-(1.d0,0.d0))
        endif
      enddo
      return
      end
