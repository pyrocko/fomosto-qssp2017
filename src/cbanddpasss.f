      subroutine cbandpass(n,f1,f2,df,fi,nf,h)
      implicit none
      integer*4 n,nf
      real*8 f1,f2,df,fi
      complex*16 h(nf)
c
      integer*4 l,k
      real*8 arg
      complex*16 cw,cw1,cw2,cdelt
c
      real*8 PI
      complex*16 ci
      data PI/3.14159265358979d0/
      data ci/(0.d0,1.d0)/
c
      if(n.le.0.or.f1.ge.f2)then
        do l=1,nf
          h(l)=(1.d0,0.d0)
        enddo
      else if(f1.le.0.d0)then
        cw2=dcmplx(-2.d0*PI*fi,2.d0*PI*f2)
        do l=1,nf
          cw=dcmplx(-2.d0*PI*fi,2.d0*PI*dble(l-1)*df)
          cdelt=cw2-cw1
          h(l)=(1.d0,0.d0)
          do k=1,n
            arg=PI*dble(2*k-1)/dble(2*n)
            h(l)=-h(l)*ci*cdelt
     &          /(cw-cdelt*cdexp(dcmplx(0.d0,arg)))
          enddo
        enddo
      else
        cw1=dcmplx(-2.d0*PI*fi,2.d0*PI*f1)
        cw2=dcmplx(-2.d0*PI*fi,2.d0*PI*f2)
        do l=1,nf
          cw=dcmplx(-2.d0*PI*fi,2.d0*PI*dble(l-1)*df)
          cdelt=cw*(cw2-cw1)
          h(l)=(1.d0,0.d0)
          do k=1,n
            arg=PI*dble(2*k-1)/dble(2*n)
            h(l)=-h(l)*ci*cdelt
     &          /(cw*cw-cw1*cw2-cdelt*cdexp(dcmplx(0.d0,arg)))
          enddo
        enddo
      endif
      return
      end
