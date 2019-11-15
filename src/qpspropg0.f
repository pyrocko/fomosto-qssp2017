      subroutine qpspropg0(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of spheroidal response (with gravity) for degree 0
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,istp,ly,ily,nly,key
      real*8 f,rr1,rr2,dlnr,h
      complex*16 y0(3),yup(3),ylw(3),coef(2,2),b(2,2)
      external qpdifmat0
c
      if(lylwa.gt.lylw)return
c
c
c===============================================================================
c
c     propagation from surface to source
c
      if(freesurf.and.lyup.eq.1)then
        yup(1)=(1.d0,0.d0)
        yup(2)=(0.d0,0.d0)
        yup(3)=(0.d0,0.d0)
      else
        call qpstart0g(lyup,2,yup)
      endif
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,3)
c
      f=dreal(comi)/PI2
c
      do ly=lyup,lys-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
          call ruku(yup,3,1,ly,0,qpdifmat0,rr1,rr2,nruku(0,ly))
        enddo
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,3)
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c     yup(3)=yup(3)
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lylwb)then
        call qpstart0g(lylw,1,ylw)
      else
        call qpstart0g(lylw,0,ylw)
      endif
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylw-1,lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call ruku(ylw,3,1,ly,0,qpdifmat0,rr1,rr2,nruku(0,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      y0(1)=y0(1)/crrup(lyr)
      y0(2)=y0(2)/crrup(lyr)**2
c     y0(3)=y0(3)
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpspropg0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
          ypsv(5,istp)=b(1,istp)*y0(3)
          ypsv(6,istp)=ypsv(5,istp)/crrup(lyr)
        enddo
      else
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
          ypsv(5,istp)=b(2,istp)*(y0(3)-ylw(3))+b(1,istp)*yup(3)
          ypsv(6,istp)=ypsv(5,istp)/crrup(lyr)
        enddo
      endif
c
      if(lylwa.le.0)return
c
c
c===============================================================================
c
c     propagation from bottom to source
c
      call qpstart0g(lylwa,1,ylw)
c
      if(lylwa.eq.lyr.and.lylwa.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylwa-1,lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call ruku(ylw,3,1,ly,0,qpdifmat0,rr1,rr2,nruku(0,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      if(lyr.gt.lys)then
        y0(1)=y0(1)/crrup(lyr)
        y0(2)=y0(2)/crrup(lyr)**2
c       y0(3)=y0(3)
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsprop0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=ypsv(i,istp)-b(1,istp)*y0(i)
          enddo
          ypsv(5,istp)=ypsv(5,istp)-b(1,istp)*y0(3)
          ypsv(6,istp)=ypsv(6,istp)-b(1,istp)*y0(3)/crrup(lyr)
        enddo
      else
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=ypsv(i,istp)-b(2,istp)*y0(i)
          enddo
c
c         add a constant to the potential below the source so that
c         it becomes continuous through the source level.
c
          ypsv(5,istp)=ypsv(5,istp)
     &            -b(2,istp)*(y0(3)-ylw(3))-b(1,istp)*yup(3)
          ypsv(6,istp)=ypsv(6,istp)
     &            -b(2,istp)*(y0(3)-ylw(3))-b(1,istp)*yup(3)/crrup(lyr)
        enddo
      endif
      return
      end