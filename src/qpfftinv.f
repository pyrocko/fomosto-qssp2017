      subroutine qpfftinv(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
      integer*4 lf,mf,ir
      real*8 f
c
      real*8,allocatable:: y0(:),y(:),dswap(:)
      complex*16,allocatable:: cy(:,:),lpf(:),cswap(:)
c
      allocate(y0(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: y0 not allocated!'
      allocate(y(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: y not allocated!'
      allocate(cy(nf,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: cy not allocated!'
      allocate(lpf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: lpf not allocated!'
      allocate(dswap(4*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: dswap not allocated!'
      allocate(cswap(2*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: cswap not allocated!'
c
      do lf=1,nf
        f=dble(lf-1)*df
        comi=dcmplx(PI2*f,PI2*fi)
        comi2=comi**2
        do ir=1,nr
          gm(lf,ir)=-gz(lf,ir)+(dcmplx(freeairgrd,0.d0)-comi2)*uz(lf,ir)
        enddo
      enddo
c
c     low-pass filter
c
      if(nlpf.gt.0)then
        if(f1corner.le.0.d0)then
          call butterworth(nlpf,f2corner,df,nf,lpf)
        else
          call bandpass(nlpf,f1corner,f2corner,df,nf,lpf)
        endif
        do lf=1,nf
          do ir=1,nr
            ue(lf,ir)=ue(lf,ir)*lpf(lf)
            un(lf,ir)=un(lf,ir)*lpf(lf)
            uz(lf,ir)=uz(lf,ir)*lpf(lf)
c
            ge(lf,ir)=ge(lf,ir)*lpf(lf)
            gn(lf,ir)=gn(lf,ir)*lpf(lf)
            gz(lf,ir)=gz(lf,ir)*lpf(lf)
c
            roe(lf,ir)=roe(lf,ir)*lpf(lf)
            ron(lf,ir)=ron(lf,ir)*lpf(lf)
            roz(lf,ir)=roz(lf,ir)*lpf(lf)
c
            uee(lf,ir)=uee(lf,ir)*lpf(lf)
            uen(lf,ir)=uen(lf,ir)*lpf(lf)
            uez(lf,ir)=uez(lf,ir)*lpf(lf)
            unn(lf,ir)=unn(lf,ir)*lpf(lf)
            unz(lf,ir)=unz(lf,ir)*lpf(lf)
            uzz(lf,ir)=uzz(lf,ir)*lpf(lf)
c
            see(lf,ir)=see(lf,ir)*lpf(lf)
            sen(lf,ir)=sen(lf,ir)*lpf(lf)
            sez(lf,ir)=sez(lf,ir)*lpf(lf)
            snn(lf,ir)=snn(lf,ir)*lpf(lf)
            snz(lf,ir)=snz(lf,ir)*lpf(lf)
            szz(lf,ir)=szz(lf,ir)*lpf(lf)
c
            gm(lf,ir)=gm(lf,ir)*lpf(lf)
          enddo
        enddo
      endif
c
      if(icmp(1).eq.1)then
c
c       output seismograms: displacement vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ue,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,dispout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,un,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,dispout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,dispout(3))
      endif
      if(icmp(2).eq.1)then
c
c       output seismograms: velocity vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ue,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,veloout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,un,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,veloout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,veloout(3))
      endif
      if(icmp(3).eq.1)then
c
c       output seismograms: acceleration vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ue,cy,cswap,dswap,
     &                 2,ntcutout,rname,y0,y,acceout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,un,cy,cswap,dswap,
     &                 2,ntcutout,rname,y0,y,acceout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uz,cy,cswap,dswap,
     &                 2,ntcutout,rname,y0,y,acceout(3))
      endif
c
      if(icmp(4).eq.1)then
c
c       output seismograms: strain tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uee,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uen,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uez,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(3))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,unn,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(4))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,unz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(5))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uzz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,strainout(6))
      endif
c
      if(icmp(5).eq.1)then
c
c       output seismograms: strain rate tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uee,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uen,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uez,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(3))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,unn,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(4))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,unz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(5))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,uzz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,strainrateout(6))
      endif
      if(icmp(6).eq.1)then
c
c       output seismograms: stress tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,see,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,sen,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,sez,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(3))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,snn,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(4))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,snz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(5))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,szz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,stressout(6))
      endif
c
      if(icmp(7).eq.1)then
c
c       output seismograms: stress rate tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,see,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,sen,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,sez,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(3))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,snn,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(4))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,snz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(5))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,szz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,stressrateout(6))
      endif
c
      if(icmp(8).eq.1)then
c
c       output seismograms: rotation vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,roe,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,rotaout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ron,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,rotaout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,roz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,rotaout(3))
      endif
c
      if(icmp(9).eq.1)then
c
c       output seismograms: rotation rate vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,roe,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,rotarateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ron,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,rotarateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,roz,cy,cswap,dswap,
     &                 1,ntcutout,rname,y0,y,rotarateout(3))
      endif
      if(icmp(10).eq.1)then
c
c       output seismograms: local gravitational acceleration vector
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,ge,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,gravout(1))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,gn,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,gravout(2))
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,gz,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,gravout(3))
      endif
c
      if(icmp(11).eq.1)then
c
c       output seismograms: gravimeter response
c
        call transfs2t(nf,nr,df,togsmin,dt,fi,tred,gm,cy,cswap,dswap,
     &                 0,ntcutout,rname,y0,y,grmout)
          
      endif
      return
      end
