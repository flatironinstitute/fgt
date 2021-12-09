c
c**********************************************************************
      subroutine g1d_directcp_vec(nd,delta,dmax,sources,ns,charge,
     $           targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(3,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),targ(ntarg),xdiff,rr,r
      real *8 rtmp,delta
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns)
c
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcg_vec(nd,delta,dmax,sources,ns,charge,
     1             targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target

c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),targ(ntarg)
      real *8 xdiff,rr,r,rtmp,delta
      real *8 pot(nd,ntarg),grad(nd,ntarg)
      real *8 dx,rinc
      real *8 charge(nd,ns)
c
      rfac=-2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               grad(ii,itarg) = grad(ii,itarg) + dx*rinc
            enddo
         endif
      enddo
      enddo
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directch_vec(nd,delta,dmax,sources,ns,charge,
     1           targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c     grad = grad
c     hess = (dxx,dxy,dyy)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c     hess(nd,3)    : Hessian is incremented
c                     using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
      
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(ns),targ(ntarg)
      real *8 xdiff,rr,r,delta
      real *8 pot(nd,ntarg),grad(nd,ntarg),hess(nd,ntarg)
      real *8 rtmp,dx,rinc,dxx
      real *8 charge(nd,ns)
c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dxx = rfac + dx*dx
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
            
               grad(ii,itarg) = grad(ii,itarg) + dx*rinc
               hess(ii,itarg) = hess(ii,itarg) + dxx*rinc
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directdp_vec(nd,delta,dmax,sources,ns,
     $           rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),rnormal(ns),targ(ntarg),xdiff
      real *8 dx,rtmp,rinc,delta,rr,r
      real *8 pot(nd,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            
            rinc = dx*rnormal(i)*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directdg_vec(nd,delta,dmax,sources,ns,
     1             rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),rnormal(ns),targ(ntarg)
      real *8 xdiff,rr,r,rtmp,rinc,dx
      real *8 rincx,delta
      real *8 pot(nd,ntarg),grad(nd,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)

            dx = rfac*xdiff
            rinc = dx*rnormal(i)
            rincx = rnormal(i)*rfac - rinc*dx
c     
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,itarg) = grad(ii,itarg) + rincx*rtmp
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directdh_vec(nd,delta,dmax,sources,ns,
     1           rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c     hess(nd)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(ns),rnormal(ns),targ(ntarg)
      real *8 xdiff,rr,r,rtmp,rinc,dx
      real *8 rincx,delta
      real *8 rincxx
      real *8 pot(nd,ntarg),grad(nd,ntarg),hess(nd,ntarg)
      real *8 dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            
            rinc = dx*rnormal(i)
            rincx = rnormal(i)*rfac - rinc*dx
            
            rincxx = -rfac*(dx*rnormal(i)+rinc) - rincx*dx
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,itarg) = grad(ii,itarg) + rincx*rtmp
               hess(ii,itarg) = hess(ii,itarg) + rincxx*rtmp
            enddo
         endif
      enddo
      enddo

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcdp_vec(nd,delta,dmax,sources,ns,
     $           charge,rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges and dipoles at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),rnormal(ns),targ(ntarg)
      real *8 xdiff,rr,r,rtmp,rinc,dx,delta
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then         
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
c
            dx = rfac*xdiff
            rinc = dx*rnormal(i)*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcdg_vec(nd,delta,dmax,sources,ns,
     1             charge,rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),rnormal(ns),targ(ntarg),delta
      real *8 xdiff,rr,r,rtmp,rinc,dx
      real *8 rincx
      real *8 pot(nd,ntarg),grad(nd,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               grad(ii,itarg) = grad(ii,itarg) + dx*rinc
            enddo
c
            rinc = -dx*rnormal(i)
            rincx = rinc*dx-rnormal(i)*rfac
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,itarg) = grad(ii,itarg) + rincx*rtmp
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcdh_vec(nd,delta,dmax,sources,ns,
     1           charge,rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd)      : gradient is incremented
c     hess(nd)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(ns),targ(ntarg),rnormal(ns)
      real *8 xdiff,rr,r,rtmp,rinc,dx
      real *8 rincx,delta
      real *8 rincxx
      real *8 dxx
      real *8 pot(nd,ntarg),grad(nd,ntarg),hess(nd,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)


c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(itarg)-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            
            dxx = rfac + dx*dx
c
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
            
               grad(ii,itarg) = grad(ii,itarg) + dx*rinc
               hess(ii,itarg) = hess(ii,itarg) + dxx*rinc
            enddo

            rinc = -dx*rnormal(i)
            rincx = rinc*dx-rnormal(i)*rfac
            rincxx = rfac*(rinc-dx*rnormal(i)) + rincx*dx
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,itarg) = grad(ii,itarg) + rincx*rtmp
               hess(ii,itarg) = hess(ii,itarg) + rincxx*rtmp
            enddo
         endif
      enddo
      enddo
      return
      end
c
c
c
