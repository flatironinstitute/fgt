c     This file contains the following subroutines
c
c     g1d_directcp_vec : direct evaluation of the potential
c                     due to sources of strength charge()
c
c     g1d_directcg_vec : direct evaluation of the potential and gradient
c                     due to sources of strength charge()
c
c     g1d_directch_vec : direct evaluation of the potential and gradient and hessian
c                     due to sources of strength charge()
c
c     g1d_directdp_vec : direct evaluation of the potential
c                     due to dipoles
c
c     g1d_directdg_vec : direct evaluation of the potential and gradient
c                     due to dipoles
c
c     g1d_directdh_vec : direct evaluation of the potential and gradient and hessian
c                     due to dipoles
c
c     g1d_directcdp_vec : direct evaluation of the potential
c                     due to charges and dipoles
c
c     g1d_directcdg_vec : direct evaluation of the potential and gradient
c                     due to charges and dipoles
c
c     g1d_directcdh_vec : direct evaluation of the potential and gradient and hessian
c                     due to chrages and dipoles
c
c
c**********************************************************************
      subroutine g1d_directcp_vec(nd,delta,eps,sources,ns,charge,
     $           targ,pot)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
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
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(ns),targ,xdiff,ydiff,rr,r
      real *8 rtmp,delta,eps,dmax
      real *8 pot(nd)
      real *8 charge(nd,ns)
c
      dmax=log(1.0d0/eps)*delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii) = pot(ii) + rtmp*charge(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcg_vec(nd,delta,eps,sources,ns,charge,
     1             targ,pot,grad)
      implicit none
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
      real *8 sources(ns),targ
      real *8 xdiff,ydiff,rr,r,rtmp,delta,eps,dmax
      real *8 pot(nd),grad(nd)
      real *8 dx,dy,rinc,rfac
      real *8 charge(nd,ns)
c
      dmax=log(1.0d0/eps)*delta
      rfac=-2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii) = grad(ii) + dx*rinc
            enddo
         endif
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
      subroutine g1d_directch_vec(nd,delta,eps,sources,ns,charge,targ,
     1           pot,grad,hess)
      implicit none
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
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(ns),targ
      real *8 xdiff,ydiff,rr,r,delta,eps,dmax
      real *8 pot(nd),grad(nd),hess(nd)
      real *8 rtmp,dx,dy,rinc,dxx,dxy,dyy,rfac
      real *8 charge(nd,ns)
c
      dmax=log(1.0d0/eps)*delta
      rfac = -2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dxx = rfac + dx*dx
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii) = grad(ii) + dx*rinc
               hess(ii) = hess(ii) + dxx*rinc
            enddo
         endif
      enddo

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directdp_vec(nd,delta,eps,sources,ns,
     $           rnormal,dipstr,targ,pot)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
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
      real *8 sources(ns),rnormal(ns),targ,xdiff,ydiff,rr,r
      real *8 dx,dy,rtmp,rinc,delta,rfac,eps,dmax
      real *8 pot(nd)
      real *8 dipstr(nd,ns)
c
      dmax=log(1.0d0/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            rinc = dx*rnormal(i)*rtmp
            do ii=1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directdg_vec(nd,delta,eps,sources,ns,
     1           rnormal,dipstr,targ,pot,grad)
      implicit none
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
      real *8 sources(ns),rnormal(ns),targ
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta,rfac
      real *8 pot(nd),grad(nd),eps,dmax
      real *8 dipstr(nd,ns)
c
      dmax=log(1.0d0/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            rinc = dx*rnormal(i)
            rincx = rnormal(i)*rfac - rinc*dx
            rinc = rinc*rtmp
            rincx = rincx*rtmp
c
            do ii = 1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
               grad(ii) = grad(ii) + rincx*dipstr(ii,i)
            enddo
         endif
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
      subroutine g1d_directdh_vec(nd,delta,eps,sources,ns,
     1           rnormal,dipstr,targ,pot,grad,hess)
      implicit none
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
      real *8 sources(ns),rnormal(ns),targ
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta
      real *8 rincxx,rincxy,rincyy,rfac,eps,dmax
      real *8 pot(nd),grad(nd),hess(nd)
      real *8 dipstr(nd,ns)

c
      dmax=log(1.0d0/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff

            rinc = dx*rnormal(i)
            rincx = rnormal(i)*rfac - rinc*dx
            rincxx = -rfac*(dx*rnormal(i)+rinc) - rincx*dx
c     
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii) = grad(ii) + rincx*rtmp
               hess(ii) = hess(ii) + rincxx*rtmp
            enddo
         endif
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcdp_vec(nd,delta,eps,sources,ns,charge,
     $           rnormal,dipstr,targ,pot)
      implicit none
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
      real *8 sources(ns),rnormal(ns),targ
      real *8 xdiff,ydiff,rr,r,rtmp,dx,dy,delta,rfac,rinc
      real *8 pot(nd),eps,dmax
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      dmax=log(1.0d0/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
         
            dx = rfac*xdiff
         
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii) = pot(ii) + rtmp*charge(ii,i)
            enddo

            rinc = dx*rnormal(i)*rtmp
            do ii = 1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g1d_directcdg_vec(nd,delta,eps,sources,ns,charge,
     1             rnormal,dipstr,targ,pot,grad)
      implicit none
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
      real *8 sources(ns),rnormal(ns),targ,delta
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,rfac
      real *8 pot(nd),grad(nd),eps,dmax
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      dmax=log(1.0d0/eps)*delta
      rfac = -2.0d0/delta
      
      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii) = grad(ii) + dx*rinc
            enddo
c
            rinc = -dx*rnormal(i)
            rincx = rinc*dx - rfac*rnormal(i) 
c
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii) = grad(ii) + rincx*rtmp
            enddo
         endif
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
      subroutine g1d_directcdh_vec(nd,delta,eps,sources,ns,charge,
     1           rnormal,dipstr,targ,pot,grad,hess)
      implicit none
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
      real *8 sources(ns),targ,rnormal(ns)
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta
      real *8 rincxx,rincxy,rincyy,dxx,dxy,dyy,rfac,eps,dmax
      real *8 pot(nd),grad(nd),hess(nd)
      real *8 charge(nd,ns),dipstr(nd,ns)


c
      dmax=log(1.0d0/eps)*delta
cccc      call prin2('dmax=*',dmax,1)
      rfac = -2.0d0/delta

      do i = 1,ns
         xdiff=targ-sources(i)
         rr=xdiff*xdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dxx = rfac + dx*dx
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii) = grad(ii) + dx*rinc
               hess(ii) = hess(ii) + dxx*rinc
            enddo

ccc         rr=xdiff*xdiff+ydiff*ydiff
ccc         rtmp = dexp(-rr/delta)
            rinc = -dx*rnormal(i)
            rincx = rinc*dx - rnormal(i)*rfac  
            rincxx = rfac*(rinc-dx*rnormal(i)) + rincx*dx
c
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii) = grad(ii) + rincx*rtmp
               hess(ii) = hess(ii) + rincxx*rtmp
            enddo
         endif
      enddo
      
      return
      end
c
c
c
