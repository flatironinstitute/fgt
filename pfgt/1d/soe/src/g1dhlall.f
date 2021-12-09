c     This file contains the following subroutines
c
c     g1dformhc_vec : computes the Hermite expansion about the center CENT
c                     due to sources of strength charge()
c
c     g1dformhd_vec : computes the Hermite expansion about the center CENT
c                     due to dipoles
c
c     g1dformhcd_vec : computes the Hermite expansion about the center CENT
c                     due to charges and dipoles
c
C     g1dhevalp_vec : evaluates the Hermite expansion 
C                         potential only.
c      
C     g1dhevalg_vec : evaluates the Hermite expansion 
C                         pot/grad
c      
C     g1dhevalh_vec : evaluates the Hermite expansion 
C                         pot/grad/Hessian
      
c     g1dformlc_vec : computes the Taylor expansion about the center CENT
c                     due to sources of strength charge()
c
c     g1dformld_vec : computes the Taylor expansion about the center CENT
c                     due to dipoles
c
c     g1dformlcd_vec : computes the Taylor expansion about the center CENT
c                     due to charges and dipoles
c
C     g1dlevalp_vec : evaluates the local Taylor expansion 
C                         potential only.
c      
C     g1dlevalg_vec : evaluates the local Taylor expansion 
C                         pot/grad
c      
C     g1dlevalh_vec : evaluates the local Taylor expansion 
C                         pot/grad/Hessian
c      
C     g1dh2sx_h2s_vec : converts the Hermite expansion to the four hybrid
C                       SOE/X expansions and four SOE expansions
c
c     g1dherm2local_vec : converts the Hermite expansion to the Taylor expansion
c      
c     g1dsoe2local_vec : converts four SOE expansions to the local Taylor expansion
c
c     g1dsx2local_vec : converts four SOE/X expansions to the local Taylor expansion
c
c     g1dh2lmat : returns the matrix converting the Hermite expansion
C                 to the local Taylor expansion (1D only)
c      
c     g1dh2smat : returns the matrix converting the Hermite expansion
C                  to the SOE+ expansions (1D only)
C     
c     g1dh2xmat : returns the matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C     
C     g1ds2lmat : returns the matrix converting the SOE+ expansion
C                 to the local Taylor expansion (1D only)
C
C     g1dx2lmat : returns the matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only)
C
c
c*********************************************************************
C
C form Hermite expansions (charge, dipole, charge & dipole)
C
C*********************************************************************
      subroutine g1dformhc_vec(nd,delta,sources,ns,charge,cent,
     1            nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     charge    = strengths of sources
C     cent      = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,nterms
      real *8 cent,sources(ns),charge(nd,ns)
      real *8 ffexp(0:nterms,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xp(:)
C
C     initialize coefficients to zero.
C
      allocate(xp(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(i) - cent)*dsq
         xp(0) = 1.0D0
         do k = 1,nterms
            tmp = dsqrt(dble(k))
            xp(k) = xp(k-1)*x/tmp
         enddo
c
         do ind=1,nd
            do j=0,nterms
               ffexp(j,ind) = ffexp(j,ind) + charge(ind,i)*xp(j)
            enddo
         enddo
      enddo
      return
      end
c
C
C
C*********************************************************************C
      subroutine g1dformhd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the dipole sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd        = vector length (parallel)
C     delta     = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center    = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center,sources(ns)
      real *8 rnormal(ns),dipstr(nd,ns)
      real *8 ffexp(0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,d1,d2
      real *8, allocatable ::  hx(:)
      real *8, allocatable ::  dxhx(:)
C
C     initialize coefficients to zero.
C
      allocate(hx(0:nterms))
      allocate(dxhx(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(i) - center)*dsq
         hx(0) = 1.0D0
         dxhx(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hx(j1)=hx(j1-1)*x/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            do j=0,nterms
               ffexp(j,ind) = ffexp(j,ind) + 
     1             dxhx(j)*rnormal(i)*dipstr(ind,i)
            enddo
         enddo
      enddo
      return
      end
C
C
C
C*********************************************************************C
      subroutine g1dformhcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the dipole sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd        = vector length (parallel)
C     delta     = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     charge    = charge strengths
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center    = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center,sources(ns)
      real *8 rnormal(ns),dipstr(nd,ns),charge(nd,ns)
      real *8 ffexp(0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,chg,ytmp,d1,d2
      real *8, allocatable ::  hx(:)
      real *8, allocatable ::  dxhx(:)
C
C     initialize coefficients to zero.
C
      allocate(hx(0:nterms))
      allocate(dxhx(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(i) - center)*dsq
         hx(0) = 1.0D0
         dxhx(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hx(j1)=hx(j1-1)*x/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            do j=0,nterms
               ffexp(j,ind) = ffexp(j,ind) + 
     1             dxhx(j)*rnormal(i)*dipstr(ind,i)+hx(j)*charge(ind,i)
            enddo
         enddo
      enddo
      return
      end
c
c
c
c
c*********************************************************************
C
C eval Hermite expansions (pot, pot+grad, pot+grad+hess)
C
C*********************************************************************
      subroutine g1dhevalp_vec(nd,delta,cent,nterms,ffexp
     1    targ,ntarg,pot)
C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location TARG.
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     POT = evaluated potential
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms))
      allocate(hexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms-1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            do k = 0,nterms
               d2=0
               do j = 0,nterms
                  d2=d2+hexpx(j)*ffexp(j,k,ind)
               enddo
               dd=dd+d2*hexpy(k)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dd
         enddo
      enddo
      
      return
      end
C
c
C***********************************************************************
      subroutine g1dhevalg_vec(nd,delta,cent,nterms,ffexp
     1           targ,ntarg,pot,grad)

C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location (XP,YP).
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     pot      = evaluated potential
C     grad     = evaluated gradient
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg),grad(nd,2,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+1))
      allocate(dhexpx(0:nterms))
      allocate(hexpy(0:nterms+1))
      allocate(dhexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            gx = 0.0d0
            gy = 0.0d0
            do k = 0,nterms
               d2=0
               g1=0
               g2=0
               do j = 0,nterms
                  d2=d2+hexpx(j)*ffexp(j,k,ind)
                  g1 = g1 + dhexpx(j)*ffexp(j,k,ind)
                  g2 = g2 + hexpx(j)*ffexp(j,k,ind)
               enddo
               dd=dd+d2*hexpy(k)
               gx=gx+g1*hexpy(k)
               gy=gy+g2*dhexpy(k)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dd
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
         enddo
      enddo
      
      return
      end
C
C
c
C***********************************************************************
      subroutine g1dhevalh_vec(nd,delta,cent,nterms,ffexp,
     1           targ,ntarg,pot,grad,hess)

C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location (XP,YP).
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     pot      = evaluated potential
C     grad     = evaluated gradient
C     hess     = evaluated Hessian
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg),grad(nd,2,ntarg),hess(nd,3,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:)
      real *8, allocatable ::  ddhexpx(:),ddhexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+2))
      allocate(dhexpx(0:nterms))
      allocate(ddhexpx(0:nterms))
      allocate(hexpy(0:nterms+2))
      allocate(dhexpy(0:nterms))
      allocate(ddhexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms+1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            ddhexpx(i)=hexpx(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
            ddhexpy(i)=hexpy(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
         enddo
c
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            gx = 0.0d0
            gy = 0.0d0
            hxx = 0.0d0
            hyy = 0.0d0
            hxy = 0.0d0
            do k = 0,nterms
               d2=0
               g1=0
               g2=0
               h11=0
               h22=0
               h12=0
               do j = 0,nterms
                  d2  = d2  +   hexpx(j)*ffexp(j,k,ind)
                  g1  = g1  +  dhexpx(j)*ffexp(j,k,ind)
                  g2  = g2  +   hexpx(j)*ffexp(j,k,ind)               
                  h11 = h11 + ddhexpx(j)*ffexp(j,k,ind)
                  h22 = h22 +   hexpx(j)*ffexp(j,k,ind)
                  h12 = h12 +  dhexpx(j)*ffexp(j,k,ind)
               enddo
               dd  = dd  + d2   *hexpy(k)
               gx  = gx  + g1   *hexpy(k)
               gy  = gy  + g2  *dhexpy(k)
               hxx = hxx + h11  *hexpy(k)
               hyy = hyy + h22*ddhexpy(k)
               hxy = hxy + h12 *dhexpy(k)
            enddo
         
            pot(ind,itarg)    = pot(ind,itarg) + dd
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
            hess(ind,1,itarg) = hess(ind,1,itarg) + hxx*dsq*dsq
            hess(ind,2,itarg) = hess(ind,2,itarg) + hxy*dsq*dsq
            hess(ind,3,itarg) = hess(ind,3,itarg) + hyy*dsq*dsq
         enddo
      enddo
      
      return
      end
C
C************************************************************************
C
C     form local subroutines (charge, dipole, charge+dipole)
C
C***********************************************************************
      subroutine g1dformlc_vec(nd,delta,sources,ns,charge,center,
     1            nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns),charge(nd,ns)
      real *8 local(0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:)
c
      allocate(hexpx(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(i) - center)*dsq
c

         facx=dexp(-x*x)

         hexpx(0)=1.0d0*facx
         hexpx(1)=2.0d0*x*facx

         do j=1,nlocal-1
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
         enddo

         do ind=1,nd
            do j2=0,nlocal
               local(j2,ind) = local(j2,ind)+charge(ind,i)*hexpx(j2)
            enddo
         enddo
      enddo
      return
      end
C
C
C
      subroutine g1dformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns),dipstr(nd,ns)
      real *8 rnormal(ns)
      real *8 local(0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:)
      real *8, allocatable :: dhexpx(:),dhexpy(:)
c
      allocate(hexpx(0:nlocal+1))
      allocate(dhexpx(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(i) - center)*dsq
c
         facx=dexp(-x*x)

         hexpx(0)=1.0d0*facx
         hexpx(1)=2.0d0*x*facx

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            r1 = rnormal(i)*dipstr(ind,i)
            do j2=0,nlocal
               local(j2,ind) = local(j2,ind)+
     1             dhexpx(j2)*r1
            enddo
         enddo
      enddo
      return
      end
C
C
C
C
C
      subroutine g1dformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns),dipstr(nd,ns)
      real *8 rnormal(ns), charge(nd,ns)
      real *8 local(0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:)
      real *8, allocatable :: dhexpx(:)
c
      allocate(hexpx(0:nlocal+1))
      allocate(dhexpx(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(i) - center)*dsq
c
         facx=dexp(-x*x)

         hexpx(0)=1.0d0*facx
         hexpx(1)=2.0d0*x*facx

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            chg=charge(ind,i)
            r1 = rnormal(i)*dipstr(ind,i)
            
            do j2=0,nlocal
               local(j2,ind) = local(j2,ind)+
     1             dhexpx(j2)*r1+chg*hexpx(j2)
            enddo
         enddo
      enddo
      return
      end
C
C
C
C
C************************************************************************
C
C     eval local subroutines (pot, pot+grad, pot+grad+hess)
C
C***********************************************************************
C
      subroutine g1dlevalp_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,nd)
      real *8 center,targ(ntarg)
      real *8 pot(nd,ntarg),xp(0:200)
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(itarg) - center)*dsq

         xp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
         enddo
c
         do ind=1,nd
            d2=0
            do j2=0,nlocal
               d2=d2+local(j2,ind)*xp(j2)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+d2
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g1dlevalg_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot,grad)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,nd)
      real *8 center,targ(ntarg)
      real *8 pot(nd,ntarg)
      real *8 grad(nd,ntarg)
      real *8, allocatable :: xp(:)
      real *8, allocatable :: xpx(:)
C
      allocate(xp(0:nlocal))
      allocate(xpx(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(itarg) - center)*dsq

         xp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
         enddo
c
         do ind=1,nd
            d2=0
            d2x=0
            do j2=0,nlocal
               d2=d2+local(j2,ind)*xp(j2)
               d2x=d2x+local(j2,ind)*xpx(j2)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+d2
            grad(ind,itarg) = grad(ind,itarg)+d2x
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g1dlevalh_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,nd)
      real *8 center,targ(ntarg)
      real *8 pot(nd,ntarg),xp(0:200)
      real *8 grad(nd,ntarg)
      real *8 hess(nd,ntarg)
      real *8, allocatable :: xpx(:)
      real *8, allocatable :: xpxx(:)
C
      allocate(xpx(0:nlocal))
      allocate(xpxx(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(itarg) - center)*dsq

         xp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
         enddo
c
         xpxx(0) = 0.0d0
         xpxx(1) = 0.0d0
         do i = 2,nlocal
            xpxx(i) = i*(i-1)*xp(i-2)*dsq*dsq
         enddo

         do ind=1,nd
            d2=0
            d2x=0
            d2xx=0
            do j2=0,nlocal
               d2=d2+local(j2,ind)*xp(j2)
               d2x=d2x+local(j2,ind)*xpx(j2)
               d2xx=d2xx+local(j2,ind)*xpxx(j2)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+d2
            grad(ind,itarg) = grad(ind,itarg)+d2x
            hess(ind,itarg) = hess(ind,itarg)+d2xx
         enddo
      enddo

      return
      end
C
C***********************************************************************
C
C 2D translations (h2l, h2sx_h2s_real, h2sx_h2s)
C
c***********************************************************************      
      subroutine g1dherm2local_vec(nd,nlocal,ntermsh,h2l,
     1    hexp,local)
C
C     This subroutine converts the Hermite expansion to 
C     local Taylor expansion
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     ntermsh       = number of terms in Hermite exp
C     h2l           = matrix converting the Hermite expansion to the local expansion
C     hexp          = Hermite expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      real *8 h2l(0:ntermsh,0:nlocal)
      real *8 hexp(0:ntermsh,nd)
      real *8 local(0:nlocal,nd)
      
      do ind=1,nd
         do k1=0,nlocal,2
            cd=0
            do j1=0,ntermsh,2
               cd=cd+h2l(j1,k1)*hexp(j1,ind)
            enddo
            local(k1,ind)=local(k1,ind)+cd
         enddo

         do k1=1,nlocal,2
            cd=0
            do j1=1,ntermsh,2
               cd=cd+h2l(j1,k1)*hexp(j1,ind)
            enddo
            local(k1,ind)=local(k1,ind)+cd
         enddo
      enddo
         
      return
      end
C
C
C
C
C
      subroutine g1dh2s_vec(nd,ntermsh,nsoe,h2s,
     1    hexp,soeall)
C
C     This subroutine converts the Hermite expansion to the four hybrid
C     SOE/X expansions  and four SOE expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nsoe          = number of exponentials in SOE 
C     h2s           = matrix converting the Hermite expansion to the SOE+ expansion (1D)
C     hexp          = Hermite expansion
C      
C     OUTPUT:
C
C     soeall        = two (unrolled) SOE expansions are incremented
c                     in the order: soep,soem
C      
      implicit real*8 (a-h,o-z)
      complex *16 h2s(nsoe,0:ntermsh)
      real *8 hexp(0:ntermsh,nd)
      complex *16 soeall(2*nsoe/2,nd)

      complex *16 pe,me,po,mo
      complex *16 px,mx,xp,xm
c
      dsq = 1.0D0/dsqrt(delta)
C
      nexp = nsoe/2
      
      do ind=1,nd
         do j1=1,nsoe/2
            pe=0
            po=0
               
            do k=0,ntermsh,2
               pe=pe+h2s(j1,k)*hexp(k,ind)
            enddo
            do k=1,ntermsh,2
               po=po+h2s(j1,k)*hexp(k,ind)
            enddo
            soeall(j1,ind) = soeall(j1,ind) + pe+po
            soeall(nexp+j1,ind) = soeall(nexp+j1,ind)+pe-po
         enddo
      enddo
         
      return
      end
C
C
C
C
C
C***********************************************************************
C
C 2D translations (soe2local, sx2local)
C
c***********************************************************************      
      subroutine g1dsoe2local_vec(nd,nlocal,nsoe,s2l,
     1    soeall,local)
C
C     This subroutine converts four SOE expansions to local
C     expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     nsoe          = number of exponentials in SOE 
C     s2l           = matrix converting SOE expansions to the local expansion
C     soeall        = SOE expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      complex *16 s2l(nsoe,0:nlocal)
      real *8 local(0:nlocal,nd)
      complex *16 soeall(2*nsoe/2,nd)
      complex *16 cp,cm,cd,yp,ym
c      
c
      nexp = nsoe/2
      do ind=1,nd
         do k=0,nlocal,2
            cd=0
            do j=1,nsoe/2
               cd=cd+s2l(j,k)*(soeall(j,ind)+soeall(nexp+j,ind))
            enddo
            local(k,ind)=local(k,ind)+dble(cd)         
         enddo
         
         do k=1,nlocal,2
            cd=0
            do j=1,nsoe/2
               cd=cd+s2l(j,k)*(soeall(j,ind)-soeall(nexp+j,ind))
            enddo
            local(k,ind)=local(k,ind)+dble(cd)         
         enddo
      enddo
         
      return
      end
C
C
C
C
C
C************************************************************************
C
C 1d translation matrices
C
C**********************************************************************
      subroutine g1dh2lmat(nlocal,ntermsh,h2l)
C
C     This subroutine returns the matrix converting the 
C     Hermite expansion to the local expansion
C     
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     ntermsh       = number of terms in the Hermite exp
C
C     OUTPUT:
C
C     h2l        = translation matrix
C
      implicit real*8 (a-h,o-z)
      real *8 h2l(0:ntermsh,0:nlocal)
      real *8, allocatable :: sqc(:,:),fac(:),hexpx(:)

      ntot = ntermsh+nlocal
      allocate(sqc(0:ntot,0:ntot),fac(0:nlocal),hexpx(0:ntot))
      
      call getsqrtbinomialcoeffs(ntot,sqc)
      
      fac(0)=1.0d0
      do i=1,nlocal
         fac(i)=fac(i-1)*sqrt(i*1.0d0)
      enddo

      x=0
      hexpx(0)=1.0d0
      hexpx(1)=2*x
      
      do i=1,ntermsh+nlocal-1
         hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
      enddo

      do i=0,ntermsh
         do j=0,nlocal
            h2l(i,j)=(-1)**j*hexpx(i+j)*sqc(i+j,i)/fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C
      subroutine g1dh2smat(ntermsh,nsoe,ws,ts,h2s)
C
C     This subroutine returns a matrix converting the Hermite expansion
C     to the SOE+ expansions
C
C     INPUT
C
C     ntermsh       = number of terms in Hermite exp
C     nsoe          = number of exponentials in SOE 
c
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     h2s           = translation matrix
C
      implicit real*8 (a-h,o-z)
      complex *16 ws(*), ts(*)
      complex *16 h2s(nsoe,0:ntermsh)
      real *8, allocatable:: fac(:)

      allocate(fac(0:ntermsh))

      fac(0)=1.0d0
      do i=1,ntermsh
         fac(i)=fac(i-1)/sqrt(dble(i))
cccc         fac(i)=fac(i-1)/i
      enddo

      do j=0,ntermsh
         do i=1,nsoe
            h2s(i,j)=ws(i)*ts(i)**j*fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C
      subroutine g1ds2lmat(nlocal,nsoe,ws,ts,s2l)
C
C     This subroutine returns the matrices converting the SOE+ expansions
C     to the local expansion
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     nsoe          = number of exponentials in SOE 
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     s2l        = translation matrix
C
      implicit real*8 (a-h,o-z)
      complex *16 ws(*), ts(*)
      complex *16 s2l(nsoe,0:nlocal)
      real *8, allocatable:: fac(:)

      allocate(fac(0:nlocal))
      
      fac(0)=1.0d0
      do i=1,nlocal
         fac(i)=fac(i-1)*i
      enddo
      
      do i=1,nsoe
         do j=0,nlocal
            s2l(i,j)=(-ts(i))**j/fac(j)
         enddo
      enddo

      return
      end
C
C
C
C
C***********************************************************************
      subroutine g1dhermzero_vec(nd,hexp,ntermsh)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector multipole expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     ntermsh :   order of multipole expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     hexp  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,ntermsh,nd,ii
      real *8 hexp(ntermsh+1,nd)
c
      do ii=1,nd
      do n=1,ntermsh+1
         hexp(n,ii)=0.0d0
      enddo
      enddo
      return
      end
c      
c      
c      
c      
C***********************************************************************
      subroutine g1dlocalzero_vec(nd,local,nlocal)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector local expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     nlocal :   order of local expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nlocal,nd,ii
      real *8 local(nlocal+1,nd)
c
      do ii=1,nd
      do n=1,nlocal+1
         local(n,ii)=0.0d0
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
      subroutine getsqrtbinomialcoeffs(n,dc)
C***********************************************************************
c
c     computes the sqrt of the binomial coefficients, from fmmcommon.f
c-----------------------------------------------------------------------
      implicit none
      integer i,j,n
      real *8 dc(0:n,0:n)
      real *8, allocatable :: d(:,:)
      allocate(d(0:n,0:n))

      do i=0,n
        do j=0,n
          d(i,j) = 0
          dc(i,j) = 0
        enddo
      enddo

      do i=0,n
        d(i,0) = 1
        dc(i,0) = 1
      enddo
      do i=1,n
        d(i,i) = 1
        dc(i,i) = 1
        do j=i+1,n
          d(j,i) =d(j-1,i)+d(j-1,i-1)
          dc(j,i) = sqrt(d(j,i))
        enddo
      enddo

      return
      end


