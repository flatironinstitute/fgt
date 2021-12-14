c     This file contains the following subroutines
c
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
c     g1dformpwc_vec : computes the PW expansion due to charges directly
c      
c     g1dformpwd_vec : computes the PW expansion due to dipoles directly
c                      
c      
c     g1dformpwcd_vec : computes the PW expansion directly due to charges and dipoles
c                      
c      
c     g1dpwevalp_vec : evaluates the PW expansion directly
c                      potential only
c      
c     g1dpwevalg_vec : evaluates the PW expansion directly
c                      potential + gradient
c      
c     g1dpwevalh_vec : evaluates the PW expansion directly
c                      potential + gradient + hessian
c 
c     g1dformpwc_fast_vec : computes the PW expansion due to charges
c                      type 1 NUFFT is called once
c      
c     g1dformpwd_fast_vec : computes the PW expansion due to dipoles
c                      type 1 NUFFT is called three times using finufft3d1many
c      
c     g1dformpwcd_fast_vec : computes the PW expansion due to charges and dipoles
c                      type 1 NUFFT is called four times using finufft3d1many
c      
c     g1dpwevalp_fast_vec : evaluates the PW expansion
c                      potential only, type 2 NUFFT is called once
c                      using finufft3d2many
c      
c     g1dpwevalg_fast_vec : evaluates the PW expansion
c                      potential + gradient
c                      4 type 2 NUFFTs are called using finufft3d2many 
c      
c     g1dpwevalh_fast_vec : evaluates the PW expansion
c                      potential + gradient + hessian
c                      10 type 2 NUFFTs are called using finufft3d2many
c      
c     get_pwnodes : returns PW weights and nodes, midpoint rule is used
c                   so the number is always an even number!
c
c     pw_translation_matrices : returns precomputed translation matrices
c                   for PW mp to loc translations
c
c     g1dshiftpw_vec : translates PW expansion via precomputed translation 
c                      matrices
c
c     g1dcopypwexp_vec : copy PW expansion
c
c
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
C     sources(3,ns) = coordinates of sources
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
            chg=charge(ind,i)
            do j3=0,nlocal
               local(j3,ind) = local(j3,ind)
     1             +chg*hexpx(j3)
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
C     sources(3,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
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
            r1 = rnormal(i)*dipstr(ind,i)
            do j3=0,nlocal
               local(j3,ind) = local(j3,ind)+
     1             dhexpx(j3)*r1
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
C     sources(3,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
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
            
            do j3=0,nlocal
               local(j3,ind) = local(j3,ind)+
     1             dhexpx(j3)*r1+hexpx(j3)*chg
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
C     targ          = target location
C     ntarg         = number of targets
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,nd)
      real *8 center,targ(ntarg)
      real *8 pot(nd,ntarg),xp(0:100)
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
            do j3=0,nlocal
               d2=d2+local(j3,ind)*xp(j3)
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
C     targ          = target location
C     ntarg         = number of targets
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
            do j3=0,nlocal
               d2=d2+local(j3,ind)*xp(j3)
               d2x=d2x+local(j3,ind)*xpx(j3)
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
C     targ          = target location
C     ntarg         = number of targets
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
      real *8 hess(nd,ntarg)
      real *8, allocatable :: xp(:)
      real *8, allocatable :: xpx(:)
      real *8, allocatable :: xpxx(:)
C
      allocate(xp(0:nlocal))
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
            do j3=0,nlocal
               d2=d2+local(j3,ind)*xp(j3)
               d2x=d2x+local(j3,ind)*xpx(j3)
               d2xx=d2xx+local(j3,ind)*xpxx(j3)
            enddo

            pot(ind,itarg) = pot(ind,itarg)+d2
            
            grad(ind,itarg) = grad(ind,itarg)+d2x

            hess(ind,itarg) = hess(ind,itarg)+d2xx
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
C     form PW exp subroutines (charge, dipole, charge+dipole)
C
C***********************************************************************
      subroutine g1dformpwc_vec(nd,delta,sources,ns,charge,center,
     1    npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
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
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns),charge(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,dsq
      complex *16 eye
      complex *16 qqx,qq1
      
      complex *16 ww1(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(i) - center)*dsq
         qqx = cdexp(-eye*ts(npw2+1)*x)
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1*ws(j1)
            qq1 = qq1*qqx
            ww1(npw-j1+1) = dconjg(ww1(j1))
         enddo
c
         do ind = 1,nd
            chg=charge(ind,i)
            do j3=1,npw/2
               pwexp(j3,ind)=pwexp(j3,ind)+
     1             chg*ww1(j3)
            enddo
         enddo
      enddo

      return
      end
C
C
C
      subroutine g1dformpwd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1    center,npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
C     the center CENT due to sources of dipoles.
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
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns)
      real *8 rnormal(ns),dipstr(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye,r1,rz1,rz2,rz3,ry1,ry2,ry3
      complex *16 qqx,qq1,ztmp
      
      complex *16 ww1(100),ww1x(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(i) - center)*dsq
c
         qqx = cdexp(-eye*ts(npw2+1)*x)
         
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1*ws(j1)
            qq1 = qq1*qqx
            
            ww1(npw-j1+1) = dconjg(ww1(j1))
         enddo

         do j1=1,npw
            ztmp = -eye*ts(j1)*dsq
            ww1x(j1) = ww1(j1)*ztmp
         enddo
c
         do ind = 1,nd
            r1 = rnormal(i)*dipstr(ind,i)
            do j1=1,npw/2
               pwexp(j1,ind)=pwexp(j1,ind)+
     1             ww1x(j1)*r1
            enddo
         enddo
      enddo

      return
      end
C
C
C
      subroutine g1dformpwcd_vec(nd,delta,sources,ns,charge,
     1    rnormal,dipstr,center,npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
C     the center CENT due to charges and dipoles.
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
C     center        = center of the expansion
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center,sources(ns),charge(nd,ns)
      real *8 rnormal(ns),dipstr(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye,r1,r2,r3,rz1,rz2,rz3,ry1,ry2,ry3
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3,ztmp
      
      complex *16 ww1(100),ww1x(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(i) - center)*dsq
c
         qqx = cdexp(-eye*ts(npw2+1)*x)
         
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1*ws(j1)
            qq1 = qq1*qqx
            ww1(npw-j1+1) = dconjg(ww1(j1))
         enddo

         do j1=1,npw
            ztmp = -eye*ts(j1)*dsq
            ww1x(j1) = ww1(j1)*ztmp
         enddo
c
         do ind = 1,nd
            chg = charge(ind,i)
            r1 = rnormal(i)*dipstr(ind,i)
            do j1=1,npw/2
               pwexp(j1,ind)=pwexp(j1,ind)+
     1             ww1x(j1)*r1+ww1(j1)*chg
            enddo
         enddo
      enddo

      return
      end
C
C
C
c*********************************************************************
C
C evaluate PW expansions (potential, pot + grad, pot + grad + hess)
C
C*********************************************************************
C
C
c
C
C
C
      subroutine g1dpwevalp_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential only
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C
      implicit none
      integer nd,npw,ntarg
      real *8 delta,center,targ(ntarg)
      real *8 pot(nd,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,dsq
      complex *16 eye
      complex *16 qqx,qq1
      
      complex *16 ww1(100)

      complex *16 c3
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(itarg) - center)*dsq
         qqx = cdexp(eye*tx(npw2+1)*x)
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1
            qq1 = qq1*qqx
            ww1(npw-j1+1) = dconjg(ww1(j1))
         enddo
c
         do ind = 1,nd
            c3=0
            do j3=1,npw/2
               c3=c3+pwexp(j3,ind)*ww1(j3)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(c3)*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g1dpwevalg_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot,grad)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential only
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit real*8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center,targ(ntarg)
      real *8 pot(nd,ntarg),grad(nd,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3
      
      complex *16 ww1(100),ww1x(100)

      complex *16 cd,g1,g2,g3,d1,d1x,d1y,d2,d2x
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(itarg) - center)*dsq
         qqx = cdexp(eye*tx(npw2+1)*x)
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1
            ww1(npw-j1+1) = dconjg(ww1(j1))
            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww1x(npw-j1+1)=dconjg(ww1x(j1))

            qq1 = qq1*qqx
         enddo
c
         do ind = 1,nd
            d2=0
            d2x=0
            do j1=1,npw/2
               d2 = d2+pwexp(j1,ind)*ww1(j1)
               d2x = d2x+pwexp(j1,ind)*ww1x(j1)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(d2)*2
            grad(ind,itarg) = grad(ind,itarg)+dreal(d2x)*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g1dpwevalh_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient + hessian
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = Hessian (or vectorized Hessians) incremented
C
      implicit real*8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center,targ(ntarg)
      real *8 pot(nd,ntarg),grad(nd,ntarg),hess(nd,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,dsq
      complex *16 eye,ztmp
      complex *16 qqx,qq1
      
      complex *16 ww1(100),ww1x(100),ww1xx(100)

      complex *16 d2,d2x,d2xx

C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(itarg) - center)*dsq
         qqx = cdexp(eye*tx(npw2+1)*x)
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            ww1(j1) = qq1
            ww1(npw-j1+1) = dconjg(ww1(j1))
            ztmp = eye*tx(j1)*dsq
            ww1x(j1)= ztmp*ww1(j1)
            ww1x(npw-j1+1)=dconjg(ww1x(j1))
            ww1xx(j1)= ztmp*ww1x(j1)
            ww1xx(npw-j1+1)=dconjg(ww1xx(j1))
            
            qq1 = qq1*qqx
         enddo
c
         do ind = 1,nd
            d2=0
            d2x=0
            d2xx=0
            do j1=1,npw/2
               d2=d2+pwexp(j1,ind)*ww1(j1)
               d2x=d2x+pwexp(j1,ind)*ww1x(j1)
               d2xx=d2xx+pwexp(j1,ind)*ww1xx(j1)                  
            enddo

            pot(ind,itarg) = pot(ind,itarg)+dreal(d2)*2
            
            grad(ind,itarg) = grad(ind,itarg)+dreal(d2x)*2

            hess(ind,itarg) = hess(ind,itarg)+dreal(d2xx)*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C                  
c      
c
c*********************************************************************
C
C form PW expansions (charge, dipole, charge & dipole) using NUFFT
C
C*********************************************************************
      subroutine g1dformpwc_fast_vec(nd,delta,eps,sources,ns,charge,
     1            cent,npw,ws,ts,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     sources       = source locations
C     ns            = number of sources
C     charge        = strengths of sources
C     cent          = center of the expansion
C     npw           = number of terms in 1D PW expansion
C     ws,ts         = 1D pw expansion weights and nodes
C     wnufft        = real *8 (nexp) weights for nufft, tensor product of 1d weights
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent,sources(ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(npw/2)
      complex *16 ffexp(npw/2,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(cj(ns,nd))
      allocate(wj(ns))
      allocate(fk(npw,nd))
c
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(j)-cent)*dsq
         wj(j) = exp(-0.5d0*eye*xj(j))
      enddo

      do ind = 1,nd
         do j=1,ns
            cj(j,ind) = charge(ind,j)*wj(j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
      call finufft1d1many(nd,ns8,xj,cj,iflag,eps,npw8,
     1       fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         do j=1,npw/2
            ffexp(j,ind)=fk(j,ind)*wnufft(j)
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g1dformpwd_fast_vec(nd,delta,eps,sources,ns,rnormal,
     1    dipstr,cent,npw,ws,ts,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources       = source locations
C     ns            = number of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
C     cent          = center of the expansion
C     npw           = number of terms in 1D PW expansion
C     ws,ts         = 1D pw expansion weights and nodes
C     wnufft        = real *8 (nexp) weights for nufft, tensor product of 1d weights
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent,sources(ns),dipstr(nd,ns)
      real *8 rnormal(ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(npw/2)
      complex *16 ffexp(npw/2,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(cj(ns,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(j) - cent)*dsq
         wj(j) = exp(-0.5d0*eye*xj(j))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,ind) = ztmp*rnormal(j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = nd
      call finufft1d1many(ntrans,ns8,xj,cj,iflag,eps,
     1       npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k3=-npw2,-1
            j=j+1
            ffexp(j,ind)=wnufft(j)*fk(k3,ind)*ts(k3)
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g1dformpwcd_fast_vec(nd,delta,eps,sources,ns,charge,
     1    rnormal,dipstr,cent,npw,ws,ts,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     sources       = source locations
C     ns            = number of sources
C     charge        = strengths of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
C     cent          = center of the expansion
C     npw           = number of terms in 1D PW expansion
C     ws,ts         = 1D pw expansion weights and nodes
C     wnufft        = real *8 weights for nufft, tensor product of 1d weights
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent,sources(ns),dipstr(nd,ns)
      real *8 rnormal(ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(npw/2)
      complex *16 ffexp(npw/2,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(cj(ns,0:1,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,0:1,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(j) - cent)*dsq
         wj(j) = exp(-0.5d0*eye*xj(j))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            cj(j,0,ind) = charge(ind,j)*wj(j)
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,1,ind) = ztmp*rnormal(j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = 2*nd
      call finufft1d1many(ntrans,ns8,xj,cj,iflag,eps,
     1       npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k3=-npw2,-1
            j=j+1
            ffexp(j,ind)=wnufft(j)*(fk(k3,0,ind)+
     2                fk(k3,1,ind)*ts(k3))
         enddo
      enddo

      return
      end
c
C
C
C      
c*********************************************************************
C
C evaluate PW expansions (potential, pot + grad, pot + grad + hess)
C using NUFFT
C
C*********************************************************************
C
C
c
C
C
C
      subroutine g1dpwevalp_fast_vec(nd,delta,eps,center,npw,ws,ts,
     1              pwexp,targ,nt,pot)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential only
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target locations
C     nt            = number of targets
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center,targ(nt)
      real *8 pot(nd,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(npw/2,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      nexp = npw/2
      
      allocate(xj(nt))
      allocate(cj(nt,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,nd))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(j) - center)*dsq
         wj(j) = exp(0.5d0*eye*xj(j))
      enddo
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,ind)=pwexp(j,ind)
         enddo
         do j=nexp+1,2*nexp
            fk(j,ind)=0
         enddo
      enddo

      nt8=nt
      npw8=npw
      iflag = 1
      call finufft1d2many(nd,nt8,xj,cj,iflag,eps,npw8,
     1    fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g1dpwevalg_fast_vec(nd,delta,eps,center,npw,ws,ts,
     1              pwexp,targ,nt,pot,grad)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target locations
C     nt            = number of targets
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,nt
      real *8 delta,center,targ(nt)
      real *8 pot(nd,nt),grad(nd,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(npw/2,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
      complex *16, allocatable :: zts(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      nexp = npw/2
      
      allocate(xj(nt))
      allocate(cj(nt,0:1,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:1,nd))
      allocate(zts(-npw/2:npw/2-1))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(j) - center)*dsq
         wj(j) = exp(0.5d0*eye*xj(j))
      enddo

      z = eye/dsqrt(delta)
      do j=-npw/2,npw/2-1
         zts(j)=z*ts(j)
      enddo
      
      do ind = 1,nd
         j=0
         do j3=-npw/2,-1
            j=j+1
            fk(j,0,ind)=pwexp(j,ind)
            fk(j,1,ind)=pwexp(j,ind)*zts(j3)
         enddo
         
         do k=0,1
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=2*nd
      iflag = 1
      call finufft1d2many(ntrans,nt8,xj,cj,iflag,eps,
     1    npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            grad(ind,j)=grad(ind,j)+dreal(cj(j,1,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g1dpwevalh_fast_vec(nd,delta,eps,center,npw,ws,ts,
     1              pwexp,targ,nt,pot,grad,hess)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient + hess
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target locations
C     nt            = number of targets
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = hessian (or vectorized hessians) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center,targ(nt)
      real *8 pot(nd,nt),grad(nd,nt),hess(nd,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(npw/2,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
      complex *16, allocatable :: cp(:)
      
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)

      nexp = npw/2
      
      allocate(xj(nt))
      allocate(cj(nt,0:2,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:2,nd))
      allocate(cp(-npw/2:npw/2-1))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(j) - center)*dsq
         wj(j) = exp(0.5d0*eye*xj(j))
      enddo

      z = eye/dsqrt(delta)
      do j=-npw/2,npw/2-1
         cp(j) = z*ts(j)
      enddo
      
      do ind = 1,nd
         j=0
         do j3=-npw/2,-1
            j=j+1
            fk(j,0,ind)=pwexp(j,ind)
            fk(j,1,ind)=pwexp(j,ind)*cp(j3)
            fk(j,2,ind)=pwexp(j,ind)*cp(j3)*cp(j3)
         enddo
         
         do k=0,2
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=3*nd
      iflag = 1
      
      call finufft1d2many(ntrans,nt8,xj,cj,iflag,eps,
     1    npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            grad(ind,j)=grad(ind,j)+dreal(cj(j,1,ind)*wj(j))*2
            hess(ind,j)=hess(ind,j)+dreal(cj(j,2,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
c*********************************************************************
C
C get plane wave approximation nodes and weights
C
C*********************************************************************
      subroutine get_pwnodes(pmax,npw,ws,ts)
C
C     Get planewave exp weights,nodes
C
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      h = pmax/npw2
      w = h/(2.0d0*dsqrt(pi))

      do j =-npw2,npw2-1
         ts(j) = (j+0.5d0)*h
         ws(j) = w*dexp(-ts(j)*ts(j)/4)
      enddo
c
      return
      end
C
C
c
C
c*********************************************************************
C
C shift PW expansions (mp to loc at the cutoff level)
C
C*********************************************************************
      subroutine nufft_weights(npw,ws,ts,
     1    nexp,wnufft)
C
C     This subroutine precomputes all translation matrices for all SOE/X
C     expansions from child to parent or vice versa.
C
c     used in mp to mp or loc to loc stage
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     npw     = number of terms in plane wave exp
C     nmax    = number of different translation lengths in the whole scheme 
C     ts      = pw nodes
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      
      real *8 wnufft(nexp)
      
      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      j=0
      do j1=-npw/2,-1
         j=j+1
         wnufft(j) = ws(j1)
      enddo
c
      return
      end
c
c
c     
c
c*********************************************************************
C
C shift PW expansions (mp to loc, mp to mp, loc to loc)
C
C*********************************************************************
      subroutine pw_translation_matrices(xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw/2,-nmax:nmax)
      
      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      do j1=1,npw/2
         ztmp = exp(eye*ts(j1)*xmin)
         wshift(j1,0)=1
         do k1=1,nmax
            wshift(j1,k1) = ztmp
            ztmp = ztmp*ztmp
            wshift(j1,-k1) = dconjg(wshift(j1,k1))
         enddo
      enddo
c
      return
      end
c
c
c     
c
      subroutine merge_split_pw_matrices(xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ws,ts    = real *8, 1d PW expansion weights and nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw/2,2,nmax)
      real *8 ts(npw)
      complex *16 ztmp
      complex *16 eye
      complex *16, allocatable:: ww(:,:)
C
      eye =dcmplx(0,1)
      
      allocate(ww(npw,nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo
      
      do k1=1,nmax
         do j=1,npw/2
c           p              
            wshift(j,1,k1) = ww(j,k1)
c           m
            wshift(j,2,k1) = conjg(ww(j,k1))
         enddo
      enddo

      return
      end
c
c
c
c      
      subroutine g1dshiftpw_vec(nd,nexp,pwexp1,
     1              pwexp2,wshift)
C
C     This subroutine converts the PW expansion (pwexp1) about
C     the center (CENT1) into an PW expansion (pwexp2) about 
C     (CENT2) using precomputed translation matrix wshift.
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     delta   = Gaussian variance
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C     wshift  = precomputed PW exp translation matrix 
C
C     OUTPUT:
C
C     pwexp2 = shifted expansion 
C
      implicit none
      integer nd,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)
      complex *16 wshift(nexp)

C
      do ind=1,nd
         do j=1,nexp
            pwexp2(j,ind) = pwexp2(j,ind)
     1          +pwexp1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g1dcopypwexp_vec(nd,nexp,pwexp1,
     1              pwexp2)
C
C     This subroutine copy one PW expansion (pwexp1) 
C     to an PW expansion (pwexp2).
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C
C     OUTPUT:
C
C     pwexp2 = copied expansion 
C
      implicit none
      integer nd,nn,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)

C
      do ind=1,nd
         do j=1,nexp
            pwexp2(j,ind) = pwexp1(j,ind)
         enddo
      enddo
c
      return
      end
c
C
c
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
cC***********************************************************************
      subroutine g1dpwzero_vec(nd,pwexp,npw)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector planewave expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     npw    :   number of terms in 1d planewave expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     pwexp  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,npw,nd,ii
      complex *16 pwexp(npw/2,nd)
c
      do ii=1,nd
         do n=1,npw/2
            pwexp(n,ii)=0.0d0
         enddo
      enddo
      
      return
      end
c      
c      
c      
c      
c
