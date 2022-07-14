c-----------------------------------------------------------------------------
c     This file contains subroutines determining the proper number terms
c     for the expansions used in the FGT.
c
c
c
c     fgthlterms - determine number of terms in Hermite/local expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c     fgtpwterms - determine number of terms in the plane wave expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c     get_pwnodes : returns PW weights and nodes, midpoint rule is used
c                   so the number is always an even number!
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine fgthlterms(ndim,bsize,delta,eps,nterms)
c      
c     Determine the number of terms in the Taylor/Hermite expansions 
c     for boxes with side length = bsize
c
c     The error formula is from Wan and Karniadakis 2006 JCP 219 pp. 7-12.
c     Strain derived a similar formula in SIAM J. Sci. Stat. Comput. 12 (5) 
c     (1991) 1131–1139, essentially with only "ndim" missing.
c
c
c
c
c     INPUT:
c     ndim : dimension of the underlying space
c     bsize: side length of the box
c     delta: Gaussian variance
c     eps: desired precision
c    
c     OUTPUT:
c     nterms: number of terms in Taylor/Hermite expansion
c      
      implicit real*8 (a-h,o-z)
      integer nterms,p0,p
      real*8 delta, bsize
c
      e1 = exp(1.0d0)

      r = bsize/sqrt(2.0d0*delta)

      p0 = int(e1*r*r) + 1

      
      rk = 1.09*(8*atan(1.0d0))**(-0.25d0)

      
      do p=p0,300
         rp = r*sqrt(e1/p)
         error=ndim*rk*((p*1.0d0)**(-0.25d0))*(rp**p)/(1-rp)
         if (error .lt. eps) then
            nterms = p
            exit
         endif
      enddo

cccc      print *, p0, r*sqrt(2.0d0), nterms
      
      return
      end 
c
c
c
c
      subroutine fgtpwterms(bsize,delta,eps,iperiod,pmax,npw)
c      
c     Determine the number of terms in the plane wave expansions 
c     for boxes with side length = bsize
c
c
c     Method: 
c
c     Let G(x)=exp(-x^2), i.e., x is already scaled by sqrt(delta).
c     Then its Fourier transform is \hat{G}(k) = exp(-k^2/4).
c      
c     By the Poisson summation formula, we have
c      
c     \sum G(x+2\pi n/h) = \frac{h}{2sqrt{\pi} \sum \hat{G}(m h) exp(i m h x)
c
c     where both sums are from -\infty to +\infty.
c
c     (1) Truncation on the Fourier series on the RHS requires that
c
c            \hat{G}(pmax) < eps.
c
c     This leads to 
c
c            pmax = 2*sqrt(log(1/eps)).
c
c     (2) We also need to make sure the aliasing error on the LHS is small, 
c     that is,
c
c            G(2\pi/h - rmax)< eps.
c
c     This leads to
c
c           h = 2 \pi / (rmax + sqrt(log(1/eps))).
c
c     INPUT:
c     bsize: side length of the box
c     delta: Gaussian variance
c     eps: desired precision
c    
c     OUTPUT:
c     npw: number of terms in the plane wave expansion
c     pmax : maximum plane wave number
c
      implicit real*8 (a-h,o-z)
      real*8 delta, bsize
c
      pi = 4.0d0*atan(1.0d0)

      if (iperiod.eq.0) then
         rmax = bsize/sqrt(delta)

         d = sqrt(log(1.0d0/eps))
cccc      print *, 'in fgtpwterms, rmax/d=',rmax/d
         pmax = 2.0d0*d

         h = 2*pi/(rmax+d)

         npw = int(pmax/h)+1
cccc         npw = int(pmax/h)

         npw = 2*npw
cccc         print *, h, pmax, pmax/h, npw
      else
         npw = 2*(int(bsize*sqrt(log(1.0d0/eps)/delta)/pi)+1)+1
cccc         npw = 2*(int(bsize*sqrt(log(1.0d0/eps)/delta)/pi))+1

      endif
      
      return
      end 
c
c
c
c
c*********************************************************************
C
C get plane wave approximation nodes and weights
C
C*********************************************************************
      subroutine get_pwnodes(itype,pmax,npw,ws,ts)
C
C     Get planewave exp weights,nodes
C
c     input parameters
c     itype : integer (in) 1: trapezoidal rule; 0: midpoint rule
c      
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      h = pmax/npw2
      w = h/(2.0d0*dsqrt(pi))

      if (itype.eq.1) then
         do j =-npw2,npw2-1
            ts(j) = j*h
            ws(j) = w*dexp(-ts(j)*ts(j)/4)
         enddo
      elseif (itype.eq.0) then
         do j =-npw2,npw2-1
            ts(j) = (j+0.5d0)*h
            ws(j) = w*dexp(-ts(j)*ts(j)/4)
         enddo
      endif
c
      return
      end
C
C
c
C

      
c*********************************************************************
C
C get plane wave approximation nodes and weights for periodic Green's function
C
C*********************************************************************
      subroutine get_periodic_pwnodes(bsize0,delta,eps,npw,ws,ts)
C
C     Get planewave exp weights,nodes for periodic Green's function
C
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:(npw-1)/2),ts(-npw/2:(npw-1)/2)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      
      h = 2*pi/bsize0
      w = sqrt(pi*delta)/bsize0

      hw=h*sqrt(delta)/2
      do j =-npw2,(npw-1)/2
         ts(j) = j*h
         ws(j) = w*dexp(-(hw*j)**2)
      enddo
c
      return
      end
C
C
c
C

      
