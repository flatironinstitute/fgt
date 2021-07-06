c     This file contains a set of subroutines that build various 1D tables
c     used in the box FGT in all dimensions.
c
c
c     mk_leg2pw: builds the table converting Legendre expansion coefficients
c                to planewave expansion coefficients
c
c     mk_pw2pot: builds the table converting planewave expansion coefficients
c                to potential values on Legendre nodes
c
c     mk_loctab_coll: builds the table converting Legendre expansion coefficients
c                in the source box to potential values on Legendre nodes in the
c                target box at the same level, three of them for each level.
c
c     mk_loctab_stob: builds the table converting Legendre expansion coefficients
c                in a small source box to potential values on Legendre nodes in a
c                large target box at the coarse level. four of them for each level.
c
c     mk_loctab_btos: builds the table converting Legendre expansion coefficients
c                in a large source box to potential values on Legendre nodes in a
c                small target box at the fine level. four of them for each level.
c
C*********************************************************************C
      subroutine mk_leg2pw_old(n,npw,nnodes,ws,ts,delta,boxdim,
     1    tab_leg2pw)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_leg2pw(n,j) = ws(j)*(D/2) * 
c              int_{-1}^1 P_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c              where D is the box dimension at current level in
c              tree hierarchy.
c
c     That is, 
c     tab_leg2pw(n,j) is the Fourier transform of P_n at a specific
c     frequency. These can be expressed in terms of half-order Bessel
c     functions. A faster scheme would be to compute the top two 
c     coefficients and then use a downward recurrence to get the rest.
c     There are also analytic formulae for these Fourier transforms.
c
c     INPUT:
c     n        dimension of coeff array
c     npw      number of plane waves
c     nnodes   number of nodes used in numerical quadrature
c     ws,ts    weights and nodes of plane wave quadrature
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_leg2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 P_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 u,v
      real *8, allocatable :: legev(:,:)
      complex *16 tab_leg2pw(n,npw),zsum,eye
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: whts(:), xnodes(:)
c
      allocate(legev(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
c
      eye = dcmplx(0.0d0,1.0d0)
      itype = 1
      call legeexps(itype,nnodes,xnodes,u,v,whts)
c
      do i = 1,nnodes
         call legepols(xnodes(i),n-1,legev(1,i))
      enddo
ccc      call prin2(' legev is *',legev,n*nnodes)
c
      fac = boxdim/(2*dsqrt(delta))
      do m = 1,n
         do j = 1,npw
            zsum = 0.0d0
            do i = 1,nnodes
               zsum = zsum +
     1         legev(m,i)*cdexp(-eye*ts(j)*fac*xnodes(i))*whts(i)
            enddo
            tab_leg2pw(m,j) = ws(j)*zsum*boxdim/2.0d0
         enddo
      enddo
      return
      end subroutine
c
c
c
c
c
c
C*********************************************************************C
      subroutine mk_leg2pw(n,npw,nnodes,ws,ts,delta,boxdim,tab_leg2pw)
C*********************************************************************C
c     Use half-order Bessel J function to construct the table.
c
c     tab_leg2pw(n,j) = ws(j)*(D/2) * 
c            int_{-1}^1 P_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c      
c     it is known that
c      
c     int_{-1}^1 P_n(x) exp(- i a x) dx = 1/i^m *sqrt(2pi/a) J_{n+1/2)(a),
c      
c     where i^2=-1, J_{n+1/2) is the Bessel J function of half integer order.
c
c      
c     we use the package TOMS 644 to evaluate J_{n+1/2}. James Bremer also has 
c     a package for the evaluation of J_{nu}, but it requires reading  
c     tables, evaluates J_{n+1/2} and Y_{n+1/2} at the same time,
c     and does the evaluate for each n one at a time. TOMS 644 evaluates
c     J functions only, and does the evaluation for a sequence of n at the
c     same time.
c      
c     Here D is the box dimension at current level in the tree hierarchy.
c
c     Since a could be large, reaching about 100 for high accuracy,
c     the quadrature scheme will need to set nnodes very large,
c     leading to an inefficient scheme.
c      
c     That is, 
c     tab_leg2pw(n,j) is the Fourier transform of P_n at a specific
c     frequency. These can be expressed in terms of half-order Bessel
c     functions. A faster scheme would be to compute the top two 
c     coefficients and then use a downward recurrence to get the rest.
c     There are also analytic formulae for these Fourier transforms.
c
c     INPUT:
c     n        dimension of coeff array
c     npw      number of plane waves
c     nnodes   number of nodes used in numerical quadrature
c     ws,ts    weights and nodes of plane wave quadrature
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_leg2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 P_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 tab_leg2pw(n,npw),eye,eyem,z
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: cyr(:),cyi(:)
c
      allocate(cyr(n))
      allocate(cyi(n))
c
      pi = 4*atan(1.0d0)
      
      eye = dcmplx(0.0d0,1.0d0)
      fac = boxdim/(2*dsqrt(delta))

      fnu=0.5d0
      kode=1
      do j = npw/2+1,npw
         dd = ws(j)*boxdim/2.0d0
         
         if (abs(ts(j)).lt.1d-12) then
            tab_leg2pw(1,j)=dd*2
            do m=2,n
               tab_leg2pw(m,j)=0
            enddo
         else
            zi=0.0d0
            zr = ts(j)*fac
            ddd = sqrt(2*pi/zr)
            z = cmplx(zr,0.0d0)
            call zbesj(zr,zi,fnu,kode,n,cyr,cyi,nz,ierr)
            eyem=1.0d0
            do m=1,n
               tab_leg2pw(m,j)=dd*dcmplx(cyr(m),cyi(m))*ddd/eyem
               tab_leg2pw(m,npw-j+1) = conjg(tab_leg2pw(m,j))
               eyem=eyem*eye
            enddo
         endif
      enddo

      return
      end subroutine
c
c
c
c
c
c
C*********************************************************************C
      subroutine mk_pw2pot(norder,npw,ts,xs,delta,boxdim,tab_pw2pot)
C*********************************************************************C
c     generates a table converting plane wave expansion to potential 
c     values
c      
c     tab_pw2pot(j,n) = exp(i ts(j) xs(n) D/(2 \sqrt{delta}))
c              where D is the box dimension at current level in
c              tree hierarchy.
c
c     INPUT:
c     norder   number of Legendre nodes
c     npw      number of plane waves
c     ts       nodes of plane wave quadrature
c     xs       Legendre nodes
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_pw2pot 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 ts(npw),xs(norder)
      complex *16 tab_pw2pot(npw,norder),eye,qqx,qq1
c
      eye = dcmplx(0.0d0,1.0d0)
c
      dsq = boxdim/2/dsqrt(delta)
C
      npw2=npw/2
      do i=1,norder
         x = xs(i)*dsq

         qqx = cdexp(eye*ts(npw2+1)*x)
         
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            tab_pw2pot(j1,i) = qq1
            qq1 = qq1*qqx
            
            tab_pw2pot(npw-j1+1,i) = dconjg(tab_pw2pot(j1,i))
         enddo

cccc         do j=1,npw
cccc            tab_pw2pot(j,i)=exp(eye*ts(j)*x)
cccc         enddo
      enddo

      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_coll_old(n,nnodes,delta,boxdim,tab_colleague)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{-D/2}^{D/2} P_n(x) exp( (\xi_j -x)^2/delta)
c              where D is the box dimension at current level in
c              tree hierarchy and \xi_j is either on 
c              [-D/2,D/2]   -> tab_colleague(n,n,0)
c              [-3D/2,-D/2] -> tab_colleague(n,n,-1)
c              [D/2,3D/2]   -> tab_colleague(n,n,1)
c              
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_colleague 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 u,v, rsum,xi,xip1,xim1
      real *8, allocatable :: legev(:,:)
      real *8 tab_colleague(n,n,-1:1)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: xnodest(:)
c
      allocate(legev(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
      allocate(xnodest(n))
c
      itype = 1
      call legeexps(itype,n,xnodest,u,v,whts)
      do i=1,n
         xnodest(i) = boxdim*xnodest(i)/2.0d0
      enddo
ccc      call prin2(' Legendre nodes are *',xnodest,n)
      call legeexps(itype,nnodes,xnodes,u,v,whts)
ccc      call prin2(' Legendre nodes are *',xnodes,nnodes)
c
      do i = 1,nnodes
         call legepols(xnodes(i),n-1,legev(1,i))
      enddo
      do i=1,nnodes
         xnodes(i) = boxdim*xnodes(i)/2.0d0
         whts(i) = whts(i)*boxdim/2.0d0
      enddo
ccc      call prin2(' legev is *',legev,n*nnodes)
c
      do m = 1,n
         do j = 1,n
            xi = xnodest(j)
            xim1 = xnodest(j) - boxdim
            xip1 = xnodest(j) + boxdim
            rsum = 0.0d0
            rsumm1 = 0.0d0
            rsump1 = 0.0d0
            do i = 1,nnodes
               dx = xi - xnodes(i)
               rsum = rsum +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xim1 - xnodes(i)
               rsumm1 = rsumm1 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xip1 - xnodes(i)
               rsump1 = rsump1 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            tab_colleague(m,j,-1) = rsumm1
            tab_colleague(m,j,0)  = rsum
            tab_colleague(m,j,1)  = rsump1
         enddo
      enddo
      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine mk_loctab_stob_old(n,nnodes,delta,boxdim,tab_stob)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{source box} P_n(x) exp( (\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-5D/4,-D/4]     -> tab_stob(n,n,1)
c              [-3D/4, D/4]     -> tab_stob(n,n,2)
c              [ -D/4,3D/4]     -> tab_stob(n,n,3)
c              [  D/4,5D/4]     -> tab_stob(n,n,4)
c              
c     Here we assume that the source box is centered at the origin.        
c              
c              
c      _____ _____ ____  
c     |     |     |    | 
c     |     |     |    | 
c     |_____|_____|____| 
c     |     |  |A |    | For target points in large box B, of
c     |     |--|--| B  | dimension D, adjacent small source box centers 
c     |_____|__|__|____| can be offset by one of -3D/4,-D/4,D/4,3D/4
c     |     |     |    | in either x, y, or z.
c     |     |     |    |  
c     |_____|_____|____|   
c                          
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   target box dimension at current level
c
c     OUTPUT:
c     tab_stob 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 u,v, rsum,xi,xi1,xi2,xi3,xi4
      real *8, allocatable :: legev(:,:)
      real *8 tab_stob(n,n,4)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: xnodest(:)
c
      allocate(legev(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
      allocate(xnodest(n))
c
c     target grid points
c
      itype = 1
      call legeexps(itype,n,xnodest,u,v,whts)
      do i=1,n
         xnodest(i) = boxdim*xnodest(i)/2.0d0
      enddo
ccc      call prin2(' Legendre nodes are *',xnodest,n)
      call legeexps(itype,nnodes,xnodes,u,v,whts)
ccc      call prin2(' Legendre nodes are *',xnodes,nnodes)
c
c     fine grid points for numerical integration on source
c
      do i = 1,nnodes
         call legepols(xnodes(i),n-1,legev(1,i))
      enddo
      do i=1,nnodes
         xnodes(i) = boxdim*xnodes(i)/4.0d0
         whts(i) = whts(i)*boxdim/4.0d0
      enddo
ccc      call prin2(' legev is *',legev,n*nnodes)
c
      do m = 1,n
         do j = 1,n
            xi = xnodest(j)
            xi1 = xnodest(j) - 3.0d0*boxdim/4.0d0
            xi2 = xnodest(j) - boxdim/4.0d0
            xi3 = xnodest(j) + boxdim/4.0d0
            xi4 = xnodest(j) + 3.0d0*boxdim/4.0d0
            
            rsum1 = 0.0d0
            rsum2 = 0.0d0
            rsum3 = 0.0d0
            rsum4 = 0.0d0

            do i = 1,nnodes
               dx = xi1 - xnodes(i)
               rsum1 = rsum1 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xi2 - xnodes(i)
               rsum2 = rsum2 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)

               dx = xi3 - xnodes(i)
               rsum3 = rsum3 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xi4 - xnodes(i)
               rsum4 = rsum4 +
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            
            tab_stob(m,j,1) = rsum1
            tab_stob(m,j,2) = rsum2
            tab_stob(m,j,3) = rsum3
            tab_stob(m,j,4) = rsum4
         enddo
      enddo
      return
      end subroutine
c
c

c
c
C*********************************************************************C
      subroutine mk_loctab_btos_old(n,nnodes,delta,boxdim,tab_btos)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{source box} P_n(x) exp( (\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-2D,-D]    -> tab_btos(n,n,1)
c              [ -D, 0]    -> tab_btos(n,n,2)
c              [  0, D]    -> tab_btos(n,n,3)
c              [  D,2D]    -> tab_btos(n,n,4)
c              
c     Here we assume that the source box is centered at the origin.        
c              
c              
c      _____ _____ ____  
c     |     |     |    | 
c     |  A  |     |    | 
c     |_____|_____|____| 
c     |     |B |  |    | For target points in small target box B, of
c     |     |--|--|    | dimension D, adjacent large source box A centers 
c     |_____|__|__|____| can be offset by one of -3D/2,-D/2,D/2,3D/2
c     |     |     |    | in either x, y, or z.
c     |     |     |    |  
c     |_____|_____|____|   
c                          
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   target box dimension at current level
c
c     OUTPUT:
c     tab_stob 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 u,v, rsum,xi,xi1,xi2,xi3,xi4
      real *8, allocatable :: legev(:,:)
      real *8 tab_btos(n,n,4)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: xnodest(:)
c
      allocate(legev(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
      allocate(xnodest(n))
c
c     target grid points
c
      itype = 1
      call legeexps(itype,n,xnodest,u,v,whts)
      do i=1,n
         xnodest(i) = boxdim*xnodest(i)/2.0d0
      enddo
ccc      call prin2(' Legendre nodes are *',xnodest,n)
      call legeexps(itype,nnodes,xnodes,u,v,whts)
ccc      call prin2(' Legendre nodes are *',xnodes,nnodes)
c
c     fine grid points for numerical integration on source
c
      do i = 1,nnodes
         call legepols(xnodes(i),n-1,legev(1,i))
      enddo
      do i=1,nnodes
         xnodes(i) = boxdim*xnodes(i)
         whts(i) = whts(i)*boxdim
      enddo
ccc      call prin2(' legev is *',legev,n*nnodes)
c
      do m = 1,n
         do j = 1,n
            xi = xnodest(j)
            xi1 = xnodest(j) - 3.0d0*boxdim/2.0d0
            xi2 = xnodest(j) - boxdim/2.0d0
            xi3 = xnodest(j) + boxdim/2.0d0
            xi4 = xnodest(j) + 3.0d0*boxdim/2.0d0
            rsum1 = 0.0d0
            rsum2 = 0.0d0
            rsum3 = 0.0d0
            rsum4 = 0.0d0
            do i = 1,nnodes
               dx = xi1 - xnodes(i)
               rsum1 = rsum1 +
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi2 - xnodes(i)
               rsum2 = rsum2 +
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi3 - xnodes(i)
               rsum3 = rsum3 +
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi4 - xnodes(i)
               rsum4 = rsum4 +
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            tab_btos(m,j,1) = rsum1
            tab_btos(m,j,2) = rsum2
            tab_btos(m,j,3) = rsum3
            tab_btos(m,j,4) = rsum4
         enddo
      enddo
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_coll(n,nnodes,delta,boxdim,tab_colleague)
C*********************************************************************C
c     This routine is an optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{-D/2}^{D/2} P_n(x*2/D) exp( -(\xi_j -x)^2/delta)
c              where D is the box dimension at current level in
c              tree hierarchy and \xi_j is either on 
c              [-D/2,D/2]   -> tab_colleague(n,n,0)
c              [-3D/2,-D/2] -> tab_colleague(n,n,-1)
c              [D/2,3D/2]   -> tab_colleague(n,n,1)
c              
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_colleague 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 xi,xip1,xim1
      real *8 tab_colleague(n,n,-1:1)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: legev(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
c
      itype = 0
      call legeexps(itype,n,xnodest,u,v,whts)

      sigma = 4*delta/boxdim**2
      lambda = sqrt(sigma)
      eps=2d-16
      dmax=sqrt(log(1/eps)*sigma)
cccc      print *, 'dmax=',dmax
c      
      if (lambda.le.0.125d0) then
         print *, 'enter small lambda zone in coll table'
         
c
c     colleague table -1, the scaled target interval is on [-3,-1], the scaled
c     source interval is on [-1,1]. 
         do j=1,n
            xim1 = xnodest(j) - 2
            dx = 1-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_colleague(m,j,-1)  = 0
               enddo
            else
               call mk_loctab_recurrence(n,xim1,lambda,dmax,boxdim,fint)
               do m=1,n
                  tab_colleague(m,j,-1)  = fint(m)
               enddo
            endif
         enddo

c     colleague table 0, the scaled target interval is on [-1,1], the scaled source
c     interval is on [-1, 1].
c      
c     by symmetry, only needs to construct the table for half of the target points.
c      
         do j=1,n/2
            xi = xnodest(j)
            call mk_loctab_recurrence(n,xi,lambda,dmax,boxdim,fint)
            do m=1,n
               tab_colleague(m,j,0)  = fint(m)
            enddo
         enddo         
      else
         nquad = 50
         allocate(legev(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c        quadrature nodes and weights on the source interval
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/2
         enddo
      
         do i = 1,nquad
            call legepols(xnodes(i),n-1,legev(1,i))
         enddo
c         
c        colleague table -1, the scaled target interval is on [-3,-1].
         do j=1,n
            xim1 = xnodest(j) - 2.0d0
            do i=1,nquad
               dx = xnodes(i) -xim1
               if (dx.lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = 1,n
            dx = 1-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_colleague(m,j,-1) = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + legev(m,i)*wexp(i,j)
                  enddo
                  tab_colleague(m,j,-1) = rsum1
               enddo
            endif
         enddo
c        colleague table 0, the scaled target interval is on [-1,1], 
c        the scaled source interval is on [-1, 1].
c        By symmetry, only needs to construct the table for half of the target points.
         do j=1,n/2
            xi = xnodest(j)
            do i=1,nquad
               dx = xnodes(i) - xi
               if (abs(dx).lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = 1,n/2
            do m = 1,n
               rsum1 = 0.0d0
               do i = 1,nquad
                  rsum1 = rsum1 + legev(m,i)*wexp(i,j)
               enddo
               tab_colleague(m,j,0)  = rsum1
            enddo
         enddo
      endif
      
c     obtain the other half of table 0 by symmetry
      do j=1,n/2
         do m=1,n,2
            tab_colleague(m,n-j+1,0)  = tab_colleague(m,j,0) 
         enddo
         do m=2,n,2
            tab_colleague(m,n-j+1,0)  = -tab_colleague(m,j,0) 
         enddo
      enddo
      
c     use symmetry to construct the table +1 for targets on the right side
      do j=1,n
         do m=1,n,2
            tab_colleague(m,j,1) = tab_colleague(m,n-j+1,-1)
         enddo
         do m=2,n,2
            tab_colleague(m,j,1) = -tab_colleague(m,n-j+1,-1)
         enddo
      enddo
      
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_stob(n,nnodes,delta,boxdim,tab_stob)
C*********************************************************************C
c     This routine is an optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-5D/4,-D/4]     -> tab_stob(n,n,1)
c              [-3D/4, D/4]     -> tab_stob(n,n,2)
c              [ -D/4,3D/4]     -> tab_stob(n,n,3)
c              [  D/4,5D/4]     -> tab_stob(n,n,4)
c              
c     Here we assume that the source box is centered at the origin.        
c              
c              
c     Small source interval to big target interval
c      
c              
c      _____ _____ ____  
c     |     |     |    | 
c     |     |     |    | 
c     |_____|_____|____| 
c     |     |  |A |    | For target points in large box B, of
c     |     |--|--| B  | dimension D, adjacent small source box centers 
c     |_____|__|__|____| can be offset by one of -3D/4,-D/4,D/4,3D/4
c     |     |     |    | in either x, y, or z.
c     |     |     |    |  
c     |_____|_____|____|   
c                          
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   target box dimension at current level
c
c     OUTPUT:
c     tab_stob 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 xi,xi1,xi2,xi3,xi4
      real *8 tab_stob(n,n,4)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: legev(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
c
      itype = 0
      call legeexps(itype,n,xnodest,u,v,whts)
c     big target interval, enlarged by a factor of 2
      do i=1,n
         xnodest(i) = xnodest(i)*2
      enddo
      
      sigma = 16*delta/boxdim**2
      lambda = sqrt(sigma)
      eps=2d-16
      dmax=sqrt(log(1/eps)*sigma)
cccc      print *, 'dmax=',dmax
c
      if (lambda.le.0.125d0) then
         print *, 'enter small lambda zone in stob table'
         
c        stob table 1, the scaled target interval is on [-5,-1], the scaled
c        source interval is on [-1,1]. 
         do j=1,n
            xi1 = xnodest(j) - 3.0d0
            dx = 2-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_stob(m,j,1) = 0
               enddo
            else
               call mk_loctab_recurrence(n,xi1,lambda,dmax,
     1             boxdim/2,fint)
               do m=1,n
                  tab_stob(m,j,1)  = fint(m)
               enddo
            endif
         enddo

c        stob table 2, the scaled target interval is on [-3,1], the scaled source
c        interval is on [-1, 1].
c      
c        split the target into two parts. 1. for targets on [-3,-1], check whether 
c        the target is far away from the source interval
         do j=1,n/2
            xi2 = xnodest(j) - 1.0d0
            dx = -xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_stob(m,j,2) = 0
               enddo
            else
               call mk_loctab_recurrence(n,xi2,lambda,dmax,
     1             boxdim/2,fint)
               do m=1,n
                  tab_stob(m,j,2)  = fint(m)
               enddo
            endif
         enddo
c        for targets on [-1,1], no need to check
         do j=n/2+1,n
            xi2 = xnodest(j) - 1.0d0       
            call mk_loctab_recurrence(n,xi2,lambda,dmax,boxdim/2,fint)
            do m=1,n
               tab_stob(m,j,2)  = fint(m)
            enddo
         enddo
      else
         nquad = 50
         allocate(legev(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c        obtain quadrature nodes and weights on the source interval
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/4
         enddo
      
         do i = 1,nquad
            call legepols(xnodes(i),n-1,legev(1,i))
         enddo
c        stob table 1, the scaled target interval is on [-5,-1]
         do j=1,n
            xi1 = xnodest(j) - 3.0d0
            do i=1,nquad
               dx = xnodes(i) - xi1
               if (dx.lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = 1,n
            dx = 2-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_stob(m,j,1) = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + legev(m,i)*wexp(i,j)
                  enddo
                  tab_stob(m,j,1) = rsum1
               enddo
            endif
         enddo
c        stob table 2, the scaled target interval is on [-3,1]
c        for targets on [-3,-1], check whether the target is far away from
c        the source interval
         do j=1,n/2
            xi2 = xnodest(j) - 1.0d0
            do i=1,nquad
               dx = xnodes(i) - xi2
               if (dx.lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = 1,n/2
            dx=-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_stob(m,j,2) = 0
               enddo
            else
               do m = 1,n
                  rsum2 = 0.0d0
                  do i = 1,nquad
                     rsum2 = rsum2 + legev(m,i)*wexp(i,j)
                  enddo
                  tab_stob(m,j,2) = rsum2
               enddo
            endif
         enddo
c        for targets on [-1,1], no need to check
         do j=n/2+1,n
            xi2 = xnodest(j) - 1.0d0
            do i=1,nquad
               dx = xnodes(i) - xi2
               if (abs(dx).lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = n/2+1,n
            do m = 1,n
               rsum1 = 0.0d0
               do i = 1,nquad
                  rsum1 = rsum1 + legev(m,i)*wexp(i,j)
               enddo
               tab_stob(m,j,2) = rsum1
            enddo
         enddo
      endif

c     use symmetry to construct the tables 3 and 4
      do j=1,n
         do m=1,n,2
            tab_stob(m,j,3) = tab_stob(m,n-j+1,2)
         enddo
         do m=2,n,2
            tab_stob(m,j,3) = -tab_stob(m,n-j+1,2)
         enddo
      enddo
      
      do j=1,n
         do m=1,n,2
            tab_stob(m,j,4) = tab_stob(m,n-j+1,1)
         enddo
         do m=2,n,2
            tab_stob(m,j,4) = -tab_stob(m,n-j+1,1)
         enddo
      enddo

      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_btos(n,nnodes,delta,boxdim,tab_btos)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-2D,-D]    -> tab_btos(n,n,1)
c              [ -D, 0]    -> tab_btos(n,n,2)
c              [  0, D]    -> tab_btos(n,n,3)
c              [  D,2D]    -> tab_btos(n,n,4)
c              
c     Here we assume that the source box is centered at the origin.        
c              
c     Big source interval to small target interval
c      
c      _____ _____ ____  
c     |     |     |    | 
c     |  A  |     |    | 
c     |_____|_____|____| 
c     |     |B |  |    | For target points in small target box B, of
c     |     |--|--|    | dimension D, adjacent large source box A centers 
c     |_____|__|__|____| can be offset by one of -3D/2,-D/2,D/2,3D/2
c     |     |     |    | in either x, y, or z.
c     |     |     |    |  
c     |_____|_____|____|   
c                          
c
c     INPUT:
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   target box dimension at current level
c
c     OUTPUT:
c     tab_btos
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 xi,xi1,xi2,xi3,xi4
      real *8 tab_btos(n,n,4)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: legev(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
c
      itype = 0
      call legeexps(itype,n,xnodest,u,v,ws)
c     small target interval, reduced by a factor of 2
      do i=1,n
         xnodest(i) = xnodest(i)/2
      enddo

      sigma = delta/boxdim**2
      lambda = sqrt(sigma)
cccc      print *, 'lambda=', lambda
      
      eps=2d-16
      dmax=sqrt(log(1/eps)*sigma)
cccc      print *, 'dmax=',dmax
c
      if (lambda.le.0.125d0) then
         print *, 'enter small lambda zone in btos table'
c        btos table 1, the scaled target interval is on [-2,-1], the scaled
c        source interval is on [-1,1]. 
         do j=1,n
            dx=0.5d0-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_btos(m,j,1) = 0
               enddo
            else
               xi1 = xnodest(j) - 1.5d0
               call mk_loctab_recurrence(n,xi1,lambda,dmax,
     1             boxdim*2,fint)
               do m=1,n
                  tab_btos(m,j,1)  = fint(m)
               enddo
            endif
         enddo
c        btos table 2, the scaled target interval is on [-1,0], the scaled source
c        interval is on [-1, 1]. 
         do j=1,n
            xi2 = xnodest(j) - 0.5d0         
            call mk_loctab_recurrence(n,xi2,lambda,dmax,boxdim*2,fint)
            do m=1,n
               tab_btos(m,j,2)  = fint(m)
            enddo
         enddo
      else
         nquad = 50
         allocate(legev(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c        quadrature nodes and weights on the source interval
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)
         do i=1,nquad
            whts(i) = whts(i)*boxdim
         enddo
         
         do i = 1,nquad
            call legepols(xnodes(i),n-1,legev(1,i))
         enddo
c        btos table 1, the scaled target interval is on [-2,-1].         
         do j=1,n
            xi1 = xnodest(j) - 1.5d0
            do i=1,nquad
               dx = xnodes(i) - xi1
               if (dx.lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo
         do j = 1,n
            dx=0.5d0-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tab_btos(m,j,1) = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + legev(m,i)*wexp(i,j)
                  enddo
                  tab_btos(m,j,1) = rsum1
               enddo
            endif
         enddo
c        btos table 2, the scaled target interval is on [-1,0]         
         do j=1,n
            xi2 = xnodest(j) - 0.5d0
            do i=1,nquad
               dx = xnodes(i) - xi2
               if (abs(dx).lt.dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*whts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo

         do j = 1,n
            do m = 1,n
               rsum1 = 0.0d0
               do i = 1,nquad
                  rsum1 = rsum1 + legev(m,i)*wexp(i,j)
               enddo
               tab_btos(m,j,2) = rsum1
            enddo
         enddo
      endif

c     use symmetry to construct tables 3 and 4
      do j=1,n
         do m=1,n,2
            tab_btos(m,j,3) = tab_btos(m,n-j+1,2)
         enddo
         do m=2,n,2
            tab_btos(m,j,3) = -tab_btos(m,n-j+1,2)
         enddo
      enddo
      
      do j=1,n
         do m=1,n,2
            tab_btos(m,j,4) = tab_btos(m,n-j+1,1)
         enddo
         do m=2,n,2
            tab_btos(m,j,4) = -tab_btos(m,n-j+1,1)
         enddo
      enddo

      return
      end subroutine
c
c
c
c
      subroutine mk_loctab_recurrence(m,targ,lambda,dmax,boxdim,fint)
      implicit real *8 (a-h,o-z)
cccc      sqrtpih=sqrt(4.0d0*atan(1.0d0))/2
      data sqrtpih/0.886226925452757940959713778283913d0/
c     calculate the values of the integral
c
c     \frac{1}{\lambda}\int_{-1}^1 P_n(x) e^{-(targ-x)^2/\lambda^2}dx
c
c     for n=0,1,...,m-1.
c      
c     Algorithm: use the five term recurrence formula
c     
c     Assumption: lambda<=1/8
c
      real *8 fint(m),lambda,lambda2
      
      lambda2=lambda*lambda
      
      a = (-1-targ)/lambda
      b = (1-targ)/lambda

      erfa = erf(a)
      erfb = erf(b)

      expa = exp(-a*a)
      expb = exp(-b*b)

      d1 = lambda*(expb - expa)
      d2 = lambda*(expb + expa)
      
      fint(1) = sqrtpih*(erfb - erfa)
      fint(2) = -d1/2 + targ*fint(1)
      fint(3) = -d2 - targ*d1
     1    +(lambda2-2.0d0/3+2*targ**2)*fint(1)
      fint(3) = 0.75d0*fint(3)
      
      fint(4) = -d1 + 2*targ*fint(3)
     1    -(0.4d0-4*lambda2)*fint(2)
      fint(4) = 5*fint(4)/8
      
      do k=5,m
         n=k-3
         fint(k)=targ*(fint(k-1)-fint(k-3)) + (n-1)*fint(k-4)/(2*n-1) 
     1       + (2*n+1)*(lambda2/2 + 1.0d0/(2*n+3)/(2*n-1))*fint(k-2) 
         fint(k) = (2*n+3)*fint(k)/(n+2)
      enddo

c     finally, multiply lambda back and also the proper weight adjustment
c     for the source box
      ww = boxdim*lambda/2
      do k=1,m
         fint(k)=fint(k)*ww
      enddo

      return
      end
c
c
c
      
