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
      subroutine mk_leg2pw(n,npw,nnodes,ws,ts,delta,boxdim,tab_leg2pw)
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

ccc         qqx = cdexp(eye*ts(npw2+1)*x)
ccc         
ccc         qq1 = qqx
ccc         qqx = qqx*qqx
ccc
ccc         do j1=npw2+1,npw
ccc            tab_pw2pot(j1,i) = qq1
ccc            qq1 = qq1*qqx
ccc            
ccc            tab_pw2pot(npw-j1+1,i) = dconjg(tab_pw2pot(j1,i))
ccc         enddo
c
         do j=1,(npw+1)/2
            tab_pw2pot(j,i)=exp(eye*ts(j)*x)
            tab_pw2pot(npw-j+1,i)=dconjg(tab_pw2pot(j,i))
         enddo
c
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
c
C*********************************************************************C
      subroutine mk_pw2pg(norder,npw,ts,xs,delta,boxdim,tab_pw2pot,
     1           tab_pw2deriv)
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
c     tab_pw2deriv 
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 ts(npw),xs(norder)
      complex *16 tab_pw2pot(npw,norder),eye,qqx,qq1
      complex *16 tab_pw2deriv(npw,norder)
c
      eye = dcmplx(0.0d0,1.0d0)
c
      dsq = boxdim/2/dsqrt(delta)
C
      npw2=npw/2
      do i=1,norder
         x = xs(i)*dsq

ccc         qqx = cdexp(eye*ts(npw2+1)*x)
ccc         
ccc         qq1 = qqx
ccc         qqx = qqx*qqx
ccc
ccc         do j1=npw2+1,npw
ccc            tab_pw2pot(j1,i) = qq1
ccc            tab_pw2deriv(j1,i) = qq1*eye*ts(np
ccc            qq1 = qq1*qqx
ccc            
ccc            tab_pw2pot(npw-j1+1,i) = dconjg(tab_pw2pot(j1,i))
ccc         enddo

         do j=1,(npw+1)/2
            tab_pw2pot(j,i)=exp(eye*ts(j)*x)
            tab_pw2deriv(j,i)= eye*ts(j)*dsq*exp(eye*ts(j)*x)*2/boxdim
            tab_pw2pot(npw-j+1,i)=dconjg(tab_pw2pot(j,i))
            tab_pw2deriv(npw-j+1,i)=dconjg(tab_pw2deriv(j,i))
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
      subroutine mk_loctab_collpg(n,nnodes,delta,boxdim,
     1           tab_colleague,tab_colleaguex)
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
      real *8 tab_colleaguex(n,n,-1:1)
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
               rsum = rsum +  (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xim1 - xnodes(i)
               rsumm1 = rsumm1 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xip1 - xnodes(i)
               rsump1 = rsump1 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            tab_colleaguex(m,j,-1) = rsumm1
            tab_colleaguex(m,j,0)  = rsum
            tab_colleaguex(m,j,1)  = rsump1
         enddo
      enddo

      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine mk_loctab_stobpg(n,nnodes,delta,boxdim,
     1           tab_stob,tab_stobx)
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
      real *8 tab_stobx(n,n,4)
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
c
c
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
               rsum1 = rsum1 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xi2 - xnodes(i)
               rsum2 = rsum2 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)

               dx = xi3 - xnodes(i)
               rsum3 = rsum3 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
               
               dx = xi4 - xnodes(i)
               rsum4 = rsum4 + (-2.0d0*dx/delta)*
     1             legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            
            tab_stobx(m,j,1) = rsum1
            tab_stobx(m,j,2) = rsum2
            tab_stobx(m,j,3) = rsum3
            tab_stobx(m,j,4) = rsum4
         enddo
      enddo

      return
      end subroutine
c
c

c
c
C*********************************************************************C
      subroutine mk_loctab_btospg(n,nnodes,delta,boxdim,
     1           tab_btos,tab_btosx)
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
      real *8 tab_btosx(n,n,4)
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
               rsum1 = rsum1 + (-2.0d0*dx/delta)*
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi2 - xnodes(i)
               rsum2 = rsum2 + (-2.0d0*dx/delta)*
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi3 - xnodes(i)
               rsum3 = rsum3 + (-2.0d0*dx/delta)*
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
               dx = xi4 - xnodes(i)
               rsum4 = rsum4 + (-2.0d0*dx/delta)*
     1         legev(m,i)*exp(-dx*dx/delta)*whts(i)
            enddo
            tab_btosx(m,j,1) = rsum1
            tab_btosx(m,j,2) = rsum2
            tab_btosx(m,j,3) = rsum3
            tab_btosx(m,j,4) = rsum4
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
      subroutine mk_loctab_stob(n,nnodes,delta,boxdim,tab_stob)
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
      subroutine mk_loctab_btos(n,nnodes,delta,boxdim,tab_btos)
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
