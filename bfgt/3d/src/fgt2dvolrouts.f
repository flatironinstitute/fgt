C*********************************************************************C
      subroutine legtrans2d(nd,n,fdat,f,texp,umat,coeff)
C*********************************************************************C
c     2d Legendre transformation of data on tensor product
c     Legendre grid
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        grid size (and expansion order)
c     fdat     tensor product function values
c              fdat(i,j,ind) = f(x_i,y_j) for ind sample of f.
c     f, texp  work arrays of length n
c     umat     1D transform matrix from legeexps
c
c     OUTPUT:
c     coeff    Legendre coefficients
c              f = sum coeff(n,m) P_n(x) P_m(y)
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fdat(n,n,nd),coeff(n,n,nd),umat(n,n)
      real *8 f(n),texp(n)
c
c     transform rows
c
      do ind = 1,nd
         do j = 1,n
            do i = 1,n
               f(i) = fdat(i,j,ind)
               texp(i) = 0.0d0
            enddo
            do k = 1,n
            do i = 1,n
               texp(i) = texp(i) + umat(i,k)*f(k)
            enddo
            enddo
            do i = 1,n
               coeff(i,j,ind) = texp(i)
            enddo
         enddo
c
c     transform columns
c
         do i = 1,n
            do j = 1,n
               f(j) = coeff(i,j,ind)
               texp(j) = 0.0d0
            enddo
            do j = 1,n
            do k = 1,n
               texp(k) = texp(k) + umat(k,j)*f(j)
            enddo
            enddo
            do j = 1,n
               coeff(i,j,ind) = texp(j)
            enddo
         enddo
      enddo
      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine leg2deval(n,coeff,xx,yy,val,w1,w2)
C*********************************************************************C
c     slow 2d Legendre series evaluation at single point.
c
c     INPUT:
c     n        grid size (and expansion order)
c     coeff    Legendre coefficients
c              f = sum coeff(n,m) P_n(x) P_m(y)
c     xx,yy    target coordinates
c     w1, w2   work arrays of length n
c
c     OUTPUT:
c     val      series value at (xx,yy)
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fdat(n,n),coeff(n,n)
      real *8 f,w1(n),w2(n)
c
      do j = 1,n
         do i = 1,n
            w1(i) = coeff(i,j)
         enddo
         call legeexev(xx,val,w1,n)
         w2(j) = val
      enddo
      call legeexev(yy,val,w2,n)
c
      return
      end subroutine
c
c
c
c
c
C*********************************************************************C
      subroutine leg2d_to_pw(nd,n,coeff,npw,ff,tab_leg2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     Legendre series coefficients (implicitly about the box center).
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       Legendre coefficients
c                 f = sum coeff(n,m) P_n(x) P_m(y)
c     npw         number of plane waves.
c                 NOTE 2D convention is pwexp(npw/2,npw)
c     ff          workspace
c     tab_leg2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,nd)
      complex *16 ff(n,npw/2),tab_leg2pw(n,npw)
      complex *16 pwexp(npw/2,npw,nd),zsum
c
      do ind = 1,nd
         do m = 1,n
            do j = 1,npw/2
               ff(m,j) = 0.0d0
               do i = 1,n
                  ff(m,j)=ff(m,j)+tab_leg2pw(i,j)*coeff(i,m,ind)
               enddo
            enddo
         enddo
c
         do k = 1,npw
         do j = 1,npw/2
               pwexp(j,k,ind) = 0.0d0
               do m = 1,n
                  pwexp(j,k,ind)=pwexp(j,k,ind)
     1                 +ff(m,j)*tab_leg2pw(m,k)
               enddo
         enddo
         enddo
c
      enddo
      return
      end subroutine
c

c
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
         call legepols(xnodes(i),n,legev(1,i))
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
      subroutine leg2d_to_potloc(nd,n,coeff,ff,pot,tab_colx,tab_coly)
C*********************************************************************C
c     This routine computes the volume Gauss transform over a 
c     single box source distribution given as a Legendre series.
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays (tab_colx, tab_coly). 
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in 2D, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x or y in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x or y.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x or y.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       Legendre coefficients
c                 f = sum coeff(n,m) P_n(x) P_m(y)
c     ff          workspace
c     tab_colx    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in x.
c     tab_coly    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in y.
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,nd),pot(n,n,nd)
      real *8 ff(n,n),tab_colx(n,n),tab_coly(n,n)
c
      do ind = 1,nd
         do m = 1,n
            do j = 1,n
               ff(m,j) = 0.0d0
               do i = 1,n
                  ff(m,j)=ff(m,j)+tab_colx(i,j)*coeff(i,m,ind)
               enddo
            enddo
         enddo
c
         do k = 1,n
         do j = 1,n
               pot(j,k,ind) = 0.0d0
               do m = 1,n
                  pot(j,k,ind)=pot(j,k,ind)
     1                 +ff(m,j)*tab_coly(m,k)
               enddo
         enddo
         enddo
c
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
      subroutine leg3d_to_potloc(nd,n,coeff,ff,ff2,pot,tabx,taby,tabz)
C*********************************************************************C
c     This routine computes the volume Gauss transform over a 
c     single box source distribution given as a Legendre series.
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays (tabx, taby,tabz). 
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in 3D, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x, y, or z in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, y, or z.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, y, or z.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       Legendre coefficients
c                 f = sum coeff(n,m,k) P_n(x) P_m(y) P_k(z)
c     ff          workspace
c     ff2          workspace
c     tabx    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in x.
c     taby    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in y.
c     tabz    precomputed table of 1D integrals
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,n,nd),pot(nd,n,n,n)
      real *8 ff(n,n,n),ff2(n,n,n),tabx(n,n),taby(n,n),tabz(n,n)
c
      do ind = 1,nd
c        transform in x
         do j3=1,n
            do j2=1,n
               do k1=1,n
                  cd=0
                  do j1=1,n
                     cd=cd+tabx(j1,k1)*coeff(j1,j2,j3,ind)
                  enddo
                  ff(k1,j2,j3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do j3=1,n
            do k2=1,n            
               do k1=1,n
                  cd=0
                  do j2=1,n
                     cd=cd+taby(j2,k2)*ff(k1,j2,j3)
                  enddo
                  ff2(k1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transfrom in z
         do k3=1,n
            do k2=1,n
               do k1=1,n
                  cd=0
                  do j3=1,n
                     cd=cd+tabz(j3,k3)*ff2(k1,k2,j3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo

      return
      end subroutine
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
         call legepols(xnodes(i),n,legev(1,i))
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
         call legepols(xnodes(i),n,legev(1,i))
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
         call legepols(xnodes(i),n,legev(1,i))
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

