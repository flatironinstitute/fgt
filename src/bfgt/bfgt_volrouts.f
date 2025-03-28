c     This file contains a set of subroutines that build various 1D tables
c     used in the box FGT in all dimensions, including polynomial values
c     to plane wave expansion, polynomial values to potential/gradient/hessian,
c     and plane wave expansion to potential/gradient/hessian.
c      
c
C*********************************************************************C
      subroutine mk_poly2pw_tables(n,ipoly,npw,nnodes,ws,ts,delta,
     1    boxdim,tab_poly2pw)
C*********************************************************************C
c     Construct the table converting polynomial values to plane-wave
c     expansion coefficients.
c
c     tab_poly2pw(n,j) = ws(j)*(D/2) * 
c            int_{-1}^1 p_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx,
c
c     where p_n is either Legendre polynomial of degree n 
c                      or Chebysheve polynomial of degree n.
c
c
c     INPUT:
c     n        dimension of coeff array
c     ipoly    polynomial type
c              0: Legendre polynomials
c              1: Chebyshev polynomials
c      
c     npw      number of plane waves
c     nnodes   number of nodes used in numerical quadrature
c     ws,ts    weights and nodes of plane wave quadrature
c     delta    Gaussian variance
c     boxdim   box dimension at current level
c
c     OUTPUT:
c     tab_poly2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 p_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 tab_poly2pw(n,npw),eye,eyem,z
      real *8 ws(npw),ts(npw)

      if (ipoly.eq.0) then
         call mk_leg2pw_tables(n,npw,nnodes,ws,ts,delta,boxdim,
     1       tab_poly2pw)
      elseif (ipoly.eq.1) then
         call mk_cheb2pw_tables(n,npw,nnodes,ws,ts,delta,boxdim,
     1       tab_poly2pw)
      endif
      
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_leg2pw_tables(n,npw,nnodes,ws,ts,delta,boxdim,
     1    tab_leg2pw)
C*********************************************************************C
c     Use spherical Bessel functions to construct the table.
c
c     tab_leg2pw(n,j) = ws(j)*(D/2) * 
c            int_{-1}^1 P_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c      
c     it is known that
c      
c     int_{-1}^1 P_n(x) exp(- i a x) dx = 1/i^n *sqrt(2pi/a) J_{n+1/2)(a)
c                                       = 2/i^n j_n(a),
c      
c     where i^2=-1, J_{n+1/2) is the Bessel J function of half integer order
c     and j_n(z)=sqrt(pi/(2z))J_{n+1/2}(z) is the spherical Bessel function
c     of order n.
c      
c     Here D is the box dimension at current level in the tree hierarchy.
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
      complex *16 tab_leg2pw(n,npw),eye,eyem,z,cd
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: xq(:),wq(:),u(:,:),v(:,:)
      complex *16, allocatable :: tabtemp(:,:),fjs(:)
      complex *16 fjder

      ifder=0
      allocate(tabtemp(n,npw))
      allocate(xq(n),wq(n),u(n,n),v(n,n))
      allocate(fjs(0:n))
c
      pi = 4*atan(1.0d0)
      
      scale=1.0d0
      eye = dcmplx(0.0d0,1.0d0)
      fac = boxdim/(2*dsqrt(delta))

      do j = npw/2+1,npw
         dd = ws(j)*boxdim
         
         if (abs(ts(j)).lt.1d-12) then
            tabtemp(1,j)=dd
            do m=2,n
               tabtemp(m,j)=0
            enddo
         else
            z = ts(j)*fac
            call besseljs3d(n,z,scale,fjs,ifder,fjder)
            eyem=1.0d0
            do m=1,n
               tabtemp(m,j)=dd*fjs(m-1)/eyem
               tabtemp(m,npw-j+1) = conjg(tabtemp(m,j))
               eyem=eyem*eye
            enddo
         endif
      enddo

      itype=2
      call legeexps(itype,n,xq,u,v,wq)

      do i=1,npw
      do j=1,n
         cd=0
         do k=1,n
            cd=cd+tabtemp(k,i)*u(k,j)
         enddo
         tab_leg2pw(j,i)=cd
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
      subroutine mk_cheb2pw_tables(n,npw,nnodes,ws,ts,delta,boxdim,
     1    tab_cheb2pw)
C*********************************************************************C
c     Use spherical Bessel functions to construct the table.
c
c     tab_cheb2pw(n,j) = ws(j)*(D/2) * 
c              int_{-1}^1 T_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c              where D is the box dimension at current level in
c              tree hierarchy.
c
c     That is, 
c     tab_cheb2pw(n,j) is the Fourier transform of T_n at a specific
c     frequency. 
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
c     tab_cheb2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 P_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 tab_cheb2pw(n,npw),zsum,eye,eyem,cd
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: chebv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      complex *16, allocatable :: legv(:)
      real *8, allocatable :: xq(:),wq(:),u(:,:),v(:,:)
      complex *16, allocatable :: tabtemp(:,:),fjs(:)
      complex *16 fjder

      allocate(tabtemp(n,npw),fjs(0:n))
      allocate(xq(n),wq(n),u(n,n),v(n,n))
      
c
      allocate(chebv(n,n))
      allocate(whts(n))
      allocate(xnodes(n))
      allocate(legv(n))
c
      itype = 2
      call legeexps(itype,n,xnodes,u,v,whts)
c
      do i = 1,n
         call chebpols(xnodes(i),n-1,chebv(1,i))
      enddo
c     v is the matrix converting Chebyshev polynomials into Legendre
c     polynomial expansions
      do j=1,n
         do i=1,n
            dd=0
            do k=1,n
               dd=dd+u(j,k)*chebv(i,k)
            enddo
            v(j,i)=dd
         enddo
      enddo
c
c
      pi = 4*atan(1.0d0)

      ifder=0
      scale=1.0d0
      eye = dcmplx(0.0d0,1.0d0)
      fac = boxdim/(2*dsqrt(delta))

      do j = npw/2+1,npw
         dd = ws(j)*boxdim
         
         if (abs(ts(j)).lt.1d-12) then
            do m=0,n-1,2
               tabtemp(m+1,j) = dd/(1-m*m)
            enddo
            do m=1,n-1,2
               tabtemp(m+1,j) = 0
            enddo
         else
            z = ts(j)*fac
            call besseljs3d(n,z,scale,fjs,ifder,fjder)
            eyem=1.0d0
c           legv contains the Fourier transform of the Legendre polynomials
            do m=1,n
               legv(m)=dd*fjs(m-1)/eyem
               eyem=eyem*eye
            enddo
c           now convert to the Fourier transform of the Chebyshev polynomials
            do i=1,n
               cd=0
               do k=1,n
                  cd=cd+v(k,i)*legv(k)
               enddo
               tabtemp(i,j) = cd
            enddo
         endif

         if (j .ne. npw-j+1) then
            do i=1,n
               tabtemp(i,npw-j+1) = conjg(tabtemp(i,j))
            enddo
         endif
      enddo
c
      itype=2
      call chebexps(itype,n,xq,u,v,wq)

      do i=1,npw
      do j=1,n
         cd=0
         do k=1,n
            cd=cd+tabtemp(k,i)*u(k,j)
         enddo
         tab_cheb2pw(j,i)=cd
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
      subroutine mk_pw2pot_tables(norder,npw,ts,xs,delta,boxdim,
     1    tab_pw2pot)
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

cccc         qqx = cdexp(eye*ts(npw2+1)*x)
cccc         qq1 = qqx
cccc         qqx = qqx*qqx
cccc         do j1=npw2+1,npw
cccc            tab_pw2pot(j1,i) = qq1
cccc            qq1 = qq1*qqx
cccc            
cccc            tab_pw2pot(npw-j1+1,i) = dconjg(tab_pw2pot(j1,i))
cccc         enddo
         do j=1,(npw+1)/2
            tab_pw2pot(j,i)=exp(eye*ts(j)*x)
            tab_pw2pot(npw-j+1,i)=dconjg(tab_pw2pot(j,i))
         enddo
      enddo

      return
      end subroutine
c
c
c
c
c
C*********************************************************************C
      subroutine mk_pw2pgh_tables(norder,npw,ts,xs,delta,boxdim,
     1    tab_pw2pot,tab_pw2potx,tab_pw2potxx)
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
c     tab_pw2pot   table for potentials
c     tab_pw2potx  table for gradients
c     tab_pw2potxx table for hessians
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 ts(npw),xs(norder)
      complex *16 tab_pw2pot(npw,norder),eye,qqx,qq1,dfac
      complex *16 tab_pw2potx(npw,norder)
      complex *16 tab_pw2potxx(npw,norder)

c
      eye = dcmplx(0.0d0,1.0d0)
c
      dsq = boxdim/2/dsqrt(delta)
      dfac = eye*dsq*2/boxdim
C
      npw2=npw/2
      do i=1,norder
         x = xs(i)*dsq

cccc         qqx = cdexp(eye*ts(npw2+1)*x)
         
cccc         qq1 = qqx
cccc         qqx = qqx*qqx

         do j1=npw2+1,npw
cccc        tab_pw2pot(j1,i) = qq1
            tab_pw2pot(j1,i) = exp(eye*ts(j1)*x)
            tab_pw2potx(j1,i)= ts(j1)*tab_pw2pot(j1,i)*dfac
            tab_pw2potxx(j1,i)= ts(j1)*tab_pw2potx(j1,i)*dfac

cccc            qq1 = qq1*qqx
            
              tab_pw2pot(npw-j1+1,i) = conjg(tab_pw2pot(j1,i))
             tab_pw2potx(npw-j1+1,i) = conjg(tab_pw2potx(j1,i))
            tab_pw2potxx(npw-j1+1,i) = conjg(tab_pw2potxx(j1,i))
         enddo
      enddo
ccc      call prin2(' in mk tab_pwpot *',tab_pw2pot,2*npw*norder)
ccc      call prin2(' in mk tab_pwpotx *',tab_pw2potx,2*npw*norder)
ccc      call prin2(' in mk tab_pwpotxx *',tab_pw2potxx,2*npw*norder)

      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_all(eps,ipoly,n,nnodes,delta,boxdim,
     1    mrefinelev,nloctab,tab_loc,tabx_loc,tabxx_loc,ind_loc)
C*********************************************************************C
c     This routine is an optimized table generator that generates
c     all 1D local interaction tables
c
c
c     We assume that we start from a level-restricted tree. Then every box 
c     may be refined mrefinelev times. The 1D tables we construct here
c     will compute all integrals of the form
c
c     \int_{B_S} exp(-(x-y)^2/delta)p_n(y) dy
c
c     for y in the source box in the original level-restricted tree, while 
c     the target point x lies on an interval obtained by mrefinelev times
c     of the target box/interval in the original level-restricted tree.
c
c     In the original 2:1 level-restricted tree, there are only three cases:
c     colleague - source box and target box are at the same level, 3 of them;
c     stob - small source box to big target box, 4 of them;
c     btos - large source box to small target box, 4 of them.
c
c     Now, each of these 11 cases will be further refined mrefinelev times,
c     leading to 11*(1+2+...+2^mrefinelev) cases.
c
c
c
c
c     (a) same level interaction
c     tab_coll(n,j,k) = 
c              int_{-D/2}^{D/2} P_n(x*2/D) exp( -(\xi_j -x)^2/delta)
c              where D is the box dimension at current level in
c              tree hierarchy and the target \xi_j is either on 
c              [-D/2,D/2]   -> tab_loc(n,n,0)
c              [-3D/2,-D/2] -> tab_loc(n,n,-4)
c              [D/2,3D/2]   -> tab_coll(n,n,4)
c
c     Here we assume that the source box is centered at the origin.        
c      
c     (b) small source to big target interaction              
c     tab_stob(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-5D/4,-D/4]     -> tab_loc(n,n,-6)
c              [-3D/4, D/4]     -> tab_loc(n,n,-2)
c              [ -D/4,3D/4]     -> tab_loc(n,n,2)
c              [  D/4,5D/4]     -> tab_loc(n,n,6)
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
c
c     (b) big source to small target interaction              
c     tab_btos(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of the TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-2D,-D]    -> tab_loc(n,n,-3)
c              [ -D, 0]    -> tab_loc(n,n,-1)
c              [  0, D]    -> tab_loc(n,n,1)
c              [  D,2D]    -> tab_loc(n,n,3)
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
c     ipoly          polynomial type
c                    0: Legendre polynomials
c                    1: Chebyshev polynomials
c     n              dimension of coeff array
c     nnodes         number of nodes used in numerical quadrature
c     delta          Gaussian variance
c     boxdim         target box dimension at current level
c     mrefinelev     maximum number of refinement
c
c     OUTPUT:
c     tab_loc        tables for computing potential values
c     tabx_loc       tables for computing gradient values
c     tabxx_loc      tables for computing hessian values
c----------------------------------------------------------------------c
      implicit none
      integer mrefinelev,nloctab,m0
      real *8 tab_loc(n,n,-nloctab:nloctab)
      real *8 tabx_loc(n,n,-nloctab:nloctab)
      real *8 tabxx_loc(n,n,-nloctab:nloctab)
      integer ind_loc(2,n+1,-nloctab:nloctab)

      real *8 eps,delta,boxdim
      real *8 scale(-nloctab:nloctab),utmp,vtmp
      real *8 cen(-nloctab:nloctab)
      real *8, allocatable :: tnodes(:), ws(:) 
      real *8, allocatable :: xnodes(:), wts(:) 
      real *8, allocatable :: umat(:,:),vmat(:,:)
      real *8, allocatable :: ptmp(:)
      real *8, allocatable :: polyv(:,:)

      real *8, allocatable :: tab(:,:)
      real *8, allocatable :: tabx(:,:)
      real *8, allocatable :: tabxx(:,:)

      integer ind(2,n+1),indx(2,n+1),indxx(2,n+1)
      
      integer itype,nquad,i,j,k,m,ipoly,n,ifzero,nnodes,nt,i1,i2
      integer i0,j0
c
      allocate(tab(n,n))
      allocate(tabx(n,n))
      allocate(tabxx(n,n))

c     target nodes
      allocate(tnodes(n),ws(n),umat(n,n),vmat(n,n))
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,n,tnodes,umat,vmat,ws)
      if (ipoly.eq.1) call chebexps(itype,n,tnodes,umat,vmat,ws)

c     quadrature nodes and weights
      nquad = 50
cccc      nquad = 100000
      allocate(xnodes(nquad),wts(nquad))
      allocate(ptmp(n),polyv(nquad,n))

      itype=1
      call legeexps(itype,nquad,xnodes,utmp,vtmp,wts)
c     polynomial values at quadrature nodes
      do i = 1,nquad
         if (ipoly.eq.0) then
            call legepols(xnodes(i),n-1,ptmp)
         elseif (ipoly.eq.1) then
            call chebpols(xnodes(i),n-1,ptmp)
         endif
         do j=1,n
            polyv(i,j)=ptmp(j)
         enddo
      enddo
c
c     scale is the ratio of the source box size to the target box size.

      m0=mrefinelev+1
      do i0=m0,-1,-1
      do j0=1,2**m0
         if (i0.gt.0) then
            k=-2**(m0-i0)-2**(m0-i0+1)*(j0-1)
         elseif (i0.eq.0) then
            k=-2**(m0-i0+1)*(j0-1)
         elseif (i0.eq.-1) then
            k=-2**m0-2**(m0+1)*(j0-1)
         endif
         if (k.lt.-nloctab) exit
         scale(k)=2.0d0**i0
         cen(k)=k/2.0d0**m0
cccc         print *, i0,j0, k, cen(k), scale(k)
c        make local tables for potential, gradient and hessian
         call mk_loctab(ipoly,n,delta,n,tnodes,boxdim,cen(k),
     1       scale(k),tab,tabx,tabxx,
     2       nquad,xnodes,wts,polyv)
c        multiply by u to obtain the tables acting on values
         call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tab,n,
     1       0.0d0,tab_loc(1,1,k),n)
         call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tabx,n,
     1       0.0d0,tabx_loc(1,1,k),n)
         call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tabxx,n,
     1       0.0d0,tabxx_loc(1,1,k),n)

c     use symmetry to construct the tables for targets on the right side
         do j=1,n
         do m=1,n
            tab_loc(m,j,-k) = tab_loc(n-m+1,n-j+1,k)
            tabx_loc(m,j,-k) = -tabx_loc(n-m+1,n-j+1,k)
            tabxx_loc(m,j,-k) = tabxx_loc(n-m+1,n-j+1,k)
         enddo
         enddo

c      
c     STEP 4: find sparse patterns of all local interaction matrices
c      
         call compute_sparse_pattern(eps,n,n,tab_loc(1,1,-k),
     1       ind,ifzero)
         call compute_sparse_pattern(eps,n,n,tabx_loc(1,1,-k),
     1       indx,ifzero)
         call compute_sparse_pattern(eps,n,n,tabxx_loc(1,1,-k),
     1       indxx,ifzero)
         do j=1,n+1
            i2=ind(2,j)
            if (i2.lt.indx(2,j)) i2=indx(2,j)
            if (i2.lt.indxx(2,j)) i2=indxx(2,j)
            ind_loc(2,j,-k)=i2

            i1=ind(1,j)
            if (indx(1,j).eq.0) then
            elseif (i1.gt.indx(1,j)) then
               i1=indx(1,j)
            endif
            if (indxx(1,j).eq.0) then
            elseif (i1.gt.indxx(1,j)) then
               i1=indxx(1,j)
            endif
            if (i2.gt.0 .and.i1.eq.0) i1=1
            ind_loc(1,j,-k)=i1
         enddo

         call compute_sparse_pattern(eps,n,n,tab_loc(1,1,k),
     1       ind,ifzero)
         call compute_sparse_pattern(eps,n,n,tabx_loc(1,1,k),
     1       indx,ifzero)
         call compute_sparse_pattern(eps,n,n,tabxx_loc(1,1,k),
     1       indxx,ifzero)
         do j=1,n+1
            i2=ind(2,j)
            if (i2.lt.indx(2,j)) i2=indx(2,j)
            if (i2.lt.indxx(2,j)) i2=indxx(2,j)
            ind_loc(2,j,k)=i2

            i1=ind(1,j)
            if (indx(1,j).eq.0) then
            elseif (i1.gt.indx(1,j)) then
               i1=indx(1,j)
            endif
            if (indxx(1,j).eq.0) then
            elseif (i1.gt.indxx(1,j)) then
               i1=indxx(1,j)
            endif
            if (i2.gt.0 .and.i1.eq.0) i1=1
            ind_loc(1,j,k)=i1
         enddo
            
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
c
C*********************************************************************C
      subroutine mk_loctab_all_old(eps,ipoly,n,nnodes,delta,boxdim,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc)
C*********************************************************************C
c     This routine is an optimized table generator that generates
c     all 1D local interaction tables
c
c     (a) same level interaction
c     tab_coll(n,j,k) = 
c              int_{-D/2}^{D/2} P_n(x*2/D) exp( -(\xi_j -x)^2/delta)
c              where D is the box dimension at current level in
c              tree hierarchy and the target \xi_j is either on 
c              [-D/2,D/2]   -> tab_loc(n,n,0)
c              [-3D/2,-D/2] -> tab_loc(n,n,-4)
c              [D/2,3D/2]   -> tab_coll(n,n,4)
c
c     Here we assume that the source box is centered at the origin.        
c      
c     (b) small source to big target interaction              
c     tab_stob(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-5D/4,-D/4]     -> tab_loc(n,n,-6)
c              [-3D/4, D/4]     -> tab_loc(n,n,-2)
c              [ -D/4,3D/4]     -> tab_loc(n,n,2)
c              [  D/4,5D/4]     -> tab_loc(n,n,6)
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
c
c     (b) big source to small target interaction              
c     tab_btos(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of the TARGET BOX 
c              at current level in tree hierarchy and 
c              \xi_j is either on 
c              [-2D,-D]    -> tab_loc(n,n,-3)
c              [ -D, 0]    -> tab_loc(n,n,-1)
c              [  0, D]    -> tab_loc(n,n,1)
c              [  D,2D]    -> tab_loc(n,n,3)
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
c     ipoly    polynomial type
c              0: Legendre polynomials
c              1: Chebyshev polynomials
c     n        dimension of coeff array
c     nnodes   number of nodes used in numerical quadrature
c     delta    Gaussian variance
c     boxdim   target box dimension at current level
c
c     OUTPUT:
c     tab_loc  tables for computing potential values
c     tabx_loc  tables for computing gradient values
c     tabxx_loc  tables for computing hessian values
c----------------------------------------------------------------------c
      implicit none
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      real *8 tabxx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)

      real *8 eps,delta,boxdim
      real *8 scale(-6:6),utmp,vtmp
      real *8 cen(-6:6)
      real *8, allocatable :: tnodes(:), ws(:) 
      real *8, allocatable :: xnodes(:), wts(:) 
      real *8, allocatable :: umat(:,:),vmat(:,:)
      real *8, allocatable :: ptmp(:)
      real *8, allocatable :: polyv(:,:)

      real *8, allocatable :: tab(:,:)
      real *8, allocatable :: tabx(:,:)
      real *8, allocatable :: tabxx(:,:)

      integer ind(2,n+1),indx(2,n+1),indxx(2,n+1)
      
      integer itype,nquad,i,j,k,m,ipoly,n,ifzero,nnodes,nt,i1,i2
c
      allocate(tab(n,n))
      allocate(tabx(n,n))
      allocate(tabxx(n,n))

c     target nodes
      allocate(tnodes(n),ws(n),umat(n,n),vmat(n,n))
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,n,tnodes,umat,vmat,ws)
      if (ipoly.eq.1) call chebexps(itype,n,tnodes,umat,vmat,ws)

c     quadrature nodes and weights
      nquad = 50
cccc      nquad = 100000
      allocate(xnodes(nquad),wts(nquad))
      allocate(ptmp(n),polyv(nquad,n))

      itype=1
      call legeexps(itype,nquad,xnodes,utmp,vtmp,wts)
c     polynomial values at quadrature nodes
      do i = 1,nquad
         if (ipoly.eq.0) then
            call legepols(xnodes(i),n-1,ptmp)
         elseif (ipoly.eq.1) then
            call chebpols(xnodes(i),n-1,ptmp)
         endif
         do j=1,n
            polyv(i,j)=ptmp(j)
         enddo
      enddo
c
c     scale is the ratio of the source box size to the target box size.
c     source and target at the same level 
      scale(-4) = 1.0d0
      scale(0) = 1.0d0
c     small source to big target 
      scale(-6) = 0.5d0
      scale(-2) = 0.5d0
c     big source to small target 
      scale(-3) = 2.0d0
      scale(-1) = 2.0d0
      
c     The scaled source interval is always [-1,1].
c     coll table -1, the scaled target interval is [-3,-1], centered at -2
c     coll table  0, the scaled target interval is [-1, 1], centered at  0
c     coll table  1, the scaled target interval is [ 1, 3], centered at  2
      cen(-4) = -2.0d0
      cen(0)  =  0.0d0
c     stob table 1, the scaled target interval is [-5,-1], centered at -3
c     stob table 2, the scaled target interval is [-3, 1], centered at -1
c     stob table 3, the scaled target interval is [-1, 3], centered at  1
c     stob table 4, the scaled target interval is [ 1, 5], centered at  3
      cen(-6) = -3.0d0
      cen(-2) = -1.0d0
c     btos table 1, the scaled target interval is [-2,-1], centered at -1.5
c     btos table 2, the scaled target interval is [-1, 0], centered at -0.5
c     btos table 3, the scaled target interval is [ 0, 1], centered at  0.5
c     btos table 4, the scaled target interval is [ 1, 2], centered at  1.5
      cen(-3) = -1.5d0
      cen(-1) = -0.5d0

      do i=-6,0
        if (i.ne.-5) then
          nt=n
          if (i.eq.0) nt=(n+1)/2
          call mk_loctab(ipoly,n,delta,nt,tnodes,boxdim,cen(i),
     1       scale(i),tab,tabx,tabxx,
     2       nquad,xnodes,wts,polyv)
          if (i.eq.0) then
c           obtain the other half of table 0 by symmetry
            do j=1,n/2
            do m=1,n,2
               tab(m,n-j+1)  = tab(m,j) 
               tabx(m,n-j+1)  = -tabx(m,j) 
               tabxx(m,n-j+1)  = tabxx(m,j) 
            enddo
            do m=2,n,2
               tab(m,n-j+1)  = -tab(m,j) 
               tabx(m,n-j+1)  = tabx(m,j) 
               tabxx(m,n-j+1)  = -tabxx(m,j) 
            enddo
            enddo
          endif
c         multiply by u to obtain the tables acting on values
          call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tab,n,
     1        0.0d0,tab_loc(1,1,i),n)
          call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tabx,n,
     1        0.0d0,tabx_loc(1,1,i),n)
          call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,tabxx,n,
     1        0.0d0,tabxx_loc(1,1,i),n)
        endif    
      enddo

c     use symmetry to construct the tables for targets on the right side
      do i=1,6
        if (i.ne.5) then
          do j=1,n
          do m=1,n
            tab_loc(m,j,i) = tab_loc(n-m+1,n-j+1,-i)
            tabx_loc(m,j,i) = -tabx_loc(n-m+1,n-j+1,-i)
            tabxx_loc(m,j,i) = tabxx_loc(n-m+1,n-j+1,-i)
          enddo
          enddo
        endif
      enddo

c      
c     STEP 4: find sparse patterns of all local interaction matrices
c      
      do k=-6,6
         if (abs(k).ne.5) then
            call compute_sparse_pattern(eps,n,n,tab_loc(1,1,k),
     1          ind,ifzero)
            call compute_sparse_pattern(eps,n,n,tabx_loc(1,1,k),
     1          indx,ifzero)
            call compute_sparse_pattern(eps,n,n,tabxx_loc(1,1,k),
     1          indxx,ifzero)
            do j=1,n+1
               i2=ind(2,j)
               if (i2.lt.indx(2,j)) i2=indx(2,j)
               if (i2.lt.indxx(2,j)) i2=indxx(2,j)
               ind_loc(2,j,k)=i2

               i1=ind(1,j)
               if (indx(1,j).eq.0) then
               elseif (i1.gt.indx(1,j)) then
                  i1=indx(1,j)
               endif
               if (indxx(1,j).eq.0) then
               elseif (i1.gt.indxx(1,j)) then
                  i1=indxx(1,j)
               endif
               if (i2.gt.0 .and.i1.eq.0) i1=1
               ind_loc(1,j,k)=i1

            enddo
c            call compute_sparse_pattern(eps,n,n,tabx_loc(1,1,k),
c     1          ind_loc(1,1,k),ifzero)
c            call prinf('ind_loc=*',ind_loc(1,1,k),2*(n+1))
         endif
      enddo

      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine mk_loctab(ipoly,n,delta,nt,tnodes,tboxdim,
     1    tboxcen,scale,tab,tabx,tabxx,
     2    nquad,xnodes,wts,polyv)
C*********************************************************************C
c     This routine is the table generator workhorse.
c
c     tab(n,j) = 
c              int_{source box} p_n(x) exp( -(targ_j -x)^2/delta) dx
c              
c     Here we assume that the source box is always centered at the origin.        
c              
c              
c     INPUT:
c     n         polynomial approximation order
c     ipoly     polynomial type
c               0: Legendre polynomial      
c               1: Chebshev polynomial
c     delta     Gaussian variance
c     nt        number of targets
c     tnodes    shifted and scaled target nodes on [-1,1]
c     tboxdim   size of the TARGET box 
c     tboxcen   center of the scaled TARGET box, relative to the 
c               scaled source box at [-1,1] 
c     scale     ratio of the source box size to the target box size
c     nquad     number of quadrature nodes on the source box
c     xnodes    quadrature nodes on [-1,1]
c     wts       quadrature weights on [-1,1]
c     polyv     [n,nquad] polynomial values at quadrature nodes
c
c     OUTPUT:
c     tab       the table for compute the potential 
c     tabx      the table for compute the first derivative of the potential 
c     tabxx     the table for compute the second derivative of the potential 
c----------------------------------------------------------------------c
      implicit none
      integer ipoly,n,nt,nquad
      real *8 delta,tboxdim,tboxcen,scale
      real *8 tab(n,nt),tabx(n,nt),tabxx(n,nt)
      real *8 tnodes(n),fint(100),fx(100),fxx(100),lambda
      real *8 xnodes(nquad),wts(nquad), polyv(nquad,n)
      integer i,j,m
      real *8 bs,bsh,reps,dmax,targ,dis,dx,sigma
      real *8 rsum,rsumx,rsumxx,dd
      real *8, allocatable :: targs(:)
      real *8, allocatable :: wexp(:,:)
c
      allocate(targs(n))
c     target position relative to the standard source box [-1,1]
      do i=1,nt
         targs(i) = tnodes(i)/scale + tboxcen
      enddo

c     size of the source box
      bs = tboxdim*scale
      
      sigma = 4*delta/bs**2
      lambda = sqrt(sigma)
cccc      print *, 'lambda=', lambda
      
      reps=1d-16
      dmax=sqrt(log(1/reps)*sigma)
cccc      print *, 'dmax=',dmax
c
c
      if (lambda.le.0.125d0) then
cccc      if (lambda.le.0.125d-5) then
         do j=1,nt
            targ = targs(j)
c           distance between the target point and the source box
            dis=0
            if (targ .lt. -1) dis = -1-targ
            if (targ .gt.  1) dis = targ-1
            if (dis.gt.dmax) then
               do m=1,n
                  tab(m,j) = 0
                  tabx(m,j) = 0
                  tabxx(m,j) = 0
               enddo
            else
               call mk_polytab_recur(ipoly,n,targ,lambda,
     1             bs,fint,fx,fxx)
               do m=1,n
                  tab(m,j) = fint(m)
                  tabx(m,j) = fx(m)
                  tabxx(m,j) = fxx(m)
               enddo
            endif
         enddo
      else
         allocate(wexp(nquad,nt))
c        compute exponentials only once
         do j=1,nt
            targ = targs(j)
            do i=1,nquad
c              distance between targ and quadrature source nodes
               dx = targ-xnodes(i)
               if (abs(dx) .lt. dmax) then
                  wexp(i,j)=exp(-dx*dx/sigma)*wts(i)
               else
                  wexp(i,j)=0
               endif
            enddo
         enddo
c     
         bsh = bs/2
         do j = 1,nt
            targ = targs(j)
c           distance between the target point and the source box
            dis=0
            if (targ .lt. -1) dis = -1-targ
            if (targ .gt.  1) dis = targ-1
            if (dis.gt.dmax) then
               do m=1,n
                  tab(m,j)  = 0
                  tabx(m,j)  = 0
                  tabxx(m,j)  = 0
               enddo
            else
               do m = 1,n
                  rsum = 0.0d0
                  rsumx = 0.0d0
                  rsumxx = 0.0d0
                  do i = 1,nquad
                     dx = targ-xnodes(i)
                     dd = -2*dx/sigma*2/bs
                     rsum = rsum + polyv(i,m)*wexp(i,j)
                     rsumx = rsumx + polyv(i,m)*wexp(i,j)*dd
                     rsumxx = rsumxx + polyv(i,m)*wexp(i,j)
     1                   *(dd*dd-2/delta)
                  enddo
                  tab(m,j) = rsum*bsh
                  tabx(m,j) = rsumx*bsh
                  tabxx(m,j) = rsumxx*bsh
               enddo
            endif
         enddo
      endif


      return
      end subroutine
c
c
c
c
      subroutine mk_polytab_recur(ipoly,m,targ,lambda,boxdim,
     1    fint,fx,fxx)
c     calculate the values of the integral
c
c     \frac{1}{\lambda}\int_{-1}^1 p_n(x) e^{-(targ-x)^2/\lambda^2}dx
c
c     for n=0,1,...,m-1.
c      
c     p is either the Legendre polynomial or the Chebyshev polynomial.
c      
c     Algorithm: use the five term recurrence formula
c     
c     Assumption: lambda<=1/8
c
      implicit none
      integer ipoly,m
      real *8 fint(*),fx(*),fxx(*),lambda,boxdim,targ

      if (ipoly.eq.0) call
     1    mk_legetab_recur(m,targ,lambda,boxdim,fint,fx,fxx)
      
      if (ipoly.eq.1) call
     1    mk_chebtab_recur(m,targ,lambda,boxdim,fint,fx,fxx)

      return
      end subroutine
c
c
c
c
      subroutine mk_legetab_recur(m,targ,lambda,boxdim,fint,fx,fxx)
      implicit real *8 (a-h,o-z)
cccc      sqrtpih=sqrt(4.0d0*atan(1.0d0))/2
      data sqrtpih/0.886226925452757940959713778283913d0/
c     calculate the values of the integral
c
c     I_n(targ)=\int_{-1}^1 P_n(x) e^{-(targ-x)^2/\lambda^2}dx
c
c     for n=0,1,...,m-1.
c      
c     Algorithm: use the five term recurrence formula to calculate I_n/lambda
c     
c     Assumption: lambda<=1/8
c
      real *8 fint(*),fx(*),fxx(*),lambda,lambda2
      real *8 fp(100),fpp(100)

c     fint(n) = I_n(targ)/lambda
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
      
      fint(4) = -d1 + 2*targ*(4.0d0*fint(3)/3-fint(1)/3)
     1    -(0.4d0-4*lambda2)*fint(2)
      fint(4) = 5.0d0*fint(4)/8
      
      do k=5,m
         n=k-3
         fint(k)=targ*(fint(k-1)-fint(k-3))
     1       + (2*n+1)*(lambda2/2 + 1.0d0/((2*n+3)*(2*n-1)))*fint(k-2) 
     2       + (n-1.0d0)*fint(k-4)/(2*n-1) 
         fint(k) = (2*n+3.0d0)*fint(k)/(n+2)
      enddo

c     fp(n) = \int_{-1}^1 G(targ-x) P'_n(x)dx
      fp(1)=0
      fp(2)=fint(1)*lambda
      do k=3,m
         fp(k)=fp(k-2)+(2*k-3)*fint(k-1)*lambda
      enddo

c     fpp(n) = \int_{-1}^1 G(targ-x) P''_n(x)dx
      fpp(1)=0
      fpp(2)=0
      do k=3,m
         fpp(k)=fpp(k-2)+(2*k-3)*fp(k-1)
      enddo

      isign=1
      do k=1,m
         px1 = (k-1)*k/2.0d0
         fxx(k)=fpp(k)-px1*(expb+isign*expa)
     1       -2*(b*expb-isign*a*expa)/lambda
         isign=-isign
      enddo

c     finally, multiply lambda back and also the proper weight adjustment
c     for the source box
      ww = boxdim*lambda/2
      do k=1,m
         fint(k)=fint(k)*ww
         fxx(k)=fxx(k)*2/boxdim
      enddo

      do k=2,m,2
         fx(k)=fp(k)-d2/lambda
      enddo
      do k=1,m,2
         fx(k)=fp(k)-d1/lambda
      enddo

cccc      call prin2('fint=*',fint,m)
cccc      call prin2('fx=*',fx,m)
cccc      call prin2('fxx=*',fxx,m)
      return
      end subroutine
c
c
c
      subroutine mk_chebtab_recur(m,targ,lambda,boxdim,fint,fx,fxx)
      implicit real *8 (a-h,o-z)
cccc      sqrtpih=sqrt(4.0d0*atan(1.0d0))/2
      data sqrtpih/0.886226925452757940959713778283913d0/
c     calculate the values of the integral
c
c     I_n(targ)=\int_{-1}^1 T_n(x) e^{-(targ-x)^2/\lambda^2}dx
c
c     for n=0,1,...,m-1.
c      
c     Algorithm: use the five term recurrence formula to calculate
c      I_n/lambda.
c     
c     Assumption: lambda<=1/8
c
      real *8 fint(*),fx(*),fxx(*),lambda,lambda2
      real *8 fp(100),fpp(100)
      
c     fint(n) = I_n(targ)/lambda
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
     1    +(lambda2-1+2*targ**2)*fint(1)

      fint(4) = -d1 + 2*targ*fint(3)
     1    -(1-4*lambda2)*fint(2)

      isign=1
      do k=5,m
         if (isign.eq.1)  dc = d2
         if (isign.eq.-1) dc = d1
         n=k-1
         dd=(n-1.0d0)/(n-3)
         fint(k)=2*targ*fint(k-1) 
     1       + 2*((n-1)*lambda2 + 1.0d0/(n-3))*fint(k-2) 
     2       + dd*(fint(k-4) - 2*targ*fint(k-3))
     3       -(1-dd)*dc
         isign=-isign
      enddo

c     fp(n) = \int_{-1}^1 G(targ-x) T'_n(x)dx
      fp(1)=0
      fp(2)=fint(1)*lambda
      fp(3)=4*fint(2)*lambda
      do k=4,m
         fp(k)=(fp(k-2)/(k-3)+2*fint(k-1)*lambda)*(k-1)
      enddo

c     fpp(n) = \int_{-1}^1 G(targ-x) T''_n(x)dx
      fpp(1)=0
      fpp(2)=0
      fpp(3)=4*fint(1)*lambda
      do k=4,m
         fpp(k)=(fpp(k-2)/(k-3)+2*fp(k-1))*(k-1)
      enddo

      isign=1
      do k=1,m
c        T'_n(1)
         tp1 = (k-1)**2
         fxx(k)=fpp(k)-tp1*(expb+isign*expa)
     1       -2*(b*expb-isign*a*expa)/lambda
c     ignore the boundary terms
c         fxx(k)=fpp(k)

         isign=-isign
      enddo

      
c     finally, multiply lambda back and also the proper weight adjustment
c     for the source box
      ww = boxdim*lambda/2
      sc=2/boxdim
      do k=1,m
         fint(k)=fint(k)*ww
         fxx(k)=fxx(k)*sc
      enddo

      do k=2,m,2
         fx(k)=fp(k)-d2/lambda
c         fx(k)=fp(k)
      enddo
      do k=1,m,2
         fx(k)=fp(k)-d1/lambda
c         fx(k)=fp(k)
      enddo


      return
      end subroutine
c
c
c
C*********************************************************************C
c      
      subroutine compute_sparse_pattern(eps,m,n,a,indc,ifzero)
c     given an mxn matrix a, determine its sparse pattern
c      
c     input:
c     eps - zero threshold
c           if all entries are < eps in magnitude, then the
c           matrix is regarded as the zero matrix. Otherwise,
c           the entries < eps * ||a||_infinity are regarded 
c           as zero.
c
c     m,n - number of rows and columns of the input matrix
c     a - mxn input matrix
c
c     output
c     indc - (2,n) start and end row indices for each column
c     ifzero - 1 if the matrix is 0 w.r.t the precision eps
c              0 if the matrix is not the zero matrix
      implicit real *8 (a-h,o-z)
      real *8 a(m,n)
      integer indc(2,n+1)

      ifzero = 0
      do i=1,n+1
         indc(1,i)=0
         indc(2,i)=-1
      enddo

      dmax = 0
      do i=1,n
         do j=1,m
            c = abs(a(j,i))
            if (dmax .lt. c) dmax=c
         enddo
      enddo

cccc      print *, 'matrix infinity norm=', dmax
      if (dmax .le. 2d-16) then
         ifzero = 1
         return
      endif

      do i=1,n
         do j=1,m
            c = abs(a(j,i))
            if (c .gt. eps*dmax/n) then
               indc(1,i) = j
               exit
            endif
         enddo
      enddo

      do i=1,n
         do j=m,1,-1
            c = abs(a(j,i))
            if (c .gt. eps*dmax/n) then
               indc(2,i) = j
               exit
            endif
         enddo
      enddo

      do i=1,n
         if (indc(1,i).gt.0) then
            indc(1,n+1)=i
            exit
         endif
      enddo
      
      do i=n,1,-1
         if (indc(2,i).gt.0) then
            indc(2,n+1)=i
            exit
         endif
      enddo
      
c      call prinf('indc=*',indc,2*(n+1))
      
      return
      end subroutine
c
c
c      
