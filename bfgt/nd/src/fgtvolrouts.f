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
      subroutine mk_poly2pw_tables(n,ipoly,npw,nnodes,ws,ts,delta,
     1    boxdim,tab_poly2pw)
C*********************************************************************C
c     Construct the table.
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
      complex *16, allocatable :: tabtemp(:,:)
      complex *16 fjs(0:n),fjder

      ifder=0
      allocate(tabtemp(n,npw))
      allocate(xq(n),wq(n),u(n,n),v(n,n))
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
      complex *16, allocatable :: tabtemp(:,:)
      complex *16 fjs(0:n),fjder

      allocate(tabtemp(n,npw))
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
c     tab_pw2pot 
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

         qqx = cdexp(eye*ts(npw2+1)*x)
         
         qq1 = qqx
         qqx = qqx*qqx

         do j1=npw2+1,npw
            tab_pw2pot(j1,i) = qq1
            tab_pw2potx(j1,i)= ts(j1)*tab_pw2pot(j1,i)*dfac
            tab_pw2potxx(j1,i)= ts(j1)*tab_pw2potx(j1,i)*dfac

            qq1 = qq1*qqx
            
              tab_pw2pot(npw-j1+1,i) = dconjg(tab_pw2pot(j1,i))
             tab_pw2potx(npw-j1+1,i) = dconjg(tab_pw2potx(j1,i))
            tab_pw2potxx(npw-j1+1,i) = dconjg(tab_pw2potxx(j1,i))
         enddo

cccc         do j=1,npw
cccc            tab_pw2pot(j,i)=exp(eye*ts(j)*x)
cccc         enddo
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
      subroutine mk_loctab_coll(eps,ipoly,n,nnodes,delta,boxdim,
     1    tab_colleague,indc_coll)
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
      integer indc_coll(2,n+1,-1:1)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: polyv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
      real *8, allocatable :: ws(:),u(:,:),v(:,:)
      real *8, allocatable :: tabtmp0(:,:),tabtmpm1(:,:)
c
      allocate(ws(n),u(n,n),v(n,n))
      allocate(tabtmp0(n,n),tabtmpm1(n,n))
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,n,xnodest,u,v,ws)
      if (ipoly.eq.1) call chebexps(itype,n,xnodest,u,v,ws)

      sigma = 4*delta/boxdim**2
      lambda = sqrt(sigma)
      reps=2d-16
      dmax=sqrt(log(1/reps)*sigma)
cccc      print *, 'dmax=',dmax
c
c
c     colleague table -1, the scaled target interval is on [-3,-1], the scaled
c     source interval is on [-1,1]. Use the same quadrature nodes for all
c     target points
      if (lambda.le.0.125d0) then
         do j=1,n
            xim1 = xnodest(j) - 2
            dx = 1-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tabtmpm1(m,j)  = 0
               enddo
            else
               if (ipoly .eq. 0) then
                  call mk_legetab_recur(n,xim1,lambda,boxdim,fint)
               elseif (ipoly .eq. 1) then
                  call mk_chebtab_recur(n,xim1,lambda,boxdim,fint)
               endif
               do m=1,n
                  tabtmpm1(m,j)  = fint(m)
               enddo
            endif
         enddo
      else
         nquad = 50
         allocate(polyv(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c
         itype=1
         call legeexps(itype,nquad,xnodes,utmp,vtmp,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/2
         enddo
      
         do i = 1,nquad
            if (ipoly.eq.0) call legepols(xnodes(i),n-1,polyv(1,i))
            if (ipoly.eq.1) call chebpols(xnodes(i),n-1,polyv(1,i))
         enddo
      
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
                  tabtmpm1(m,j) = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + polyv(m,i)*wexp(i,j)
                  enddo
                  tabtmpm1(m,j) = rsum1
               enddo
            endif
         enddo
      endif
c     colleague table 0, the scaled target interval is on [-1,1], the scaled source
c     interval is on [-1, 1]. use different quadrature nodes for each target.
c      
c     by symmetry, only needs to construct the table for half of the target points.
c      
      if (lambda.le.0.125d0) then
         do j=1,n/2
            xi = xnodest(j)
            if (ipoly.eq.0) then
               call mk_legetab_recur(n,xi,lambda,boxdim,fint)
            elseif (ipoly.eq.1) then
               call mk_chebtab_recur(n,xi,lambda,boxdim,fint)
            endif
            do m=1,n
               tabtmp0(m,j)  = fint(m)
            enddo
         enddo
      else
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
                  rsum1 = rsum1 + polyv(m,i)*wexp(i,j)
               enddo
               tabtmp0(m,j)  = rsum1
            enddo
         enddo
      endif

c     obtain the other half of table 0 by symmetry
      do j=1,n/2
         do m=1,n,2
            tabtmp0(m,n-j+1)  = tabtmp0(m,j) 
         enddo
         do m=2,n,2
            tabtmp0(m,n-j+1)  = -tabtmp0(m,j) 
         enddo
      enddo

c     multiply by u to obtain the tables acting on values
      do j=1,n
      do m=1,n
         dd=0
         do k=1,n
            dd=dd+tabtmp0(k,j)*u(k,m)
         enddo
         tab_colleague(m,j,0)=dd
      enddo
      enddo
      
      do j=1,n
      do m=1,n
         dd=0
         do k=1,n
            dd=dd+tabtmpm1(k,j)*u(k,m)
         enddo
         tab_colleague(m,j,-1)=dd
      enddo
      enddo
      
c     use symmetry to construct the table +1 for targets on the right side
      do j=1,n
         do m=1,n
            tab_colleague(m,j,1) = tab_colleague(n-m+1,n-j+1,-1)
         enddo
      enddo

      do k=-1,1
         call compute_sparse_pattern(eps,n,n,tab_colleague(1,1,k),
     1       indc_coll(1,1,k),ifzero)
      enddo
      
      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine mk_loctab_stob(eps,ipoly,n,nnodes,delta,boxdim,
     1    tab_stob,indc_stob)
C*********************************************************************C
c     This routine is an optimized table generator.
c
c     tab_stob(n,j,k) = 
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
      integer indc_stob(2,n+1,4)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: polyv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
      real *8, allocatable :: ws(:),u(:,:),v(:,:)
      real *8, allocatable :: tabtmp1(:,:),tabtmp2(:,:)
c
      allocate(ws(n),u(n,n),v(n,n))
      allocate(tabtmp1(n,n),tabtmp2(n,n))
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,n,xnodest,u,v,ws)
      if (ipoly.eq.1) call chebexps(itype,n,xnodest,u,v,ws)
      do i=1,n
         xnodest(i) = xnodest(i)*2
      enddo
      
      sigma = 16*delta/boxdim**2
      lambda = sqrt(sigma)
      reps=2d-16
      dmax=sqrt(log(1/reps)*sigma)
cccc      print *, 'dmax=',dmax
c
c     stob table 1, the scaled target interval is on [-5,-1], the scaled
c     source interval is on [-1,1]. Use the same quadrature nodes for all
c     target points
      if (lambda.le.0.125d0) then
         do j=1,n
            xi1 = xnodest(j) - 3.0d0
            dx = 2-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tabtmp1(m,j) = 0
               enddo
            else
               if (ipoly.eq.0) then
                  call mk_legetab_recur(n,xi1,lambda,
     1                boxdim/2,fint)
               elseif (ipoly.eq.1) then
                  call mk_chebtab_recur(n,xi1,lambda,
     1                boxdim/2,fint)
               endif
               do m=1,n
                  tabtmp1(m,j)  = fint(m)
               enddo
            endif
         enddo

         do j=1,n/2
            xi2 = xnodest(j) - 1.0d0
            dx=-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tabtmp2(m,j) = 0
               enddo
            else
               if (ipoly.eq.0) then
                  call mk_legetab_recur(n,xi2,lambda,
     1                boxdim/2,fint)
               elseif (ipoly.eq.1) then
                  call mk_chebtab_recur(n,xi2,lambda,
     1                boxdim/2,fint)
               endif
               do m=1,n
                  tabtmp2(m,j)  = fint(m)
               enddo
            endif
         enddo               
      else
         nquad = 50
         allocate(polyv(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c
         itype=1
         call legeexps(itype,nquad,xnodes,utmp,vtmp,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/4
         enddo

         if (ipoly.eq.0) then
            do i = 1,nquad
               call legepols(xnodes(i),n-1,polyv(1,i))
            enddo
         elseif (ipoly.eq.1) then
            do i = 1,nquad
               call chebpols(xnodes(i),n-1,polyv(1,i))
            enddo
         endif
      
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
                  tabtmp1(m,j) = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + polyv(m,i)*wexp(i,j)
                  enddo
                  tabtmp1(m,j) = rsum1
               enddo
            endif
         enddo
c     stob table 2, the scaled target interval is on [-3,1], the scaled source
c     interval is on [-1, 1].
c      
c     split the target into two parts. 1. for targets on [-3,-1], use previous 
c     quadrature
c      
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
                  tabtmp2(m,j) = 0
               enddo
            else
               do m = 1,n
                  rsum2 = 0.0d0
                  do i = 1,nquad
                     rsum2 = rsum2 + polyv(m,i)*wexp(i,j)
                  enddo
                  tabtmp2(m,j) = rsum2
               enddo
            endif
         enddo
      endif
      
c      
c      
c     Use different quadrature nodes for each target on [-1,1]
c     to ensure accuracy when dmax is small.
c
      if (lambda.le.0.125d0) then
         do j=n/2+1,n
            xi2 = xnodest(j) - 1.0d0
            if (ipoly.eq.0) then
               call mk_legetab_recur(n,xi2,lambda,boxdim/2,fint)
            elseif (ipoly.eq.1) then
               call mk_chebtab_recur(n,xi2,lambda,boxdim/2,fint)
            endif
            do m=1,n
               tabtmp2(m,j)  = fint(m)
            enddo
         enddo
      else
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
                  rsum1 = rsum1 + polyv(m,i)*wexp(i,j)
               enddo
               tabtmp2(m,j) = rsum1
            enddo
         enddo
      endif

c     multiply by u to obtain the tables acting on values
      do j=1,n
      do m=1,n
         dd=0
         do k=1,n
            dd=dd+tabtmp1(k,j)*u(k,m)
         enddo
         tab_stob(m,j,1)=dd
      enddo
      enddo
      
      do j=1,n
      do m=1,n
         dd=0
         do k=1,n
            dd=dd+tabtmp2(k,j)*u(k,m)
         enddo
         tab_stob(m,j,2)=dd
      enddo
      enddo
      
      
      
c     use symmetry to construct the tables 3 and 4
      do j=1,n
         do m=1,n
            tab_stob(m,j,3) = tab_stob(n-m+1,n-j+1,2)
         enddo
      enddo
      
      do j=1,n
         do m=1,n
            tab_stob(m,j,4) = tab_stob(n-m+1,n-j+1,1)
         enddo
      enddo

      do k=1,4
         call compute_sparse_pattern(eps,n,n,tab_stob(1,1,k),
     1       indc_stob(1,1,k),ifzero)
      enddo
      

      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine mk_loctab_btos(eps,ipoly,n,nnodes,delta,boxdim,
     1    tab_btos,indc_btos)
C*********************************************************************C
c     This routine is an optimized table generator.
c
c     tab_btos(n,j,k) = 
c              int_{source box} P_n(x) exp( -(\xi_j -x)^2/delta)
c              where boxdim is the box dimension of the TARGET BOX 
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
c     tab_btos
c----------------------------------------------------------------------c
      implicit none
      real *8 tab_btos(n,n,4)
      integer indc_btos(2,n+1,4)

      real *8 eps,delta,boxdim
      real *8 btosscale,utmp,vtmp
      real *8 btoscen(2)
      real *8, allocatable :: tnodes(:), ws(:) 
      real *8, allocatable :: xnodes(:), wts(:) 
      real *8, allocatable :: umat(:,:),vmat(:,:)
      real *8, allocatable :: ptmp(:)
      real *8, allocatable :: polyv(:,:)
      real *8, allocatable :: btos(:,:)
      real *8, allocatable :: btosx(:,:)
      real *8, allocatable :: btosxx(:,:)
      integer itype,nquad,i,j,k,m,ipoly,n,ifzero,nnodes
c
      allocate(btos(n,n))
      allocate(btosx(n,n))
      allocate(btosxx(n,n))

c     target nodes
      allocate(tnodes(n),ws(n),umat(n,n),vmat(n,n))
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,n,tnodes,umat,vmat,ws)
      if (ipoly.eq.1) call chebexps(itype,n,tnodes,umat,vmat,ws)

c     quadrature nodes and weights
      nquad = 50
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
c     scale is the ratio of the source box size to the target box size.
c     big source to small target 
      btosscale = 2.0d0
c     The scaled source interval is always [-1,1].
c     btos table 1, the scaled target interval is [-2,-1], centered at -1.5
c     btos table 2, the scaled target interval is [-1, 0], centered at -0.5
c     btos table 3, the scaled target interval is [0,1],   centered at  0.5
c     btos table 4, the scaled target interval is [1, 2],  centered at  1.5
      btoscen(1) = -1.5d0
      btoscen(2) = -0.5d0

      do i=1,2
         call mk_loctab(ipoly,n,delta,n,tnodes,boxdim,btoscen(i),
     1       btosscale,btos,btosx,btosxx,
     2       nquad,xnodes,wts,polyv)
c        multiply by u to obtain the tables acting on values
         call dgemm_f77('t','n',n,n,n,1.0d0,umat,n,btos,n,
     1       0.0d0,tab_btos(1,1,i),n)
      enddo
      
c     use symmetry to construct tables 3 and 4
      do j=1,n
         do m=1,n
            tab_btos(m,j,3) = tab_btos(n-m+1,n-j+1,2)
            tab_btos(m,j,4) = tab_btos(n-m+1,n-j+1,1)
         enddo
      enddo

      do k=1,4
         call compute_sparse_pattern(eps,n,n,tab_btos(1,1,k),
     1       indc_btos(1,1,k),ifzero)
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
      real *8 tnodes(n),fint(100),lambda
      real *8 xnodes(nquad),wts(nquad), polyv(nquad,n)
      integer i,j,m
      real *8 bs,bsh,reps,dmax,targ,dis,dx,sigma
      real *8 rsum
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
         do j=1,nt
            targ = targs(j)
c           distance between the target point and the source box
            dis=0
            if (targ .lt. -1) dis = -1-targ
            if (targ .gt.  1) dis = targ-1
            if (dis.gt.dmax) then
               do m=1,n
                  tab(m,j) = 0
               enddo
            else
               call mk_polytab_recur(ipoly,n,targ,lambda,
     1             bs,fint)
               do m=1,n
                  tab(m,j) = fint(m)
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
               enddo
            else
               do m = 1,n
                  rsum = 0.0d0
                  do i = 1,nquad
                     rsum = rsum + polyv(i,m)*wexp(i,j)
                  enddo
                  tab(m,j) = rsum*bsh
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
      subroutine mk_polytab_recur(ipoly,m,targ,lambda,boxdim,fint)
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
      real *8 fint(m),lambda,boxdim,targ

      if (ipoly.eq.0) call
     1    mk_legetab_recur(m,targ,lambda,boxdim,fint)
      
      if (ipoly.eq.1) call
     1    mk_chebtab_recur(m,targ,lambda,boxdim,fint)

      return
      end subroutine
c
c
c
c
      subroutine mk_legetab_recur(m,targ,lambda,boxdim,fint)
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

c     finally, multiply lambda back and also the proper weight adjustment
c     for the source box
      ww = boxdim*lambda/2
      do k=1,m
         fint(k)=fint(k)*ww
      enddo

      return
      end subroutine
c
c
c
      subroutine mk_chebtab_recur(m,targ,lambda,boxdim,fint)
      implicit real *8 (a-h,o-z)
cccc      sqrtpih=sqrt(4.0d0*atan(1.0d0))/2
      data sqrtpih/0.886226925452757940959713778283913d0/
c     calculate the values of the integral
c
c     \frac{1}{\lambda}\int_{-1}^1 T_n(x) e^{-(targ-x)^2/\lambda^2}dx
c
c     for n=0,1,...,m-1.
c      
c     Algorithm: use the five term recurrence formula in Shravan's CFGT paper.
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

c     finally, multiply lambda back and also the proper weight adjustment
c     for the source box
      ww = boxdim*lambda/2
      do k=1,m
         fint(k)=fint(k)*ww
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
      
cccc      call prinf('indc=*',indc,2*(n+1))
      
      return
      end subroutine
      
