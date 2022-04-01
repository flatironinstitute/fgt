c     This file contains a set of subroutines that build various 1D tables
c     used in the box FGT in all dimensions.
c
c
c     mk_cheb2pw: builds the table converting Chebyshev expansion coefficients
c                to planewave expansion coefficients
c
c     mk_pw2pot: builds the table converting planewave expansion coefficients
c                to potential values on Chebyshev nodes
c
c     mk_loctab_coll: builds the table converting Chebyshev expansion coefficients
c                in the source box to potential values on Chebyshev nodes in the
c                target box at the same level, three of them for each level.
c
c     mk_loctab_stob: builds the table converting Chebyshev expansion coefficients
c                in a small source box to potential values on Chebyshev nodes in a
c                large target box at the coarse level. four of them for each level.
c
c     mk_loctab_btos: builds the table converting Chebyshev expansion coefficients
c                in a large source box to potential values on Chebyshev nodes in a
c                small target box at the fine level. four of them for each level.
c
C*********************************************************************C
      subroutine mk_cheb2pw_old(n,npw,nnodes,ws,ts,delta,boxdim,
     1    tab_cheb2pw)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_cheb2pw(n,j) = ws(j)*(D/2) * 
c              int_{-1}^1 P_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c              where D is the box dimension at current level in
c              tree hierarchy.
c
c     That is, 
c     tab_cheb2pw(n,j) is the Fourier transform of P_n at a specific
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
c     tab_cheb2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 P_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 u,v
      real *8, allocatable :: chebv(:,:)
      complex *16 tab_cheb2pw(n,npw),zsum,eye
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: whts(:), xnodes(:)
      complex *16, allocatable :: wexp(:,:)
      
c
      allocate(wexp(nnodes,(npw+1)/2))
      allocate(chebv(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
c
      eye = dcmplx(0.0d0,1.0d0)
      itype = 1
      call legeexps(itype,nnodes,xnodes,u,v,whts)
c
      do i = 1,nnodes
         call chebpols(xnodes(i),n-1,chebv(1,i))
      enddo

      fac = boxdim/(2*dsqrt(delta))
      do j = 1,(npw+1)/2
cccc         print *, abs(ts(j)*fac)
         do i=1,nnodes
            wexp(i,j)=whts(i)*cdexp(-eye*ts(j)*fac*xnodes(i))
         enddo
      enddo
ccc      call prin2(' chebv is *',chebv,n*nnodes)
c
      do j = 1,(npw+1)/2
         dd = ws(j)*boxdim/2.0d0
         do m = 1,n
            zsum = 0.0d0
            do i = 1,nnodes
               zsum = zsum + chebv(m,i)*wexp(i,j)
            enddo
            tab_cheb2pw(m,j) = dd*zsum
            tab_cheb2pw(m,npw-j+1) = conjg(tab_cheb2pw(m,j))
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
      subroutine mk_cheb2pw(n,npw,nnodes,ws,ts,delta,boxdim,
     1    tab_cheb2pw)
C*********************************************************************C
c     Use half-order Bessel J function to construct the table.
c
c     tab_cheb2pw(n,j) = ws(j)*(D/2) * 
c              int_{-1}^1 P_n(x) exp(- i ts(j)Dx/(2 \sqrt{delta})) dx
c              where D is the box dimension at current level in
c              tree hierarchy.
c
c     That is, 
c     tab_cheb2pw(n,j) is the Fourier transform of P_n at a specific
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
c     tab_cheb2pw  (n,j) entry is:  ws(j) * int_{-1}^{1}  
c                 P_n(x) exp(- i ts(j) D x/(2 \sqrt(delta)  dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 tab_cheb2pw(n,npw),zsum,eye,eyem,cd
      real *8 ws(npw),ts(npw)
      real *8, allocatable :: cyr(:),cyi(:),chebv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      complex *16, allocatable :: legv(:)
      real *8, allocatable :: xq(:),wq(:),u(:,:),v(:,:)
      complex *16, allocatable :: tabtemp(:,:)

      allocate(tabtemp(n,npw))
      allocate(xq(n),wq(n),u(n,n),v(n,n))
      
c
      allocate(chebv(n,n))
      allocate(whts(n))
      allocate(xnodes(n))
      allocate(legv(n))
      allocate(cyr(n))
      allocate(cyi(n))
      
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
      
      eye = dcmplx(0.0d0,1.0d0)
      fac = boxdim/(2*dsqrt(delta))

      fnu=0.5d0
      kode=1
      do j = npw/2+1,npw
         dd = ws(j)*boxdim/2.0d0
         
         if (abs(ts(j)).lt.1d-12) then
            do m=0,n-1,2
               tabtemp(m+1,j) = dd*2/(1-m*m)
            enddo
            do m=1,n-1,2
               tabtemp(m+1,j) = 0
            enddo
         else
            zi=0.0d0
            zr = ts(j)*fac
            ddd = sqrt(2*pi/zr)
            call zbesj(zr,zi,fnu,kode,n,cyr,cyi,nz,ierr)
            eyem=1.0d0
c           legv contains the Fourier transform of the Legendre polynomials
            do m=1,n
               legv(m)=dd*dcmplx(cyr(m),cyi(m))*ddd/eyem
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
c     norder   number of Chebyshev nodes
c     npw      number of plane waves
c     ts       nodes of plane wave quadrature
c     xs       Chebyshev nodes
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
cccc         
cccc         qq1 = qqx
cccc         qqx = qqx*qqx
cccc
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
C*********************************************************************C
      subroutine mk_loctab_coll_old(n,nnodes,delta,boxdim,tab_colleague)
C*********************************************************************C
c     This routine is a correct but not optimized table generator.
c
c     tab_colleague(n,j,k) = 
c              int_{-D/2}^{D/2} P_n(D x/2) exp( -(\xi_j -x)^2/delta)
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
      real *8, allocatable :: chebv(:,:)
      real *8 tab_colleague(n,n,-1:1)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: xnodest(:),wexp(:,:)
      real *8, allocatable :: wexpp1(:,:),wexpm1(:,:)
c
      allocate(chebv(n,nnodes))
      allocate(whts(nnodes))
      allocate(xnodes(nnodes))
      allocate(xnodest(n))
      allocate(wexp(nnodes,n))
      allocate(wexpp1(nnodes,n))
      allocate(wexpm1(nnodes,n))
c
      itype = 1
      call chebexps(itype,n,xnodest,u,v,whts)
      do i=1,n
         xnodest(i) = boxdim*xnodest(i)/2.0d0
      enddo
ccc      call prin2(' Chebyshev nodes are *',xnodest,n)
      call legeexps(itype,nnodes,xnodes,u,v,whts)
ccc      call prin2(' Chebyshev nodes are *',xnodes,nnodes)
c
      do i = 1,nnodes
         call chebpols(xnodes(i),n-1,chebv(1,i))
      enddo
      do i=1,nnodes
         xnodes(i) = boxdim*xnodes(i)/2.0d0
         whts(i) = whts(i)*boxdim/2.0d0
      enddo
cccc      call prin2('mk_loctab_coll, D^2/4delta=*',(boxdim/2)**2/delta,1)

      do j=1,n
         xi = xnodest(j)
         xim1 = xnodest(j) - boxdim
         xip1 = xnodest(j) + boxdim
         do i=1,nnodes
            dx=xi-xnodes(i)
            wexp(i,j)=exp(-dx*dx/delta)*whts(i)
               
            dx = xim1 - xnodes(i)
            wexpm1(i,j)=exp(-dx*dx/delta)*whts(i)
            
            dx = xip1 - xnodes(i)
            wexpp1(i,j)=exp(-dx*dx/delta)*whts(i)
         enddo
      enddo
            
ccc      call prin2(' chebv is *',chebv,n*nnodes)
c
      do j = 1,n
      do m = 1,n
         rsum = 0.0d0
         rsumm1 = 0.0d0
         rsump1 = 0.0d0
         do i = 1,nnodes
            rsum = rsum + chebv(m,i)*wexp(i,j)
            rsumm1 = rsumm1 + chebv(m,i)*wexpm1(i,j)
            rsump1 = rsump1 + chebv(m,i)*wexpp1(i,j)
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
      subroutine mk_loctab_coll(eps,n,nnodes,delta,boxdim,
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
      real *8, allocatable :: chebv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
      real *8, allocatable :: ws(:),u(:,:),v(:,:)
      real *8, allocatable :: tabtmp0(:,:),tabtmpm1(:,:)
c
      allocate(ws(n),u(n,n),v(n,n))
      allocate(tabtmp0(n,n),tabtmpm1(n,n))
c
      itype = 2
      call chebexps(itype,n,xnodest,u,v,ws)

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
               call mk_loctab_recurrence(n,xim1,lambda,dmax,boxdim,fint)
               do m=1,n
                  tabtmpm1(m,j)  = fint(m)
               enddo
            endif
         enddo
      else
         nquad = 50
         allocate(chebv(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/2
         enddo
      
         do i = 1,nquad
            call chebpols(xnodes(i),n-1,chebv(1,i))
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
                     rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
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
            call mk_loctab_recurrence(n,xi,lambda,dmax,boxdim,fint)
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
                  rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
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
      subroutine mk_loctab_stob(eps,n,nnodes,delta,boxdim,
     1    tab_stob,indc_stob)
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
      real *8, allocatable :: chebv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
      real *8, allocatable :: ws(:),u(:,:),v(:,:)
      real *8, allocatable :: tabtmp1(:,:),tabtmp2(:,:)
c
      allocate(ws(n),u(n,n),v(n,n))
      allocate(tabtmp1(n,n),tabtmp2(n,n))
c
      itype = 2
      call chebexps(itype,n,xnodest,u,v,ws)
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
               call mk_loctab_recurrence(n,xi1,lambda,dmax,
     1             boxdim/2,fint)
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
               call mk_loctab_recurrence(n,xi2,lambda,dmax,
     1             boxdim/2,fint)
               do m=1,n
                  tabtmp2(m,j)  = fint(m)
               enddo
            endif
         enddo               
      else
         nquad = 50
         allocate(chebv(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
c
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)

         do i=1,nquad
            whts(i) = whts(i)*boxdim/4
         enddo
      
         do i = 1,nquad
            call chebpols(xnodes(i),n-1,chebv(1,i))
         enddo
      
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
                     rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
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
                     rsum2 = rsum2 + chebv(m,i)*wexp(i,j)
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
            call mk_loctab_recurrence(n,xi2,lambda,dmax,boxdim/2,fint)
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
                  rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
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
      subroutine mk_loctab_btos(eps,n,nnodes,delta,boxdim,
     1    tab_btos,indc_btos)
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
      integer indc_btos(2,n+1,4)
      real *8 xnodest(100),fint(100),lambda
      real *8, allocatable :: chebv(:,:)
      real *8, allocatable :: whts(:), xnodes(:)
      real *8, allocatable :: wexp(:,:)
      real *8, allocatable :: ws(:),u(:,:),v(:,:)
      real *8, allocatable :: tabtmp1(:,:),tabtmp2(:,:)
c
      allocate(ws(n),u(n,n),v(n,n))
      allocate(tabtmp1(n,n),tabtmp2(n,n))
c
      itype = 2
      call chebexps(itype,n,xnodest,u,v,ws)
      do i=1,n
         xnodest(i) = xnodest(i)/2
      enddo

      sigma = delta/boxdim**2
      lambda = sqrt(sigma)
cccc      print *, 'lambda=', lambda
      
      reps=2d-16
      dmax=sqrt(log(1/reps)*sigma)
cccc      print *, 'dmax=',dmax
c
c     btos table 1, the scaled target interval is on [-2,-1], the scaled
c     source interval is on [-1,1]. Use the same quadrature nodes for all
c     target points
c
      if (lambda.le.0.125d0) then
         do j=1,n
            dx=0.5d0-xnodest(j)
            if (dx.gt.dmax) then
               do m=1,n
                  tabtmp1(m,j) = 0
               enddo
            else
               xi1 = xnodest(j) - 1.5d0
               call mk_loctab_recurrence(n,xi1,lambda,dmax,
     1             boxdim*2,fint)
               do m=1,n
                  tabtmp1(m,j)  = fint(m)
               enddo
            endif
         enddo
      else
         nquad = 50
         allocate(chebv(n,nquad))
         allocate(whts(nquad))
         allocate(xnodes(nquad))
         allocate(wexp(nquad,n))
         itype=1
         call legeexps(itype,nquad,xnodes,u,v,whts)
         do i=1,nquad
            whts(i) = whts(i)*boxdim
         enddo
         
         do i = 1,nquad
            call chebpols(xnodes(i),n-1,chebv(1,i))
         enddo
         
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
                  tabtmp1(m,j)  = 0
               enddo
            else
               do m = 1,n
                  rsum1 = 0.0d0
                  do i = 1,nquad
                     rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
                  enddo
                  tabtmp1(m,j) = rsum1
               enddo
            endif
         enddo
      endif
      
c     btos table 2, the scaled target interval is on [-1,0], the scaled source
c     interval is on [-1, 1]. Use different quadrature nodes for each target
c     to ensure accuracy when dmax is small.

      if (lambda.le.0.125d0) then
         do j=1,n
            xi2 = xnodest(j) - 0.5d0         
            call mk_loctab_recurrence(n,xi2,lambda,dmax,boxdim*2,fint)
            do m=1,n
               tabtmp2(m,j)  = fint(m)
            enddo
         enddo
      else
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
                  rsum1 = rsum1 + chebv(m,i)*wexp(i,j)
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
         tab_btos(m,j,1)=dd
      enddo
      enddo
      
      do j=1,n
      do m=1,n
         dd=0
         do k=1,n
            dd=dd+tabtmp2(k,j)*u(k,m)
         enddo
         tab_btos(m,j,2)=dd
      enddo
      enddo
      
      
      
c     use symmetry to construct tables 3 and 4
      do j=1,n
         do m=1,n
            tab_btos(m,j,3) = tab_btos(n-m+1,n-j+1,2)
         enddo
      enddo
      
      do j=1,n
         do m=1,n
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
c
      subroutine mk_loctab_1col(n,targ,sigma,dmax,boxdim,fint)
      implicit real *8 (a-h,o-z)
c     calculate the values of the integral
c
c     \int_{-1}^1 T_m(x) e^{-(targ-x)^2/sigma}dx
c
c     for m=0,1,...,n-1.
c     Algorithm: 1. find the effective range [a,b] of the source interval;
c     2. if b<a, return 0 for all integrals;
c     3. if b>a, use Gauss-Legendre quadrature to calculate the integral.
c
c
c
      real *8 fint(n)
      real *8 xs(100),ws(100),wexp(100)
      real *8, allocatable :: chebv(:,:)

      a=targ-dmax
      if (a.lt.-1) a=-1

      b=targ+dmax
      if (b.gt.1) b=1

      if (b.le.a) then
         do i=1,n
            fint(i)=0
         enddo
         return
      endif

      if (b-a.lt.0.2d0) then
         nquad = 40
      elseif (b-a.lt.1.0d0) then
         nquad = 40
      else
         nquad = 44
      endif

      allocate(chebv(n,nquad))
      
      itype=1
      call legeexps(itype,nquad,xs,u,v,ws)

      do i=1,nquad
         xs(i)= xs(i)*(b-a)/2+(b+a)/2
         ws(i)= ws(i)*(b-a)/2
      enddo

      do i=1,nquad
         call chebpols(xs(i),n-1,chebv(1,i))
      enddo

      do i=1,nquad
         wexp(i)=exp(-(targ-xs(i))**2/sigma)*ws(i)*boxdim/2
      enddo

      do m=1,n
         dd=0
         do i=1,nquad
            dd=dd+chebv(m,i)*wexp(i)
         enddo
         fint(m)=dd
      enddo

      return
      end
c
c
c
      subroutine mk_loctab_recurrence(m,targ,lambda,dmax,boxdim,fint)
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
      end
c
c
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
      end
