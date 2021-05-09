c
      implicit real *8 (a-h,o-z)
      integer *4 n,n1,i,j,nm1
      real *8 umat(100*100),vmat(100*100)
      real *8 endinter(100)
      real *8 a,b,y(12000),u,v,u1,v1,thresh,pi,h,ealg,salg
      real *8 fprime(12000),f2(12000),relerr(12000),dint,rint,stot
      real *8 fdiff(12000)
      real *8 w1(12000)
      real *8 w2(12000)
      real *8 targ(2),center(2),center2(2),centers(2),centerb(2)
      real *8 fint(12000)
      real *8 texp(12000),tdiff(12000),tj,xtes,sum,val,fx,fpx,fa,fb
      real *8 delta,whts(12000),work(10000),y2(12000)
      real *8, allocatable :: f(:,:)
      real *8, allocatable :: fsmall(:,:)
      real *8, allocatable :: fbig(:,:)
      real *8, allocatable :: coeff(:,:)
      real *8, allocatable :: coeffs(:,:)
      real *8, allocatable :: coeffb(:,:)
      real *8, allocatable :: tab_colleague(:,:,:)
      real *8, allocatable :: tab_stob(:,:,:)
      real *8, allocatable :: tab_btos(:,:,:)
      real *8, allocatable :: ff(:,:)
      real *8, allocatable :: potdir(:,:)
      real *8, allocatable :: pot(:,:)
c
c----
c
      call prini(6,13)
      pi = 4*datan(1.0D0)
C
c---- set interval of calculation [a,b].
c
      a = 0.3d0
      b = 1.0d0
c
c---- set order of approximation
c
      print *, ' enter order of approximation '
      read *, n
      print *, ' number of nodes is ',n
      allocate(f(n,n))
      allocate(fsmall(n,n))
      allocate(fbig(n,n))
      allocate(ff(n,n))
      allocate(tab_colleague(n,n,-1:1))
      allocate(tab_stob(n,n,4))
      allocate(tab_btos(n,n,4))
      allocate(coeff(n,n))
      allocate(coeffs(n,n))
      allocate(coeffb(n,n))
      allocate(potdir(n,n))
      allocate(pot(n,n))
c
c---- compute N Legendre nodes :
c
      itype =  2
      call legeexps(itype,n,y,umat,vmat,whts)
      call prin2(' Legendre nodes are *',y,n)
c
C     compute a function F
C
      boxdim = b-a
c
c     center is center of targ box
c     center2 is center of colleague source box
c     centers is center of small adjacent source box
c
      center(1) = (b+a)/2.0d0
      center(2) = (b+a)/2.0d0
      center2(1) = center(1) - boxdim
      center2(2) = center(2) - boxdim
      centers(1) = center(1) - boxdim/4
      centers(2) = center(2) - 3*boxdim/4
      centerb(1) = center(1) - boxdim/2
      centerb(2) = center(2) - 3*boxdim/2
      delta = 0.8d0
      do i=1,n
      do j=1,n
        xx = (b-a)*y(i)/2.0d0 + center2(1)
        yy = (b-a)*y(j)/2.0d0 + center2(2)
	call rfun(a,b,xx,yy,fx)
	f(i,j)= fx
        xx = (b-a)*y(i)/4.0d0 + centers(1)
        yy = (b-a)*y(j)/4.0d0 + centers(2)
	call rfun(a,b,xx,yy,fx)
	fsmall(i,j)= fx
        xx = (b-a)*y(i)/1.0d0 + centerb(1)
        yy = (b-a)*y(j)/1.0d0 + centerb(2)
	call rfun(a,b,xx,yy,fx)
	fbig(i,j)= fx
      enddo
      enddo
      call prin2(' function values are *',f,n*n)
      call prin2(' function values small are *',fsmall,n*n)
c
c     compute forward and inverse transform as test 
c
      nd = 1
      call legtrans2d(nd,n,f,w1,w2,umat,coeff)
      call legtrans2d(nd,n,fsmall,w1,w2,umat,coeffs)
      call legtrans2d(nd,n,fbig,w1,w2,umat,coeffb)
ccc      call prin2(' Legendre expansion for f is *',coeff,n*n)

      nnodes = 100
      itype =  1
      call legeexps(itype,nnodes,y2,umat,vmat,whts)
      do i=1,nnodes
         whts(i) = whts(i)*(b-a)/2.0d0
      enddo
c
      call mk_loctab_coll(n,nnodes,delta,boxdim,tab_colleague)
ccc      call prin2(' tab_colleagues is *',tab_colleague,n*n*3)
      t0 = second()
      do ii = 1,1
         call leg2d_to_potloc(nd,n,coeff,ff,pot,
     1        tab_colleague(1,1,1), 
     1        tab_colleague(1,1,1))
      enddo
      t1 = second()
      call prin2(' estimate for n*n*10000 is *',9*(t1-t0),1)
      call prin2(' pot = *',pot,n*n)
c
      etot = 0
      stot = 0
      do ii = 1,n
      do jj = 1,n
         targ(1) = (b-a)*y(ii)/2.0d0 + center(1)
         targ(2) = (b-a)*y(jj)/2.0d0 + center(2)
         potdir(ii,jj) = 0.0d0
         do i=1,nnodes
         do j=1,nnodes
           xx = (b-a)*y2(i)/2.0d0 + center2(1)
           yy = (b-a)*y2(j)/2.0d0 + center2(2)
	   call rfun(a,b,xx,yy,fx)
           dx = targ(1) - xx
           dy = targ(2) - yy
           potdir(ii,jj) = potdir(ii,jj)+whts(i)*whts(j)*fx*
     1         exp(-(dx*dx+dy*dy)/delta)
         enddo
         enddo
         etot = etot + (potdir(ii,jj) - pot(ii,jj))**2
         stot = stot + potdir(ii,jj)**2
      enddo
      enddo
      call prin2(' potdir = *',potdir,n*n)
      call prin2(' error = *',sqrt(etot/stot),1)
c
c
c
      nnodes = 100
      itype =  1
      call legeexps(itype,nnodes,y2,umat,vmat,whts)
      do i=1,nnodes
         whts(i) = whts(i)*(b-a)/4.0d0
      enddo
c
      call mk_loctab_stob(n,nnodes,delta,boxdim,tab_stob)
      t0 = second()
      do ii = 1,1
         call leg2d_to_potloc(nd,n,coeffs,ff,pot,
     1        tab_stob(1,1,3), 
     1        tab_stob(1,1,4))
      enddo
      t1 = second()
      call prin2(' estimate for n*n*10000 is *',9*(t1-t0),1)
      call prin2(' pot = *',pot,n*n)
c
      etot = 0
      stot = 0
      do ii = 1,n
      do jj = 1,n
         targ(1) = (b-a)*y(ii)/2.0d0 + center(1)
         targ(2) = (b-a)*y(jj)/2.0d0 + center(2)
         potdir(ii,jj) = 0.0d0
         do i=1,nnodes
         do j=1,nnodes
           xx = (b-a)*y2(i)/4.0d0 + centers(1)
           yy = (b-a)*y2(j)/4.0d0 + centers(2)
	   call rfun(a,b,xx,yy,fx)
           dx = targ(1) - xx
           dy = targ(2) - yy
           potdir(ii,jj) = potdir(ii,jj)+whts(i)*whts(j)*fx*
     1         exp(-(dx*dx+dy*dy)/delta)
         enddo
         enddo
         etot = etot + (potdir(ii,jj) - pot(ii,jj))**2
         stot = stot + potdir(ii,jj)**2
      enddo
      enddo
      call prin2(' potdir = *',potdir,n*n)
      call prin2(' error = *',sqrt(etot/stot),1)
c
c
      nnodes = 100
      itype =  1
      call legeexps(itype,nnodes,y2,umat,vmat,whts)
      do i=1,nnodes
         whts(i) = whts(i)*(b-a)
      enddo
c
      call mk_loctab_btos(n,nnodes,delta,boxdim,tab_btos)
      t0 = second()
      do ii = 1,1
         call leg2d_to_potloc(nd,n,coeffb,ff,pot,
     1        tab_btos(1,1,3), 
     1        tab_btos(1,1,4))
      enddo
      t1 = second()
      call prin2(' estimate for n*n*10000 is *',9*(t1-t0),1)
      call prin2(' pot = *',pot,n*n)
c
      etot = 0
      stot = 0
      do ii = 1,n
      do jj = 1,n
         targ(1) = (b-a)*y(ii)/2.0d0 + center(1)
         targ(2) = (b-a)*y(jj)/2.0d0 + center(2)
         potdir(ii,jj) = 0.0d0
         do i=1,nnodes
         do j=1,nnodes
           xx = (b-a)*y2(i)/1.0d0 + centerb(1)
           yy = (b-a)*y2(j)/1.0d0 + centerb(2)
	   call rfun(a,b,xx,yy,fx)
           dx = targ(1) - xx
           dy = targ(2) - yy
           potdir(ii,jj) = potdir(ii,jj)+whts(i)*whts(j)*fx*
     1         exp(-(dx*dx+dy*dy)/delta)
         enddo
         enddo
         etot = etot + (potdir(ii,jj) - pot(ii,jj))**2
         stot = stot + potdir(ii,jj)**2
      enddo
      enddo
      call prin2(' potdir = *',potdir,n*n)
      call prin2(' error = *',sqrt(etot/stot),1)
c

c
      stop
      end
c
      subroutine rfun(a,b,x,y,fval)
      implicit real *8 (a-h,o-z)
      real *8 x,y,temp,fval
      integer *4 n,i
C
c---- compute a function fval(x,y)
c
      rm = 0.51d0
      fval = dcos(x)*dsin(y)
ccc      fval = 1.0d0
ccc      fval = x
      return
      end
c
