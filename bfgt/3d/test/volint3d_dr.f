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
      real *8 targ(3),center(3)
      real *8 fint(12000)
      real *8 texp(12000),tdiff(12000),tj,xtes,sum,val,fx,fpx,fa,fb
      real *8 delta,whts(12000),work(10000),y2(12000)
      real *8, allocatable :: f(:,:,:)
      real *8, allocatable :: coeff(:,:,:)
      real *8  ws(100),ts(100)
      complex *16, allocatable :: tab_leg2pw(:,:)
      complex *16, allocatable :: ff(:,:,:)
      complex *16, allocatable :: ff2(:,:,:)
      complex *16, allocatable :: pwexp(:,:,:)
c
c----
c
      call prini(6,13)
      pi = 4*datan(1.0D0)
C
c---- set interval of calculation [a,b].
c
      a = 0.3d0
      b = 1.7d0
      a = 0.3d0
      b = 1.0d0
c
c---- set order of approximation
c
      print *, ' enter order of approximation '
      read *, n
      print *, ' number of nodes is ',n
      allocate(f(n,n,n))
      allocate(coeff(n,n,n))
c
c---- compute N Legendre nodes :
c
      itype =  2
      call legeexps(itype,n,y,umat,vmat,whts)
      call prin2(' Legendre nodes are *',y,n)
      do i=1,n
         whts(i) = whts(i)*(b-a)/2.0d0
      enddo
c
C     compute a function F
C
      center(1) = (b+a)/2.0d0
      center(2) = (b+a)/2.0d0
      center(3) = (b+a)/2.0d0
      targ(1) = 2.0d0
      targ(2) = 2.1d0
      targ(3) = 2.3d0
      targ(1) = center(1)
      targ(2) = center(2)
      targ(3) = center(3)
      delta = 8d-1
      rsum = 0.0d0
      do i=1,n
      do j=1,n
      do k=1,n
        xx = (b-a)*y(i)/2.0d0 + center(1)
        yy = (b-a)*y(j)/2.0d0 + center(2)
        zz = (b-a)*y(k)/2.0d0 + center(3)
	call rfun(a,b,xx,yy,zz,fx)
	f(i,j,k)= fx
        dx = targ(1) - xx
        dy = targ(2) - yy
        dz = targ(3) - zz
        rsum = rsum + whts(i)*whts(j)*whts(k)*f(i,j,k)*
     1         exp(-(dx*dx+dy*dy+dz*dz)/delta)
      enddo
      enddo
      enddo
ccc      call prin2(' function values are *',f,n*n)
      call prin2(' rsum = *',rsum,1)
c
c     compute forward and inverse transform as test 
c
      nd = 1
ccc      call legtrans3d(nd,n,f,w1,w2,umat,coeff)
      call legval2coefs2_3d(nd,n,f,coeff,umat)
ccc      call prin2(' Legendre expansion for f is *',coeff,n*n)
      pmax = 12.0d0
      npw = 40
      allocate(tab_leg2pw(n,npw))
      allocate(ff(npw,n,n))
      allocate(ff2(npw,npw,n))
      allocate(pwexp(npw,npw,npw/2))

      boxdim = b-a
      ntarg = 1
      call get_pwnodes(pmax,npw,ws,ts)
      call prin2(' ws is *',ws,npw)
      call prin2(' ts is *',ts,npw)
      nnodes = 100
      call mk_leg2pw(n,npw,nnodes,ws,ts,delta,boxdim,tab_leg2pw)
cccc      call prin2(' tab_leg2pw is *',tab_leg2pw,n*npw*2)
      t0 = second()
      do ii = 1,1
         call leg3d_to_pw(nd,n,coeff,npw,ff,ff2,tab_leg2pw,pwexp)
      enddo
      t1 = second()
      call prin2(' time for n*n*n*1000 is *',t1-t0,1)
cccc      call prin2(' pwexp is *',pwexp,npw*npw*npw)
      call g3dpwevalp_vec(nd,delta,center,npw,ws,ts,
     1              pwexp,targ,ntarg,pot)
      call prin2(' pot = *',pot,1)
      call prin2(' error = *',abs(pot-rsum)/abs(rsum),1)
c
c
      stop
      end
c
      subroutine rfun(a,b,x,y,z,fval)
      implicit real *8 (a-h,o-z)
      real *8 x,y,temp,fval
      integer *4 n,i
C
c---- compute a function fval(x,y)
c
      rm = 0.51d0
      fval = dcos(x)*dsin(y)*dcos(dsin(z))
ccc      fval = 1.0d0
ccc      fval = x
      return
      end
c
      subroutine get_pwnodes0(pmax,npw,ws,ts)
C
C     Get planewave exp weights,nodes
C
      implicit real *8 (a-h,o-z)
      real *8 ws(npw),ts(npw)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      h = pmax/npw2
      w = h/(2.0d0*dsqrt(pi))

      do j =1,npw2
         ts(j) = (j-0.5d0)*h
         ws(j) = w*dexp(-ts(j)*ts(j)/4)
         ts(j+npw2)=-ts(j)
         ws(j+npw2)=ws(j)
      enddo
c
      return
      end
C

