c
c
c
c
c
      implicit real *8 (a-h,o-z)
      real *8 epsvals(5),deltas(20),pps(20,5)
      integer nsrcs(20)
      integer dim
c	  
      call prini(6,13)

      pi = 4*atan(1.0d0)
      
      epsvals(1) = 1d-3
      epsvals(2) = 1d-6
      epsvals(3) = 1d-9
      epsvals(4) = 1d-12
c
      do i=1,10
         deltas(i)=10.0d0**(-i)
      enddo

      do i=1,5
         nsrcs(i)=10**5*2**i
      enddo
      do i=6,10
         nsrcs(i)=nsrcs(5)
      enddo
      do i=1,10
         nsrcs(i)=10**6
      enddo
c     nd - number of different densities
      nd=1
c     dim - dimension of the problem
      dim=3
c     ifpgh = 1: potential; 2: pot+grad; 3: pot+grad+hess
      ifpgh=2
c     ifpghtarg: flag for arbitrary targets
      ifpghtarg=3
c     ntarg: number of extra targets
      ntarg=0
c     nsrc: number of source points

c     
c
c     test all parameters
c
      iw = 70
      iw2 = 80

      ifuniform=0
      if (ifuniform.eq.1) then
         open(iw, file='error3du.txt', position='append')
         open(iw2, file='timing3du.txt', position='append')
      else
         open(iw, file='error3dn.txt', position='append')
         open(iw2, file='timing3dn.txt', position='append')
      endif
      
      do i=2,2
         eps = epsvals(i)
         do j=1,4
            delta=deltas(j)
            nsrc=nsrcs(j)
            ntarg=nsrc
c
 2000       format(2x, i8, 2x)
c     write(iw,2000,advance="no") nint(1/ratio)

            call testfgt3d(nd,dim,eps,delta,nsrc,ntarg,ifpgh,ifpghtarg,
     1          ifuniform,pps(j,i),rerr)

 4800       format(2x,'&',2x,D8.2,1x,'&',2x,D8.2,1x,'&',2x,D8.2,1x,'\\')
            write(iw,4800) eps, delta, rerr
c     
c            if (isnan(rerr) .or. rerr.gt.eps) then               
c               call prin2('eps=*',eps,1)
c               call prin2('delta=*',delta,1)
c            endif
c     
         enddo  
      enddo

      do j=1,4
         write(iw2,*) deltas(j), pps(j,1), pps(j,2), pps(j,3)
      enddo
      
      close(iw)
      close(iw2)

      end
c
c
c
c
      subroutine testfgt3d(nd,dim,eps,delta,nsrc,ntarg,ifpgh,ifpghtarg,
     1    ifuniform,pps,errps)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: rnormal(:,:)
      real *8, allocatable :: charges(:,:),dipstr(:,:)
      real *8, allocatable :: pot(:,:),grad(:,:,:),hess(:,:,:)
      real *8, allocatable :: pottarg(:,:),gradtarg(:,:,:),
     1    hesstarg(:,:,:)
      real *8, allocatable :: potex(:,:),gradex(:,:,:),hessex(:,:,:)
      real *8, allocatable :: pottargex(:,:),gradtargex(:,:,:),
     1                             hesstargex(:,:,:)

      integer dim
      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      call prin2(' delta = *',delta,1)
c      call prinf_long(' nsrc = *',nsrc,1)
c
      allocate(sources(dim,nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(rnormal(dim,nsrc))
      allocate(targ(dim,ntarg))
      nhess=dim*(dim+1)/2
      allocate(pot(nd,nsrc),grad(nd,dim,nsrc),hess(nd,nhess,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,nhess,ntarg))
      allocate(hesstarg(nd,nhess,ntarg))

      rin = 0.3d0
      rwig = 0.12d0
      nwig = 6
      
      do i=1,nsrc
         if (ifuniform.eq.0) then
c        nonuniform source distribution
            theta = hkrand(0)*pi
            rr=rin+rwig*cos(nwig*theta)
            ct=cos(theta)
            st=sin(theta)

            phi = hkrand(0)*2*pi
            cp=cos(phi)
            sp=sin(phi)

            if (dim.eq.3) then
               sources(1,i) = rr*st*cp+0.5d0
               sources(2,i) = rr*st*sp+0.5d0
               sources(3,i) = rr*ct+0.5d0
            elseif (dim.eq.2) then
               sources(1,i) = rr*cp+0.5d0
               sources(2,i) = rr*sp+0.5d0
            endif
         else
c        uniform source distribution
            do k=1,dim
               sources(k,i) = hkrand(0)
            enddo
         endif

         do k=1,dim
            rnormal(k,i) = hkrand(0)
         enddo
         do ind = 1,nd
            charges(ind,i) = hkrand(0) 
            dipstr(ind,i) = hkrand(0)
         enddo
      enddo

      do ind=1,nd
         dipstr(ind,1) = 0.0d0 
         dipstr(ind,2) = 0.0d0 
         charges(ind,1) = 0.0d0 
         charges(ind,2) = 0.0d0 
         charges(ind,3) = 1.0d0
      enddo

      do k=1,dim
         sources(k,1) = 0.0d0
      enddo

      do k=1,dim
         sources(k,2) = 1.0d0
      enddo

      do k=1,dim
         sources(k,3) = 0.05d0
      enddo
 
      do i=1,ntarg
         do k=1,dim
            targ(k,i) = hkrand(0)
         enddo
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nd,nts),gradex(nd,dim,nts),hessex(nd,nhess,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,dim,ntt))
      allocate(hesstargex(nd,nhess,ntt))
c
      call dzero(potex,nts*nd)
      call dzero(gradex,dim*nts*nd)
      call dzero(hessex,nhess*nts*nd)
      
      call dzero(pottargex,ntt*nd)
      call dzero(gradtargex,dim*ntt*nd)
      call dzero(hesstargex,nhess*ntt*nd)
      
      ifcharge = 0
      ifdipole = 1
      call prinf(' ifcharge is *',ifcharge,1)
      call prinf(' ifdipole is *',ifdipole,1)
      call prinf(' ifpgh is *',ifpgh,1)
      call prinf(' ifpghtarg is *',ifpghtarg,1)

      iper = 0
      call cpu_time(t1)
      call fgt3d(nd,delta,eps,nsrc,sources,ifcharge,charges,
     1    ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2    ntarg,targ,ifpghtarg,pottarg,gradtarg,
     3    hesstarg)
      call cpu_time(t2)
      pps=(nsrc*ifpgh+ntarg*ifpghtarg+0.0d0)/(t2-t1)
      call prin2('points per sec=*',pps,1)

cccc      call prin2('pot=*',pot,nd*nts)
cccc      call prin2('pottarg=*',pottarg,nd*ntt)
cccc      call prin2('gradtarg=*',gradtarg,nd*3*ntt)
cccc      call prin2('hesstarg=*',hesstarg,nd*6*ntt)

      dmax = log(1.0d2/eps)*delta
      call cpu_time(t1)
      call fgtpart_direct(nd,dim,delta,dmax,1,nsrc,1,nts,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    sources,ifpgh,potex,gradex,hessex)
cccc      call prin2('potex=*',potex,nts*nd)
      call cpu_time(t2)
      print *, 'direct eval time = ', t2-t1
      call fgtpart_direct(nd,dim,delta,dmax,1,nsrc,1,ntt,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      call prin2('pottargex=*',pottargex,ntt*nd)
cccc      call prin2('gradtargex=*',gradtargex,3*ntt*nd)
cccc      call prin2('hesstargex=*',hesstargex,6*ntt*nd)

      if (ifpgh .gt. 0) call derr(potex,pot,nts*nd,errps)
      if (ifpgh .gt. 1) call derr(gradex,grad,dim*nts*nd,errgs)
      if (ifpgh .gt. 2) call derr(hessex,hess,nhess*nts*nd,errhs)
      
      if (ifpghtarg .gt. 0) call derr(pottargex,pottarg,ntt*nd,errpt)
      if (ifpghtarg .gt. 1)
     1    call derr(gradtargex,gradtarg,dim*ntt*nd,errgt)
      if (ifpghtarg .gt. 2) 
     1    call derr(hesstargex,hesstarg,nhess*ntt*nd,errht)
cc      errgs = 0
cc      errgt = 0

cc      errhs = 0
cc      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht,ifpgh,ifpghtarg)

      return
      end
c-----------------------------------------------------     
      subroutine dzero(vec,n)
      implicit real *8 (a-h,o-z)
      real *8 vec(*)

      do i=1,n
         vec(i) = 0
      enddo

      return
      end
c------------------------------------
      subroutine derr(vec1,vec2,n,erra)
      implicit real *8 (a-h,o-z)
      real *8 vec1(*),vec2(*)

      ra = 0
      erra = 0
      do i=1,n
         ra = ra + vec1(i)**2
         erra = erra + (vec1(i)-vec2(i))**2
      enddo

      if (sqrt(ra)/n .lt. 1d-10) then
         call prin2('vector norm =*', sqrt(ra)/n,1)
         call prin2('switch to absolute error*',a,0)
         erra = sqrt(erra)/n
      else
         erra = sqrt(erra/ra)
      endif
ccc      erra = sqrt(erra)/n

      return
      end
c----------------------------------
      subroutine errprint(errps,errgs,errhs,errpt,errgt,errht,
     1    ifpgh,ifpghtarg)
      implicit real *8 (a-h,o-z)
 1100 format(3(2x,e11.5))


      write(6,*) 'ifpgh is ', ifpgh
      write(6,*) 'ifpghtarg is ', ifpghtarg
      if (ifpgh .gt. 0) write(6,*) 'error in sources'
      if (ifpgh .gt. 0) call prin2('pot err =*', errps,1)
      if (ifpgh .gt. 1) call prin2('grad err=*', errgs,1)
      if (ifpgh .gt. 2) call prin2('hess err=*', errhs,1)
      write(6,*) 
      if (ifpghtarg .gt. 0) write(6,* ) 'error in targets'
      if (ifpghtarg .gt. 0) call prin2('pot err =*', errpt,1)
      if (ifpghtarg .gt. 1) call prin2('grad err=*', errgt,1)
      if (ifpghtarg .gt. 2) call prin2('hess err=*', errht,1)
      write(6,*)
      write(6,*)'==================='

      return
      end
      
