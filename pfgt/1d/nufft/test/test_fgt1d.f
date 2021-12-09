c
c
c
c
c
c
c
c
c
c
c
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:),targ(:)
      real *8, allocatable :: rnormal(:)
      real *8, allocatable :: charges(:,:),dipstr(:,:)
      real *8, allocatable :: pot(:,:),grad(:,:),hess(:,:)
      real *8, allocatable :: pottarg(:,:),gradtarg(:,:),
     1    hesstarg(:,:)
      real *8, allocatable :: potex(:,:),gradex(:,:),hessex(:,:)
      real *8, allocatable :: pottargex(:,:),gradtargex(:,:),
     1                             hesstargex(:,:)

      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nsrc =  2 000 000
cccc      nsrc =  1 000 000
      ntarg = nsrc
      nd = 1
      delta = 5.0d-4
      bb = 1.0d0/(2.0d0**6)
      n = 1
      delta = bb*bb/(1.5*1.5)*2.0d0**n
cccc      delta = 0.1*delta
      delta=1.0d-8
      
      call prin2(' delta = *',delta,1)
      call prinf_long(' nsrc = *',nsrc,1)
c
      allocate(sources(nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(rnormal(nsrc))
      allocate(targ(ntarg))
      allocate(pot(nd,nsrc),grad(nd,nsrc),hess(nd,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,ntarg))
      allocate(hesstarg(nd,ntarg))

      rin = 0.3d0
      rwig = 0.15d0
      nwig = 10
      
      do i=1,nsrc
c        nonuniform source distribution
         theta = hkrand(0)*pi
         phi = hkrand(0)*2*pi
         sources(i) = (rin+rwig*cos(theta))*sin(theta)*cos(phi)+0.5d0
         
c        uniform source distribution
         sources(i) = hkrand(0)
         
         rnormal(i) = hkrand(0)
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
      
      sources(1) = 0.0d0
      sources(2) = 1.0d0
      sources(3) = 0.05d0
 
      do i=1,ntarg
         targ(i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nd,nts),gradex(nd,nts),hessex(nd,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,ntt))
      allocate(hesstargex(nd,ntt))

      eps = 1d-10
c
      call dzero(potex,nts*nd)
      call dzero(gradex,nts*nd)
      call dzero(hessex,nts*nd)
      
      call dzero(pottargex,ntt*nd)
      call dzero(gradtargex,ntt*nd)
      call dzero(hesstargex,ntt*nd)
      
      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 0
      call prinf(' ifcharge is *',ifcharge,1)
      call prinf(' ifdipole is *',ifdipole,1)
      call prinf(' ifpgh is *',ifpgh,1)
      call prinf(' ifpghtarg is *',ifpghtarg,1)

      iper = 0
      call fgt1d(nd,delta,eps,nsrc,sources,ifcharge,charges,
     1    ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2    ntarg,targ,ifpghtarg,pottarg,gradtarg,
     3    hesstarg)

cccc      call prin2('pot=*',pot,nd*nts)
cccc      call prin2('pottarg=*',pottarg,nd*ntt)
cccc      call prin2('gradtarg=*',gradtarg,nd*3*ntt)
cccc      call prin2('hesstarg=*',hesstarg,nd*6*ntt)

      dmax = 1.0d30
      call cpu_time(t1)
      call fgt1dpart_direct_vec(nd,delta,dmax,1,nsrc,1,nts,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    sources,ifpgh,potex,gradex,hessex)
cccc      call prin2('potex=*',potex,nts*nd)
      call cpu_time(t2)
      print *, 'direct eval time = ', t2-t1
      call fgt1dpart_direct_vec(nd,delta,dmax,1,nsrc,1,ntt,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      call prin2('pottargex=*',pottargex,ntt*nd)
cccc      call prin2('gradtargex=*',gradtargex,3*ntt*nd)
cccc      call prin2('hesstargex=*',hesstargex,6*ntt(nd)

      if (ifpgh .gt. 0) call derr(potex,pot,nts*nd,errps)
      if (ifpgh .gt. 1) call derr(gradex,grad,nts*nd,errgs)
      if (ifpgh .gt. 2) call derr(hessex,hess,nts*nd,errhs)
      
      if (ifpghtarg .gt. 0) call derr(pottargex,pottarg,ntt*nd,errpt)
      if (ifpghtarg .gt. 1)
     1    call derr(gradtargex,gradtarg,ntt*nd,errgt)
      if (ifpghtarg .gt. 2) 
     1    call derr(hesstargex,hesstarg,ntt*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht,ifpgh,ifpghtarg)

      stop
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
      
