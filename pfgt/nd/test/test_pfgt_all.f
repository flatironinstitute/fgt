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
         deltas(i)=10.0d0**(-i+2)
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
c     iperiod -> 0: free space; 1: periodic in all dimensions
      iperiod=0
c     ntarg: number of extra targets
      ntarg=0
c     nsrc: number of source points

c     whether to have monopole sources -> 1: yes; 0: no
      ifcharge=1
c     whether to have dipole sources -> 1: yes; 0: no
      ifdipole=0
c     evaluation flag for sources -> 1: pot; 2: pot+grad; 3: pot+grad+hess
      ifpgh=1
c     evaluation flag for targets -> 1: pot; 2: pot+grad; 3: pot+grad+hess
      ifpghtarg=2
c
c     
c
c     test all parameters
c
      iw = 70
      iw2 = 80

      ifuniform=1
      if (iperiod.eq.1 .and. dim.eq.2) then
         open(iw, file='errorp2d.txt', position='append')
         open(iw2, file='timingp2d.txt', position='append')
      elseif (iperiod.eq.0 .and. dim.eq.2) then
         open(iw, file='errorf2d.txt', position='append')
         open(iw2, file='timingf2d.txt', position='append')
      elseif (iperiod.eq.1 .and. dim.eq.3) then
         open(iw, file='errorp3d.txt', position='append')
         open(iw2, file='timingp3d.txt', position='append')
      elseif (iperiod.eq.0 .and. dim.eq.3) then
         open(iw, file='errorf3d.txt', position='append')
         open(iw2, file='timingf3d.txt', position='append')
      endif
      
      do i=1,3
         eps = epsvals(i)
         do j=1,10
            delta=deltas(j)
            nsrc=nsrcs(j)
            ntarg=nsrc
c
 2000       format(2x, i8, 2x)
c     write(iw,2000,advance="no") nint(1/ratio)

            call testpfgt(nd,dim,eps,delta,iperiod,nsrc,ntarg,ifcharge,
     1          ifdipole,ifpgh,ifpghtarg,ifuniform,pps(j,i),rerr)

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

      do j=1,10
         write(iw2,*) deltas(j), pps(j,1), pps(j,2), pps(j,3)
      enddo
      
      close(iw)
      close(iw2)

      end
c
c
c
c
      subroutine testpfgt(nd,dim,eps,delta,iperiod,nsrc,ntarg,
     1    ifcharge,ifdipole,ifpgh,ifpghtarg,ifuniform,pps,rerr)
      implicit real *8 (a-h,o-z)
      integer dim
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: sim(:,:)
      real *8, allocatable :: rnormal(:,:)
      real *8, allocatable :: charges(:,:),dipstr(:,:)
      real *8, allocatable :: charge1(:,:)
      real *8, allocatable :: pot(:,:),grad(:,:,:),hess(:,:,:)
      real *8, allocatable :: pottarg(:,:),gradtarg(:,:,:),
     1    hesstarg(:,:,:)
      real *8, allocatable :: potex(:,:),gradex(:,:,:),hessex(:,:,:)
      real *8, allocatable :: pottargex(:,:),gradtargex(:,:,:),
     1                             hesstargex(:,:,:)

      real *8 cen0(dim),bs0,shifts(dim)
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
      allocate(pottarg(nd,ntarg),gradtarg(nd,dim,ntarg))
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
            elseif (dim.eq.1) then
               sources(1,i) = (cos(i*pi/(nsrc+1))+1)/2
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

      if (iperiod.eq.1) then
         bs0=1.0d0
         do k=1,dim
            cen0(k)=0.5d0
         enddo

         nx=12
         h=bs0/(nx+1)

         if (dim.eq.1) then
            targ(1,1) = -bs0/2+cen0(1)
            ntthalf=1
            targ(1,2) = bs0/2+cen0(1)
            ntt=2
         elseif (dim.eq.2) then
            ii=0
            do kk=1,dim
            do i=1,nx
               ii=ii+1
               do k=1,dim
                  if (k.eq.kk) then
                     targ(k,ii)=-bs0/2+i*h+cen0(k)
                  else
                     targ(k,ii)=-bs0/2+cen0(k)
                  endif
               enddo
            enddo
            enddo
            ntthalf=ii
            do kk=1,dim
               do i=1,nx
                  ii=ii+1
                  do k=1,dim
                     if (k.eq.kk) then
                        targ(k,ii)=-bs0/2+i*h+cen0(k)
                     else
                        targ(k,ii)=bs0/2+cen0(k)
                     endif
                  enddo
               enddo
            enddo
            ntt=ii
         elseif (dim.eq.3) then
            ii=0
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)=-bs0/2+j*h+cen0(1)
               targ(2,ii)=-bs0/2+i*h+cen0(2)
               targ(3,ii)=-bs0/2+cen0(3)
            enddo
            enddo
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)=-bs0/2+j*h+cen0(1)
               targ(2,ii)=-bs0/2+cen0(2)
               targ(3,ii)=-bs0/2+i*h+cen0(3)
            enddo
            enddo
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)=-bs0/2+cen0(1)
               targ(2,ii)=-bs0/2+j*h+cen0(2)
               targ(3,ii)=-bs0/2+i*h+cen0(3)
            enddo
            enddo
            ntthalf=ii
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)=-bs0/2+j*h+cen0(1)
               targ(2,ii)=-bs0/2+i*h+cen0(2)
               targ(3,ii)= bs0/2+cen0(3)
            enddo
            enddo
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)=-bs0/2+j*h+cen0(1)
               targ(2,ii)= bs0/2+cen0(2)
               targ(3,ii)=-bs0/2+i*h+cen0(3)
            enddo
            enddo
            do i=1,nx
            do j=1,nx      
               ii=ii+1
               targ(1,ii)= bs0/2+cen0(1)
               targ(2,ii)=-bs0/2+j*h+cen0(2)
               targ(3,ii)=-bs0/2+i*h+cen0(3)
            enddo
            enddo
         endif
      endif

cccc      call prinf('ntt=*',ntt,1)
cccc      call prin2('targ=*',targ,dim*ntt)
      
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
      
      call prinf(' ifcharge is *',ifcharge,1)
      call prinf(' ifdipole is *',ifdipole,1)
      call prinf(' ifpgh is *',ifpgh,1)
      call prinf(' ifpghtarg is *',ifpghtarg,1)

      
      call cpu_time(t1)
      call pfgt(nd,dim,delta,eps,iperiod,bs0,cen0,
     1    nsrc,sources,
     2    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    ifpgh,pot,grad,hess,ntarg,targ,
     3    ifpghtarg,pottarg,gradtarg,hesstarg)
      call cpu_time(t2)
      pps=(nsrc*ifpgh+ntarg*ifpghtarg+0.0d0)/(t2-t1)
      call prin2('time in pfgt=*',t2-t1,1)
      call prin2('points per sec=*',pps,1)

cccc      call prin2('pot=*',pot,nd*nts)
cccc      call prin2('pottarg=*',pottarg,nd*ntt)
cccc      call prin2('gradtarg=*',gradtarg,nd*dim*ntt)
cccc      call prin2('hesstarg=*',hesstarg,nd*6*ntt)

      dmax = log(1.0d2/eps)*delta
      call cpu_time(t1)
c     try to compute the true potential for periodic conditions by including 
c     the nearest neighbors of the original computational box
c      
c     this will give the true potential only if the cutoff level is >=0, i.e.,
c     delta is not too large.
      if (iperiod.eq.1) then
         n0=1
      else
         n0=0
      endif
      
      ncell=(2*n0+1)**dim
      nsim=nsrc*ncell
      allocate(sim(dim,nsim))
      if (dim.eq.1) then
         kk=0
         do ii=-n0,n0
         do i=1,nsrc
            kk=kk+1
            sim(1,kk)=sources(1,i)+jj*bs0
            sim(2,kk)=sources(2,i)+ii*bs0
         enddo
         enddo
      elseif (dim.eq.2) then
         kk=0
         do ii=-n0,n0
         do jj=-n0,n0
         do i=1,nsrc
            kk=kk+1
            sim(1,kk)=sources(1,i)+jj*bs0
            sim(2,kk)=sources(2,i)+ii*bs0
         enddo
         enddo
         enddo
      elseif (dim.eq.3) then
         mm=0
         do ii=-n0,n0
         do jj=-n0,n0
         do kk=-n0,n0
         do i=1,nsrc
            mm=mm+1
            sim(1,mm)=sources(1,i)+kk*bs0
            sim(2,mm)=sources(2,i)+jj*bs0
            sim(3,mm)=sources(3,i)+ii*bs0
         enddo
         enddo
         enddo
         enddo
      endif
      
      if (ifcharge.eq.1) then
         allocate(charge1(nd,nsim))
         k=0
         do nc=1,ncell
         do i=1,nsrc
            k=k+1
            do j=1,nd
               charge1(j,k)=charges(j,i)
            enddo
         enddo
         enddo
      endif
      
      iperiod0=0
      call pfgt_direct(nd,dim,delta,dmax,iperiod0,shifts,
c     1    1,nsrc,1,nts,sources,
c     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     1    1,nsim,1,nts,sim,
     1    ifcharge,charge1,ifdipole,rnormal,dipstr,
     2    sources,ifpgh,potex,gradex,hessex)
cccc      call prin2('potex=*',potex,nts*nd)
      call cpu_time(t2)
      print *, 'direct eval time = ', t2-t1
      call pfgt_direct(nd,dim,delta,dmax,iperiod0,shifts,
     1    1,nsim,1,ntt,sim, 
     1    ifcharge,charge1,ifdipole,rnormal,dipstr,
     2    targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      call prin2('pottargex=*',pottargex,ntt*nd)
cccc      call prin2('gradtargex=*',gradtargex,3*ntt*nd)
cccc      call prin2('hesstargex=*',hesstargex,6*ntt*nd)

      if (ifpgh .gt. 0)
     1    call derr(potex,pot,nts*nd,errps,pnorm,errpa)
      if (ifpgh .gt. 1)
     1    call derr(gradex,grad,dim*nts*nd,errgs,gnorm,errga)
      if (ifpgh .gt. 2)
     1    call derr(hessex,hess,nhess*nts*nd,errhs,hnorm,errha)

      if (iperiod.eq.0) then
         if (ifpghtarg .gt. 0)
     1       call derr(pottargex,pottarg,ntt*nd,errpt,tmp1,tmp2)
         if (ifpghtarg .gt. 1)
     1       call derr(gradtargex,gradtarg,dim*ntt*nd,errgt,tmp1,tmp2)
         if (ifpghtarg .gt. 2) 
     1       call derr(hesstargex,hesstarg,nhess*ntt*nd,errht,
     2       tmp1,tmp2)
      else
         if (ifpghtarg .gt. 0) call derr(pottarg,pottarg(1,ntthalf+1),
     1       ntthalf*nd,errpt,tmp1,tmp2)
         if (ifpghtarg .gt. 1)
     1       call derr(gradtarg,gradtarg(1,1,ntthalf+1),
     2       dim*ntthalf*nd,errgt,tmp1,tmp2)
         if (ifpghtarg .gt. 2) 
     1       call derr(hesstarg,hesstarg(1,1,ntthalf+1),
     2       nhess*ntthalf*nd,errht,tmp1,tmp2)
      endif
cc      errgs = 0
cc      errgt = 0

cc      errhs = 0
cc      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht,ifpgh,ifpghtarg)

      if (iperiod.eq.1) then
         rerr=errpt/pnorm
         if (rerr.gt.eps*10) print *, pnorm,errpa,gnorm,errga
      else
         rerr=errps
      endif
      
      return
      end
c
c
c
c
      subroutine dzero(vec,n)
      implicit real *8 (a-h,o-z)
      real *8 vec(*)

      do i=1,n
         vec(i) = 0
      enddo

      return
      end
c
c
c
c
      subroutine derr(vec1,vec2,n,relerr,rnorm1,abserr)
      implicit real *8 (a-h,o-z)
      real *8 vec1(*),vec2(*)

      ra = 0
      erra = 0
      do i=1,n
         ra = ra + vec1(i)**2
         erra = erra + (vec1(i)-vec2(i))**2
      enddo

      rnorm1=sqrt(ra)
      abserr=sqrt(erra)
      relerr=abserr/rnorm1

      return
      end
c
c
c
c
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
c
c
c
c
      subroutine print_tree2d_matlab(ndim,itree,ltree,nboxes,centers,
     1   boxsize,nlevels,iptr,ns,src,nt,targ,fname1,fname2,fname3)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(12)
c            pointer to various arrays inside itree
c          ns - integer
c            number of sources
c          src - real *8 (2,ns)
c            xy coorindates of source locations
c          nt - integer
c            number of targets
c          targ - real *8 (2,nt)
c            xy cooridnates of target locations
c          fname1 - character *
c            file name to which tree info is to be written
c          fname1 - character *
c            file name to which source points are to be written
c          fname3 - character *
c            file name to which target points are to be written
c 
c          output
c          files with name fname1, fname2, fname3,
c            which contains the tree info, source points, target points
c            file can be plotted using the matlab script
c            tree_plot.m

      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(12),ns,nt
      real *8 centers(ndim,nboxes),boxsize(0:nlevels)
      real *8 x(5),y(5),src(ndim,ns),targ(ndim,nt)
      character (len=*) fname1,fname2,fname3

      open(unit=33,file=trim(fname1))
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(4)+ibox-1).eq.0) then
           ilev = itree(iptr(2)+ibox-1)
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           if (ndim.eq.2) then
              y1 = centers(2,ibox) - bs/2
              y2 = centers(2,ibox) + bs/2
           endif

           if (ndim.eq.2) then
              write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
           else
              write(33,1111) x1,x2
           endif
         endif
      enddo
      close(33)

 2222 format(2(2x,e11.5))

      open(unit=33,file=trim(fname2))
      if (ns .gt. 0) then
         do i=1,ns
            write(33,2222) src(1,i),src(2,i)
         enddo
      endif
      
      close(33)
      open(unit=33,file=trim(fname3))
      if (nt .gt. 0) then
         do i=1,nt
            write(33,2222) targ(1,i),targ(2,i)
         enddo
      endif

      close(33)

      return
      end
c
c
c
c
