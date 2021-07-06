      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: xref(:,:)
      real *8 xyztmp(3),rintl(0:200),umat,vmat,wts
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zpars

      real *8, allocatable :: pot(:,:,:),potex(:,:,:)
      complex *16 ima,zz,ztmp,zk

      real *8 alpha,beta
      character *9 fname1,fname3
      character *8 fname2
      real *8 src(2),targ(2)

      character *1 type
      data ima/(0.0d0,1.0d0)/

      external fgaussn,fgauss1
      logical flag

      call prini(6,13)
      zk = ima
      done = 1
      pi = atan(done)*4
c
c      initialize function parameters
c
      delta = 4d-4
      delta = 1d-1/5120*(1-1/sqrt(5.0d0))/2
      
      boxlen = 1.0d0
      
      rsig = 1.0d0/4000.0d0
      rsig = 0.00025d0
      rsig = 1.0d-4
      
      nd = 1
      dpars(1) = 0.75d0
      dpars(2) = 0.5d0

      dpars(3) = rsig
      dpars(4) = 1/pi/rsig
      
      dpars(5) = 0.25d0
      dpars(6) = 0.5d0

      dpars(7) = rsig
      dpars(8) = -0.5d0/pi/rsig
      
      dpars(9) = 0.678d0
      dpars(10) = 0.4d0

      dpars(11) = rsig/4.5
      dpars(12) = 1.0d0/pi/rsig
      
      dpars(13) = 0.412d0
      dpars(14) = 0.8d0

      dpars(15) = rsig/1.2
      dpars(16) = 1/pi/dpars(15)
      
      dpars(17) = 0.12d0
      dpars(18) = 0.45d0

      dpars(19) = rsig/3.3
      dpars(20) = 0.5d0/pi/rsig
      
      norder = 16
      iptype = 1
      eta = 1.0d0

      npbox = norder*norder

      eps = 1.5d-13
      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1   fgaussn,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)

      allocate(fvals(nd,npbox,nboxes),centers(2,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,fgaussn,nd,
     1  dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,fvals,
     2  centers,boxsize)

      fname1 = 'tree.data'
      fname2 = 'src.data'
      fname3 = 'targ.data'
      
      ns=0
      nt=0
      call print_tree_matlab(itree,ltree,nboxes,centers,boxsize,nlevels,
     1   iptr,ns,src,nt,targ,fname1,fname2,fname3)
      

      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1   (nboxes*norder**2+0.0d0)/(t2-t1),1)
c
c
c       convert values to coefs
c
      
cccc      npols = norder*(norder+1)*(norder+2)/6
      npols = norder*norder

      allocate(pot(nd,npbox,nboxes))

      do i=1,nboxes
        do j=1,npbox
           do ind=1,nd
              pot(ind,j,i) = 0
           enddo
        enddo
      enddo

      type = 'f'
      
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
      call bfgt2d(nd,delta,eps,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
     2   pot,timeinfo,tprecomp)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time taken in fgt=*',t2-t1,1)

      nlfbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) nlfbox = nlfbox+1
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo
      
      call prin2('speed in pps with precomputation included=*',
     1    (npbox*nlfbox+0.0d0)/(t2-t1),1)
      call prin2('speed in pps with precomputation excluded=*',
     1    (npbox*nlfbox+0.0d0)/d,1)

      erra = 0.0d0
      ra = 0.0d0

      allocate(potex(nd,npbox,nboxes))

      itype = 0
      allocate(xref(2,npbox))
      call legetens_exps_2d(itype,norder,type,xref,umat,1,vmat,1,wts)

      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              targ(1)=centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
              targ(2)=centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0

              call exact(nd,delta,targ,dpars,potex(1,j,ibox))

              do ind=1,nd
                 erra = erra + (pot(ind,j,ibox)-potex(ind,j,ibox))**2
                 ra = ra + potex(ind,j,ibox)**2
              enddo
            enddo
          endif
        enddo
      enddo


      erra = sqrt(erra/ra)
      call prin2('relative l2 error=*',erra,1)
cccc      call prin2('ra=*',ra,1)

      stop
      end
c
c
c
c 
      subroutine fgaussn(nd,xy,dpars,zpars,ipars,f)
c
c       compute three gaussians, their
c       centers are given in dpars(1:3*nd), and their 
c       variances in dpars(3*nd+1:4*nd)
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(*),f(nd),xy(2)

      k=2

      do ind=1,nd
         f(ind)=0
         do i=1,k
            idp = (i-1)*4
            rr = (xy(1)+0.5d0 - dpars(idp+1))**2 + 
     1          (xy(2)+0.5d0 - dpars(idp+2))**2
            sigma = dpars(idp+3)
            f(ind) = f(ind)+dpars(idp+4)*exp(-rr/sigma)
         enddo
      enddo

      return
      end

c
c
c
c 
      subroutine exact(nd,delta,targ,dpars,pot)

      implicit real*8 (a-h,o-z)
      real*8 targ(2),pot(nd)
      real*8 gf(2),dpars(*)
c
      one=1.0d0
      pi=4*atan(one)

c-----------------------
      do ind=1,nd
        pot(ind)=0.0d0
      enddo
c-----------------------
      ng=2
      do ind=1,nd
         do i=1,ng
            idp = (i-1)*4
            sigma = dpars(idp+3)
            dc = sigma
            d = delta
         
            do k=1,2
               c=dpars(idp+k)
               x=targ(k)+0.5d0
               gf(k)=sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     1             *(-erf((-dc-d+x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc))
     2             +erf((x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)))
     3             /dsqrt(((dc+d)/d/dc))
            enddo
            pot(ind)=pot(ind)+dpars(idp+4)*gf(1)*gf(2)
         enddo
      enddo

      return
      end
c      
c      
c      
c      
      subroutine print_tree_matlab(itree,ltree,nboxes,centers,boxsize,
     1   nlevels,iptr,ns,src,nt,targ,fname1,fname2,fname3)
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
      real *8 centers(2,nboxes),boxsize(0:nlevels),src(2,ns),targ(2,nt)
      real *8 x(5),y(5)
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

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
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
c
c
c
c
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
ccc      

      return
      end
c----------------------------------











