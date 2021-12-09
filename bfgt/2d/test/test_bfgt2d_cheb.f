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

      real *8 alpha,beta,targ(2)

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
      boxlen = 1.0d0
      
      rsig = 1.0d0/4000.0d0
cc      rsig = 0.005d0

      nd = 1
      dpars(1) = 0.2d0
      dpars(2) = 0.1d0

      dpars(3) = rsig
      
      dpars(4) = 0.312d0
      dpars(5) = 0.5d0

      dpars(6) = rsig/2.1

      dpars(7) = 0.678d0
      dpars(8) = 0.4d0

      dpars(9) = rsig/4.5
      
      dpars(10) = 0.412d0
      dpars(11) = 0.8d0

      dpars(12) = rsig/1.2
      
      dpars(13) = 0.12d0
      dpars(14) = 0.45d0

      dpars(15) = rsig/3.3
      
      norder = 16
      iptype = 0
      eta = 1.0d0

      npbox = norder*norder

      eps = 0.5d-10
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
      call chebtens_exps_2d(itype,norder,type,xref,umat,1,vmat,1,wts)

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

      k=5

      do ind=1,nd
         f(ind)=0
         do i=1,5
            idp = (i-1)*3
            rr = (xy(1)+0.5d0 - dpars(idp+1))**2 + 
     1          (xy(2)+0.5d0 - dpars(idp+2))**2
            sigma = dpars(idp+3)
            f(ind) = f(ind)+exp(-rr/sigma)
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
      ng=5
      do ind=1,nd
         do i=1,ng
            idp = (i-1)*3
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
            pot(ind)=pot(ind)+gf(1)*gf(2)
         enddo
      enddo

      return
      end







