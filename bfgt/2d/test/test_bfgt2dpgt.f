      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: xref(:,:)
      real *8 xyztmp(3),rintl(0:200),umat,vmat,wts
      real *8 targs(2,1000000)
      real *8 pote(1000000)
      real *8 grade(2,1000000)
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zpars

      real *8, allocatable :: pot(:,:,:), potex(:,:)
      real *8, allocatable :: grad(:,:,:,:), gradex(:,:,:)
      real *8, allocatable :: potexe(:)
      real *8, allocatable :: gradexe(:,:)
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
      delta = 1d-4
ccc      delta = 1.0d-4
      boxlen = 1.0d0
      
      rsig = 1.0d0/4000.0d0
ccc      rsig = 0.005d0

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
ccc      call prin2('centers =*',centers(1,1),2*nboxes)
c
c
c       convert values to coefs
c
      
cccc      npols = norder*(norder+1)*(norder+2)/6
      npols = norder*norder

      allocate(pot(nd,npbox,nboxes))
      allocate(grad(nd,npbox,2,nboxes))

      do i=1,nboxes
        do j=1,npbox
           do ind=1,nd
              pot(ind,j,i) = 0
              grad(ind,j,1,i) = 0
              grad(ind,j,2,i) = 0
           enddo
        enddo
      enddo

      ntarg = 1000000
      do i = 1,ntarg
         targs(1,i) = rand()-0.51d0
         targs(2,i) = rand()-0.52d0
      enddo
      targs(1,1) = 0.32767D-1
      targs(2,1) = -.28104D0
      targs(1,1) = dpars(1) - 0.52d0
      targs(2,1) = dpars(2) - 0.51d0
      call prin2('targ is =*',targs(1,1),2*ntarg)
      type = 'f'
      
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
      call bfgt2dpgt(nd,delta,eps,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
     2   pot,grad,ntarg,targs,pote,grade,timeinfo,tprecomp)
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
ccc      call prin2('pote is =*',pote,ntarg)
ccc      call prin2('grade is =*',grade,2*ntarg)

      erra = 0.0d0
      ra = 0.0d0
      errb = 0.0d0
      rb = 0.0d0

      allocate(potex(npbox,nboxes))
      allocate(gradex(2,npbox,nboxes))
      allocate(potexe(ntarg))
      allocate(gradexe(2,ntarg))

      itype = 0
      allocate(xref(2,npbox))
      call legetens_exps_2d(itype,norder,type,xref,umat,1,vmat,1,wts)

      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              targ(1)=centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
              targ(2)=centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0

              call exact(nd,delta,targ,dpars,potex(j,ibox),
     1           gradex(1,j,ibox))
              do ind=1,nd
                 erra = erra + (pot(1,j,ibox)-potex(j,ibox))**2
                 ra = ra + potex(j,ibox)**2
                 errb = errb + 
     1                  (grad(1,j,1,ibox)-gradex(1,j,ibox))**2
                 rb = rb + gradex(1,j,ibox)**2
                 errb = errb + 
     1                  (grad(1,j,2,ibox)-gradex(2,j,ibox))**2
                 rb = rb + gradex(2,j,ibox)**2
              enddo
            enddo
          endif
        enddo
      enddo


      erra = sqrt(erra/ra)
      errb = sqrt(errb/rb)
      call prin2('relative l2 error=*',erra,1)
      call prin2('relative grad l2 error=*',errb,1)
cccc      call prin2('ra=*',ra,1)
      erra = 0.0d0
      ra = 0.0d0
      errb = 0.0d0
      rb = 0.0d0

      do ii = 1,ntarg
         call exact(nd,delta,targs(1,ii),dpars,potexe(ii),
     1           gradexe(1,ii))
ccc      call prin2('potex at targ is =*',potex(1,1),1)
ccc      call prin2('gradex at targ is =*',gradex(1,1,1),2)
ccc      call prin2('pote is =*',pote(ii),1)
ccc      call prin2('grade is =*',grade(1,ii),2)
         erra = erra + (pote(ii)-potexe(ii))**2
         ra = ra + potexe(ii)**2
         errb = errb + 
     1          (grade(1,ii)-gradexe(1,ii))**2
         rb = rb + gradexe(1,ii)**2
         errb = errb + 
     1          (grade(2,ii)-gradexe(2,ii))**2
         rb = rb + gradexe(2,ii)**2
      enddo
ccc      call prin2('potexe is =*',potexe,ntarg)
ccc      call prin2('gradexe is =*',gradexe,2*ntarg)
      erra = sqrt(erra)
      errb = sqrt(errb)
      call prin2(' targ ra =*',sqrt(ra),1)
      call prin2('relative targ l2 error=*',erra,1)
      call prin2('relative targ grad l2 error=*',errb,1)
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

      one=1.0d0
      pi=4*atan(one)
c
      dt = 0.1d0
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
      subroutine exact(nd,delta,targ,dpars,pot,grad)

      implicit real*8 (a-h,o-z)
      real*8 targ(2),pot(nd),grad(nd,2)
      real*8 gf(2),gfp(2),dpars(*)
c
      one=1.0d0
      pi=4*atan(one)
      dt = 0.1d0

c-----------------------
      do ind=1,nd
        pot(ind)=0.0d0
        grad(ind,1)=0.0d0
        grad(ind,2)=0.0d0
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
               arg1 = (-dc-d+x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)
               darg = 1.0d0/d/dsqrt((dc+d)/d/dc)
               arg2 = (x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)
               gfp(k)=gf(k)*(-2.0d0*(x-c)/(dc+d)) +
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2             *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3             darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc))
            enddo
            pot(ind)=pot(ind)+gf(1)*gf(2)
            grad(ind,1)=grad(ind,1)+gfp(1)*gf(2)
            grad(ind,2)=grad(ind,2)+gf(1)*gfp(2)
         enddo
      enddo
      return
      end







