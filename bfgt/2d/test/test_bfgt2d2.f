      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer ifpgh,ifpghtarg
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: xref(:,:)
      real *8 xyztmp(3),rintl(0:200),wts
      real *8 targs(2,1000 000)
      real *8 pote(1000 000)
      real *8 grade(2,1000 000)
      real *8 hesse(3,1000 000)
c
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zpars

      real *8, allocatable :: pot(:,:,:), potex(:,:,:)
      real *8, allocatable :: grad(:,:,:,:), gradex(:,:,:,:)
      real *8, allocatable :: hess(:,:,:,:), hessex(:,:,:,:)

      real *8, allocatable :: coefsp(:,:,:)
      real *8, allocatable :: coefsg(:,:,:,:)
      real *8, allocatable :: coefsh(:,:,:,:)

      real *8, allocatable :: adiff(:,:)
      
      real *8, allocatable :: potexe(:,:)
      real *8, allocatable :: gradexe(:,:,:)
      real *8, allocatable :: hessexe(:,:,:)
      complex *16 ima,zz,ztmp,zk

      real *8 xs(100),ws(100),umat(2000),vmat(2000)
      real *8 ainte(2000),endinter(1000),work(10000)
      real *8 polin(100),polout(100)
      
      real *8 alpha,beta
      character *9 fname1,fname3
      character *8 fname2
      real *8 src(2),targ(2)

      character *1 type
      data ima/(0.0d0,1.0d0)/

      external fgaussn
cccc      logical flag

      call prini(6,13)
cccc      call prini_off()
      zk = ima
      done = 1
      pi = atan(done)*4
c
c      initialize function parameters
c
      delta = 1d-1/5120*(1-1/sqrt(5.0d0))/2
      delta = 1d-4
      
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

      ntarg = 1 000 000
      do i=1,ntarg
         do j=1,2
            targs(j,i) = hkrand(0)-0.5d0
         enddo
      enddo
      
      eps = 1.5d-10
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
      allocate(grad(nd,2,npbox,nboxes))
      allocate(hess(nd,3,npbox,nboxes))

      allocate(potexe(nd,ntarg))
      allocate(gradexe(nd,2,ntarg))
      allocate(hessexe(nd,3,ntarg))

      do i=1,nboxes
        do j=1,npbox
           do ind=1,nd
              pot(ind,j,i) = 0
              grad(ind,1,j,i) = 0
              grad(ind,2,j,i) = 0
              hess(ind,1,j,i) = 0
              hess(ind,2,j,i) = 0
              hess(ind,3,j,i) = 0
           enddo
        enddo
      enddo

      type = 'f'
      ifpgh= 3
      ifpghtarg=1
      
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
c      call bfgt2dpght(nd,delta,eps,nboxes,nlevels,ltree,itree,
c     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
c     2   ifpgh,pot,grad,hess,ntarg,targs,
c     3   ifpghtarg,pote,grade,hesse,timeinfo,tprecomp)
      call bfgt2d(nd,delta,eps,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
     2   pot,timeinfo,tprecomp)

      call cpu_time(t2) 
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

      allocate(potex(nd,npbox,nboxes))
      allocate(gradex(nd,2,npbox,nboxes))
      allocate(hessex(nd,3,npbox,nboxes))

      itype = 0
      allocate(xref(2,npbox))
      call legetens_exps_2d(itype,norder,type,xref,umat,1,vmat,1,wts)

      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              targ(1)=centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
              targ(2)=centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0

ccc             call exact_old(nd,delta,targ,dpars,potex(nd,j,ibox))
             call exact(nd,delta,targ,dpars,potex(nd,j,ibox),
     1           gradex(nd,1,j,ibox),hessex(nd,1,j,ibox))
            enddo
          endif
        enddo
      enddo

      call treedata_derror(nd,nlevels,itree,iptr,
     1    npbox,potex,pot,abserrp,rnormp,nleaf)
      errp = abserrp/rnormp
      call prin2('relative pot l2 error=*',errp,1)

      do j=1,ntarg
         call exact(nd,delta,targs(1,j),dpars,potexe(nd,j),
     1           gradexe(nd,1,j),hessexe(nd,1,j))
      enddo

cccc  call prin2('ra=*',ra,1)

c
c       convert values to coefs
c
      itype=2
      call legeexps(itype,norder,xs,umat,vmat,ws)

      allocate(coefsp(nd,npbox,nboxes))
      allocate(coefsg(nd,2,npbox,nboxes))
      allocate(coefsh(nd,3,npbox,nboxes))

      itype=0
      call cpu_time(t1) 
C     $     t1 = omp_get_wtime()
      call treedata_trans2d(nd,itype,nlevels,itree,iptr,boxsize,
     1    norder,pot,coefsp,umat)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time on treedata_trans2d=*',t2-t1,1)
      
      call prin2('speed in pps=*',
     1    (npbox*nlfbox+0.0d0)/(t2-t1),1)

c     construct 1D Chebyshev differentiation matrix
cccc      itype=2
cccc      call chebinmt(norder,ainte,adiff,xs,ws,endinter,
cccc     1    itype,work)
      allocate(adiff(norder,norder))
      do 2600 i=1,norder
c 
         do 2200 j=1,norder+3
            polin(j)=0
 2200    continue
c 
         polin(i)=1
         call legediff(polin,norder+1,polout)
c 
         polout(norder)=0
         do 2400 j=1,norder
            adiff(j,i)=polout(j)
 2400    continue
c 
 2600 continue
      
c     compute the expansion coefficients of gradient and hessian
      call cpu_time(t1) 
C     $     t1 = omp_get_wtime()
      call treedata_coefs_p_to_gh2d(nd,nlevels,itree,iptr,
     1    boxsize,norder,coefsp,coefsg,coefsh,adiff)
c     evaluate the gradient
      itype=0
      call treedata_trans2d(nd*2,itype,nlevels,itree,iptr,boxsize,
     1    norder,coefsg,grad,vmat)
c     evaluate the hessian
      call treedata_trans2d(nd*3,itype,nlevels,itree,iptr,boxsize,
     1    norder,coefsh,hess,vmat)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time on calculating grad and hess=*',t2-t1,1)
      
      call prin2('speed in pps=*',
     1    (npbox*nlfbox+0.0d0)/(t2-t1),1)

      call treedata_derror(nd*2,nlevels,itree,iptr,
     1    npbox,gradex,grad,abserrg,rnormg,nleaf)
      call treedata_derror(nd*3,nlevels,itree,iptr,
     1    npbox,hessex,hess,abserrh,rnormh,nleaf)
      errg = abserrg/rnormg
      errh = abserrh/rnormh
      call prin2('relative grad l2 error=*',errg,1)
      call prin2('relative hess l2 error=*',errh,1)

c
c
      itype=1
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
c     evaluate the potential, gradient and hessian at targets 
      call treedata_evalpgh2d(nd,itype,itree,ltree,nboxes,nlevels,
     1    iptr,centers,boxsize,norder,coefsp,
     2    ntarg,targs,pote,grade,hesse)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time on extra targ pot eval=*',t2-t1,1)
      call prin2('speed in pps=*',(ntarg+0.0d0)/(t2-t1),1)

c     compute relative error
      call derr(potexe,pote,nd*ntarg,errpe)
      call derr(gradexe,grade,nd*2*ntarg,errge)
      call derr(hessexe,hesse,nd*3*ntarg,errhe)

      call prin2('relative pottarg l2 error=*',errpe,1)
      call prin2('relative gradtarg l2 error=*',errge,1)
      call prin2('relative hesstarg l2 error=*',errhe,1)
      
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
c
c 
      subroutine exact(nd,delta,targ,dpars,pot,grad,hess)

      implicit real*8 (a-h,o-z)
      real*8 targ(2),pot(nd),grad(nd,2),hess(nd,3)
      real*8 gf(2),gfp(2),dpars(*)
      real*8 gfpp(2)
c
      one=1.0d0
      pi=4*atan(one)

c-----------------------
      do ind=1,nd
        pot(ind)=0.0d0
        grad(ind,1)=0.0d0
        grad(ind,2)=0.0d0
        hess(ind,1)=0.0d0
        hess(ind,2)=0.0d0
        hess(ind,3)=0.0d0
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
c
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
               gfpp(k)=gfp(k)*(-2.0d0*(x-c)/(dc+d)) +
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2             *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3             darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc))
               gfpp(k)=gfpp(k)+gf(k)*(-2.0d0/(dc+d)) +
     1            ( (-2.0d0*(x-c)/(dc+d))*
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2             *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3             darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc)))+
     1            (sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2             *darg*(2*arg1*dexp(-arg1*arg1) - 
     1                   2*arg2*dexp(-arg2*arg2))
     2             *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3             darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc)))
            enddo
            pot(ind)=pot(ind)+dpars(idp+4)*gf(1)*gf(2)
            grad(ind,1)=grad(ind,1)+dpars(idp+4)*gfp(1)*gf(2)
            grad(ind,2)=grad(ind,2)+dpars(idp+4)*gf(1)*gfp(2)
            hess(ind,1)=hess(ind,1)+dpars(idp+4)*gfpp(1)*gf(2)
            hess(ind,2)=hess(ind,2)+dpars(idp+4)*gfp(1)*gfp(2)
            hess(ind,3)=hess(ind,3)+dpars(idp+4)*gf(1)*gfpp(2)
         enddo
      enddo
      return
      end
c
c 
      subroutine exact_old(nd,delta,targ,dpars,pot)

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











