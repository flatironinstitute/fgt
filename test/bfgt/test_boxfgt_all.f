      implicit real *8 (a-h,o-z)
      real *8 epsvals(5),deltas(20),pps(20,5)
      integer ntot(20)
c	  
      call prini(6,13)

      pi = 4*atan(1.0d0)
      
      epsvals(1) = 1d-3
      epsvals(2) = 1d-6
      epsvals(3) = 1d-9
      epsvals(4) = 1d-12
c
      do i=1,20
         deltas(i)=10.0d0**(-i+1)
      enddo

c     ipoly=0: Legendre; 1: Chebyshev
      ipoly=0
c     iperiod=0: free-space; 1: doubly periodic
      iperiod=1
c     polynomial expansion order in each dimension
      norder=16
c     nd - number of different densities
      nd=1
c     ndim - dimension of the underlying space
      ndim=2
c     ifpgh = 1: potential; 2: pot+grad; 3: pot+grad+hess
      ifpgh=1
c     ifpghtarg: flag for arbitrary targets
      ifpghtarg=0
c     ntarg: number of extra targets
      ntarg=0
c     
c
c     test all parameters
c
      iw = 70
      iw2 = 80
      if (iperiod.eq.1 .and. ndim.eq.2) then
         open(iw, file='errorp2d.txt', position='append')
         open(iw2, file='timingp2d.txt', position='append')
      elseif (iperiod.eq.0 .and. ndim.eq.2) then
         open(iw, file='errorf2d.txt', position='append')
         open(iw2, file='timingf2d.txt', position='append')
      elseif (iperiod.eq.1 .and. ndim.eq.3) then
         open(iw, file='errorp3d.txt', position='append')
         open(iw2, file='timingp3d.txt', position='append')
      elseif (iperiod.eq.0 .and. ndim.eq.3) then
         open(iw, file='errorf3d.txt', position='append')
         open(iw2, file='timingf3d.txt', position='append')
      endif
      
c      write(iw,*) 'ipoly=',ipoly,' iperiod=',iperiod,' norder =',norder
c      write(iw,*) 'ndim=',ndim,' ifpgh=',ifpgh
      
      do i=3,3
         eps = epsvals(i)
c     figure 6.1
         do j=4,4
c     figures 6.2 and 6.3
c         do j=1,11
            delta=deltas(j)
            if (ndim.eq.2) then
c     figure 6.1 2d
               kdelta=2**(j+2)
               delta=1.0d0/kdelta**2/pi**2
            elseif (ndim.eq.3) then
c     figure 6.1 3d
               kdelta=j+1
            endif
c
 2000       format(2x, i8, 2x)
c     write(iw,2000,advance="no") nint(1/ratio)

c     figure 6.1
            call testboxfgt(nd,ndim,norder,ipoly,kdelta,
c     figures 6.2 and 6.3
c            call testboxfgt(nd,ndim,norder,ipoly,delta,
     1          iperiod,eps,ntarg,ifpgh,ifpghtarg,rerr,pps(j,i),
     2          ntot(j))

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
c     figure 6.1
      do j=4,4
         write(iw2,*) ntot(j), pps(j,1), pps(j,2), pps(j,3)
c     figures 6.2 and 6.3
c     do j=1,11
c         write(iw2,*) deltas(j), pps(j,1), pps(j,2), pps(j,3)
      enddo
      
      close(iw)
      close(iw2)

      end
c
c
c
c     figure 6.1
      subroutine testboxfgt(nd,ndim,norder,ipoly,kdelta,iperiod,eps,
c     figures 6.2 and 6.3
c      subroutine testboxfgt(nd,ndim,norder,ipoly,delta,iperiod,eps,
     1    ntarg,ifpgh,ifpghtarg,errp,pps,ntot)
      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(8),ltree,ipars(100)
      integer ifpgh,ifpghtarg
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: xref(:,:)
      real *8 rintl(0:200),wts
c
      real *8 timeinfo(100)
      complex *16 zpars(10)

      real *8, allocatable :: pot(:,:,:), potex(:,:,:)
      real *8, allocatable :: grad(:,:,:,:), gradex(:,:,:,:)
      real *8, allocatable :: hess(:,:,:,:), hessex(:,:,:,:)

      real *8, allocatable :: coefs(:,:,:)
      real *8, allocatable :: coefsg(:,:,:,:)
      real *8, allocatable :: adiff(:,:)
      
      real *8, allocatable :: targs(:,:)

      real *8, allocatable :: pote(:,:)
      real *8, allocatable :: grade(:,:,:)
      real *8, allocatable :: hesse(:,:,:)

      real *8, allocatable :: potexe(:,:)
      real *8, allocatable :: gradexe(:,:,:)
      real *8, allocatable :: hessexe(:,:,:)

      complex *16 ima,zz,ztmp,zk

      real *8 xs(100),ws(100),vmat(2000)
      real *8 vpmat(2000),vppmat(2000)
      real *8 ainte(2000),endinter(1000),work(10000)
      real *8 polin(100),polout(100)

      real *8, allocatable :: umat(:,:),umat_nd(:,:,:)
      real *8 alpha,beta,src(ndim,1)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      real *8, allocatable :: targ(:)

      character *1 type
      data ima/(0.0d0,1.0d0)/

      external rhsfun,uexact
      done = 1
      pi = atan(done)*4

      ifnewtree=0
      
c     figure 6.1 2d
      delta=1.0d0/kdelta**2/pi**2
c     figure 6.1 3d
c      delta=1.0d-3
      
      ipars(1) = ndim
c     for exact solution
      dpars(51)=delta

      ipars(5)=iperiod
      
c     p in L^p norm - 0: L^infty norm; 1: L^1 norm; 2: L^2 norm
      iptype = 2
cccc      if (iperiod.eq.1) iptype = 2
cccc      eta = 1.0d0
      eta = 0.0d0
      
      if (iperiod.eq.0) then
c        number of gaussians in the rhs function
c     figure 6.3 3d
         ng = 5
c     figure 6.1 3d
         ng = kdelta
         ipars(2) = ng
      elseif (iperiod.eq.1) then
         do i=1,ndim
cccc  ipars(i+1)=4*i
c     figure 6.1 2d
            ipars(i+1)=kdelta
c            if (ndim.eq.2) then
c     figure 6.2
c               ipars(i+1)=32
c            elseif (ndim.eq.3) then
c     figure 6.3
c               ipars(i+1)=4
c            endif
               
         enddo
      endif
c     number of points per box
      npbox = norder**ndim

      type = 'f'

      ipars(10) = max(ifpgh,ifpghtarg)
      
cccc  call prini_off()

      zk = ipars(2)*7.05d0
c
c     initialize function parameters
c
      boxlen = 1.0d0
c     gaussian variance of the input data
      rsig = 1.0d0/4000.0d0
      rsig = 0.00025d0
      if (ndim.eq.2) then
         rsig = 1.0d-5
      elseif (ndim.eq.3) then
c     figure 6.3
         rsig = 4.0d-3
c     figure 6.1
         rsig = 1d-2
         rsig = rsig/kdelta
      endif
      
c     proper normalization of the input data
      rsign=(rsig*delta)**(ndim/2.0d0)
      
c     first gaussian
c     centers
      dpars(1) = 0.2d0
      dpars(2) = 0.1d0
      dpars(3) = 0.44d0
c     variance
      dpars(4) = rsig
c     strength
      dpars(5) = 1/pi/rsign

c     second gaussian
      dpars(6) = 0.3d0
      dpars(7) = 0.5d0
      dpars(8) = 0.25d0

      dpars(9) = rsig/2
      dpars(10) = -0.5d0/pi/rsign

      
c     third gaussian
      dpars(11) = 0.680
      dpars(12) = 0.4d0
      dpars(13) = 0.53d0

      dpars(14) = rsig/3
      dpars(15) = 1.0d0/pi/rsign
      
c     fourth gaussian
      dpars(16) = 0.41d0
      dpars(17) = 0.8d0
      dpars(18) = 0.67d0

      dpars(19) = rsig/4
      dpars(20) = 1/pi/rsign

c     fifth gaussian
      dpars(21) = 0.12d0
      dpars(22) = 0.45d0
      dpars(23) = 0.33d0

      dpars(24) = rsig/5
      dpars(25) = 1/pi/rsign


      nhess = ndim*(ndim+1)/2
      allocate(targs(ndim,ntarg),pote(nd,ntarg))
      allocate(grade(nd,ndim,ntarg),hesse(nd,nhess,ntarg))

      do i=1,ntarg
         do j=1,ndim
            targs(j,i) = hkrand(0)-0.5d0
         enddo
      enddo

      call prin2('delta=*',delta,1)

      epstree=1.0d-10
cccc      epstree=eps
      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,
     2    nboxes,nlevels,ltree,rintl)
      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)
      
      allocate(fvals(nd,npbox,nboxes),centers(ndim,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)
      call prinf('laddr=*',itree,2*(nlevels+1))
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c     compute the number of leaf boxes
      nlfbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            nlfbox = nlfbox+1
          endif
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)

      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1   (nboxes*npbox+0.0d0)/(t2-t1),1)

c     plot the tree
      ifplot=0
      if (ifplot.eq.1) then
         fname1 = 'trgtree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         ns=0
         nt=0
         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
      endif
      
c     allocate memory and initialization

      allocate(pot(nd,npbox,nboxes))
      allocate(grad(nd,ndim,npbox,nboxes))
      allocate(hess(nd,nhess,npbox,nboxes))

      do i=1,nboxes
      do j=1,npbox
      do ind=1,nd
         pot(ind,j,i) = 0
         if (ifpgh.ge.2) then
            do k=1,ndim
               grad(ind,k,j,i) = 0
            enddo
         endif
         if (ifpgh.ge.3) then
            do k=1,nhess
               hess(ind,k,j,i) = 0
            enddo
         endif
      enddo
      enddo
      enddo

      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      

      call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,pote,grade,hesse,timeinfo)

      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time taken in fgt=*',t2-t1,1)
      ntot=npbox*nlfbox+ntarg
      pps=(npbox*nlfbox+ntarg+0.0d0)/(t2-t1)
      call prinf('ntotal=*',ntot,1)
      call prin2('speed in pps=*',pps,1)
      
      allocate(potex(nd,npbox,nboxes))
      allocate(gradex(nd,ndim,npbox,nboxes))
      allocate(hessex(nd,nhess,npbox,nboxes))

c     compute exact solutions on tensor grid
      allocate(xref(ndim,npbox))
      itype = 0
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
      
      allocate(targ(ndim))
      do ilevel=1,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo

               call uexact(nd,targ,dpars,zpars,ipars,
     1             potex(1,j,ibox),gradex(1,1,j,ibox),
     2             hessex(1,1,j,ibox))
cccc               call rhsfun(nd,targ,dpars,zpars,ipars,fval)

cccc               write(34,*) targ(1), fval, potex(1,j,ibox), pot(1,j,ibox)
cccc               write(34,*) targ(1), potex(1,j,ibox), pot(1,j,ibox)
             enddo
          endif
        enddo
      enddo

      if (ifpgh.ge.1) then
         call treedata_derror(nd,nlevels,itree,iptr,
     1       npbox,potex,pot,abserrp,rnormp,nleaf)
         call treedata_lpnorm(ndim,iptype,ipoly,nd,nlevels,itree,
     1    iptr,boxsize,norder,npbox,fvals,rnormf,nleaf)
c     example 1a
         if (rnormp.lt.1d-6) then
            errp = abserrp
         else
            errp = abserrp/rnormp
         endif
         call prin2('absolute pot l2 error=*',abserrp,1)
         call prin2('pot l2 norm=*',rnormp,1)
         call prin2('relative pot l2 error=*',errp,1)
      endif

      if (ifpgh.ge.2) then
         call treedata_derror(nd*ndim,nlevels,itree,iptr,
     1       npbox,gradex,grad,abserrg,rnormg,nleaf)
         errg = abserrg/rnormg
         call prin2('relative grad l2 error=*',errg,1)
      endif

      if (ifpgh.ge.3) then
         call treedata_derror(nd*nhess,nlevels,itree,iptr,
     1       npbox,hessex,hess,abserrh,rnormh,nleaf)
         errh = abserrh/rnormh
         call prin2('relative hess l2 error=*',errh,1)
      endif
c
c     compute exact solutions on arbitrary targets      
      allocate(potexe(nd,ntarg))
      allocate(gradexe(nd,ndim,ntarg))
      allocate(hessexe(nd,nhess,ntarg))

      do j=1,ntarg
         call uexact(nd,targs(1,j),dpars,zpars,ipars,
     1       potexe(1,j),gradexe(1,1,j),hessexe(1,1,j))
      enddo
c
c     compute relative error
      if (ifpghtarg.ge.1) then
         call derr(potexe,pote,nd*ntarg,errpe)
         call prin2('relative pottarg l2 error=*',errpe,1)
      endif

      if (ifpghtarg.ge.2) then
         call derr(gradexe,grade,nd*ndim*ntarg,errge)
         call prin2('relative gradtarg l2 error=*',errge,1)
      endif

      if (ifpghtarg.ge.3) then
         call derr(hessexe,hesse,nd*nhess*ntarg,errhe)
         call prin2('relative hesstarg l2 error=*',errhe,1)
      endif

      return
      end
c
c
c
      subroutine rhsfun(nd,xyz,dpars,zpars,ipars,f)
c     right-hand-side function
c     for free space problem:  consisting of several gaussians, their
c       centers are given in dpars(1:3), their 
c       variances in dpars(4), and their strength in dpars(5)
c
c     for periodic conditions: simple sin/cos functions
c     
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars
      real *8 dpars(*),f(nd),xyz(*)
c
      iperiod=ipars(5)

      if (iperiod.eq.0) then
         call fgaussn(nd,xyz,dpars,zpars,ipars,f)
      elseif (iperiod.eq.1) then
         call fperiodic(nd,xyz,dpars,zpars,ipars,f)
      endif
c     
c     
      return
      end
c
c
c
c
      subroutine uexact(nd,targ,dpars,zpars,ipars,
     1    pot,grad,hess)
c     exact solution for the given rhs function
      implicit real*8 (a-h,o-z)
      real*8 targ(*),pot(nd),grad(nd,*),hess(nd,*)
      real*8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)
      

      iperiod=ipars(5)

      if (iperiod.eq.0) then
         call ugexact(nd,targ,dpars,zpars,ipars,pot,grad,hess)
      elseif (iperiod.eq.1) then
         call upexact(nd,targ,dpars,zpars,ipars,pot,grad,hess)
      endif
c     
      
      return
      end
c
c
c
c     
      subroutine fgaussn(nd,xyz,dpars,zpars,ipars,f)
c     right-hand-side function
c       consisting of several gaussians, their
c       centers are given in dpars(1:3), their 
c       variances in dpars(4), and their strength in dpars(5)
c
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars
      real *8 dpars(*),f(nd),xyz(*)
c     number of Gaussians
      ng=ipars(2)
      ndim=ipars(1)
      
      do ind=1,nd
         f(ind)=0
         do i=1,ng
            idp = (i-1)*5
            rr=0
            do k=1,ndim
               rr = rr + ( xyz(k) - (dpars(idp+k)-0.5d0) )**2  
            enddo
            sigma = dpars(idp+4)
            f(ind) = f(ind)+dpars(idp+5)*exp(-rr/sigma)
         enddo
      enddo

      return
      end
c
c
c
c
      subroutine ugexact(nd,targ,dpars,zpars,ipars,
     1    pot,grad,hess)
      implicit real*8 (a-h,o-z)
      real*8 targ(*),pot(nd),grad(nd,*),hess(nd,*)
      real*8 gf(10),gfp(10),dpars(*)
      real*8 gfpp(10)
      integer ipars(*)
      complex *16 zpars(*)
      real *8 pi
      data pi/3.14159265358979323846264338327950288419716939937510d0/
c

      ndim=ipars(1)
      ng=ipars(2)
      delta=dpars(51)

      ifpgh=ipars(10)

c-----------------------
      do ind=1,nd
         pot(ind)=0.0d0
         if (ifpgh.ge.2) then
            do k=1,ndim
               grad(ind,k)=0.0d0
            enddo
         endif
         if (ifpgh.ge.3) then
            do k=1,ndim*(ndim+1)/2
               hess(ind,k)=0.0d0
            enddo
         endif
      enddo
c-----------------------
      do ind=1,nd
         do i=1,ng
            idp = (i-1)*5
            sigma = dpars(idp+4)
            dc = sigma
            d = delta
         
            do k=1,ndim
               c=dpars(idp+k)
               x=targ(k)+0.5d0
c
               gf(k)=sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     1             *(-erf((-dc-d+x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc))
     2             +erf((x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)))
     3             /dsqrt(((dc+d)/d/dc))

               if (ifpgh.ge.2) then
                  arg1 = (-dc-d+x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)
                  darg = 1.0d0/d/dsqrt((dc+d)/d/dc)
                  arg2 = (x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)

                  gfp(k)=gf(k)*(-2.0d0*(x-c)/(dc+d)) +
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2                *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3                darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc))
               endif

               if (ifpgh.ge.3) then
                  gfpp(k)=gfp(k)*(-2.0d0*(x-c)/(dc+d)) +
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2                *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3                darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc))

                  gfpp(k)=gfpp(k)+gf(k)*(-2.0d0/(dc+d)) +
     1                ( (-2.0d0*(x-c)/(dc+d))*
     1                sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2                *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3                darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc)))+
     1                (sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     2                *darg*(2*arg1*dexp(-arg1*arg1) - 
     1                2*arg2*dexp(-arg2*arg2))
     2                *(-dexp(-arg1*arg1) +dexp(-arg2*arg2))*
     3                darg*(2.0d0/sqrt(pi))/dsqrt(((dc+d)/d/dc)))
               endif
            enddo

            str=dpars(idp+5)
            if (ndim.eq.1) then
               pot(ind)=pot(ind)+str*gf(1)
               if (ifpgh.ge.2) grad(ind,1)=grad(ind,1)+str*gfp(1)
               if (ifpgh.ge.3) hess(ind,1)=hess(ind,1)+str*gfpp(1)
            elseif (ndim.eq.2) then
               pot(ind)=pot(ind)+str*gf(1)*gf(2)

               if (ifpgh.ge.2) then
                  grad(ind,1)=grad(ind,1)+str*gfp(1)*gf(2)
                  grad(ind,2)=grad(ind,2)+str*gf(1)*gfp(2)
               endif

               if (ifpgh.ge.3) then
                  hess(ind,1)=hess(ind,1)+str*gfpp(1)*gf(2)
                  hess(ind,2)=hess(ind,2)+str*gfp(1)*gfp(2)
                  hess(ind,3)=hess(ind,3)+str*gf(1)*gfpp(2)
               endif
            elseif (ndim.eq.3) then
               pot(ind)=pot(ind)+str*gf(1)*gf(2)*gf(3)

               if (ifpgh.ge.2) then
                  grad(ind,1)=grad(ind,1)+str*gfp(1)*gf(2)*gf(3)       
                  grad(ind,2)=grad(ind,2)+str*gf(1)*gfp(2)*gf(3)
                  grad(ind,3)=grad(ind,3)+str*gf(1)*gf(2)*gfp(3)
               endif

               if (ifpgh.ge.3) then
                  hess(ind,1)=hess(ind,1)+str*gfpp(1)*gf(2)*gf(3)
                  hess(ind,2)=hess(ind,2)+str*gf(1)*gfpp(2)*gf(3)
                  hess(ind,3)=hess(ind,3)+str*gf(1)*gf(2)*gfpp(3)

                  hess(ind,4)=hess(ind,4)+str*gfp(1)*gfp(2)*gf(3)
                  hess(ind,5)=hess(ind,5)+str*gfp(1)*gf(2)*gfp(3)
                  hess(ind,6)=hess(ind,6)+str*gf(1)*gfp(2)*gfp(3)
               endif
            endif
         enddo
      enddo
      return
      end
c
c 
c
c 
      subroutine fperiodic(nd,xyz,dpars,zpars,ipars,f)
c     periodic right-hand-side function
c     
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars
      real *8 dpars(*),f(nd),xyz(*)
      integer nxyz(10)
      real *8 pi
      data pi/3.14159265358979323846264338327950288419716939937510D0/
c
      ndim=ipars(1)
      pi2=pi*2

      delta=dpars(51)
      c0=delta**(-ndim/2.0d0)
      
      do ind=1,nd
c         f(ind)=1.0d0
         f(ind)=c0
         do i=2,ndim,2
            f(ind)=f(ind)*cos(pi2*xyz(i)*ipars(1+i))
         enddo
         do i=1,ndim,2
            f(ind)=f(ind)*sin(pi2*xyz(i)*ipars(1+i))
         enddo
      enddo

c      do ind=1,nd
c         f(ind)=1.0d0
c         do i=1,ndim
c            f(ind)=f(ind)*sin(pi2*xyz(i)*ipars(1+i))
c         enddo
c      enddo

c     
c     
      return
      end
c
c
c
      subroutine upexact(nd,targ,dpars,zpars,ipars,
     1    pot,grad,hess)
c     exact solution for periodic rhs
      implicit real*8 (a-h,o-z)
      real*8 targ(*),pot(nd),grad(nd,*),hess(nd,*)
      real*8 dpars(*)
      real*8 cvals(10),svals(10)
      integer ipars(*),nxyz(10)
      complex *16 zpars(*)
      real *8 pi
      data pi/3.14159265358979323846264338327950288419716939937510d0/
c
      pisq=pi**2
      pi2=2*pi
      pi22=pi2**2
      
      ndim=ipars(1)
      delta=dpars(51)

      do i=1,ndim
         nxyz(i)=ipars(1+i)
      enddo

      ifpgh=ipars(10)
      
c-----------------------
      do ind=1,nd
         pot(ind)=0.0d0
         if (ifpgh.ge.2) then
            do k=1,ndim
               grad(ind,k)=0.0d0
            enddo
         endif
         if (ifpgh.ge.3) then
            do k=1,ndim*(ndim+1)/2
               hess(ind,k)=0.0d0
            enddo
         endif
      enddo

c      c0=(pi*delta)**(ndim/2.0d0)
      c0=pi**(ndim/2.0d0)

      xi2=0.0d0
      do i=1,ndim
         xi2=xi2+nxyz(i)**2
      enddo

      c1=c0*exp(-xi2*pisq*delta)

      do i=1,ndim
         svals(i)= sin(pi2*targ(i)*nxyz(i))
         cvals(i)= cos(pi2*targ(i)*nxyz(i))
      enddo
c-----------------------
      do ind=1,nd
         pot(ind)=c1
         do i=2,ndim,2
            pot(ind)=pot(ind)*cvals(i)
         enddo
         do i=1,ndim,2
c         do i=1,ndim
            pot(ind)=pot(ind)*svals(i)
         enddo

         if (ifpgh.ge.2) then
            do i=1,ndim
               grad(ind,i)=pi2*nxyz(i)*c1
               do j=2,ndim,2
                  if (j.ne.i) then
                     grad(ind,i)=grad(ind,i)*cvals(j)
                  elseif (j.eq.i) then
                     grad(ind,i)=-grad(ind,i)*svals(j)
                  endif
               enddo
               do j=1,ndim,2
c               do j=1,ndim
                  if (j.ne.i) then
                     grad(ind,i)=grad(ind,i)*svals(j)
                  elseif (j.eq.i) then
                     grad(ind,i)=grad(ind,i)*cvals(j)
                  endif
               enddo
            enddo
         endif

         if (ifpgh.ge.3) then
            if (ndim.eq.1) then
               hess(ind,1)=-pi22*nxyz(1)**2*pot(ind)

            elseif (ndim.eq.2) then
               cd = -pi22*pot(ind)
               dxy = -pi22*c1*cvals(1)*svals(2)
c               dxy = pi22*c1*cvals(1)*cvals(2)
               hess(ind,1)=cd*nxyz(1)**2
               hess(ind,2)=dxy*nxyz(1)*nxyz(2)
               hess(ind,3)=cd*nxyz(2)**2

            elseif (ndim.eq.3) then
               cd = -pi22*pot(ind)
               dxy = -pi22*c1*cvals(1)*svals(2)*svals(3)
               dxz =  pi22*c1*cvals(1)*cvals(2)*cvals(3)
               dyz = -pi22*c1*svals(1)*svals(2)*cvals(3)

               hess(ind,1)=cd*nxyz(1)**2
               hess(ind,2)=cd*nxyz(2)**2
               hess(ind,3)=cd*nxyz(3)**2
               
               hess(ind,4)=dxy*nxyz(1)*nxyz(2)
               hess(ind,5)=dxz*nxyz(1)*nxyz(3)
               hess(ind,6)=dyz*nxyz(2)*nxyz(3)
            endif
         endif
         
      enddo
      
      return
      end
c
c 
c
c     
      subroutine exact_pot(nd,targ,dpars,zpars,ipars,
     1    pot)
      implicit real*8 (a-h,o-z)
      real*8 targ(*),pot(nd)
      real*8 gf(10),dpars(*)
      integer ipars(*)
c
      one=1.0d0
      pi=4*atan(one)

      ndim=ipars(1)
      ng=ipars(2)
      delta=dpars(51)
      
c-----------------------
      do ind=1,nd
         pot(ind)=0.0d0
      enddo
c-----------------------
      do ind=1,nd
         do i=1,ng
            idp = (i-1)*5
            sigma = dpars(idp+4)
            dc = sigma
            d = delta
         
            do k=1,ndim
               c=dpars(idp+k)
               x=targ(k)+0.5d0
c
               gf(k)=sqrt(pi)/2.0d0*dexp(-(x-c)**2/(dc+d))
     1             *(-erf((-dc-d+x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc))
     2             +erf((x*dc+d*c)/d/dc/dsqrt((dc+d)/d/dc)))
     3             /dsqrt(((dc+d)/d/dc))
            enddo

            str=dpars(idp+5)
            if (ndim.eq.1) then
               pot(ind)=pot(ind)+str*gf(1)
            elseif (ndim.eq.2) then
               pot(ind)=pot(ind)+str*gf(1)*gf(2)
            elseif (ndim.eq.3) then
               pot(ind)=pot(ind)+str*gf(1)*gf(2)*gf(3)
            endif
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
      subroutine print_tree_matlab(ndim,itree,ltree,nboxes,centers,
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
      integer itree(ltree),ltree,nboxes,nlevels,iptr(*),ns,nt
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
      subroutine derr(vec1,vec2,n,erra)
      implicit real *8 (a-h,o-z)
      real *8 vec1(*),vec2(*)

      ra = 0
      erra = 0
      
      do i=1,n
         ra = ra + vec1(i)**2
c         if (ra .lt. abs(vec1(i))) ra=abs(vec1(i))
c         if (erra.lt.abs(vec1(i)-vec2(i))) then
c            erra=abs(vec1(i)-vec2(i))
c         endif
         erra = erra + (vec1(i)-vec2(i))**2
      enddo

c      erra=erra/ra

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











