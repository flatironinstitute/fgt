cc Copyright (C) 2020-2021: Leslie Greengard, Shidong Jiang, Manas Rachh
cc Contact: lgreengard@flatironinstitute.org
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

      subroutine fgt3d(nd,delta,eps,ns,sources,ifcharge,charge,
     1            ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of FGTs (same source and target locations, 
c                   different charge, dipole strengths)
c   delta         : Gaussian variance 
c   eps           : precision requested
c   ns            : number of sources
c   sources(3,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(3,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(3,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,3,*)    : gradients at the source locations
c   hess(nd,6,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient at the target locations
c   hesstarg(nd,6,*): hessian at the target locations
c
c   Note: hessians are in the order xx, yy, zz, xy, xz, yz
c      
      implicit none
c
cc      calling sequence variables
c 
      integer nd,dim
      real *8 eps,delta
      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 sources(3,ns),targ(3,nt)
      real *8 rnormal(3,ns)
      real *8 charge(nd,*),dipstr(nd,*)

      real *8 pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      real *8 pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer iptr(8)
      integer iper,nlmin,nlmax,ifunif
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer idivflag,nlevels,nboxes,ndiv,ndiv0,npwlevel
      integer ltree

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: rnormalsort(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: chargesort(:,:),dipstrsort(:,:)
      real *8, allocatable :: potsort(:,:),gradsort(:,:,:),
     1                             hesssort(:,:,:)
      real *8, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1                              hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer *8 lmptot

c
cc      temporary variables
c
      integer i,ilev,lmptmp,id,k,nhess
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      integer ibox,istart,iend
      real *8 omp_get_wtime,pps
      real *8 time1,time2,pi,done,pmax,bs0,cen0(3),bsize,pweps

      done = 1
      pi = atan(done)*4.0d0
      
      ifprint = 1

      dim=3
      call pts_tree_boxsize0(dim,delta,eps,sources,ns,targ,nt,
     1    npwlevel,bs0,cen0)
cccc      write(6,*) ' npwlevel',npwlevel

C
C     if box is too small as compared with delta, there is no need
C     to do anything - simply form local, then eval it.
c
      if (npwlevel .lt. -2) then
         call gndhlterms(dim,bs0,delta,eps,nlocal0)
         nlocal0 = nlocal0 + max(ifpgh,ifpghtarg) -1 
         write(6,*) ' nlocal0',nlocal0
      
         call cpu_time(time1)
C$      time1=omp_get_wtime()
         call fgt3d_large_delta(nd,delta,eps,ns,sources,
     1       ifcharge,charge,ifdipole,rnormal,dipstr,iper,
     2       nlocal0,cen0,
     3       ifpgh,pot,grad,hess,nt,targ,
     4       ifpghtarg,pottarg,gradtarg,hesstarg)
         call cpu_time(time2)
C$        time2=omp_get_wtime()
         if( ifprint .eq. 1 ) call prin2('time in fgt main=*',
     1       time2-time1,1)
         return
      endif     
c
c     need to fix the dependence of ndiv on eps
c   
c
      idivflag =0
c     ndiv0 is the maximum number of points per box above the cutoff level
c     it determines the speed of the algorithm when delta goes to zero.
c      if (eps.le.1d-8 .and. delta.le.1d-5) then
c         ndiv0=1000
c      else
c         ndiv0 = 250
c      endif
      ndiv0=100
c     ndiv is the maximum number of points per box at or below the cutoff level
c     it's determined by numerical experiments on finding the crossover point
c     between direct evaluation and the fast scheme.
      ndiv = ndiv0
c
      ifunif = 0
      iper = 0

c
c     call the tree memory management
c     code to determine number of boxes,
c     number of levels and length of tree
c
c     1. determine the maximum level - it seems that the best strategy for point FGT
c     is to stop refinement as long as the Hermite expansion length is less
c     than 20 for high precision calculation

cccc      npwlevel=npwlevel-1
      if (npwlevel .ge. 0) nlmax=npwlevel+1

cccc      write(6,*) ' nlmax',nlmax

      
      nlmin = 0
      call pts_tree_mem(dim,sources,ns,targ,nt,idivflag,
     1    ndiv,nlmin,nlmax,ifunif,iper,
     2    ndiv0,npwlevel,bs0,cen0,
     3    nlevels,nboxes,ltree) 

c 
cccc      write(6,*) ' nlevels',nlevels
cccc      write(6,*) ' nboxes',nboxes

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(dim,nboxes))
c
c     call the tree code
c
      call pts_tree_build(dim,sources,ns,targ,nt,idivflag,ndiv,
     1    nlmin,nlmax,ifunif,iper,nlevels,nboxes,
     2    ndiv0,npwlevel,bs0,cen0,
     3    ltree,itree,iptr,tcenters,boxsize)
cccc      write(6,*) ' boxsize',boxsize(0)
cccc      write(6,*) ' cutoff length', boxsize(0)/2**npwlevel/sqrt(delta)
cccc      write(6,*) ' nperbox',ns*4/(3*nboxes)

      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(dim,ns,sources,itree,ltree,nboxes,nlevels,
     1    iptr,tcenters,isrc,isrcse)
cccc      call prinf('isrcse=*',isrcse,20)

      call pts_tree_sort(dim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,itarg,itargse)
cccc      call prinf('itargse=*',itargse,nboxes)

      allocate(sourcesort(dim,ns))
      allocate(targsort(dim,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,ns),dipstrsort(nd,1))
        allocate(rnormalsort(dim,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns))
        allocate(rnormalsort(dim,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns))
        allocate(rnormalsort(dim,ns))
      endif

      nhess = dim*(dim+1)/2
      if(ifpgh.eq.1) then
         allocate(potsort(nd,ns),gradsort(nd,dim,1),
     1    hesssort(nd,nhess,1))
      else if(ifpgh.eq.2) then
         allocate(potsort(nd,ns),gradsort(nd,dim,ns),
     1       hesssort(nd,nhess,1))
      else if(ifpgh.eq.3) then
         allocate(potsort(nd,ns),gradsort(nd,dim,ns),
     1       hesssort(nd,nhess,ns))
      else
         allocate(potsort(nd,1),gradsort(nd,dim,1),
     1       hesssort(nd,nhess,1))
      endif
c      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,dim,1),
     1     hesstargsort(nd,nhess,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,dim,nt),
     1      hesstargsort(nd,nhess,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,dim,nt),
     1     hesstargsort(nd,nhess,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,dim,1),
     1     hesstargsort(nd,nhess,1))
      endif
c
c     initialize potentials,hessians,gradients
c
      if(ifpgh.eq.1) then
        do i=1,ns
          do id=1,nd
            potsort(id,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ns
          do id=1,nd
             potsort(id,i) = 0
             do k=1,dim
                gradsort(id,k,i) = 0
             enddo
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do id=1,nd
            potsort(id,i) = 0
            do k=1,dim
               gradsort(id,k,i) = 0
            enddo
            do k=1,nhess
               hesssort(id,k,i) = 0
            enddo
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,nt
          do id=1,nd
            pottargsort(id,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,nt
          do id=1,nd
            pottargsort(id,i) = 0
            do k=1,dim
               gradtargsort(id,k,i) = 0
            enddo
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do id=1,nd
            pottargsort(id,i) = 0
            do k=1,dim
               gradtargsort(id,k,i) = 0
            enddo
            do k=1,nhess
               hesstargsort(id,k,i) = 0
            enddo
          enddo
        enddo
      endif
c
c     compute the length of plane wave expansion
      
      bsize = 2*bs0/(2.0d0**(max(npwlevel,0)))
      pweps = eps
      if (nadd .gt. 2) pweps=pweps/10
      call gndpwterms(bsize,delta,pweps,pmax,npw)
      
cccc      call prinf(' nlocal =*',nlocal,nlevels+1)
cccc      call prinf(' ntermmax =*',ntermmax,1)
cccc      call prin2(' pmax =*',pmax,1)
      if (ifprint.eq.1) call prinf(' npw =*',npw,1)
c
c
c     reorder sources
c
      call dreorderf(dim,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1    call dreorderf(nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(nd,ns,dipstr,dipstrsort,isrc)
         call dreorderf(dim,ns,rnormal,rnormalsort,isrc)
      endif
c
cc     reorder targets
c
      call dreorderf(dim,nt,targ,targsort,itarg)
c
c     call the main FGT routine
c

c     Memory allocation is complete. 
c     Call main FGT routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call fgt3dmain(nd,dim,delta,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,rnormalsort,dipstrsort,
     $   nt,targsort,
     $   itree,ltree,iptr,nlevels,npwlevel,ndiv,
     $   nboxes,iper,boxsize,tcenters,itree(iptr(1)),
     $   isrcse,itargse,pmax,npw,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) then
         call prin2('time in fgt main=*',time2-time1,1)
         pps=(ns*ifpgh+nt*ifpghtarg+0.0d0)/(time2-time1)
         call prin2('points per sec=*',pps,1)
      endif

c
c     resort the output arrays in input order
c
      if(ifpgh.eq.1) then
        call dreorderi(nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderiv(nd,dim,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderiv(nd,dim,ns,gradsort,grad,isrc)
        call dreorderiv(nd,nhess,ns,hesssort,hess,isrc)
      endif

      if(ifpghtarg.eq.1) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,dim,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,dim,nt,gradtargsort,gradtarg,itarg)
        call dreorderiv(nd,nhess,nt,hesstargsort,hesstarg,itarg)
      endif

      return
      end
c
c
c
c
c
      subroutine fgt3dmain(nd,dim,delta,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,rnormalsort,dipstrsort,
     $     ntarget,targetsort,
     $     itree,ltree,iptr,nlevels,npwlevel,ndiv,
     $     nboxes,iper,boxsize,centers,laddr,
     $     isrcse,itargse,pmax,npw,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg)
c
c
c   the FGT in R^3: evaluate all pairwise particle
c   interactions
c   and interactions with targets using PW expansions.
c
c
c   \phi(x_i) = \sum_{j\ne i} charge_j e^{-|x_i-x_j|^2/delta)
c   + dipstr_j 2*(x_i - x_j)\cdot rnormal/delta * e^{-|x_i-x_j|^2/delta)
c
c
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FGT precision requested
c
c   nsource:     integer:  number of sources
c   sourcesort: real *8 (3,ns):  source locations
c
c   ifcharge:  charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c   chargesort: complex *16 (nsource): charge strengths
c
c   ifdipole:  dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c   dipstrsort: complex *16 (nsource): dipole strengths
c   ntarget: integer:  number of targets
c   targetsort: real *8 (3,ntarget):  target locations
c   iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree2d.f
c
c   ltree    in: integer *8
c            length of tree
c
c    iptr in: integer(8)
c             iptr is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     npwlevel in: integer
c             cutoff level at which the X expansion is valid
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c     iper    in: integer
c             flag for periodic implementation
c
c     centers in: real *8(3,nboxes)
c                 array containing the centers of all the boxes
c
c     isrcse in: integer(2,nboxes)
c               starting and ending location of sources in ibox
c                in sorted list of sources
c
c     itargse in: integer(2,nboxes)
c               starting and ending location of targets in ibox
c                in sorted list of sources
c
c     pmax    in:  cutoff limit in the planewave expansion
c     npw     in:  length of planewave expansions
c
c     ifpgh  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at sources.
c             ifpgh = 1, only potentials will be evaluated
c             ifpgh = 2, potentials/gradients will be evaluated
c             ifpgh = 3, potentials/gradients/hessians will be evaluated
c
c     ifpghtarg  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at targets.
c             ifpghtarg = 1, only potentials will be evaluated
c             ifpghtarg = 2, potentials/gradients will be evaluated
c             ifpghtarg = 3, potentials/gradients/hessians will be evaluated
c
c   OUTPUT
c
c   pot: potential at the source locations
c   grad: gradient at the source locations
c   hess: gradient at the source locations
c  
c   pottarg: potential at the target locations
c   gradtarg: gradient at the target locations
c   hesstarg: gradient at the target locations
c------------------------------------------------------------------

      implicit none

      integer nd

c     our fortran-header, always needed
      include '/home/shidong/finufft/include/finufft.fh'
      
      integer dim
      integer iper
      integer nsource,ntarget
      integer ndiv,nlevels,npwlevel
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      real *8 eps,delta
      real *8 sourcesort(dim,nsource)
      real *8 rnormalsort(dim,nsource)
      real *8 chargesort(nd,*)
      real *8 dipstrsort(nd,*)
      real *8 targetsort(dim,ntarget)
      real *8 pot(nd,*)
      real *8 grad(nd,dim,*)
      real *8 hess(nd,dim*(dim+1)/2,*)
      real *8 pottarg(nd,*)
      real *8 gradtarg(nd,dim,*)
      real *8 hesstarg(nd,dim*(dim+1)/2,*)
c
      real *8 pmax
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(dim,*)
c
      integer laddr(2,0:nlevels)
      integer npw
      integer iptr(8)
      integer ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      real *8 boxsize(0:nlevels)
c
c     temp variables
      integer i,j,k,l,idim,ip,isx,ii,j1,j2
      integer ibox,jbox,ilev,klev,npts,nptssrc,nptstarg,nptsj,nptstmp
      integer nchild,ncoll,nb
      integer isep, mnbors

      integer mnlist1,mnlist2,mnlist3,mnlist4,mnlistpw
      integer nlist1
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist1s(:)

      integer, allocatable :: nlistpw(:), listpw(:,:)
      
c
      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt
      integer ier
      integer *8 lmptot
      
      integer ifprint

      integer nexp,nmax,ncutlevbox,nhess,mc,ind,iperiod
      integer ncdir,ncddir,nddir,npdir,ngdir,nhdir
      
      integer nn,jx,jy,jz,nlevstart,nlevend

      real *8 d,time1,time2,omp_get_wtime
      real *8 dx,dy,dz,bs0
      real *8 tt1,tt2,xmin,xmin2,t1,t2,dt,dtt
      real *8 dmax,hpw
      
      integer, allocatable :: ifhung(:),ifpwexp(:),iaddr(:,:)
      integer ndirect
      
      real *8, allocatable :: rmlexp(:)
      real *8, allocatable :: ws(:),ts(:)
      real *8, allocatable :: wnufftcd(:,:),wnufftgh(:,:)

      real *8, allocatable :: targtmp(:,:),pottmp(:,:)
      real *8, allocatable :: gradtmp(:,:,:),hesstmp(:,:,:)
      
      complex *16, allocatable :: wpwshift(:,:)

c     this is what you use as the "opaque" ptr to ptr to finufft_plan...
      integer *8 fftplan,fftplantarg
c     this is how you create the options struct in fortran...
      type(finufft_opts) opts
c     or this is if you want default opts, make a null pointer...
      type(finufft_opts), pointer :: defopts => null()
      integer ttype,ntrans,ntranstarg,ncd,npgh,npghtarg,iflag
      integer *8 npw8
      integer *8, allocatable :: n_modes(:)
      
      double precision pi
      complex *16 eye

C     
      eye = dcmplx(0,1)
      iperiod=0
c      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c
      ifprint=1
      pi = 4*atan(1.0d0)

      bs0 = boxsize(0)
      mc = 2**dim
      mnbors=3**dim
      nhess=dim*(dim+1)/2

      if (npwlevel.ge.0 .and. npwlevel.le.nlevels) then
         ncutlevbox=laddr(2,npwlevel)-laddr(1,npwlevel)+1
cccc         write(6,*) ' n per box on the cutoff level',nsource/ncutlevbox      
      endif
c
c     compute list info
c
      call gt3d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1    iptr,centers,
     2    boxsize,iper,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call gt3d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,
     2    itree(iptr(4)),itree(iptr(5)),
     1    boxsize,laddr,
     2    mnlistpw,nlistpw,listpw)
cccccc
      nlevend=nlevels
      if (npwlevel.le.nlevels) nlevend=npwlevel

c
c
c     check whether we need to create and evaluate the PW expansion for boxes 
c     at the cutoff level
      allocate(ifpwexp(nboxes))
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
         do i=1,nboxes
            ifpwexp(i)=0
         enddo

         do ilev=npwlevel,npwlevel
            do ibox=itree(2*ilev+1),itree(2*ilev+2)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
               if (nlistpw(ibox).gt.0 .or. npts.gt.ndiv) then
                  ifpwexp(ibox)=1
               endif
            enddo
         enddo
      endif
      
c
c     allocate memory need by multipole, local expansions at all levels
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c
c       ... allocate iaddr and temporary arrays
c
      allocate(iaddr(2,nboxes))
c     
c     irmlexp is pointer for workspace need by various expansions.
c
      call gndmpalloc(nd,dim,itree,iaddr,
     1    nlevels,npwlevel,ifpwexp,lmptot,npw)
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
c
      allocate(rmlexp(lmptot),stat=ier)
      do i=1,lmptot
         rmlexp(i)=0
      enddo
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      ndirect = 0
      do ilev = 0,nlevend
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
c     Check if the current box is a nonempty leaf box            
            if(nchild.eq.0.and.npts.gt.0.and.npts.le.ndiv) then
               ifhung(ibox) = 1
               ndirect = ndirect + 1
            endif
         enddo
      enddo
      if (ifprint.eq.1) call prinf('# of direct evaluation boxes=*',
     1    ndirect,1)
      
      do i=1,10
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo

c     compute list info
c
      mnbors = 27
      isep = 1
      call compute_mnlist1(dim,nboxes,nlevels,itree(iptr(1)),
     1  centers,boxsize,itree(iptr(3)),itree(iptr(4)),
     2  itree(iptr(5)),isep,itree(iptr(6)),
     2  itree(iptr(7)),iper,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1s(nboxes))

c     modified list1 for direct evaluation
      call compute_modified_list1(nlevels,npwlevel,
     1  nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),
     2  itree(iptr(5)),isep,itree(iptr(6)),mnbors,
     3  itree(iptr(7)),iper,nlist1s,mnlist1,list1)

c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800
c      
c     get planewave nodes and weights
      allocate(ws(-npw/2:npw/2-1))
      allocate(ts(-npw/2:npw/2-1))
      call get_pwnodes(pmax,npw,ws,ts)
      hpw = ts(1)

c     compute translation matrices for PW expansions
      xmin  = boxsize(npwlevel)/sqrt(delta)

      nmax = 1

      nexp = npw**dim

      allocate(wpwshift(nexp,(2*nmax+1)**dim))
      call gnd_mk_translation_matrices(dim,xmin,npw,ts,nmax,
     1    wpwshift)

      allocate(wnufftcd(nexp,dim+1),wnufftgh(nexp,dim+nhess))
      call nufft_weights(dim,npw,ws,ts,nexp,wnufftcd,wnufftgh)

c     xmin is used in shiftpw subroutines to
c     determine the right translation matrices
c      
      xmin  = boxsize(npwlevel)
c
c
      ncdir = ndiv
      nddir = 750
      ncddir = 1000

      npdir = ndiv
      ngdir = 900
      nhdir = 1200

      if (ifprint.eq.1) call prinf('laddr=*',laddr,2*(nlevels+1))
      
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, form multipole pw expansions at the cutoff level
c       
c     mandatory parameters to FINUFFT guru interface... (ttype = trans type)
      ttype = 1
      ncd=0
      if (ifcharge.eq.1) ncd=ncd+1
      if (ifdipole.eq.1) ncd=ncd+dim
      ntrans = nd*ncd
      iflag = -1
      npw8 = npw
      allocate(n_modes(dim))
      do k=1,dim
         n_modes(k) = npw8
      enddo
      
      call finufft_default_opts(opts)
      opts%fftw = 0
cccc      opts%upsampfac = 1.1d0

      call finufft_makeplan(ttype,dim,n_modes,iflag,ntrans,
     $     eps,fftplan,opts,ier)

      do 1100 ilev = npwlevel,npwlevel
ccc         nb=0
ccc         dt=0
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1 
c           Check if current box needs to form pw exp         
            if(npts.gt.ndiv) then
ccc   nb=nb+1
ccc   call cpu_time(t1)
c     form the pw expansion
               call gnd_formpw(nd,dim,delta,eps,sourcesort(1,istart),
     1             npts,ifcharge,chargesort(1,istart),ifdipole,
     2             rnormalsort(1,istart),dipstrsort(1,istart),
     3             centers(1,ibox),hpw,nexp,wnufftcd,
     4             rmlexp(iaddr(1,ibox)),fftplan)
cccc  call cpu_time(t2)
cccc  dt=dt+t2-t1
c     copy the multipole PW exp into local PW exp
c     for self interaction 
               call gndcopypwexp_vec(nd,nexp,rmlexp(iaddr(1,ibox)),
     1             rmlexp(iaddr(2,ibox)))
            endif
         enddo
C$OMP END PARALLEL DO 
 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f6.2)
cccc  write(6,111) ilev,nb,dt
c     end of ilev do loop
 1100 continue

      call finufft_destroy(fftplan,ier)
      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1


      
      if(ifprint.ge.1)
     $    call prinf('=== Step 2 (mp to loc) ===*',i,0)
c      ... step 2, convert multipole pw expansions into local
c       pw expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = npwlevel,npwlevel
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          npts = 0
          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = npts + iend-istart+1
          endif
c
          if(npts.gt.0) then
c
c           shift PW expansions
c
            do j=1,nlistpw(ibox)
              jbox=listpw(j,ibox)
              jstart = isrcse(1,jbox)
              jend = isrcse(2,jbox)
              nptsj = jend-jstart+1
              if (nptsj .gt. ndiv) then
cccc              if (nptsj .gt. ndiv.or. ifpwexp(jbox).eq.1) then
                 call gnd_find_pwshift_ind(dim,iperiod,centers(1,ibox),
     1               centers(1,jbox),bs0,xmin,nmax,ind)
                 call gnd_shiftpw(nd,nexp,rmlexp(iaddr(1,jbox)),
     1               rmlexp(iaddr(2,ibox)),wpwshift(1,ind))
              endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1300 continue
c

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2) = time2-time1

cccc      call prin2('timeinfo2=*',time2-time1,1)

      if(ifprint.ge.1)
     $    call prinf('=== step 3 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local pw expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()

c     mandatory parameters to FINUFFT guru interface... (ttype = trans type)
      ttype = 2
      if (ifpgh.eq.1) npgh=1
      if (ifpgh.eq.2) npgh=1+dim
      if (ifpgh.eq.3) npgh=1+dim+nhess
      ntrans = nd*npgh

      if (ifpghtarg.ne.ifpgh .and. ifpghtarg.gt.0) then
         npghtarg=0
         if (ifpghtarg.eq.1) npghtarg=1
         if (ifpghtarg.eq.2) npghtarg=1+dim
         if (ifpghtarg.eq.3) npghtarg=1+dim+nhess
         ntranstarg = nd*npghtarg
      endif
      iflag = 1
      
c     use default options
      call finufft_makeplan(ttype,dim,n_modes,iflag,ntrans,
     $     eps,fftplan,opts,ier)

      if (ifpghtarg.ne.ifpgh .and. ifpghtarg.gt.0) then
         call finufft_makeplan(ttype,dim,n_modes,iflag,ntranstarg,
     $       eps,fftplantarg,opts,ier)
      endif
      
      do 1500 ilev = npwlevel,npwlevel
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
         call cpu_time(t1)
         nb=0
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if (ifpwexp(ibox).eq.1) then
               istartt = itargse(1,ibox) 
               iendt = itargse(2,ibox)
               nptstarg = iendt-istartt + 1
               
               istarts = isrcse(1,ibox)
               iends = isrcse(2,ibox)
               nptssrc = iends-istarts+1

               nb=nb+1
               if (ifpghtarg.ne.ifpgh) then
c                 evaluate local expansion at targets
                  if (nptstarg.gt.0 .and. ifpghtarg.gt.0) then
                     call gnd_pweval(nd,dim,delta,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targetsort(1,istartt),nptstarg,ifpghtarg,
     3                   pottarg(1,istartt),gradtarg(1,1,istartt),
     4                   hesstarg(1,1,istartt),fftplantarg)
                  endif
c     
c                 evaluate local expansion at sources
                  if (nptssrc.gt.0) then
                     call gnd_pweval(nd,dim,delta,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   sourcesort(1,istarts),nptssrc,ifpgh,
     3                   pot(1,istarts),grad(1,1,istarts),
     4                   hess(1,1,istarts),fftplan)
                  endif
               else
c                 evaluate local expansion at targets
                  if (nptstarg.gt.0 .and. nptssrc.eq.0) then
                     call gnd_pweval(nd,dim,delta,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targetsort(1,istartt),nptstarg,ifpghtarg,
     3                   pottarg(1,istartt),gradtarg(1,1,istartt),
     4                   hesstarg(1,1,istartt),fftplan)
                  endif
c     
c                 evaluate local expansion at sources
                  if (nptssrc.gt.0 .and. nptstarg.eq.0) then
                     call gnd_pweval(nd,dim,delta,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   sourcesort(1,istarts),nptssrc,ifpgh,
     3                   pot(1,istarts),grad(1,1,istarts),
     4                   hess(1,1,istarts),fftplan)
                  endif
                  
c     evaluate local expansion at sources and targets together using one NUFFT
                  if (nptssrc.gt.0 .and. nptstarg.gt.0) then
                     nptstmp=nptssrc+nptstarg

                     allocate(targtmp(dim,nptstmp))
                     j=0
                     do i=istarts,iends
                        j=j+1
                        do k=1,dim
                           targtmp(k,j)=sourcesort(k,i)
                        enddo
                     enddo
                     do i=istartt,iendt
                        j=j+1
                        do k=1,dim
                           targtmp(k,j)=targetsort(k,i)
                        enddo
                     enddo

                     allocate(pottmp(nd,nptstmp))
                     do i=1,nptstmp
                        do ind=1,nd
                           pottmp(ind,i)=0
                        enddo
                     enddo
                     
                     if (ifpgh.ge.2) then
                        allocate(gradtmp(nd,dim,nptstmp))
                        do i=1,nptstmp
                        do k=1,dim
                        do ind=1,nd
                           gradtmp(ind,k,i)=0
                        enddo
                        enddo
                        enddo
                     else
                        allocate(gradtmp(nd,dim,1))
                     endif
                     
                     if (ifpgh.eq.3) then
                        allocate(hesstmp(nd,nhess,nptstmp))
                        do i=1,nptstmp
                        do k=1,nhess
                        do ind=1,nd
                           hesstmp(ind,k,i)=0
                        enddo
                        enddo
                        enddo
                     else
                        allocate(hesstmp(nd,nhess,1))
                     endif
                     
                     call gnd_pweval(nd,dim,delta,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targtmp,nptstmp,ifpgh,
     3                   pottmp,gradtmp,hesstmp,fftplan)
                     j=0
                     do i=istarts,iends
                        j=j+1
                        do ind=1,nd
                           pot(ind,i)=pot(ind,i)+pottmp(ind,j)
                        enddo
                     enddo
                     do i=istartt,iendt
                        j=j+1
                        do ind=1,nd
                           pottarg(ind,i)=pottarg(ind,i)+pottmp(ind,j)
                        enddo
                     enddo

                     if (ifpgh.ge.2) then
                        j=0
                        do i=istarts,iends
                           j=j+1
                           do k=1,dim
                           do ind=1,nd
                              grad(ind,k,i)=grad(ind,k,i)
     1                            +gradtmp(ind,k,j)
                           enddo
                           enddo
                        enddo
                        do i=istartt,iendt
                           j=j+1
                           do k=1,dim
                           do ind=1,nd
                              gradtarg(ind,k,i)=gradtarg(ind,k,i)
     1                            +gradtmp(ind,k,j)
                           enddo
                           enddo
                        enddo
                     endif
                     
                     if (ifpgh.eq.3) then
                        j=0
                        do i=istarts,iends
                           j=j+1
                           do k=1,nhess
                           do ind=1,nd
                              hess(ind,k,i)=hess(ind,k,i)
     1                            +hesstmp(ind,k,j)
                           enddo
                           enddo
                        enddo
                        do i=istartt,iendt
                           j=j+1
                           do k=1,nhess
                           do ind=1,nd
                              hesstarg(ind,k,i)=hesstarg(ind,k,i)
     1                            +hesstmp(ind,k,j)
                           enddo
                           enddo
                        enddo
                     endif
                     deallocate(targtmp,pottmp,gradtmp,hesstmp)
                  endif
               endif
            endif
         enddo
         call cpu_time(t2)
 222     format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f6.2)
         write(6,222) ilev,nb,t2-t1
C     $OMP END PARALLEL DO        
 1500 continue

      call finufft_destroy(fftplan,ier)
      if (ifpghtarg.ne.ifpgh) call finufft_destroy(fftplantarg,ier)
      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(3) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 4 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
      dmax = log(1.0d0/eps)*delta
      
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(istarts,iends,nptssrc,nptstarg)
C$OMP$SCHEDULE(DYNAMIC)  
         do jbox = laddr(1,ilev),laddr(2,ilev)
c        jbox is the source box            
            if (ifhung(jbox) .eq. 1 ) then
               
               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)
               
               nlist1 = nlist1s(jbox)
               do i=1,nlist1
cccc              ibox is the target box
                  ibox = list1(i,jbox)

                  istarts = isrcse(1,ibox)
                  iends = isrcse(2,ibox)
                  nptssrc = iends-istarts + 1
                  
                  istartt = itargse(1,ibox)
                  iendt = itargse(2,ibox)
                  nptstarg = iendt-istartt + 1
                  if (nptstarg.gt.0) then
                     call fgtpart_direct(nd,dim,delta,dmax,
     1                   jstart,jend,istartt,iendt,sourcesort,
     2                   ifcharge,chargesort,
     3                   ifdipole,rnormalsort,dipstrsort,targetsort,
     4                   ifpghtarg,pottarg,gradtarg,hesstarg)
                  endif
                  if (nptssrc.gt.0) then
                     call fgtpart_direct(nd,dim,delta,dmax,
     1                   jstart,jend,istarts,iends,sourcesort,
     2                   ifcharge,chargesort,
     3                   ifdipole,rnormalsort,dipstrsort,sourcesort,
     4                   ifpgh,pot,grad,hess)
                  endif
                  
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(4) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,4)
      d = 0
      do i = 1,4
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
cccc      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
c
c
c
c------------------------------------------------------------------     
      subroutine fgtpart_direct(nd,dim,delta,dmax,istart,iend,
     $    jstart,jend,source,ifcharge,charge,
     2    ifdipole,rnormal,dipstr,
     $    targ,ifpgh,pot,grad,hess)
c--------------------------------------------------------------------
c     This subroutine adds the contribution due to sources
c     istart to iend in the source array to the fields at targets
c     jstart to jend in the target array.
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
c
c     delta        in: Gaussian variance
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: real *8(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                 dipole strengths at the source locations
c
c     targ        in: real *8(3,nt)
c                 target locations
c
c     ifpgh        in: Integer
c                  Flag for computing the potential/gradient/hessian.
c                  ifpgh = 1, only potential is computed
c                  ifpgh = 2, potential/gradient are computed
c                  ifpgh = 3, potential/gradient/hessian are computed
c
c------------------------------------------------------------
c     OUTPUT
c
c     pot          potential incremented at targets
c     grad         gradients incremented at targets
c     hess         Hessians  incremented at targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i,ntarg
        integer ifcharge,ifdipole
        integer nd
        integer dim
        integer ifpgh
c
        real *8 source(3,*)
        real *8 rnormal(3,*)
        real *8 charge(nd,*),dipstr(nd,*)
        real *8 targ(3,*),delta,eps,dmax
        real *8 pot(nd,*)
        real *8 grad(nd,3,*)
        real *8 hess(nd,6,*)
c
        
        ns = iend - istart + 1
        ntarg = jend-jstart+1
        
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             call gnd_directcp(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directcg(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then 
             call gnd_directch(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call gnd_directdp(nd,dim,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directdg(nd,dim,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call gnd_directdh(nd,dim,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif 

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call gnd_directcdp(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directcdg(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call gnd_directcdh(nd,dim,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart),
     3           hess(1,1,jstart))
          endif
        endif


c
        return
        end
c
c
c
c
c------------------------------------------------------------------    
      subroutine gndmpalloc(nd,dim,laddr,iaddr,
     1    nlevels,npwlevel,ifpwexp,lmptot,npw)
c     This subroutine determines the size of the array
c     to be allocated for multipole/local expansions
c
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     npwlevel    in: Integer
c                 cutoff level where the plane wave expansion is
c                 valid at or below ilev = npwlevel
c
c     npw         in: Integer
c                 Number of terms in the plane wave expansion
c
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer dim
      integer nlevels,npwlevel,npw,nd
      integer iaddr(2,*), laddr(2,0:nlevels), ifpwexp(*)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart,itmp2
c
      istart = 1
      if (npwlevel .gt. nlevels) then
         lmptot=0
         return
      endif
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      nn = npw**dim
      nn = nn*2*nd

      itmp=0
      do i = nlevstart,nlevstart
cccc         print *, 'nboxes at npwlevel=',laddr(2,i)-laddr(1,i)+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole PW expansion         
c
           if (ifpwexp(ibox).eq.1) then
              iaddr(1,ibox) = istart + itmp*nn
              itmp = itmp+1
           endif
         enddo
C$OMP END PARALLEL DO         
         istart = istart + itmp*nn
      enddo
c
      itmp2=0
      do i = nlevstart,nlevstart

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local PW expansion         
c
           if (ifpwexp(ibox).eq.1) then
              iaddr(2,ibox) = istart + itmp2*nn
              itmp2 = itmp2+1
           endif
         enddo
C$OMP END PARALLEL DO         
         istart = istart + itmp2*nn
      enddo
            
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
c
c
c
      subroutine fgt3d_large_delta(nd,delta,eps,ns,sources,
     1    ifcharge,charge,ifdipole,rnormal,dipstr,iper,
     2    nlocal,center,
     3    ifpgh,pot,grad,hess,nt,targ,
     4    ifpghtarg,pottarg,gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of FGTs (same source and target locations, 
c                   different charge, dipole strengths)
c   delta         : Gaussian variance 
c   eps           : precision requested
c   ns            : number of sources
c   sources(3,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(3,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   nlocal        : number of terms in the local expansion
c   center(3)     : the center of the bounding box
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(3,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,3,*)    : gradients at the source locations
c   hess(nd,6,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient at the target locations
c   hesstarg(nd,6,*): hessian at the target locations
c
      implicit none
c
cc      calling sequence variables
c 
      integer nd,ns,nt,i,j,k,ind
      integer iper,ifcharge,ifdipole
      integer ifpgh,ifpghtarg,nlocal

      real *8 omp_get_wtime
      real *8 eps,delta,bgf,bsize
      real *8 sources(3,ns),targ(3,nt)
      real *8 rnormal(3,ns)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 center(3)

      real *8 pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      real *8 pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)


      real *8, allocatable :: local(:,:,:,:)

      allocate(local(0:nlocal,0:nlocal,0:nlocal,nd))
      do ind=1,nd
         do i=0,nlocal
            do j=0,nlocal
               do k=0,nlocal
                  local(k,j,i,ind)=0
               enddo
            enddo
         enddo
      enddo
      
c
c     form local
c
      if (ifcharge.eq.1 .and. ifdipole.eq.0) then
         call gndformlc_vec(nd,delta,sources,ns,charge,center,
     1       nlocal,local)
      elseif (ifcharge.eq.0 .and. ifdipole.eq.1) then
         call  gndformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1       center,nlocal,local)
      elseif (ifcharge.eq.1 .and. ifdipole.eq.1) then
         call gndformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1       dipstr,center,nlocal,local)
      endif

c
c     local eval
c
c     targets
c      
      if(ifpghtarg.eq.1) then
         call gndlevalp_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg)
      endif
      if(ifpghtarg.eq.2) then
         call gndlevalg_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg)
      endif
      if(ifpghtarg.eq.3) then
         call gndlevalh_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg,hesstarg)
      endif
c
c     sources
c      
      if(ifpgh.eq.1) then
         call gndlevalp_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot)
      endif
      if(ifpgh.eq.2) then
         call gndlevalg_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad)
      endif
      if(ifpgh.eq.3) then
         call gndlevalh_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad,hess)
      endif
      
      return
      end
c
c
c
c
