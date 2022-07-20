ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

      subroutine pfgt(nd,dim,delta,eps,iperiod,bs0,cen0,
     1    ns,sources,
     1    ifcharge,charge,ifdipole,rnormal,dipstr,
     2    ifpgh,pot,grad,hess,nt,targ,
     3    ifpghtarg,pottarg,gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of FGTs (same source and target locations, 
c                   different charge, dipole strengths)
c   dim           : dimension of the space
c   delta         : Gaussian variance 
c   eps           : precision requested
c   ns            : number of sources
c   sources(dim,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(dim,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(dim,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,dim,*)    : gradients at the source locations
c   hess(nd,dim*(dim+1)/2,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,dim,*): gradient at the target locations
c   hesstarg(nd,dim*(dim+1)/2,*): hessian at the target locations
c
c   Note: hessians are in the order xx, yy, zz, xy, xz, yz for 3D
c     and xx, xy, yy for 2D
c      
      implicit none
c
cc      calling sequence variables
c 
      integer nd,dim
      real *8 eps,delta
      integer ns,nt,iperiod
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 sources(dim,ns),targ(dim,nt)
      real *8 rnormal(dim,ns)
      real *8 charge(nd,*),dipstr(nd,*)

      real *8 pot(nd,*),grad(nd,dim,*),hess(nd,dim*(dim+1)/2,*)
      real *8 pottarg(nd,*),gradtarg(nd,dim,*)
      real *8 hesstarg(nd,dim*(dim+1)/2,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer iptr(8)
      integer nlmin,nlmax,ifunif
      real *8, allocatable :: centers(:,:),boxsize(:)
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
cc      temporary variables
c
      integer i,ilev,lmptmp,id,k,nhess
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      integer ibox,istart,iend,ifplot
      real *8 omp_get_wtime,pps
      real *8 time1,time2,pmax,bs0,cen0(dim),bsize,pweps
      character *9 fname1
      character *8 fname2
      character *9 fname3
      
      ifprint = 1

      call pts_tree_boxsize0(dim,delta,eps,iperiod,sources,ns,targ,nt,
     1    npwlevel,bs0,cen0)
      if (ifprint.eq.1) call prinf('npwlevel=*',npwlevel,1)
C
c     divide on sources
c      
      idivflag =0
c     ndiv is the maximum number of points per box at or above the cutoff level
c     it's determined by numerical experiments on finding the crossover point
c     between direct evaluation and the fast scheme.
c      
c     need to fix the dependence of ndiv on eps
      if (dim.eq.1) then
         ndiv = 20
      elseif (dim.eq.2) then
         ndiv = 80
      elseif (dim.eq.3) then
         ndiv = 250
      endif
c
      ifunif = 0

c
c     call the tree memory management code to determine number of boxes,
c     number of levels and length of tree
c
      if (npwlevel .ge. 0) then
         nlmax=npwlevel+1
      else
         nlmax=0
      endif
      
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c     find the memory requirements for the tree
      nlmin = 0
      call pts_tree_mem(dim,sources,ns,targ,nt,idivflag,
     1    ndiv,nlmin,nlmax,ifunif,iperiod,
     2    npwlevel,bs0,cen0,
     3    nlevels,nboxes,ltree) 
c 
      if (ifprint.eq.1) call prinf('nlevels=*',nlevels,1)
c     memory allocation for the tree
      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(centers(dim,nboxes))
c
c     build the actual tree
c
      call pts_tree_build(dim,sources,ns,targ,nt,idivflag,ndiv,
     1    nlmin,nlmax,ifunif,iperiod,nlevels,nboxes,
     2    npwlevel,bs0,cen0,
     3    ltree,itree,iptr,centers,boxsize)
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      if( ifprint .eq. 1 ) then
         call prin2('time in tree building=*',time2-time1,1)
         pps=(ns*ifpgh+nt*ifpghtarg+0.0d0)/(time2-time1)
         call prin2('points per sec=*',pps,1)
      endif

c     plot the tree
      ifplot=0
      if (ifplot.eq.1) then
         fname1 = 'tree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         call print_tree2d_matlab(dim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,sources,nt,targ,
     2       fname1,fname2,fname3)
      endif
      
      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes))

      call cpu_time(time1)
C$    time1=omp_get_wtime()
c     sort source points to the tree
      call pts_tree_sort(dim,ns,sources,itree,ltree,nboxes,nlevels,
     1    iptr,centers,isrc,isrcse)
cccc      call prinf('isrcse=*',isrcse,20)

c     sort target points to the tree
      call pts_tree_sort(dim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   centers,itarg,itargse)
cccc      call prinf('itargse=*',itargse,nboxes)
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      if( ifprint .eq. 1 ) then
         call prin2('time in pts_tree_sort=*',time2-time1,1)
         pps=(ns*ifpgh+nt*ifpghtarg+0.0d0)/(time2-time1)
         call prin2('points per sec=*',pps,1)
      endif

c     allocate memory for sorted quantities
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
      npw=0
      if (npwlevel.le.nlevels) then
c         if (npwlevel.le.1.and.iperiod.eq.1) then
         if (npwlevel.le.1) then
            bsize=2*boxsize(max(npwlevel,0))
            call fgtpwterms(bsize,delta,eps,iperiod,pmax,npw)
            if (iperiod.eq.1) npw=npw+1
         else
            bsize=2*boxsize(npwlevel)
            call fgtpwterms(bsize,delta,eps,0,pmax,npw)
         endif
      endif
      if (npw.lt.4) npw=4
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
C$    time1=omp_get_wtime()
      call pfgtmain(nd,dim,delta,eps,iperiod,
     $   ifcharge,chargesort,ifdipole,rnormalsort,dipstrsort,
     $   ns,sourcesort,nt,targsort,
     $   nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     $   npwlevel,ndiv,pmax,npw,
     $   isrcse,itargse,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,hesstargsort)
      call cpu_time(time2)
C$    time2=omp_get_wtime()
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
      subroutine pfgtmain(nd,dim,delta,eps,iperiod,
     $    ifcharge,chargesort,ifdipole,rnormalsort,dipstrsort,
     $    nsource,sourcesort,ntarget,targetsort,
     $    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     $    npwlevel,ndiv,pmax,npw,
     $    isrcse,itargse,
     $    ifpgh,pot,grad,hess,
     $    ifpghtarg,pottarg,gradtarg,hesstarg)
c
c
c   the FGT in R^dim: evaluate all pairwise particle
c   interactions
c   and interactions with targets using PW expansions.
c
c
c   \phi(x_i) = \sum_{j\ne i} charge_j e^{-|x_i-x_j|^2/delta)
c   + dipstr_j 2*((x_i - x_j)\cdot rnormal(:,j))/delta * e^{-|x_i-x_j|^2/delta)
c
c     Note that the dipole term is the gradient w.r.t. the source point.
c      
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   dim:  dimension of the space
c
c   eps:  FGT precision requested
c
c   iperiod in: integer
c             flag for periodic implementation
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
c   nsource:     integer:  number of sources
c   sourcesort: real *8 (dim,ns):  source locations
c
c   ntarget: integer:  number of targets
c   targetsort: real *8 (dim,ntarget):  target locations
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree_nd.f
c
c   ltree    in: integer 
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
c             cutoff level at which the plane-wave expansion is valid
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c
c     centers in: real *8(dim,nboxes)
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
      implicit real *8 (a-h,o-z)
      integer nd,dim,iperiod,nsource,ntarget

c     our fortran-header, always needed
      include 'finufft.fh'

      integer ndiv,nlevels,npwlevel,ncutoff
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      real *8 eps,delta
      real *8 sourcesort(dim,nsource)
      real *8 rnormalsort(dim,*)
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
      real *8 timeinfo(20)
      real *8 centers(dim,*)
c
      integer npw
      integer iptr(8)
      integer ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      real *8 boxsize(0:nlevels)
c
c     temp variables
      integer, allocatable :: nlist1(:), list1(:,:)
      integer, allocatable :: nlistpw(:), listpw(:,:)
c
      integer *8 lmptot
      
      integer ifprint,itype

      integer, allocatable :: ifpwexp(:),iaddr(:,:)
      
      real *8, allocatable :: rmlexp(:)
      real *8, allocatable :: ws(:),ts(:)
      real *8, allocatable :: wnufftcd(:,:),wnufftgh(:,:)

      real *8, allocatable :: targtmp(:,:),pottmp(:,:)
      real *8, allocatable :: gradtmp(:,:,:),hesstmp(:,:,:)
      real *8 shifts(dim)
      
      complex *16, allocatable :: wpwshift(:,:)

c     this is what you use as the "opaque" ptr to ptr to finufft_plan...
      integer *8, allocatable :: fftplan(:),fftplantarg(:)
c     this is how you create the options struct in fortran...
      type(finufft_opts) opts
c     or this is if you want default opts, make a null pointer...
      type(finufft_opts), pointer :: defopts => null()
      integer ttype,ntrans,ntranstarg,ncd,npgh,npghtarg,iflag
      integer *8 npw8
      integer *8, allocatable :: n_modes(:)
      integer nthd,ithd
      integer omp_get_max_threads,omp_get_thread_num

      nthd = 1
C$    nthd = omp_get_max_threads()
      allocate(fftplan(nthd),fftplantarg(nthd))
c      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c
      ifprint=1

      ncutoff=max(npwlevel,0)
      
      bs0 = boxsize(0)
      mc = 2**dim
      mnbors=3**dim
      nhess=dim*(dim+1)/2

c     get planewave nodes and weights
      allocate(ws(-npw/2:npw/2-1))
      allocate(ts(-npw/2:npw/2-1))
      iperiod0=iperiod
      delta0=delta
      if ((npwlevel.le.1).and.(iperiod.eq.1)) then
         call get_periodic_pwnodes(bs0,delta,eps,npw,ws,ts)
         iperiod0=0
         delta0=1.0d0
c        in this case, periodic Gaussian kernel is used. Thus, the colleagues
c        should be computed as the free-space case to avoid double counting
         call computecoll(dim,nlevels,nboxes,itree(iptr(1)),
     1       boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2       itree(iptr(5)),iperiod0,itree(iptr(6)),itree(iptr(7)))
      else
c     itype = 1 is the trapezoidal rule, which is convenient for the use of
c     NUFFTs
         itype=1
         call get_pwnodes(itype,pmax,npw,ws,ts)
      endif
c     hpw is the step size in the Fourier space, needed in NUFFT calls
      hpw = ts(1)
c
c     compute list info
c
      call gnd_compute_mnlistpw(dim,nboxes,nlevels,ltree,itree,
     1    iptr,centers,boxsize,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call pfgt_compute_listpw(dim,npwlevel,nboxes,nlevels,
     1    ltree,itree,iptr,centers,boxsize,itree(iptr(1)),
     3    mnlistpw,nlistpw,listpw)      
cccc      if (ifprint.eq.1) call prinf('listpw=*',listpw,nboxes*mnlistpw)

      nlevend=nlevels
      if (npwlevel.le.nlevels) nlevend=npwlevel
c
c     check whether we need to create and evaluate planewave expansions 
c     for boxes
      allocate(ifpwexp(nboxes))
      call pfgt_find_pwexp_boxes(dim,npwlevel,nboxes,
     1    nlevels,ltree,itree,iptr,ifpwexp)
cccc      if (ifprint.eq.1) call prinf('ifpwexp=*',ifpwexp,nboxes)
c
c     compute list info
c
      isep = 1
      call compute_mnlist1(dim,nboxes,nlevels,itree(iptr(1)),
     1  centers,boxsize,itree(iptr(3)),itree(iptr(4)),
     2  itree(iptr(5)),isep,itree(iptr(6)),
     2  itree(iptr(7)),iperiod,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
c     modified list1 for direct evaluation
c     list1 of a childless source box ibox at ilev<=npwlevel
c     contains all childless target boxes that are neighbors of ibox
c     at or above npwlevel
      call pfgt_compute_modified_list1(dim,npwlevel,ifpwexp,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,iperiod0,
     2    mnlist1,nlist1,list1)

c
c     allocate memory need by multipole, local expansions at all levels
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c       ... allocate iaddr and temporary arrays
c
      allocate(iaddr(2,nboxes))
c     
      call pfgt_mpalloc(nd,dim,itree,iaddr,
     1    nlevels,npwlevel,ifpwexp,lmptot,npw)
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
c
      allocate(rmlexp(lmptot),stat=ier)
      do i=1,lmptot
         rmlexp(i)=0
      enddo
      
      do i=1,10
         timeinfo(i)=0
      enddo

c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800
c      
c     get planewave nodes and weights


c     compute translation matrices for PW expansions
      xmin  = boxsize(ncutoff)/sqrt(delta0)
c     translation only at the cutoff level
      nmax = 1
c     number of plane-wave modes
      nexp = npw**dim

      allocate(wpwshift(nexp,(2*nmax+1)**dim))
      call gnd_mk_full_translation_matrices(dim,xmin,npw,ts,nmax,
     1    wpwshift)

      allocate(wnufftcd(nexp,dim+1),wnufftgh(nexp,dim+nhess))
      call nufft_weights(dim,npw,ws,ts,nexp,wnufftcd,wnufftgh)

c     xmin is used in shiftpw subroutines to
c     determine the right translation matrices
c      
      xmin  = boxsize(ncutoff)
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()

c     
c       ... step 1, form multipole pw expansions at the cutoff level
c       
c     mandatory parameters to FINUFFT guru interface... (ttype = trans type)
c     type-1 for form multipole
      ttype = 1
c     charge -> 1 transform; dipole -> dim transforms; charge+dipole -> dim+1 transforms
      ncd=0
      if (ifcharge.eq.1) ncd=ncd+1
      if (ifdipole.eq.1) ncd=ncd+dim
c     ntrans is the total number of transforms
      ntrans = nd*ncd
c     minus sign for sources
      iflag = -1
c     integer *8 for finufft calls
      npw8 = npw
c     number of Fourier modes in each dimension
      allocate(n_modes(dim))
      do k=1,dim
         n_modes(k) = npw8
      enddo
      
      call finufft_default_opts(opts)
c     change the FFTW plan to MEASURE, which slows down makeplan, but speeds up
c     subsequent FFTW calls.
      opts%nthreads=1

      print *, "nthd=",nthd 
      do ithd=1,nthd
        call finufft_makeplan(ttype,dim,n_modes,iflag,ntrans,
     $      eps,fftplan(ithd),opts,ier)
      enddo


      do 1100 ilev = ncutoff,ncutoff
ccc         nb=0
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,ithd)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            ithd = 0
C$          ithd = omp_get_thread_num()
            ithd = ithd+1

c           Check if current box needs to form pw exp         
            if(npts.gt.ndiv) then
ccc               nb=nb+1
c              form the pw expansion
               call gnd_formpw(nd,dim,delta0,eps,sourcesort(1,istart),
     1             npts,ifcharge,chargesort(1,istart),ifdipole,
     2             rnormalsort(1,istart),dipstrsort(1,istart),
     3             centers(1,ibox),hpw,nexp,wnufftcd,
     4             rmlexp(iaddr(1,ibox)),fftplan(ithd))

c              copy the multipole PW exp into local PW exp
c              for self interaction 
               call gnd_copy_pwexp(nd,nexp,rmlexp(iaddr(1,ibox)),
     1             rmlexp(iaddr(2,ibox)))
            endif
         enddo
C$OMP END PARALLEL DO 
 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f6.2)
ccc         write(6,111) ilev,nb
c     end of ilev do loop
 1100 continue

      do i=1,nthd
        call finufft_destroy(fftplan(i),ier)
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1


      
      if(ifprint.ge.1)
     $    call prinf('=== Step 2 (mp to loc) ===*',i,0)
c      ... step 2, convert multipole pw expansions into local
c       pw expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = ncutoff,ncutoff
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
c           ibox is the target box
c           shift PW expansions
            do j=1,nlistpw(ibox)
               jbox=listpw(j,ibox)
c              jbox is the source box
               call gnd_find_pwshift_ind(dim,iperiod0,centers(1,ibox),
     1             centers(1,jbox),bs0,xmin,nmax,ind)
               call gnd_shiftpw(nd,nexp,rmlexp(iaddr(1,jbox)),
     1             rmlexp(iaddr(2,ibox)),wpwshift(1,ind))
            enddo
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

      do ithd=1,nthd
        call finufft_makeplan(ttype,dim,n_modes,iflag,ntrans,
     $       eps,fftplan(ithd),opts,ier)

        if (ifpghtarg.ne.ifpgh .and. ifpghtarg.gt.0) then
           call finufft_makeplan(ttype,dim,n_modes,iflag,ntranstarg,
     $         eps,fftplantarg(ithd),opts,ier)
        endif
      enddo
      
      do 1500 ilev = ncutoff,ncutoff
ccc         nb=0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istartt,iendt,istarts,iends,nptssrc,nptstarg)
C$OMP$PRIVATE(nptstmp,i,j,k,ind,targtmp,pottmp,gradtmp,hesstmp,ithd)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            if (ifpwexp(ibox).eq.1) then
               istartt = itargse(1,ibox) 
               iendt = itargse(2,ibox)
               nptstarg = iendt-istartt + 1
              
               istarts = isrcse(1,ibox)
               iends = isrcse(2,ibox)
               nptssrc = iends-istarts+1
               ithd = 0
C$             ithd = omp_get_thread_num()
               ithd = ithd+1
               
ccc               nb=nb+1
               if (ifpghtarg.ne.ifpgh) then
c                 evaluate local expansion at targets
                  if (nptstarg.gt.0 .and. ifpghtarg.gt.0) then
                     call gnd_pweval(nd,dim,delta0,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targetsort(1,istartt),nptstarg,ifpghtarg,
     3                   pottarg(1,istartt),gradtarg(1,1,istartt),
     4                   hesstarg(1,1,istartt),fftplantarg(ithd))
                  endif
c     
c                 evaluate local expansion at sources
                  if (nptssrc.gt.0) then
                     call gnd_pweval(nd,dim,delta0,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   sourcesort(1,istarts),nptssrc,ifpgh,
     3                   pot(1,istarts),grad(1,1,istarts),
     4                   hess(1,1,istarts),fftplan(ithd))
                  endif
               else
                  if (nthd.gt.1) then
c                 evaluate local expansion at targets
                  if (nptstarg.gt.0 .and. nptssrc.eq.0) then
                     call gnd_pweval(nd,dim,delta0,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targetsort(1,istartt),nptstarg,ifpghtarg,
     3                   pottarg(1,istartt),gradtarg(1,1,istartt),
     4                   hesstarg(1,1,istartt),fftplan(ithd))
                  endif
c     
c                 evaluate local expansion at sources
                  if (nptssrc.gt.0 .and. nptstarg.eq.0) then
                     call gnd_pweval(nd,dim,delta0,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   sourcesort(1,istarts),nptssrc,ifpgh,
     3                   pot(1,istarts),grad(1,1,istarts),
     4                   hess(1,1,istarts),fftplan(ithd))
                  endif

                  else
c                 evaluate local expansion at sources and targets 
c                 together using one NUFFT
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
                     
                     call gnd_pweval(nd,dim,delta0,eps,centers(1,ibox),
     1                   hpw,nexp,wnufftgh,rmlexp(iaddr(2,ibox)),
     2                   targtmp,nptstmp,ifpgh,
     3                   pottmp,gradtmp,hesstmp,fftplan(ithd))
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
            endif
         enddo
C$OMP END PARALLEL DO        
 222     format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f6.2)
ccc         write(6,222) ilev,nb,t2-t1
 1500 continue

      do i=1,nthd
        call finufft_destroy(fftplan(i),ier)
        if (ifpghtarg.ne.ifpgh .and. ifpghtarg.gt.0) then
           call finufft_destroy(fftplantarg(i),ier)
        endif
      enddo
      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(3) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 4 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      dmax = log(1.0d0/eps)*delta
      
ccc      nb=0
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,jstartt,jendt,jstarts,jends)
C$OMP$PRIVATE(ns,n1,nptssrc,nptstarg,shifts,npts)
C$OMP$SCHEDULE(DYNAMIC)  
         do jbox = itree(2*ilev+1),itree(2*ilev+2)
c        jbox is the target box            
            jstarts = isrcse(1,jbox)
            jends = isrcse(2,jbox)
            nptssrc = jends-jstarts+1

            jstartt = itargse(1,jbox)
            jendt = itargse(2,jbox)
            nptstarg = jendt-jstartt + 1

            npts = nptssrc+nptstarg
            
            n1 = nlist1(jbox)
ccc            if (n1.gt.0) nb=nb+1
            if (npts.gt.0 .and. n1.gt.0) then
            do i=1,n1
cccc           ibox is the source box
               ibox = list1(i,jbox)
               if (iperiod0.eq.1) call pfgt_find_local_shift(dim,
     1             centers(1,jbox),centers(1,ibox),bs0,shifts)
               
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               ns = iend-istart + 1
                  
               if (nptstarg.gt.0) then
                  call pfgt_direct(nd,dim,delta,dmax,iperiod0,shifts,
     1                istart,iend,jstartt,jendt,sourcesort,
     2                ifcharge,chargesort,
     3                ifdipole,rnormalsort,dipstrsort,targetsort,
     4                ifpghtarg,pottarg,gradtarg,hesstarg)
               endif
               if (nptssrc.gt.0) then
                  call pfgt_direct(nd,dim,delta,dmax,iperiod0,shifts,
     1                istart,iend,jstarts,jends,sourcesort,
     2                ifcharge,chargesort,
     3                ifdipole,rnormalsort,dipstrsort,sourcesort,
     4                ifpgh,pot,grad,hess)
               endif
                  
            enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue
ccc      print *, '# of direct source boxes=', nb
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(4) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,4)
      d = 0
      do i = 1,4
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
c
c
c
c
c------------------------------------------------------------------     
      subroutine pfgt_direct(nd,dim,delta,dmax,iperiod,shifts,
     $    istart,iend,jstart,jend,source,ifcharge,charge,
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
c     dim          in: integer
c                  dimension of the space
c
c     delta        in: Gaussian variance
c
c     dmax         in: maximum distance squared at which the Gaussian kernel is regarded as 0
c
c     iperiod      in: 0: free space; 1: periodic
c
c     shifts       in: real *8 (dim) the source center shifts when iperiod=1
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
c     source       in: real *8(dim,ns)
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
c     rnormal        in: real *8(dim,ns)
c                 dipole directions at the source locations
c     dipstr        in: real *8(nd,ns)
c                 dipole strengths at the source locations
c
c     targ        in: real *8(dim,nt)
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
      integer nd
      integer dim,iperiod
      integer istart,iend,jstart,jend,ns,ntarg
      integer ifcharge,ifdipole
      integer ifpgh
      integer i,j,k
c
      real *8 source(dim,*)
      real *8 rnormal(dim,*)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 targ(dim,*),delta,eps,dmax
      real *8 pot(nd,*)
      real *8 grad(nd,dim,*)
      real *8 hess(nd,dim*(dim+1)/2,*)
      real *8 shifts(*)
      real *8, allocatable:: sim(:,:)
c
        
      ns = iend - istart + 1
      ntarg = jend-jstart+1

      if (iperiod.eq.0) then
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

      else
        allocate(sim(dim,ns))
        do i=1,ns
           do k=1,dim
              sim(k,i)=source(k,istart+i-1)+shifts(k)
           enddo
        enddo
        
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             call gnd_directcp(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directcg(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then 
             call gnd_directch(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call gnd_directdp(nd,dim,delta,dmax,sim,ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directdg(nd,dim,delta,dmax,sim,ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call gnd_directdh(nd,dim,delta,dmax,sim,ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif 

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call gnd_directcdp(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call gnd_directcdg(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call gnd_directcdh(nd,dim,delta,dmax,sim,ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart),
     3           hess(1,1,jstart))
          endif
        endif
        deallocate(sim)
      endif
c
      return
      end
c
c
c
c
c------------------------------------------------------------------    
      subroutine pfgt_mpalloc(nd,dim,laddr,iaddr,
     1    nlevels,npwlevel,ifpwexp,lmptot,npw)
c     This subroutine determines the size of the array
c     to be allocated for multipole/local expansions
c
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     dim         in: integer
c                 dimension of the space
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     npwlevel    in: Integer
c                 cutoff level where the plane wave expansion is valid
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
c     the factor 2 is the (complex *16)/(real *8) ratio
      nn = nn*2*nd

      itmp=0
      do i = nlevstart,nlevstart
cccc C$OMP PARALLEL DO DEFAULT(SHARED)
cccc C$OMP$PRIVATE(ibox) REDUCTION(+:itmp)        
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole PW expansion         
c
           if (ifpwexp(ibox).eq.1) then
              iaddr(1,ibox) = istart + itmp*nn
              itmp = itmp+1
           endif
         enddo
cccc C$OMP END PARALLEL DO         
         istart = istart + itmp*nn
      enddo
c
      itmp2=0
      do i = nlevstart,nlevstart
cccc C$OMP PARALLEL DO DEFAULT(SHARED)
cccc C$OMP$PRIVATE(ibox) REDUCTION(+:itmp2)        
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local PW expansion         
c
           if (ifpwexp(ibox).eq.1) then
              iaddr(2,ibox) = istart + itmp2*nn
              itmp2 = itmp2+1
           endif
         enddo
cccc C$OMP END PARALLEL DO         
         istart = istart + itmp2*nn
      enddo
            
      lmptot = istart

      return
      end
c
