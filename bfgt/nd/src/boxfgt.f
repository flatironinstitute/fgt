c     
c    $Date$
c    $Revision$

      subroutine boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,pote,grade,hesse,timeinfo)
c
c
c       This code compute the Gauss transform for a collection of functions
c       defined on a tensor product grid of each leaf node in an adaptive tree
c 
c       input
c         nd - integer
c            number of right hand sides
c         ndim - integer
c            dimension of the underlying space
c         delta - double precision
c            Gaussian variance 
c         eps - double precision
c            precision requested
c         ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c         iperiod - integer
c            0: free space
c            1: periodic
c         norder - integer
c           order of expansions for input function value array
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c         fvals - double precision (nd,npbox,nboxes)
c           function values tabulated on a tensor grid in each leaf node
c         ifpgh   : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c         ifpghtarg: flag for computing pottarg/gradtarg/hesstarg
c                    ifpghtarg = 1, only potential is computed at targets
c                    ifpghtarg = 2, potential and gradient are 
c                    computed at targets
c                    ifpghtarg = 3, potential, gradient, and hessian are 
c                    computed at targets
c
c       input and output
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         centers - double precision (ndim,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     output:
c         pot - double precision (nd,npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes of the new tree
c         grad - double precision (nd,ndim,npbox,nboxes)
c            gradient of the volume potential on the tree structure 
c         hess - double precision (nd,ndim*(ndim+1)/2,npbox,nboxes)
c            hessian of the volume potential on the tree structure 
c            in 2d, the order is xx, xy, yy 
c            in 3d, the order is xx, yy, zz, xy, xz, yz
c         pote - double precision (nd,ntarg)
c            volume potential at targets
c         grade - double precision (nd,ndim,ntarg)
c            gradient of the volume potential at targets
c         hesse - double precision (nd,ndim*(ndim+1)/2,ntarg)
c            hessian of the volume potential at targets
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer nboxes,nlevels,ntarg,ifpgh,ifpghtarg
      integer iptr(8),ltree
      integer itree(ltree),norder,npbox
      real *8 targs(ndim,ntarg)
      real *8 fvals(nd,npbox,nboxes)

      real *8 pot(nd,npbox,nboxes)
      real *8 grad(nd,ndim,npbox,*)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,*)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,*)
      real *8 hesse(nd,ndim*(ndim+1)/2,*)

      real *8 centers(ndim,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(*)
c
cc     additional fgt variables
c
      integer *8 lmptot
      integer, allocatable :: iaddr(:,:)
      real *8, allocatable :: rmlexp(:)

      real *8, allocatable :: pot2(:,:,:)
      real *8, allocatable :: grad2(:,:,:,:)
      real *8, allocatable :: hess2(:,:,:,:)
c
cc      temporary variables
c
      integer npwlevel
      integer i,ilev,lmptmp,idim,ndim
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      real *8 dcutoff
      real *8 omp_get_wtime
      real *8 time1,time2,pi,done,pmax,bs0,bsize,pweps

      ifprint = 0

c     cutoff length      
      dcutoff = sqrt(delta*log(1.0d0/eps))
      if (ifprint.eq.1) call prin2(' dcutoff=*',dcutoff,1)
      
c     find the cutoff level
      npwlevel = nlevels+1
      do i=nlevels,0,-1
         if (boxsize(i).ge. dcutoff) then
            npwlevel=i
            exit
         endif
      enddo
c
      if (boxsize(nlevels).gt.dcutoff) npwlevel=nlevels+1
      if (boxsize(0) .le. dcutoff) npwlevel=0
      if(ifprint.eq.1) call prinf(' npwlevel =*',npwlevel,1)
      if(ifprint.eq.1) call prinf(' nlevels =*',nlevels,1)
c
c
c     compute the length of plane wave expansion
      npw=0
      if (npwlevel.ge.0.and.npwlevel.le.nlevels) then
c         if (npwlevel.le.1.and.iperiod.eq.1) then
         if (npwlevel.le.1) then
            bsize=boxsize(0)
            call fgtpwterms(bsize,delta,eps,iperiod,pmax,npw)
         else
            bsize=2*boxsize(npwlevel)
            call fgtpwterms(bsize,delta,eps,0,pmax,npw)
         endif
      endif
      if (ifprint.eq.1) call prinf(' npw =*',npw,1)
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
c
c     allocate memory need by multipole, local expansions at all levels
c     
c     iaddr is pointer for workspace need by various expansions.
c
      if (npwlevel.le.1.and.iperiod.eq.1) then
         call gnd_mpalloc(nd,ndim,itree,iaddr,
     1       nboxes,nlevels,0,lmptot,npw)
      else
         call gnd_mpalloc(nd,ndim,itree,iaddr,
     1       nboxes,nlevels,npwlevel,lmptot,npw)
      endif
      
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
      
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      allocate(rmlexp(lmptot),stat=ier)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in allocating rmlexp=*',
     1   time2-time1,1)
c

      nhess=ndim*(ndim+1)/2
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      if (ifpgh.ge.ifpghtarg) then
         call bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2       iaddr,rmlexp,npwlevel,pmax,npw,
     3       ifpgh,pot,grad,hess,ifnewtree,
     4       ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
      elseif (ifpghtarg.eq.3 .and. ifpgh.eq.1) then
         allocate(grad2(nd,ndim,npbox,nboxes))
         allocate(hess2(nd,nhess,npbox,nboxes))
         do i=1,nboxes
         do j=1,npbox
         do k=1,ndim
         do ind=1,nd
            grad2(ind,k,j,i) = 0
         enddo
         enddo
         enddo
         enddo
      
         do i=1,nboxes
         do j=1,npbox
         do k=1,nhess
         do ind=1,nd
            hess2(ind,k,j,i) = 0
         enddo
         enddo
         enddo
         enddo

         ifpgh2=ifpghtarg

         call bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2       iaddr,rmlexp,npwlevel,pmax,npw,
     3       ifpgh2,pot,grad2,hess2,ifnewtree,
     4       ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
      elseif (ifpghtarg.eq.3 .and. ifpgh.eq.2) then
         allocate(hess2(nd,nhess,npbox,nboxes))
         do i=1,nboxes
         do j=1,npbox
         do k=1,nhess
         do ind=1,nd
            hess2(ind,k,j,i) = 0
         enddo
         enddo
         enddo
         enddo
         ifpgh2=ifpghtarg
         call bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2       iaddr,rmlexp,npwlevel,pmax,npw,
     3       ifpgh2,pot,grad,hess2,ifnewtree,
     4       ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
      elseif (ifpghtarg.eq.2 .and. ifpgh.eq.1) then
         allocate(grad2(nd,ndim,npbox,nboxes))
         do i=1,nboxes
         do j=1,npbox
         do k=1,ndim
         do ind=1,nd
            grad2(ind,k,j,i) = 0
         enddo
         enddo
         enddo
         enddo
         ifpgh2=ifpghtarg
         call bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2       iaddr,rmlexp,npwlevel,pmax,npw,
     3       ifpgh2,pot,grad2,hess,ifnewtree,
     4       ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
      endif
      
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in box fgt main=*',
     1   time2-time1,1)
      
      return
      end
c
c
c
c
      subroutine bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2    iaddr,rmlexp,npwlevel,pmax,npw,
     3    ifpgh,pot,grad,hess,ifnewtree,
     4    ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
c
c
c       This code compute the Gauss transform for a collection of functions
c       defined on a tensor product grid of each leaf node in an adaptive tree
c 
c       input
c         nd - integer
c            number of right hand sides
c         ndim - integer
c            dimension of the underlying space
c         delta - double precision
c            Gaussian variance 
c         eps - double precision
c            precision requested
c         ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c         iperiod - integer
c            0: free space
c            1: doubly periodic
c         norder - integer
c           order of expansions for input function value array
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c         fvals - double precision (nd,npbox,nboxes)
c           function values tabulated on a tensor grid in each leaf node
c         iaddr - (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c         rmlexp - double precision, stores multipole and local PW expansions
c                  for each box below the cutoff level npwlevel
c         npwlevel - integer
c             cutoff level at which the PW expansion is valid
c         pmax - double precision,  cutoff limit in the planewave expansion
c         npw  - integer,  length of planewave expansions
      
c         ifpgh   : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c         ifpghtarg: flag for computing pottarg/gradtarg/hesstarg
c                    ifpghtarg = 1, only potential is computed at targets
c                    ifpghtarg = 2, potential and gradient are 
c                    computed at targets
c                    ifpghtarg = 3, potential, gradient, and hessian are 
c                    computed at targets
c
c       input and output
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         centers - double precision (ndim,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     output:
c         pot - double precision (nd,npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes of the new tree
c         grad - double precision (nd,ndim,npbox,nboxes)
c            gradient of the volume potential on the tree structure 
c         hess - double precision (nd,ndim*(ndim+1)/2,npbox,nboxes)
c            hessian of the volume potential on the tree structure 
c            in 2d, the order is xx, xy, yy 
c            in 3d, the order is xx, yy, zz, xy, xz, yz
c         pote - double precision (nd,ntarg)
c            volume potential at targets
c         grade - double precision (nd,ndim,ntarg)
c            gradient of the volume potential at targets
c         hesse - double precision (nd,ndim*(ndim+1)/2,ntarg)
c            hessian of the volume potential at targets
c
      implicit real *8 (a-h,o-z)
      integer nd,ndim
      real *8 delta,eps
      integer nboxes,nlevels,ntarg
      integer iptr(8),ltree
      integer itree(ltree),npbox
      real *8 fvals(nd,npbox,nboxes)
      real *8 targs(ndim,ntarg)

      real *8 pot(nd,npbox,nboxes)
      real *8 grad(nd,ndim,npbox,*)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,*)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,*)
      real *8 hesse(nd,ndim*(ndim+1)/2,*)

      real *8 boxsize(0:nlevels),centers(ndim,nboxes)
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(*)

c
cc        temporary variables
c
c     direction interaction list
      integer, allocatable :: nlist1(:),list1(:,:)
c     plane wave interaction list
      integer, allocatable :: nlistpw(:), listpw(:,:)
c     box flag array
      integer, allocatable :: ifpwexp(:)

      integer ndirect
      integer ixyz(ndim)

      real *8, allocatable :: rmask(:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: umat_nd(:,:,:)

      real *8 ws(200),ts(200)

c     plane wave m2m,l2l, and m2l translation operators
      complex *16, allocatable :: wpwshift(:,:)
      complex *16, allocatable :: wpwmsshift(:,:,:)

c     1d plane wave related tables
      complex *16, allocatable :: tab_poly2pw(:,:,:)
      complex *16, allocatable :: tab_pw2pot(:,:,:)
      complex *16, allocatable :: tab_pw2der(:,:,:)
      complex *16, allocatable :: tab_pw2dxx(:,:,:)

c     for extra targets
      real *8, allocatable :: coefsp(:,:,:)
      real *8, allocatable :: coefsg(:,:,:,:)
      real *8, allocatable :: coefsh(:,:,:,:)

c     for asymptotic calculations
      real *8, allocatable :: fcoefs(:,:,:)
      real *8, allocatable :: gvals(:,:,:,:)
      real *8, allocatable :: hvals(:,:,:,:)

c     for refinement and coarsening
      integer, allocatable :: ifrefine(:),irefinelev(:),imaxrefinelev(:)
      integer nblock(0:nlevels),ilevstart(0:nlevels+1)
      integer isgn(ndim,2**ndim)
      integer, allocatable :: nboxid(:)
      integer, allocatable :: ifdelete(:)
      
c     1d direct evaluation tables
      real *8, allocatable :: tab_loc(:,:,:,:)
      real *8, allocatable :: tabx_loc(:,:,:,:)
      real *8, allocatable :: tabxx_loc(:,:,:,:)
      integer, allocatable :: ind_loc(:,:,:,:)

      
      ifprint = 1

      done = 1
      pi = atan(done)*4

      bs0 = boxsize(0)
      mc = 2**ndim
      mnbors=3**ndim
      nhess=ndim*(ndim+1)/2
c
c
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 0 (precomputation) =====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c      
c     get planewave nodes and weights
      iperiod0=iperiod
      npwlevel0=npwlevel
      delta0=delta
      if (npwlevel.le.1.and.iperiod.eq.1) then
c     Method: in this case, the periodic Green's function admits an efficient
c     seprable Fourier series expansion, which is used in the entire computational
c     box. 1. compute planewave expansion coefficients for each leaf box; 2. merge
c     pw exp. all the way to the root box - this is a valid expansion in the entire
c     box; 3. split the pw exp. into its child boxes; 4. evaluate the pw exp. in each
c     leaf box.
c
c     Remark 1: (a) Since periodic Green's function is used, there is no need to change
c     colleague list to accommodate periodicity. Thus we need to recompute the colleagues
c     using iperiod=0. npwlevel is set to 0 internally since the plane wave expansions are
c     valid in the entire computational box. (b) delta should be set to 1.0d0
c     in mk_poly2pw_tables, mk_pw2pgh_tables and gnd_mk_merge_split_pw_matrices.
c     Hence the use of delta0 in those routines.
c
c     Remark 2: when npwlevel >= 2 and iperiod == 1, periodic Green's function will require
c     too many terms of Fourier series. So we will compute G_P[f] using free-space Green's
c     function, but tile up the whole space via identical copies of sources. However, since
c     npwlevel >=2, the Gaussian kernel is sharply peaked. Thus, we only need to wrap the 
c     original computational box by small boxes at level 2, modify the colleague list, list 1, 
c     listpw, and make corresponding changes to accommodate these modifications when 
c     calculating plane wave and local direct interactions. This is done in 
c     gnd_find_pwshift_ind and gnd_find_loctab_ind.
c         
         call get_periodic_pwnodes(bs0,delta,eps,npw,ws,ts)
         iperiod=0
         npwlevel=0
         nboxes0=itree(2*nlevels+2)
         delta0=1.0d0
         call computecoll(ndim,nlevels,nboxes0,itree(iptr(1)),
     1       boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2       itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
      else
         call get_pwnodes(pmax,npw,ws,ts)
      endif
      if(ifprint.eq.1) call prinf('npw=*',npw,1)
      
      allocate(xq(norder),umat(norder,norder),
     1    vmat(norder,norder),wts(norder))
     
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c      
      allocate(umat_nd(norder,norder,ndim))
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel

c     compute the number of leaf boxes
      nleafbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) nleafbox = nleafbox+1
        enddo
      enddo
      if(ifprint.eq.1) call prinf('nleafbox=*',nleafbox,1)
      
c     estimate number of direct evaluation boxes
      ndirect=0
      do ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
c           Check if the current box is a leaf box            
            if(nchild.eq.0) then
               ndirect = ndirect+1
            endif
         enddo
      enddo
      if(ifprint.eq.1) call 
     1    prinf('number of direct evaluation source boxes=*',ndirect,1)

c      nasym=3
c      sigma=(delta/boxsize(nlevels)**2)**nasym/8
c      if(ifprint.eq.1) call 
c     1    prin2('estimated asymptotic expansion error=*',sigma,1)
      
cccc      if (nleafbox.eq.ndirect .and. sigma.le.2*eps) goto 3000
cccc      goto 3000
      

      do i=1,20
         timeinfo(i)=0
      enddo
c
c
c       compute list info
c
c     check whether we need to create and evaluate planewave expansions 
c     for boxes
      allocate(ifpwexp(nboxes))
      call gnd_find_pwexp_boxes(ndim,npwlevel,nboxes,
     1    nlevels,ltree,itree,iptr,iperiod,ifpwexp)

      isep = 1
      call compute_mnlist1(ndim,nboxes,nlevels,itree(iptr(1)),centers,
     1    boxsize,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     2    isep,itree(iptr(6)),itree(iptr(7)),iperiod,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))

c     modified list1 for direct evaluation
c     list1 of a childless source box ibox at ilev<=npwlevel
c     contains all childless target boxes that are neighbors of ibox
c     at or above npwlevel
      call gnd_compute_modified_list1(ndim,npwlevel,ifpwexp,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,iperiod,
     2    mnlist1,nlist1,list1)

c     compute the tables converting Legendre polynomial expansion to potential
c     values, used in direct evaluation

      if (ifnewtree.eq.1) then
         mrefinelev=3
      else
         mrefinelev=0
      endif
      
      nloctab=2**(mrefinelev+1)*(mrefinelev+3)

      allocate(  tab_loc(norder,norder,-nloctab:nloctab,0:nlevels))
      allocate( tabx_loc(norder,norder,-nloctab:nloctab,0:nlevels))
      allocate(tabxx_loc(norder,norder,-nloctab:nloctab,0:nlevels))
      allocate(ind_loc(2,norder+1,-nloctab:nloctab,0:nlevels))

      nnodes=50
      do ilev = 0,nlevels
         call mk_loctab_all(eps,ipoly,norder,nnodes,delta,boxsize(ilev),
     1       mrefinelev,nloctab,tab_loc(1,1,-nloctab,ilev),
     2       tabx_loc(1,1,-nloctab,ilev),tabxx_loc(1,1,-nloctab,ilev),
     3       ind_loc(1,1,-nloctab,ilev))
      enddo

      
c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800


c     compute translation matrices for PW expansions
      nlevstart = max(npwlevel,0)

c     compute the tables converting Legendre polynomial expansion
c     to planewave expansion
      nnodes = 100
      allocate(tab_poly2pw(norder,npw,0:nlevels))

      do ilev=nlevstart,nlevels
         call mk_poly2pw_tables(norder,ipoly,npw,nnodes,ws,ts,delta0,
     1       boxsize(ilev),tab_poly2pw(1,1,ilev))
      enddo

c     compute the tables converting planewave expansions to potential values
      allocate(tab_pw2pot(npw,norder,0:nlevels))
      allocate(tab_pw2der(npw,norder,0:nlevels))
      allocate(tab_pw2dxx(npw,norder,0:nlevels))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pgh_tables(norder,npw,ts,xq,delta0,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev),
     1       tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev))  
      enddo

      nexp=(npw+1)/2
      do i=1,ndim-1
         nexp = nexp*npw
      enddo
      
      nmax = nlevels-max(npwlevel,0)      
      allocate(wpwmsshift(nexp,mc,nmax))
      xmin2 = boxsize(nlevels)/sqrt(delta0)/2
      call gnd_mk_merge_split_pw_matrices(ndim,xmin2,npw,ts,nmax,
     1    wpwmsshift)

      nmax = 1
      allocate(wpwshift(nexp,(2*nmax+1)**ndim))
      xmin  = boxsize(nlevstart)/sqrt(delta)
      call gnd_mk_pw_translation_matrices(ndim,xmin,npw,ts,nmax,
     1    wpwshift)
      
c     xmin is used in shiftpw subroutines to
c     determine the right translation matrices
c
      xmin  = boxsize(nlevstart)
      xmin2 = boxsize(nlevels)/2
      
c
c     compute list info
c
      call gnd_compute_mnlistpw(ndim,nboxes,nlevels,ltree,itree,
     1    iptr,centers,boxsize,iper,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call gnd_compute_listpw(ndim,npwlevel,nboxes,nlevels,
     1    ltree,itree,iptr,centers,boxsize,itree(iptr(1)),
     3    mnlistpw,nlistpw,listpw)      
c
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1) = time2-time1
      
      if (ifprint.eq.1) 
     1    call prin2('fgt precomputation time =*', time2-time1,1)
      
      if (ifprint.eq.1)
     1    call prinf('laddr=*',itree(iptr(1)),2*(nlevels+1))




      



      
c
c        step 1: convert function values to planewave expansions
c
    
      if(ifprint.eq.1) 
     1   call prinf("=== STEP 1 (values -> mp pwexps) ===*",i,0)
      
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do 1100 ilev = nlevels,nlevstart,-1
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=itree(2*ilev+1),itree(2*ilev+2)
          if (ifpwexp(ibox).eq.1) then
            nchild = itree(iptr(4)+ibox-1)
c           Check if current box is a leaf box            
            if(nchild.eq.0) then
c             form PW expansion directly
              call gnd_tens_prod_to_pw(ndim,nd,norder,fvals(1,1,ibox),
     1             npw,tab_poly2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
            endif
          endif
        enddo
C$OMP END PARALLEL DO
 1100 continue
      
      call cpu_time(time2)
C$       time2 = omp_get_wtime()
      timeinfo(2) = time2-time1









      
      if(ifprint .ge. 1)
     1    call prinf('=== STEP 2 (merge mp pwexps) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do 1200 ilev=nlevels-1,nlevstart,-1
         klev = nlevels - ilev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild .gt. 0) then
              call gnd_pwzero(nd,rmlexp(iaddr(1,ibox)),nexp)
           endif
           do i=1,nchild
              jbox = itree(iptr(5)+mc*(ibox-1)+i-1)
              call gnd_find_msshift_ind(ndim,centers(1,ibox),
     1            centers(1,jbox),k)
              call gnd_shiftpw(nd,nexp,rmlexp(iaddr(1,jbox)),
     1            rmlexp(iaddr(1,ibox)),wpwmsshift(1,k,klev))
           enddo
        enddo
C$OMP END PARALLEL DO    
 1200 continue
     
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1





      


      
      if(ifprint.ge.1)
     1    call prinf('=== STEP 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = nlevstart,nlevstart
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            call gnd_copy_pwexp(nd,nexp,rmlexp(iaddr(1,ibox)),
     1          rmlexp(iaddr(2,ibox)))
c           shift PW expansions
            do j=1,nlistpw(ibox)
               jbox=listpw(j,ibox)
               call gnd_find_pwshift_ind(ndim,iperiod,centers(1,ibox),
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
      timeinfo(4) = time2-time1

cccc      call prin2('timeinfo4=*',time2-time1,1)








      
      if(ifprint.ge.1)
     1    call prinf('=== STEP 4 (split loc pwexps) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      do 1400 ilev = nlevstart,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild.gt.0) then
              do i=1,nchild
                 jbox = itree(iptr(5)+mc*(ibox-1)+i-1)
                 call gnd_find_msshift_ind(ndim,centers(1,jbox),
     1               centers(1,ibox),k)
                 call gnd_shiftpw_loc(nd,nexp,rmlexp(iaddr(2,ibox)),
     1               rmlexp(iaddr(2,jbox)),wpwmsshift(1,k,nlevels-ilev))
              enddo
           endif
        enddo
C$OMP END PARALLEL DO        
 1400 continue
      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1








      
      if(ifprint.ge.1)
     1    call prinf('=== STEP 5 (eval loc pwexps) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      do 1500 ilev = nlevstart,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           if (ifpwexp(ibox).eq.1) then
             nchild = itree(iptr(4)+ibox-1)
             if(nchild.eq.0) then
               call gnd_pw2pgh(ndim,nd,norder,npw,
     1            rmlexp(iaddr(2,ibox)),tab_pw2pot(1,1,ilev),
     2            tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev),ifpgh,
     3            pot(1,1,ibox),grad(1,1,1,ibox),hess(1,1,1,ibox))
             endif
           endif
        enddo
C$OMP END PARALLEL DO        
 1500 continue
      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(6) = time2 - time1




      




      
 1800 continue

      if(ifprint .ge. 1)
     1     call prinf('=== STEP 6 (direct interactions) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nl1,bs,ixyz,jlev)
C$OMP$SCHEDULE(DYNAMIC)  
         bs = boxsize(ilev)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
c        ibox is the source box            
           nl1 = nlist1(ibox)
           do j=1,nl1
cccc       jbox is the target box
             jbox = list1(j,ibox)
             jlev = itree(iptr(2)+jbox-1)
             call gnd_find_loctab_ind(ndim,iperiod,
     1           centers(1,jbox),centers(1,ibox),bs,bs0,mrefinelev,ixyz)
cccc             call prinf('ixyz=*',ixyz,ndim)
             call gnd_tens_prod_to_pghloc(ndim,nd,norder,
     1         fvals(1,1,ibox),ifpgh,pot(1,1,jbox),
     2         grad(1,1,1,jbox),hess(1,1,1,jbox),nloctab,
     3         tab_loc(1,1,-nloctab,jlev),tabx_loc(1,1,-nloctab,jlev),
     4         tabxx_loc(1,1,-nloctab,jlev),ind_loc(1,1,-nloctab,jlev),
     5         ixyz)
            enddo
         enddo
C$OMP END PARALLEL DO         
 2000 continue
c
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(7) = time2-time1
      allocate(coefsp(nd,npbox,nboxes))
      if (ifnewtree .eq.0) then
         goto 5000
      else
         goto 4000
      endif
ccc      goto 3000








      
 3000 continue
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (asymptotic regime) =====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()

c     evaluate potential
c      call treedata_eval_pot_nd_asym(ndim,nd,delta,ipoly,nasym,
c     1    nlevels,itree,iptr,boxsize,norder,fvals,pot)

      itype=0
      iffast=1
      if (iffast.eq.0) goto 3200 
c     use spectral differentiation to evaluate gradient and hessian
      if (ifpgh.eq.2) then
         allocate(fcoefs(nd,npbox,nboxes))
         call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,pot,fcoefs,umat_nd)
         call treedata_evalg_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1       boxsize,norder,fcoefs,grad)
      endif
      
      if (ifpgh.eq.3) then
         allocate(fcoefs(nd,npbox,nboxes))
         call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,pot,fcoefs,umat_nd)
         call treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1       boxsize,norder,fcoefs,grad,hess)
      endif
      goto 3500
      
 3200 continue

c     slow version of asymptotics
c     also use asymptotic expansion to evaluate gradient and hessian
      if (ifpgh.eq.2) then
         allocate(fcoefs(nd,npbox,nboxes))
         allocate(gvals(nd,ndim,npbox,nboxes))

         call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,fvals,fcoefs,umat_nd)
         call treedata_evalg_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1       boxsize,norder,fcoefs,gvals)
         ng=nd*ndim
         call treedata_eval_pot_nd_asym(ndim,ng,delta,ipoly,nasym,
     1       nlevels,itree,iptr,boxsize,norder,gvals,grad)
      endif
      
      if (ifpgh.eq.3) then
         allocate(fcoefs(nd,npbox,nboxes))
         allocate(gvals(nd,ndim,npbox,nboxes))
         allocate(hvals(nd,nhess,npbox,nboxes))
         
         call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,fvals,fcoefs,umat_nd)
         call treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1       boxsize,norder,fcoefs,gvals,hvals)
         ng=nd*ndim
         nh=nd*nhess
         call treedata_eval_pot_nd_asym(ndim,ng,delta,ipoly,nasym,
     1       nlevels,itree,iptr,boxsize,norder,gvals,grad)
         call treedata_eval_pot_nd_asym(ndim,nh,delta,ipoly,nasym,
     1       nlevels,itree,iptr,boxsize,norder,hvals,hess)
      endif

 3500 continue
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(8) = time2-time1
      







      
 4000 continue
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (refine if needed) ===*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      mc=2**ndim
      call get_child_box_sign(ndim,isgn)

      allocate(ifrefine(nboxes))
      allocate(irefinelev(nboxes))
      allocate(imaxrefinelev(nboxes))

      iptype=2
      eta=1.0d0

      allocate(rmask(npbox))
      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0

c     compute the L_{iptype} norm of the pot by going throught each leaf box
      call treedata_lpnorm(ndim,iptype,ipoly,nd,nlevels,itree,
     1    iptr,boxsize,norder,npbox,pot,rnorm,nleaf)
      call prin2('lp norm of the potential=*',rnorm,1)
      
      rsc = rsc*rnorm
c     scaled desired precision for error
      reps = eps*rsc
      
      nboxes0 = itree(2*nlevels+2)
      call prinf('before refinement, nboxes0=*',nboxes0,1)


      do ibox=1,nboxes
         ifrefine(ibox)=0
         irefinelev(ibox)=0
         imaxrefinelev(ibox)=0
      enddo
      maxrefinelev=0

c     convert potential to its polynomial expansion coefficients
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,pot,coefsp,umat_nd)

      nrefine=0
      nmaxrefine=0
c     no need to check leaf boxes at the finest level, hence nlevels-1
      do ilev = 0,nlevels-1
         sc = boxsize(ilev)/2
         rscale=sc**eta
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            if (itree(iptr(4)+ibox-1).eq.0) then
c              call fun_err on each leaf box
               call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
     
               erra = erra/rsum
               if(erra.gt.reps) then
                  ifrefine(ibox)=1
                  nrefine=nrefine+1
c                 empirical formula determining the maximum level of refinement
                  ii = ceiling(log(erra/(reps))/log(2.0d0)/norder)
                  imaxrefinelev(ibox)=ii
                  if (ii.gt.maxrefinelev) then
                     maxrefinelev=ii
                     nmaxrefine=1
                  elseif (ii.eq.maxrefinelev) then
                     nmaxrefine=nmaxrefine+1
                  endif
               endif
            endif
         enddo
      enddo
      call prinf('nrefine=*',nrefine,1)
      call prinf('nmaxrefine=*',nmaxrefine,1)
      call prinf('maxrefinelev=*',maxrefinelev,1)

cccc      call prinf('ifrefine=*',ifrefine,nboxes)
      do i=1,nboxes
         if (ifrefine(i).gt.0) then
cccc            print *, i, imaxrefinelev(i)
         endif
      enddo

      do i=0,nlevels
         nblock(i)=0
      enddo

c     recalculate pot, grad, hessian if necessary
      nnewboxes=0
      do ibox = 1,nboxes
         if (ifrefine(ibox).eq.1 .and.
     1       irefinelev(ibox).le.maxrefinelev) then
            irefinelev(ibox)=irefinelev(ibox)+1
            ilev=itree(iptr(2)+ibox-1)
            bsh = boxsize(ilev)/4
            rscale=bsh**eta
c           nchild of ibox
            itree(iptr(4)+ibox-1)=mc
            
            do j=1,mc
               nnewboxes=nnewboxes+1
               jbox = nboxes0+nnewboxes

c              copy list1,ifpwexp info
               ifpwexp(jbox)=ifpwexp(ibox)
               nl1=nlist1(ibox)
               nlist1(jbox)=nl1
               do i=1,nl1
                  list1(i,jbox)=list1(i,ibox)
               enddo
               
               irefinelev(jbox)=irefinelev(ibox)+1
cccc               print *, ibox, jbox
               do k=1,ndim
                  centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
               enddo
c              parent
               itree(iptr(3)+jbox-1)=ibox
c              nchild of jbox
               itree(iptr(4)+jbox-1)=0
c              jbox's children
               do l=1,mc
                  itree(iptr(5)+mc*(jbox-1)+l-1) = -1
               enddo
c              ibox's children
               itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c              jbox's level
               jlev=ilev+1
               itree(iptr(2)+jbox-1) = jlev
               nblock(jlev)=nblock(jlev)+1

c              evaluate plane wave expansion
               if (ifpwexp(jbox).eq.1) then
                  call gnd_find_msshift_ind(ndim,centers(1,jbox),
     1                centers(1,ibox),k)
                  call gnd_shiftpw_loc(nd,nexp,rmlexp(iaddr(2,ibox)),
     1                rmlexp(iaddr(2,jbox)),
     2                wpwmsshift(1,k,nlevels-ilev))
                  call gnd_pw2pgh(ndim,nd,norder,npw,
     1                rmlexp(iaddr(2,jbox)),tab_pw2pot(1,1,jlev),
     2                tab_pw2der(1,1,jlev),tab_pw2dxx(1,1,jlev),
     3                ifpgh,pot(1,1,jbox),grad(1,1,1,jbox),
     4                hess(1,1,1,jbox))
               endif
c              evaluate local interactions
               do k=1,nl1
c              Note that jbox is the target box, kbox is the source box
c              list1 can be used without any change since it's symmetric
c              w.r.t source and target boxes
                  kbox = list1(k,jbox)
                  klev = itree(iptr(2)+kbox-1)
                  bs = boxsize(klev)
                  call gnd_find_loctab_ind(ndim,iperiod,
     1                centers(1,jbox),centers(1,kbox),
     2                bs,bs0,mrefinelev,ixyz)
                  call gnd_tens_prod_to_pghloc(ndim,nd,norder,
     1                fvals(1,1,kbox),ifpgh,pot(1,1,jbox),
     2                grad(1,1,1,jbox),hess(1,1,1,jbox),
     3                nloctab,tab_loc(1,1,-nloctab,jlev),
     4                tabx_loc(1,1,-nloctab,jlev),
     5                tabxx_loc(1,1,-nloctab,jlev),
     6                ind_loc(1,1,-nloctab,jlev),ixyz)
               enddo
               
               call ortho_trans_nd(ndim,nd,0,norder,pot(1,1,jbox),
     1             coefsp(1,1,jbox),umat_nd)
               call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1             iptype,rscale,erra)
               
               erra = erra/rsum
cccc               print *, 'erra=',erra, ibox, jbox
               if(erra.gt.reps) then
                  ifrefine(jbox)=1
c                 the following lines are trying to prevent excessive refinement
                  nrefinelev=1
                  kbox=jbox
                  do k=1,maxrefinelev
                     kdad=itree(iptr(3)+kbox-1)
                     if (kdad.le.nboxes0) exit
                     kbox=kdad
                  enddo
                  
                  if (k.gt.imaxrefinelev(kdad)) then
cccc                     print *, k, kdad, imaxrefinelev(kdad)
                     ifrefine(jbox)=0
                  endif
               endif
                  
            enddo
         endif
      enddo

      allocate(nboxid(nboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))

      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
      if(ifprint .ge. 1) call prinf('nnewboxes=*',nnewboxes,1)

c     call tree_reorg
      if (nnewboxes.eq.0) goto 4600
      call fgt_vol_tree_reorg(ndim,nboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,
     2    centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3    itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4    ifpgh,pot,coefsp,grad,hess)

 4600 continue
      nboxes1=itree(2*nlevels+2)
      if(ifprint .ge. 1) then
         call prinf('after refinement, nlevels=*', nlevels,1)
         call prinf('and nboxes=*', nboxes1,1)
      endif
      
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(9) = time2-time1
      call prinf('laddr=*', itree(iptr(1)),2*(nlevels+1))
      if (nnewboxes.eq.0) goto 5000
      if (nnewboxes.lt.0.1d0*nboxes0) goto 4800      
      
      




      

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 9 (coarsen if possible) ===*',i,0)
c
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      allocate(ifdelete(nboxes))
      call vol_tree_coarsen(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,pot,coefsp,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
      if(ifprint .ge. 1) call prinf('ndelboxes=*', ndelboxes,1)
      if(ifprint .ge. 1) call prinf('nblock=*', nblock,nlevels+1)
      if (ndelboxes.gt.0) then
         call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,pot,coefsp,grad,hess)
      endif
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(10) = time2-time1
      if(ifprint.ge.1)
     1    call prinf('after coarsening, nlevels=*', nlevels,1)
      nboxes0=itree(2*nlevels+2)
      if(ifprint.ge.1)
     1    call prinf('and nboxes0=*', nboxes0,1)
cccc      call prinf('laddr=*', itree(iptr(1)),2*(nlevels+1))
cccc      goto 5000      


      





 4800 continue
      if(ifprint .ge. 1)
     1    call prinf('=== STEP 10 (2:1 rebalance) ===*',i,0)
c
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      nboxes0=itree(2*nlevels+2)
      call computecoll(ndim,nlevels,nboxes0,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))

      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,pot,coefsp,grad,hess,
     2    nboxes,nlevels,centers,boxsize,nboxes0,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(11) = time2-time1
      nboxes0=itree(2*nlevels+2)
      if(ifprint .ge. 1)
     1    call prinf('after 2:1 rebalancing, nboxes=*', nboxes0,1)
cccc      call prinf('laddr=*', itree(iptr(1)),2*(nlevels+1))
      







      
 5000 continue
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 11 (extra targets) ===*',i,0)
c
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      allocate(coefsg(nd,ndim,npbox,nboxes))
      allocate(coefsh(nd,nhess,npbox,nboxes))

c     compute expansion coefficients for potential, gradient, and hessian

      iffast=0

      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,pot,coefsp,umat_nd)

      if (iffast.eq.1) goto 5200
      
      if (ifpghtarg.ge.1) then
         call treedata_evalt_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,coefsp,
     2       ntarg,targs,pote)
      endif
      if (ifpghtarg.ge.2) then
         call treedata_trans_nd(ndim,nd*ndim,itype,nlevels,itree,
     1       iptr,boxsize,norder,grad,coefsg,umat_nd)
         call treedata_evalt_nd(ndim,nd*ndim,ipoly,norder,nboxes,
     1       nlevels,ltree,itree,iptr,centers,boxsize,coefsg,
     2       ntarg,targs,grade)
      endif
      if (ifpghtarg.ge.3) then
         call treedata_trans_nd(ndim,nd*nhess,itype,nlevels,itree,
     1       iptr,boxsize,norder,hess,coefsh,umat_nd)
         call treedata_evalt_nd(ndim,nd*nhess,ipoly,norder,nboxes,
     1       nlevels,ltree,itree,iptr,centers,boxsize,coefsh,
     2       ntarg,targs,hesse)
      endif
      goto 5500
      
 5200 continue
c     evaluate the potential, gradient and hessian at targets
c     using coefficients for potential only
      if (ifpghtarg.eq.1) then
         call treedata_evalt_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,coefsp,
     2       ntarg,targs,pote)
      elseif (ifpghtarg.eq.2) then
         call treedata_evalpgt_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,coefsp,
     2       ntarg,targs,pote,grade)
      elseif (ifpghtarg.eq.3) then
         call treedata_evalpght_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,coefsp,
     2       ntarg,targs,pote,grade,hesse)
      endif

 5500 continue
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(12) = time2-time1
      





      npwlevel=npwlevel0
      iperiod=iperiod0
      
      if(ifprint.eq.1) then 
         call prin2('timeinfo=*',timeinfo,11)
         d= 0
         do i = 1,11
            d = d + timeinfo(i)
         enddo
         call prin2('time on tensor grid=*',d,1)
         call prin2('tensor grid speed in pps=*',
     1       (nleafbox*npbox*nd+0.0d0)/d,1)
         d=timeinfo(12)
         call prin2('time on extra targets=*',d,1)
         call prin2('extra targets speed in pps=*',
     1       (ntarg*nd+0.0d0)/d,1)
      ENDIF
      
      return
      end
c
c
c
c
c
c------------------------------------------------------------------    
      subroutine gnd_mpalloc(nd,ndim,laddr,iaddr,
     1    nboxes,nlevels,npwlevel,lmptot,npw)
c     This subroutine determines the size of the array
c     to be allocated for multipole/local expansions
c
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     ndim        in: integer
c                 dimension of the underlying space
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
      integer nd,ndim,nlevels,npwlevel,npw,nboxes
      integer iaddr(2,*), laddr(2,0:nlevels)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart,jbox
c
      istart = 1
      if (npwlevel .eq. nlevels) then
         lmptot=0
         return
      endif
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      nn = (npw+1)/2
      do i=1,ndim-1
         nn = nn*npw
      enddo
      nn = nn*2*nd

      do i = nlevstart,nlevels

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole PW expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i = nlevstart,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local PW expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(2,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo

c     for possible refinement
      do jbox=laddr(2,nlevels)+1,nboxes
         iaddr(2,jbox) = istart + nn
         istart = istart+nn
      enddo
      
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
c
c
c
