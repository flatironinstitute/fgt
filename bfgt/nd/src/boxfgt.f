c     
c    $Date$
c    $Revision$

      subroutine boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2    ifpgh,pot,grad,hess,ntarg,targs,ifpghtarg,pote,grade,hesse,
     3    timeinfo)
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
      real *8 grad(nd,ndim,npbox,nboxes)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,nboxes)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,ntarg)
      real *8 hesse(nd,ndim*(ndim+1)/2,ntarg)

      real *8 centers(ndim,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(*)
c
cc     additional fgt variables
c
      integer *8 lmptot
      integer, allocatable :: iaddr(:,:)
      real *8, allocatable :: rmlexp(:)

c
cc      temporary variables
c
      integer npwlevel
      integer i,ilev,lmptmp,idim,ndim
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      real *8 dcutoff
      real *8 omp_get_wtime
      real *8 time1,time2,pi,done,pmax,bs0,bsize,pweps

      ifprint = 1

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
     1       nlevels,0,lmptot,npw)
      else
         call gnd_mpalloc(nd,ndim,itree,iaddr,
     1       nlevels,npwlevel,lmptot,npw)
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
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      if (npwlevel.le.1.and.iperiod.eq.1) then
         call periodic_bfgtmain_large_delta(
     1       nd,ndim,delta,eps,ipoly,norder,npbox,
     2       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     3       iaddr,rmlexp,npw,
     4       ifpgh,pot,grad,hess,
     5       ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
      else
         call bfgtmain(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2       iaddr,rmlexp,npwlevel,pmax,npw,
     3       ifpgh,pot,grad,hess,
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
     3    ifpgh,pot,grad,hess,
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
      real *8 grad(nd,ndim,npbox,nboxes)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,nboxes)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,ntarg)
      real *8 hesse(nd,ndim*(ndim+1)/2,ntarg)

      real *8 boxsize(0:nlevels),centers(ndim,nboxes)
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(*)

c
cc        temporary variables
c
      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlistpw(:), listpw(:,:)
      integer, allocatable :: ifhung(:),iflocal(:)

      integer ndirect
      integer ixyz(ndim)
      
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: umat_nd(:,:,:)

      real *8 ws(100),ts(100)
      
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
      call prinf('nleafbox=*',nleafbox,1)
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      ndirect=0
      do ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
c           Check if the current box is a nonempty leaf box            
            if(nchild.eq.0) then
               ifhung(ibox) = 1
               ndirect = ndirect+1
            endif
         enddo
      enddo

      nasym=3
      sigma=(delta/boxsize(nlevels)**2)**nasym/8
      if(ifprint.eq.1) call 
     1    prin2('estimated asymptotic expansion error=*',sigma,1)
      
cccc      if (nleafbox.eq.ndirect .and. sigma.le.2*eps) goto 3000
cccc      goto 3000
      
      if(ifprint.eq.1) call 
     1    prinf('number of direct evaluation source boxes=*',ndirect,1)

      do i=1,10
         timeinfo(i)=0
      enddo
c
c
c       compute list info
c
      isep = 1
      call compute_mnlist1(ndim,nboxes,nlevels,itree(iptr(1)),centers,
     1    boxsize,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     2    isep,itree(iptr(6)),itree(iptr(7)),iperiod,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))

c     modified list1 for direct evaluation
      call compute_modified_list1(ndim,npwlevel,nboxes,nlevels,
     1    ltree,itree,iptr,centers,boxsize,iperiod,mnlist1,nlist1,list1)

c     compute the tables converting Legendre polynomial expansion to potential
c     values, used in direct evaluation
      allocate(tab_loc(norder,norder,-6:6,0:nlevels))
      allocate(tabx_loc(norder,norder,-6:6,0:nlevels))
      allocate(tabxx_loc(norder,norder,-6:6,0:nlevels))
      allocate(ind_loc(2,norder+1,-6:6,0:nlevels))

      nnodes=50
      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_all(eps,ipoly,norder,nnodes,delta,
     1       boxsize(ilev),tab_loc(1,1,-6,ilev),
     2       tabx_loc(1,1,-6,ilev),tabxx_loc(1,1,-6,ilev),
     3       ind_loc(1,1,-6,ilev))
      enddo

c
c     check whether we need to create and evaluate planewave expansions 
c     for boxes at the cutoff level
      allocate(iflocal(nboxes))
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
         do i=1,nboxes
            iflocal(i)=0
         enddo

         do ilev=npwlevel,npwlevel
            do ibox=itree(2*ilev+1),itree(2*ilev+2)
               nchild = itree(iptr(4)+ibox-1)
               if (nchild .eq. 0) then
                  ncoll = itree(iptr(6)+ibox-1)
                  do j=1,ncoll
                     jbox = itree(iptr(7) + (ibox-1)*mnbors+j-1)
                     nchild = itree(iptr(4)+jbox-1)
                     jlev = itree(iptr(2)+jbox-1)
                     if (nchild .gt. 0 .and. jlev.eq.ilev) then
                        iflocal(ibox)=1
                        exit
                     endif
                  enddo
               endif
            enddo
         enddo
      endif
      
c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800

c      
c     get planewave nodes and weights
      call get_pwnodes(pmax,npw,ws,ts)

c     compute translation matrices for PW expansions
      nlevstart = max(npwlevel,0)

c     compute the tables converting Legendre polynomial expansion
c     to planewave expansion
      nnodes = 100
      allocate(tab_poly2pw(norder,npw,0:nlevels))

      do ilev=nlevstart,nlevels
         call mk_poly2pw_tables(norder,ipoly,npw,nnodes,ws,ts,delta,
     1       boxsize(ilev),tab_poly2pw(1,1,ilev))
      enddo

c     compute the tables converting planewave expansions to potential values
      allocate(tab_pw2pot(npw,norder,0:nlevels))
      allocate(tab_pw2der(npw,norder,0:nlevels))
      allocate(tab_pw2dxx(npw,norder,0:nlevels))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pgh_tables(norder,npw,ts,xq,delta,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev),
     1       tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev))  
      enddo

      nexp=(npw+1)/2
      do i=1,ndim-1
         nexp = nexp*npw
      enddo
      
      nmax = nlevels-max(npwlevel,0)      
      allocate(wpwmsshift(nexp,mc,nmax))
      xmin2 = boxsize(nlevels)/sqrt(delta)/2
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
      call gtnd_compute_mnlistpw(ndim,nboxes,nlevels,ltree,itree,
     1    iptr,centers,boxsize,iper,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call gtnd_compute_listpw(ndim,npwlevel,nboxes,nlevels,
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
c
c        step 1: convert function values to planewave expansions
c
    
      if(ifprint.eq.1) 
     1   call prinf("=== STEP 1 (values -> mp) ====*",i,0)
      
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do 1100 ilev = nlevels,nlevstart,-1
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=itree(2*ilev+1),itree(2*ilev+2)
          if (ilev .ne. npwlevel .or. iflocal(ibox).eq.1) then
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
     1      call prinf('=== STEP 2 (merge mp) ====*',i,0)
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
c
      
      if(ifprint.ge.1)
     1    call prinf('=== Step 3 (mp to loc) ===*',i,0)
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
c          shift PW expansions
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
     1    call prinf('=== Step 4 (split loc) ===*',i,0)

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
     1           centers(1,ibox),k)
             call gnd_shiftpw_loc(nd,nexp,rmlexp(iaddr(2,ibox)),
     1           rmlexp(iaddr(2,jbox)),wpwmsshift(1,k,nlevels-ilev))
          enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1400 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1
      
      if(ifprint.ge.1)
     1    call prinf('=== step 5 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      do 1500 ilev = nlevstart,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           if (ilev .eq. npwlevel .and. iflocal(ibox).eq.0) then
c            do nothing here
           else    
             nchild = itree(iptr(4)+ibox-1)
             if(nchild.eq.0) then
               call gnd_pw2pgh(ndim,nd,norder,npw,
     1            rmlexp(iaddr(2,ibox)),tab_pw2pot(1,1,ilev),
     2            tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev),ifpgh,
     3            pot(1,1,ibox),grad(1,1,1,ibox),hess(1,1,1,ibox))
             endif
           endif
        enddo
cccc 222    format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f5.2)         
cccc        write(6,222) ilev,nb,dt
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(6) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     1     call prinf('=== STEP 6 (direct) =====*',i,0)
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
            if (ifhung(ibox) .eq. 1) then
               
               nl1 = nlist1(ibox)
               do j=1,nl1
cccc              jbox is the target box
                  jbox = list1(j,ibox)
                  jlev = itree(iptr(2)+jbox-1)
                  if ((ilev .ne. jlev) .or.
     1                ((iflocal(ibox).ne.1) .or. (jbox.ne.ibox))) then
                     call gnd_find_loctab_ind(ndim,iperiod,
     1                   centers(1,jbox),centers(1,ibox),
     2                   bs,bs0,ixyz)
                     call gnd_tens_prod_to_pghloc(ndim,nd,norder,
     1                   fvals(1,1,ibox),ifpgh,pot(1,1,jbox),
     2                   grad(1,1,1,jbox),hess(1,1,1,jbox),
     3                   tab_loc(1,1,-6,jlev),tabx_loc(1,1,-6,jlev),
     4                   tabxx_loc(1,1,-6,jlev),ind_loc(1,1,-6,jlev),
     3                   ixyz)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue
c
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(7) = time2-time1
      goto 4000

 3000 continue
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (asymptotic regime) =====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()

c     evaluate potential
      call treedata_eval_pot_nd_asym(ndim,nd,delta,ipoly,nasym,
     1    nlevels,itree,iptr,boxsize,norder,fvals,pot)

      itype=0
      iffast=1
      if (iffast.eq.1) then
c     use spectral differentiation to evaluate gradient and hessian
         if (ifpgh.eq.2) then
            allocate(fcoefs(nd,npbox,nboxes))
            ng=nd*ndim
            call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1          iptr,boxsize,norder,pot,fcoefs,umat_nd)
            call treedata_evalg_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1          boxsize,norder,fcoefs,grad)
         endif
      
         if (ifpgh.eq.3) then
            allocate(fcoefs(nd,npbox,nboxes))
            ng=nd*ndim
            nh=nd*nhess
            call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1          iptr,boxsize,norder,pot,fcoefs,umat_nd)
            call treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1          boxsize,norder,fcoefs,grad,hess)
         endif
      elseif (iffast.eq.0) then
c     also use asymptotic expansion to evaluate gradient and hessian
         if (ifpgh.eq.2) then
            allocate(fcoefs(nd,npbox,nboxes))
            allocate(gvals(nd,ndim,npbox,nboxes))
            ng=nd*ndim
            call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1          iptr,boxsize,norder,fvals,fcoefs,umat_nd)
            call treedata_evalg_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1          boxsize,norder,fcoefs,gvals)

            call treedata_eval_pot_nd_asym(ndim,ng,delta,ipoly,nasym,
     1          nlevels,itree,iptr,boxsize,norder,gvals,grad)
         endif
      
         if (ifpgh.eq.3) then
            allocate(fcoefs(nd,npbox,nboxes))
            allocate(gvals(nd,ndim,npbox,nboxes))
            allocate(hvals(nd,nhess,npbox,nboxes))
            ng=nd*ndim
            nh=nd*nhess
            call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1          iptr,boxsize,norder,fvals,fcoefs,umat_nd)
            
            call treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1          boxsize,norder,fcoefs,gvals,hvals)
            call treedata_eval_pot_nd_asym(ndim,ng,delta,ipoly,nasym,
     1          nlevels,itree,iptr,boxsize,norder,gvals,grad)
            call treedata_eval_pot_nd_asym(ndim,nh,delta,ipoly,nasym,
     1          nlevels,itree,iptr,boxsize,norder,hvals,hess)
      endif
      endif
      
      
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(8) = time2-time1

 4000 continue
      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (extra targets) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      allocate(coefsp(nd,npbox,nboxes))
      allocate(coefsg(nd,ndim,npbox,nboxes))
      allocate(coefsh(nd,nhess,npbox,nboxes))

c     compute expansion coefficients for potential, gradient, and hessian
      itype=0
      if (ifpghtarg.ge.1) then
         call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,pot,coefsp,umat_nd)
      endif

      iffast=1
      if (iffast.eq.0) then
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
      endif

      if (iffast.eq.1) then
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
      endif
      
      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(9) = time2-time1
      
      
      if(ifprint.eq.1) then 
         call prin2('timeinfo=*',timeinfo,9)
         d= 0
         do i = 1,8
            d = d + timeinfo(i)
         enddo
         call prin2('time on tensor grid=*',d,1)
         call prin2('tensor grid speed in pps=*',
     1       (nleafbox*npbox*nd+0.0d0)/d,1)
         d=timeinfo(9)
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
      subroutine periodic_bfgtmain_large_delta(
     1    nd,ndim,delta,eps,ipoly,norder,npbox,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     3    iaddr,rmlexp,npw,
     4    ifpgh,pot,grad,hess,
     5    ntarg,targs,ifpghtarg,pote,grade,hesse,timeinfo)
c
c     This code computes doubly periodic Gauss transform
c     to a collection of right hand sides when delta is very large and
c     doubly periodic conditions are imposed.
c      
c     method: in this case, the periodic Green's function admits an efficient
c     seprable Fourier series expansion, which is used in the entire computational
c     box. 1. compute planewave expansion coefficients for each leaf box; 2. merge
c     pw exp. all the way to the root box - this is a valid expansion in the entire
c     box; 3. split the pw exp. into its child boxes; 4. evaluate the pw exp. in each
c     leaf box.
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
      real *8 grad(nd,ndim,npbox,nboxes)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,nboxes)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,ntarg)
      real *8 hesse(nd,ndim*(ndim+1)/2,ntarg)

      real *8 boxsize(0:nlevels),centers(ndim,nboxes)
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(*)
c
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)

      real *8 ws(200),ts(200)
      
      complex *16, allocatable :: wpwmsshift(:,:,:)
      
      complex *16, allocatable :: tab_poly2pw(:,:,:)
      complex *16, allocatable :: tab_pw2pot(:,:,:)
      complex *16, allocatable :: tab_pw2der(:,:,:)
      complex *16, allocatable :: tab_pw2dxx(:,:,:)
      
      complex *16, allocatable :: ff(:,:)

      ifprint = 0

      mc = 2**ndim
      mnbors = 3**ndim
      
      done = 1
      pi = atan(done)*4
c
      allocate(xq(norder),umat(norder,norder),
     1   vmat(norder,norder),wts(norder))
     
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)

c
c     get planewave nodes and weights
      call get_periodic_pwnodes(boxsize(0),delta,eps,npw,ws,ts)

c     compute translation matrices for PW expansions
      nlevstart = max(npwlevel,0)

c     compute the tables converting function values at tensor product grid
c     to planewave expansion
      nnodes = 100
      allocate(tab_poly2pw(norder,npw,0:nlevels))
      allocate(ff(npw,norder))

      do ilev=nlevstart,nlevels
         call mk_poly2pw_tables(norder,ipoly,npw,nnodes,ws,ts,1.0d0,
     1       boxsize(ilev),tab_poly2pw(1,1,ilev))
      enddo
c     compute the tables converting planewave expansions to potential values
      allocate(tab_pw2pot(npw,norder,0:nlevels))
      allocate(tab_pw2der(npw,norder,0:nlevels))
      allocate(tab_pw2dxx(npw,norder,0:nlevels))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pgh_tables(norder,npw,ts,xq,delta,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev),
     1       tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev))  
      enddo
      
      nexp=(npw+1)/2
      do i=1,ndim-1
         nexp = nexp*npw
      enddo
      
      nmax = nlevels
      allocate(wpwmsshift(nexp,mc,nmax))
      xmin2 = boxsize(nlevels)/2
      call gnd_mk_merge_split_pw_matrices(ndim,xmin2,npw,ts,nmax,
     1    wpwmsshift)

      if(ifprint.eq.1)
     1    call prinf('laddr=*',itree(iptr(1)),2*(nlevels+1))
      
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
c
c        step 1: convert function values to planewave expansions
c
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (values -> mp) ====*",i,0)
      
      do 1100 ilev = nlevels,0,-1
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
c           Check if current box is a leaf box            
            if(nchild.eq.0) then
c              form PW expansion directly
               call gnd_tens_prod_to_pw(ndim,nd,norder,fvals(1,1,ibox),
     1             npw,tab_poly2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
            endif
         enddo
C$OMP END PARALLEL DO
cccc 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f5.2)         
cccc         write(6,111) ilev,nb,dt
 1100 continue
      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()
      timeinfo(1) = time2-time1
      
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do 1200 ilev=nlevels-1,0,-1
        klev = nlevels - ilev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,dx,dy,dz,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild .gt. 0) then
              call gnd_pwzero(nd,rmlexp(iaddr(1,ibox)),nexp)
           endif
           
          do i=1,nchild
            jbox = itree(iptr(5)+mc*(ibox-1)+i-1)
            call gnd_find_msshift_ind(ndim,centers(1,ibox),
     1          centers(1,jbox),k)
            call gnd_shiftpw(nd,nexp,rmlexp(iaddr(1,jbox)),
     1          rmlexp(iaddr(1,ibox)),wpwmsshift(1,k,klev))
          enddo
        enddo
C$OMP END PARALLEL DO    
 1200 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1
c
      
      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      ibox=1
      call gnd_copy_pwexp(nd,nexp,rmlexp(iaddr(1,ibox)),
     1    rmlexp(iaddr(2,ibox)))
c
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3) = time2-time1

cccc      call prin2('timeinfo4=*',time2-time1,1)

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      do 1400 ilev = 0,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4)+ibox-1)
          do i=1,nchild
             jbox = itree(iptr(5)+mc*(ibox-1)+i-1)
             call gnd_find_msshift_ind(ndim,centers(1,jbox),
     1           centers(1,ibox),k)
             call gnd_shiftpw_loc(nd,nexp,rmlexp(iaddr(2,ibox)),
     1           rmlexp(iaddr(2,jbox)),wpwmsshift(1,k,nlevels-ilev))
          enddo
        enddo
C$OMP END PARALLEL DO        
 1400 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1
      
      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      do 1500 ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
               call gnd_pw2pgh(ndim,nd,norder,npw,
     1            rmlexp(iaddr(2,ibox)),tab_pw2pot(1,1,ilev),
     2            tab_pw2der(1,1,ilev),tab_pw2dxx(1,1,ilev),ifpgh,
     3            pot(1,1,ibox),grad(1,1,1,ibox),hess(1,1,1,ibox))
            endif
         enddo
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(5) = time2 - time1

      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,5)
      d = 0
      do i = 1,5
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
      subroutine gnd_mpalloc(nd,ndim,laddr,iaddr,
     1    nlevels,npwlevel,lmptot,npw)
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
      integer nd,ndim,nlevels,npwlevel,npw
      integer iaddr(2,*), laddr(2,0:nlevels)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart
c
      istart = 1
      if (npwlevel .eq. nlevels) then
         lmptot=0
         return
      endif
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      nn = 1
      do i=1,ndim-1
         nn = nn*npw
      enddo
      nn = nn*((npw+1)/2)*2*nd

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
            
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
c
c
c
