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
     3       iaddr,rmlexp,npw,pot,timeinfo)
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
      
      real *8, allocatable :: fcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)

      real *8 timelev(0:200)
      real *8 ws(100),ts(100)
      
      complex *16, allocatable :: wpwshift(:,:,:)
      complex *16, allocatable :: wpwmsshift(:,:,:)
      
      complex *16, allocatable :: tab_poly2pw(:,:,:)
      complex *16, allocatable :: tab_pw2pot(:,:,:)
      complex *16, allocatable :: tab_pw2potx(:,:,:)
      complex *16, allocatable :: tab_pw2potxx(:,:,:)
      complex *16, allocatable :: tab_poly2pw2(:,:,:)

      complex *16, allocatable :: ff(:,:)
      complex *16, allocatable :: gg(:,:)
      complex *16, allocatable :: gg2(:,:)
      complex *16, allocatable :: gg3(:,:)

      real *8, allocatable :: hh(:,:)
      real *8, allocatable :: hh2(:,:)
      real *8, allocatable :: hh3(:,:)
      real *8, allocatable :: tabf(:,:)
      real *8, allocatable :: tabfx(:,:)
      real *8, allocatable :: tabfxx(:,:)
      real *8, allocatable :: lgcoefs(:,:,:)
      real *8, allocatable :: px(:)
      real *8, allocatable :: py(:)
      real *8, allocatable :: pxd(:)
      real *8, allocatable :: pyd(:)

      real *8, allocatable :: tab_coll(:,:,:,:)
      real *8, allocatable :: tab_stob(:,:,:,:)
      real *8, allocatable :: tab_btos(:,:,:,:)

      integer, allocatable :: ind_coll(:,:,:,:)
      integer, allocatable :: ind_stob(:,:,:,:)
      integer, allocatable :: ind_btos(:,:,:,:)
      
      ifprint = 1

      done = 1
      pi = atan(done)*4

      bs0 = boxsize(0)
c
c       compute coefs
c
      allocate(fcoefs(npbox,nd,nboxes),xq(norder),umat(norder,norder),
     1   vmat(norder,norder),wts(norder))
     
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 0 (precomputation) =====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do ilev = 0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             do ind=1,nd
                do i=1,npbox
                   fcoefs(i,ind,ibox)=fvals(ind,i,ibox)
c                   if (abs(fvals(ind,i,ibox)).lt.1d-16) 
c     1                 fcoefs(i,ind,ibox)=0
                enddo
             enddo
c     1        call tens_prod_trans_nd(nd,norder,fvals(1,1,ibox),
c     2        fcoefs(1,1,ibox),umat)
          endif
        enddo
      enddo
c      
      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      ndirect=0
      do ilev = 0,nlevend
C
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
      if(ifprint.eq.1) call 
     1    prinf('number of direct evaluation source boxes=*',ndirect,1)

      do i=1,10
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo

c
c
c       compute list info
c
      isep = 1
      call compute_mnlist1(ndim,nlevels,nboxes,itree(iptr(1)),boxsize,
     1    centers,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     2    isep,itree(iptr(6)),itree(iptr(7)),iperiod,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))

c     modified list1 for direct evaluation
      call compute_modified_list1(ndim,nlevels,npwlevel,
     1  nboxes,itree,ltree,iptr,
     2  centers,boxsize,
     3  iperiod,mnlist1,nlist1,list1)

c     compute the tables converting Legendre polynomial expansion to potential
c     values, used in direct evaluation
      allocate(tab_coll(norder,norder,-1:1,0:nlevels))
      allocate(tab_stob(norder,norder,4,0:nlevels))
      allocate(tab_btos(norder,norder,4,0:nlevels))

      allocate(ind_coll(2,norder+1,-1:1,0:nlevels))
      allocate(ind_stob(2,norder+1,4,0:nlevels))
      allocate(ind_btos(2,norder+1,4,0:nlevels))
      
      allocate(hh(norder,norder))
      allocate(hh2(norder,norder))
      allocate(hh3(norder,norder))

      allocate(tabf(norder,norder))
      allocate(tabfx(norder,norder))
      allocate(tabfxx(norder,norder))

      allocate(lgcoefs(norder,norder,nd))
      allocate(px(norder+1))
      allocate(py(norder+1))
      allocate(pxd(norder+1))
      allocate(pyd(norder+1))



      nnodes=50
      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_coll(eps,ipoly,norder,nnodes,delta,
     1       boxsize(ilev),tab_coll(1,1,-1,ilev),ind_coll(1,1,-1,ilev))
      enddo

      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_stob(eps,ipoly,norder,nnodes,delta,
     1       boxsize(ilev),tab_stob(1,1,1,ilev),ind_stob(1,1,1,ilev))
      enddo

      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_btos(eps,ipoly,norder,nnodes,delta,
     1       boxsize(ilev),tab_btos(1,1,1,ilev),ind_btos(1,1,1,ilev))
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
                     jbox = itree(iptr(7) + (ibox-1)*9+j-1)
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
      allocate(tab_poly2pw2(norder,npw,0:nlevels))
      allocate(ff(norder,npw/2))

      do ilev=nlevstart,nlevels
         call mk_poly2pw_tables(norder,ipoly,npw,nnodes,ws,ts,delta,
     1       boxsize(ilev),tab_poly2pw(1,1,ilev))
      enddo

c     compute the tables converting planewave expansions to potential values
      allocate(tab_pw2pot(npw,norder,0:nlevels))
      allocate(tab_pw2potx(npw,norder,0:nlevels))
      allocate(tab_pw2potxx(npw,norder,0:nlevels))
      allocate(gg(npw/2,norder))
      allocate(gg2(npw/2,norder))
      allocate(gg3(npw/2,norder))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pgh_tables(norder,npw,ts,xq,delta,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev),
     1       tab_pw2potx(1,1,ilev),tab_pw2potxx(1,1,ilev))  
      enddo

      nexp=(npw+1)/2
      do i=1,ndim-1
         nexp = nexp*npw
      enddo
      
      nmax = 1
      allocate(wpwshift(nexp,-nmax:nmax,-nmax:nmax))
      xmin  = boxsize(nlevstart)/sqrt(delta)
      call mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1    wpwshift)

      nmax = nlevels-max(npwlevel,0)      
      allocate(wpwmsshift(nexp,4,nmax))
      xmin2 = boxsize(nlevels)/sqrt(delta)/2
      call mk_merge_split_pw_matrices(xmin2,npw,ts,nmax,wpwmsshift)
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
c        step 1: convert coeffs to planewave expansions
c
    
      if(ifprint.eq.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0)
      
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do 1100 ilev = nlevels,nlevstart,-1
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            if (ilev .ne. npwlevel .or. iflocal(ibox).eq.1) then
               nchild = itree(iptr(4)+ibox-1)
c              Check if current box is a leaf box            
               if(nchild.eq.0) then
c                 form PW expansion directly
                  call gnd_tens_prod_to_pw(nd,norder,fcoefs(1,1,ibox),
     1                npw,ff,
     2                tab_poly2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
               endif
            endif
         enddo
C$OMP END PARALLEL DO
c     end of ilev do loop
 1100 continue
      
      call cpu_time(time2)
C$       time2 = omp_get_wtime()
      timeinfo(2) = time2-time1

      
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do 1200 ilev=nlevels-1,nlevstart,-1
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
              jbox = itree(iptr(5)+4*(ibox-1)+i-1)

              dx= (centers(1,ibox) - centers(1,jbox))
              dy= (centers(2,ibox) - centers(2,jbox))
              if (dx.gt.0 .and. dy.gt.0) then
                 k=1
              elseif (dx.gt.0 .and. dy.lt.0) then
                 k=2
              elseif (dx.lt.0 .and. dy.gt.0) then
                 k=3
              elseif (dx.lt.0 .and. dy.lt.0) then
                 k=4
              endif

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
     $    call prinf('=== Step 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = nlevstart,nlevstart
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,j)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           call gnd_copypwexp(nd,nexp,rmlexp(iaddr(1,ibox)),
     1         rmlexp(iaddr(2,ibox)))
c          shift PW expansions
           do j=1,nlistpw(ibox)
              jbox=listpw(j,ibox)
              jx= nint((centers(1,ibox) - centers(1,jbox))/xmin)
              jy= nint((centers(2,ibox) - centers(2,jbox))/xmin)

              if (iperiod .eq. 1) then
                 jxp1=nint((centers(1,ibox) - centers(1,jbox)-bs0)/xmin)
                 jxm1=nint((centers(1,ibox) - centers(1,jbox)+bs0)/xmin)
                 jyp1=nint((centers(2,ibox) - centers(2,jbox)-bs0)/xmin)
                 jym1=nint((centers(2,ibox) - centers(2,jbox)+bs0)/xmin)
cccc                 print *, jx,jxp1,jxm1,jy,jyp1,jym1
                 if (abs(jx).gt.abs(jxp1)) jx=jxp1
                 if (abs(jx).gt.abs(jxm1)) jx=jxm1
                 if (abs(jy).gt.abs(jyp1)) jy=jyp1
                 if (abs(jy).gt.abs(jym1)) jy=jym1
               endif
              

               call gnd_shiftpw(nd,nexp,rmlexp(iaddr(1,jbox)),
     1             rmlexp(iaddr(2,ibox)),wpwshift(1,jx,jy))
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
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      do 1400 ilev = nlevstart,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,dx,dy,dz,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild.gt.0) then
          do i=1,nchild
             jbox = itree(iptr(5)+4*(ibox-1)+i-1)
             dx= centers(1,jbox) - centers(1,ibox)
             dy= centers(2,jbox) - centers(2,ibox)
             if (dx.gt.0 .and. dy.gt.0) then
                k=1
             elseif (dx.gt.0 .and. dy.lt.0) then
                k=2
             elseif (dx.lt.0 .and. dy.gt.0) then
                k=3
             elseif (dx.lt.0 .and. dy.lt.0) then
                k=4
             endif
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
     $    call prinf('=== step 5 (eval loc) ===*',i,0)

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
                if (ifpgh.eq.1) then
                   call gnd_pw2pot(nd,norder,npw,
     1              rmlexp(iaddr(2,ibox)),gg,
     2              tab_pw2pot(1,1,ilev),
     3              pot(1,1,ibox))
                endif
                if (ifpgh.eq.2) then
                   call gnd_pw2pg(nd,norder,npw,
     1              rmlexp(iaddr(2,ibox)),gg,gg2,
     2              tab_pw2pot(1,1,ilev),tab_pw2potx(1,1,ilev),
     3              pot(1,1,ibox),grad(1,1,1,ibox))
                endif
                if (ifpgh.eq.3) then
                   call gnd_pw2pgh(nd,norder,npw,
     1              rmlexp(iaddr(2,ibox)),gg,gg2,gg3,
     2              tab_pw2pot(1,1,ilev),tab_pw2potx(1,1,ilev),
     2              tab_pw2potxx(1,1,ilev),
     3              pot(1,1,ibox),grad(1,1,1,ibox),hess(1,1,1,ibox))
                endif
             endif
           endif
ccc    end of ibox loop        
        enddo
 222    format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f5.2)         
cccc        write(6,222) ilev,nb,dt
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(6) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nl1,bs,ix,iy,iz,jlev)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
c        ibox is the source box            
            if (ifhung(ibox) .eq. 1) then
               
               nl1 = nlist1(ibox)
               do j=1,nl1
cccc              jbox is the target box
                  jbox = list1(j,ibox)
                  jlev = itree(iptr(2)+jbox-1)
                  bs = boxsize(jlev)
c                 colleague                  
                  if (ilev .eq. jlev) then
                     if (iflocal(ibox).eq. 1 .and. jbox.eq.ibox) then
                     else
                     ix = (centers(1,jbox)-centers(1,ibox))/bs
                     iy = (centers(2,jbox)-centers(2,ibox))/bs
                     if (iperiod .eq. 1) then
                        ixp1=ix-2**jlev
                        ixm1=ix+2**jlev
                        iyp1=iy-2**jlev
                        iym1=iy+2**jlev
                        if (abs(ix).gt.abs(ixp1)) ix=ixp1
                        if (abs(ix).gt.abs(ixm1)) ix=ixm1
                        if (abs(iy).gt.abs(iyp1)) iy=iyp1
                        if (abs(iy).gt.abs(iym1)) iy=iym1
                     endif

cccc                     print *, 'coll jlev=',jlev,'ix=',ix,'iy=',iy
                     call gnd_tens_prod_to_potloc2(nd,norder,
     1                   fcoefs(1,1,ibox),hh,pot(1,1,jbox),
     2                   tab_coll(1,1,ix,jlev),
     3                   tab_coll(1,1,iy,jlev),
     4                   ind_coll(1,1,ix,jlev),
     5                   ind_coll(1,1,iy,jlev))
                     
                     endif
c                 big source box to small target box                     
                  elseif (ilev .eq. jlev-1) then
                     dx = (centers(1,jbox)-centers(1,ibox))/bs
                     dy = (centers(2,jbox)-centers(2,ibox))/bs
cccc                     print *, 'ilev=',ilev,'dx=',dx,'dy=',dy
                     if (iperiod .eq. 1) then
                        dxp1=dx-bs0/bs
                        dxm1=dx+bs0/bs
                        dyp1=dy-bs0/bs
                        dym1=dy+bs0/bs
                        if (abs(dx).gt.abs(dxp1)) dx=dxp1
                        if (abs(dx).gt.abs(dxm1)) dx=dxm1
                        if (abs(dy).gt.abs(dyp1)) dy=dyp1
                        if (abs(dy).gt.abs(dym1)) dy=dym1
                     endif
                     ix=dx+2.55d0
                     iy=dy+2.55d0
cccc                     print *, 'btos jlev=',jlev,'ix=',ix,'iy=',iy

                     call gnd_tens_prod_to_potloc2(nd,norder,
     1                   fcoefs(1,1,ibox),hh,pot(1,1,jbox),
     2                   tab_btos(1,1,ix,jlev),
     3                   tab_btos(1,1,iy,jlev),
     4                   ind_btos(1,1,ix,jlev),
     5                   ind_btos(1,1,iy,jlev))
c                 small source box to big target box 
                  elseif (ilev .eq. jlev+1) then
                     dx = (centers(1,jbox)-centers(1,ibox))/bs
                     dy = (centers(2,jbox)-centers(2,ibox))/bs

                     if (iperiod .eq. 1) then
                        dxp1=dx-bs0/bs
                        dxm1=dx+bs0/bs
                        dyp1=dy-bs0/bs
                        dym1=dy+bs0/bs
                        if (abs(dx).gt.abs(dxp1)) dx=dxp1
                        if (abs(dx).gt.abs(dxm1)) dx=dxm1
                        if (abs(dy).gt.abs(dyp1)) dy=dyp1
                        if (abs(dy).gt.abs(dym1)) dy=dym1
                     endif
                     ix=dx*2+2.55d0
                     iy=dy*2+2.55d0
cccc                     print *, 'stob jlev=',jlev,'ix=',ix,'iy=',iy
                     call gnd_tens_prod_to_potloc2(nd,norder,
     1                   fcoefs(1,1,ibox),hh,pot(1,1,jbox),
     2                   tab_stob(1,1,ix,jlev),
     3                   tab_stob(1,1,iy,jlev),
     4                   ind_stob(1,1,ix,jlev),
     5                   ind_stob(1,1,iy,jlev))
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue
c
      call mk_poly_tables(norder,xq,tabf,tabfx,tabfxx)

      
      do 3000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
c        ibox is the source box            
            if (ifhung(ibox) .eq. 1) then
               call tens_prod_trans_nd(nd,norder,pot(1,1,ibox),
     1             lgcoefs,umat)
               if (ifpgh.eq.2) then
                  call ortho_to_g(nd,norder,lgcoefs,hh,hh2,
     1                grad(1,1,1,ibox),
     1                tabf,tabfx,boxsize(ilev))
                endif
               if (ifpgh.eq.3) then
                  call ortho_to_gh(nd,norder,lgcoefs,hh,hh2,hh3,
     1                grad(1,1,1,ibox),hess(1,1,1,ibox),
     1                tabf,tabfx,tabfxx,boxsize(ilev))
                endif
            endif
         enddo
C$OMP END PARALLEL DO         
 3000 continue
 3001 continue

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(7) = time2-time1

      if(ifprint.eq.1) call prin2('timeinfo=*',timeinfo,6)

      d= 0
      do i = 1,7
         d = d + timeinfo(i)
      enddo

      if(ifprint.eq.1) call prin2('sum(timeinfo)=*',d,1)
cccc      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
c
c
c
      subroutine periodic_bfgtmain_large_delta(
     1    nd,ndim,delta,eps,ipoly,norder,npbox,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     3    iaddr,rmlexp,npw,pot,timeinfo)
c
c     This code computes doubly periodic Gauss transform
c     to a collection of right hand sides when delta is very large
c 
c     method: in this case, the periodic Green's function admits an efficient
c     seprable Fourier series expansion, which is used in the entire computational
c     box. 1. compute planewave expansion coefficients for each leaf box; 2. merge
c     pw exp. all the way to the root box - this is a valid expansion in the entire
c     box; 3. split the pw exp. into its child boxes; 4. evaluate the pw exp. in each
c     leaf box.
c 
c 
c       input
c
c         nd - integer,   number of box FGTs with the same tree
c         delta - double precision, Gaussian variance 
c         eps - double precision
c            tolerance requested
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
c         norder - integer
c           order of expansions for input coefficients array
c         ncbox - integer
c           number of coefficients of expansions of functions
c           in each of the boxes
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fvals - double complex (npbox,nboxes)
c            function values tabulated on the tree
c         centers - double precision (2,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
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
c
c
c         npwlevel - integer
c             cutoff level at which the PW expansion is valid
c         npw  - integer,  length of planewave expansions
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c
c     output:
c         pot - double precision (nd,npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      integer nd
      real *8 delta,eps
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree),ncbox,npbox
      real *8 fvals(nd,npbox,nboxes)
      real *8 pot(nd,npbox,nboxes)
      real *8 boxsize(0:nlevels),centers(2,nboxes)
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 timeinfo(*)
c
      real *8, allocatable :: fcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)

      real *8 timelev(0:200)
      real *8 ws(200),ts(200)
      
      complex *16, allocatable :: wpwmsshift(:,:,:)
      
      complex *16, allocatable :: tab_poly2pw(:,:,:),tab_pw2pot(:,:,:)
      complex *16, allocatable :: ff(:,:)
      complex *16, allocatable :: gg(:,:)

      ifprint = 0

      done = 1
      pi = atan(done)*4
c
c       compute coefs
c
      allocate(fcoefs(npbox,nd,nboxes),xq(norder),umat(norder,norder),
     1   vmat(norder,norder),wts(norder))
     
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)

      do ilev = 0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             do ind=1,nd
                do i=1,npbox
                   fcoefs(i,ind,ibox)=fvals(ind,i,ibox)
                   if (abs(fvals(ind,i,ibox)).lt.1d-16) 
     1                 fcoefs(i,ind,ibox)=0
                enddo
             enddo
          endif
        enddo
      enddo

      do i=1,6
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo

c
c     get planewave nodes and weights
      call get_periodic_pwnodes(boxsize(0),delta,eps,npw,ws,ts)

c     compute the tables converting Chebyshev polynomial expansion
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
      allocate(gg(norder,(npw+1)/2))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pot_tables(norder,npw,ts,xq,1.0d0,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev))      
      enddo
      
      nexp=1
      do i=1,ndim-1
         nexp = nexp*npw
      enddo
      nexp=nexp*((npw+1)/2)
      
      nmax = nlevels
      allocate(wpwmsshift(nexp,4,nmax))
      xmin2 = boxsize(nlevels)/2
      call mk_merge_split_pw_matrices(xmin2,npw,ts,nmax,wpwmsshift)

      if(ifprint.eq.1)
     1    call prinf('laddr=*',itree(iptr(1)),2*(nlevels+1))
      
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
c
c        step 1: convert coeffs to planewave expansions
c
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0)
      
      do 1100 ilev = nlevels,0,-1
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
c           Check if current box is a leaf box            
            if(nchild.eq.0) then
c              form PW expansion directly
               call gnd_tens_prod_to_pw(nd,norder,fcoefs(1,1,ibox),npw,
     1             ff,tab_poly2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
            endif
         enddo
C$OMP END PARALLEL DO
 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f5.2)         
cccc         write(6,111) ilev,nb,dt
c     end of ilev do loop
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
            jbox = itree(iptr(5)+4*(ibox-1)+i-1)

            dx= (centers(1,ibox) - centers(1,jbox))
            dy= (centers(2,ibox) - centers(2,jbox))
            if (dx.gt.0 .and. dy.gt.0) then
               k=1
            elseif (dx.gt.0 .and. dy.lt.0) then
               k=2
            elseif (dx.lt.0 .and. dy.gt.0) then
               k=3
            elseif (dx.lt.0 .and. dy.lt.0) then
               k=4
            endif
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
      call gnd_copypwexp(nd,nexp,rmlexp(iaddr(1,ibox)),
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
C$OMP$PRIVATE(ibox,jbox,i,nchild,dx,dy,dz,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4)+ibox-1)
          do i=1,nchild
             jbox = itree(iptr(5)+4*(ibox-1)+i-1)
             dx= centers(1,jbox) - centers(1,ibox)
             dy= centers(2,jbox) - centers(2,ibox)
             if (dx.gt.0 .and. dy.gt.0) then
                k=1
             elseif (dx.gt.0 .and. dy.lt.0) then
                k=2
             elseif (dx.lt.0 .and. dy.gt.0) then
                k=3
             elseif (dx.lt.0 .and. dy.lt.0) then
                k=4
             endif
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
               call gnd_pw2pot(nd,norder,npw,rmlexp(iaddr(2,ibox)),
     1             gg,tab_pw2pot(1,1,ilev),pot(1,1,ibox))
            endif
ccc    end of ibox loop        
         enddo
 222     format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f5.2)         
cccc        write(6,222) ilev,nb,dt
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
cccc      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

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
