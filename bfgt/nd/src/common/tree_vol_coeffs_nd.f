c
c
c    adaptive tree in ndim=1,2,3 dimensions with ipoly=0 for Legendre polynomial
c    and ipoly=1 for Chebyshev polynomials
c
c    Created by Shidong Jiang on 03/31/2022
c
c
c
c    generate level restricted tree based on resolving a 
c    function to desired precision
c
c
c    The function handle is of the form
c    call fun(nd,xyz,dpars,zpars,ipars,f)
c
c      where xyz is the location in (-L/2,L/2)^d
c
c    dpars are a collection of real parameters, zpars are complex
c    parameters and ipars are integer parameters, and f is a 
c    real array of size nd
c     
c
c    For the Helmholtz/Maxwell tree, the boxes are refined until 
c     Re(zk)*boxsize<5, beyond which the function resolution criterion
c    kicks in. 
c
c    A function is said to be resolved if it's interpolant at the 8
c    children nodes, agrees with the function values at those nodes
c    upto the user specified tolerance. 
c    The error is scaled by h**(eta)
c    where eta is user specified and h is the boxsize. If there is 
c    any confusion, the user should seet \eta to 0
c    Let \tilde{f} denote the interpolant of f, then
c    the refinement criterion is 
c      \int_{B_{j} |\tilde{f}-f|^{p} *h^{\eta} < 
c        \varepsilon V_{j}^{1/p}/V_{0}^{1/p}/(\int_{B_{0}}|f|^{p})^{1/p}
c    This implies that
c       \int_{B_{0}} |\tilde{f}-f|^{p} =
c          \sum_{j} \int_{B_{j}} |\tilde{f}-f|^{p}
c          \leq \sum_{j} \eps^{p}*h^{\eta p} 
c                  V_{j}/V_{0}/(\int_{B_{0}} |f|^p)
c      If \eta = 0,
c          \leq \eps^{p}/(\int_{B_{0}} |f|^{p})
c
c    i.e., this strategy guarantees that the interpolated function
c      approximates the function with relative lp accuracy of \eps
c      
c    This code has 2 main user callable routines
c      make_vol_tree_mem -> Returns the memory requirements, 
c          tree length, number of boxes, number of levels
c      make_vol_tree -> Makes the actual tree, returns centers of boxes,
c          colleague info, function values on leaf boxes
c       
c          
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c




      subroutine vol_tree_mem(ndim,ipoly,eps,zk,boxlen,norder,iptype,
     1    eta,fun,nd,dpars,zpars,ipars,nboxes,nlevels,ltree,rintl)
c
c      get memory requirements for the tree
c
c      input parameters:
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer
c           order of discretization
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c 
c        output:
c           nlevels - integer
c             number of levels
c           nboxes - integer
c             number of boxes
c           nlboxes - integer
c             number of leaf boxes
c           ltree - integer
c             length of tree
c           rintl(0:nlevels) - real *8
c             lp norm to scale the functions by
c             (on input rintl should be of size(0:200)
c
c      
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*),delta
      real *8, allocatable :: fvals(:,:,:)
      complex *16 zpars(*),zk
      integer nd,ipars(*),iptype
      integer ltree
      integer ndim,ipoly,nlevels,nboxes,norder

      external fun

      integer, allocatable :: laddr(:,:),ilevel(:),iparent(:),nchild(:)
      integer, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      real *8, allocatable :: centers(:,:)
      integer, allocatable :: nbors(:,:),nnbors(:),ifrefine(:)

      integer, allocatable :: ilevel2(:),iparent2(:),nchild2(:),
     1    ichild2(:,:)
      real *8, allocatable :: centers2(:,:),fvals2(:,:,:)

      integer nbmax,nlmax,npbox,npc,mc,mnbors
      real *8, allocatable :: grid(:,:),qwts(:)
      real *8, allocatable:: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8 xy(ndim)
      real *8, allocatable :: wts2(:),xref2(:,:)
      real *8 rintl(0:200)
      real *8 rint
      real *8, allocatable :: rintbs(:),rintbs2(:)
      integer i,itype,j,npols

      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      real *8, allocatable :: boxsize(:)
      integer, allocatable :: irefinebox(:)

      real *8 rsc,ra,utmp,vtmp
      integer nbloc,nbctr,nbadd,irefine,ilev,ifirstbox,ilastbox
      integer nbtot,iii,idim,iper,isep,nrefine

      
      nbmax = 100000
      nlmax = 200

      allocate(boxsize(0:nlmax))

      mc=2**ndim
      mnbors=3**ndim
      npbox = norder**ndim
      
      allocate(laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax))
      allocate(nchild(nbmax),ichild(mc,nbmax))

      allocate(fvals(nd,npbox,nbmax),centers(ndim,nbmax))

      allocate(rintbs(nbmax))

c
c      set tree info for level 0
c
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparent(1) = -1
      nchild(1) = 0
      do i=1,mc
        ichild(i,1) = -1
      enddo

      do i=1,ndim
         centers(i,1) = 0
      enddo
c
c

      allocate(grid(ndim,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif

      npols = norder**ndim
      allocate(wts2(npbox),xref2(ndim,npbox))

      itype = 1
      call polytens_exps_nd(ndim,ipoly,itype,norder,'f',xref2,
     1    utmp,1,vtmp,1,wts2)
c
c       compute fvals at the grid
c

      boxsize(0) =boxlen

      rint = 0
      rintbs(1) = 0


c
c   note extra factor of 4 sincee wts2 are on [-1,1]^2 
c   as opposed to [-1/2,1/2]^2
c
      rsc = boxlen**2/mc

      do i=1,npbox
        do j=1,ndim
           xy(j) = grid(j,i)*boxlen
        enddo
        call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,1))
        if(iptype.eq.0) then
          do idim=1,nd
            if(abs(fvals(idim,i,1)).gt.rintbs(1)) rintbs(1) = 
     1          abs(fvals(idim,i,1))
          enddo
        endif

        if(iptype.eq.1) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + abs(fvals(idim,i,1))*wts2(i)*rsc
          enddo
        endif

        if(iptype.eq.2) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + fvals(idim,i,1)**2*wts2(i)*rsc
          enddo
        endif
      enddo

      if(iptype.eq.0.or.iptype.eq.1) rint = rintbs(1)
      if(iptype.eq.2) rint = sqrt(rintbs(1))

      rintl(0) = rint

      nbctr = 1



     

      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,npbox,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
        if(iptype.eq.0) rsc = 1.0d0

        rsc = rsc*rint

        print *, ilev, rint,rsc

        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,boxsize(ilev),nbmax,
     2      ifirstbox,nbloc,rsc,irefinebox,irefine)
     

c
c
c          figure out if current set of boxes is sufficient
c

        nbadd = 0 
        do i=1,nbloc
          if(irefinebox(i).eq.1) nbadd = nbadd+mc
        enddo

        nbtot = nbctr+nbadd

c
c         if current memory is not sufficient reallocate
c
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(ndim,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(mc,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call vol_tree_copy(ndim,nd,nbctr,npbox,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)
          call dcopy_f77(nbctr,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(ndim,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(mc,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call vol_tree_copy(ndim,nd,nbctr,npbox,centers2,ilevel2,
     1        iparent2,nchild2,ichild2,fvals2,centers,ilevel,iparent,
     2        nchild,ichild,fvals)
          call dcopy_f77(nbctr,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif


        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          laddr(1,ilev+1) = nbctr+1

          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nbmax,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,ilevel,iparent,nchild,ichild)


          rsc = boxsize(ilev+1)**ndim/mc
          call update_rints(ndim,nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts2,rsc,rintbs,rint)
          
          rintl(ilev+1) = rint

          
          laddr(2,ilev+1) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      nboxes = nbctr
      nlevels = ilev

      if(nlevels.ge.2) then

        nbtot = 2*mc*nboxes
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(ndim,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(mc,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call vol_tree_copy(ndim,nd,nboxes,npbox,centers,ilevel,
     1       iparent,nchild,ichild,fvals,centers2,ilevel2,
     2       iparent2,nchild2,ichild2,fvals2)
          call dcopy_f77(nboxes,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(ndim,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(mc,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call vol_tree_copy(ndim,nd,nboxes,npbox,centers2,ilevel2,
     1        iparent2,nchild2,ichild2,fvals2,centers,ilevel,
     2        iparent,nchild,ichild,fvals)
          call dcopy_f77(nboxes,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif

        allocate(nnbors(nbmax))
        allocate(nbors(mnbors,nbmax))

        do i=1,nboxes
          nnbors(i) = 0
          do j=1,mnbors
            nbors(j,i) = -1
          enddo
        enddo

        iper = 0

        call computecoll(ndim,nlevels,nboxes,laddr,boxsize,centers,
     1        iparent,nchild,ichild,iper,nnbors,nbors)

        if(nlevels.ge.2) then
          call vol_tree_fix_lr(ndim,fun,nd,dpars,zpars,ipars,
     1      norder,npbox,fvals,grid,centers,nlevels,nboxes,boxsize,
     2      nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors)
        endif
      endif

c     double the number of boxes for possible refinement
c     added by Shidong Jiang 04/13/2022
      nboxes = nboxes*2
      
      ltree = (4+mc+mnbors)*nboxes + 2*(nlevels+1)

      return
      end
c
c
c
c
c

      subroutine vol_tree_build(ndim,ipoly,eps,zk,boxlen,norder,iptype,
     1    eta,fun,nd,dpars,zpars,ipars,rintl,nboxes,nlevels,ltree,itree,
     2    iptr,centers,boxsize,fvals)
c
c      compute the tree
c
c      input parameters:
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer
c           order of discretization
c        iptype - integer
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c        nlevels - integer
c          number of levels
c        nboxes - integer
c          number of boxes
c        ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c        rintl - real *8 (0:nlevels)
c          estimate of lp norm for scaling the errors
c          at various levels. 
c          We require the estimate at each level to make sure
c          that the memory estimate code is consitent
c          with the build code else there could be potential
c          memory issues 
c         
c
c      output:
c        itree - integer(ltree)
c          tree info
c        iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c        fvals - double precision (nd,norder**2,nboxes)
c          function values at discretization nodes
c        centers - double precision (2,nboxes)
c          xyz coordinates of box centers in the oct tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      complex *16 zk,zpars(*)
      integer ndim,ipoly,nd,ipars(*),iptype
      integer nlevels,nboxes,norder
      integer iptr(8),ltree
      integer itree(ltree),ier
      real *8 fvals(nd,norder**ndim,nboxes),centers(ndim,nboxes)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:),ifrefine(:)
      real *8 boxsize(0:nlevels)
c      real *8 xq(norder),wts(norder),umat(norder,norder)
c      real *8 vmat(norder,norder)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8 rintl(0:nlevels)
      real *8 xy(ndim)

      integer i,ilev,irefine,itype,nbmax,nlmax,npbox,npc,ii
      integer ifirstbox,ilastbox,nbctr,nbloc
      real *8 rsc

      real *8 ra
      integer j,nboxes0,npols,iper,mc,mnbors,nrefine

      external fun

      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
c
      mc=2**ndim
      mnbors=3**ndim
      
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + mc*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + mnbors*nboxes



      boxsize(0) = boxlen

      do i=1,ndim
         centers(i,1) = 0
      enddo

c
c      set tree info for level 0
c
      itree(1) = 1
      itree(2) = 1
      itree(iptr(2)) = 0
      itree(iptr(3)) = -1
      itree(iptr(4)) = 0
      do i=1,mc
        itree(iptr(5)+i-1) = -1
      enddo
c
c
      npbox = norder**ndim
      allocate(grid(ndim,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
      
      npols = norder**ndim

      do i=1,npbox
         do j=1,ndim
            xy(j) = grid(j,i)*boxlen
         enddo
        call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,1))
      enddo


c
c       Reset nlevels, nboxes
c
      nbctr = 1

      do ilev=0,nlevels-1
        irefine = 0

        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,npols,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
        if(iptype.eq.0) rsc = 1.0d0
        rsc = rsc*rintl(ilev)
        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,boxsize(ilev),nboxes,
     2      ifirstbox,nbloc,rsc,irefinebox,irefine)
        

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,itree(iptr(2)),itree(iptr(3)),
     3      itree(iptr(4)),itree(iptr(5)))
          
          itree(2*ilev+4) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      nboxes0 = nbctr
      nlevels = ilev


      do i=1,nboxes0
         itree(iptr(6)+i-1) = 0
         do j=1,mnbors
            itree(iptr(7)+mnbors*(i-1)+j-1) = -1
         enddo
      enddo

      iper = 0
      call computecoll(ndim,nlevels,nboxes0,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iper,itree(iptr(6)),itree(iptr(7)))

      if(nlevels.ge.2) then
         call vol_tree_fix_lr(ndim,fun,nd,dpars,zpars,ipars,
     1       norder,npbox,fvals,grid,centers,nlevels,nboxes0,boxsize,
     2       nboxes,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       itree(iptr(6)),itree(iptr(7)))
      endif

      call prinf('nboxes in tree_build=*',nboxes,1)

cccc      call prinf('nboxes0=*',nboxes0,1)
cccc      call prinf('nlevels=*',nlevels,1)
      return
      end
c
c
c      
c
c
      subroutine vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,
     1  zk,norder,npbox,fvals,npols,umat,
     2  boxsize,nboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
      implicit none
      integer ndim,nd,npc,npbox,norder,iptype,npols
      integer nboxes,nbloc
      real *8 eta,eps,fvals(nd,npbox,nboxes)
      real *8 umat(norder,norder)
      real *8 rsc
      real *8, allocatable :: fcoefs(:,:,:),rmask(:)
      integer, allocatable :: iind2p(:,:)
      real *8 alpha,beta,boxsize,rsum
      integer irefinebox(nbloc)
      
      complex *16 zk
      integer ifirstbox

      integer irefine

      integer i,j,k,l,ibox,ifunif,i1
      real *8 rscale2,erra,bs,bs2

      character *1 transa,transb

      external fun

      ifunif = 0

      transa = 'n'
      transb = 't'


      allocate(fcoefs(nd,npols,nbloc))
      allocate(iind2p(ndim,npols))
      call polytens_ind2pow(ndim,norder-1,'f',iind2p)

      
      allocate(rmask(npols))

      rsum = 0
      do i=1,npols
        rmask(i) = 0.0d0
        i1=0
        do k=1,ndim
           i1=i1+iind2p(k,i)
        enddo
        if(i1.eq.norder-1) then
          rmask(i) = 1.0d0
          rsum = rsum + 1
        endif
      enddo

      if(iptype.eq.2) rsum = sqrt(rsum)
      if(iptype.eq.0) rsum = 1
      



      alpha = 1
      beta = 0

      bs = boxsize/4.0d0
      bs2 = 2*bs
      rscale2 = bs2**eta

      if(real(zk)*boxsize.gt.5) then
        do i=1,nbloc
          irefinebox(i) = 1
        enddo
        goto 1000
      endif


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,erra)
      do i=1,nbloc
        irefinebox(i) = 0

        ibox = ifirstbox + i-1

        call tens_prod_trans_nd(ndim,nd,norder,fvals(1,1,ibox),
     1      fcoefs(1,1,i),umat)
        
        call fun_err(nd,npols,fcoefs(1,1,i),rmask,
     1     iptype,rscale2,erra)
     
        erra = erra/rsum

        
        if(erra.gt.eps*rsc) then
          irefinebox(i) = 1
        endif
      enddo
C$OMP END PARALLEL DO     

 1000 continue
      irefine = maxval(irefinebox(1:nbloc))


c
c       make tree uniform
c

      if(ifunif.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nbloc
          irefinebox(i) = irefine
        enddo
C$OMP END PARALLEL DO      
      endif

      return
      end
c
c
c
c
c


      subroutine vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer ndim,nd,npbox
      integer ipars(*)
      real *8 dpars(*)
      complex *16 zpars(*)
      real *8 fvals(nd,npbox,nboxes)
      integer nboxes,nbloc,nbctr,nlctr
      real *8 grid(ndim,npbox),bs,xyz(ndim)
      real *8 centers(ndim,nboxes)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(2**ndim,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl,k,mc

      real *8 bsh
      integer isgn(ndim,2**ndim)
      external fun
      
      call get_child_box_sign(ndim,isgn)

      mc=2**ndim

      allocate(isum(nbloc))
      call cumsum(nbloc,irefinebox,isum)

      bsh = bs/2

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,l,xyz,k)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*mc
          nchild(ibox) = mc
          do j=1,mc
            jbox = nbl+j
            do k=1,ndim
               centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo

            do l=1,npbox
              do k=1,ndim
                 xyz(k) = centers(k,jbox) + grid(k,l)*bs
              enddo
              call fun(nd,xyz,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,mc
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      nbctr = nbctr + isum(nbloc)*mc


      return
      end
c
c
c
c
c
      subroutine fun_err(nd,n,fcoefs,rmask,iptype,rscale,erra)
c       this subroutine estimates the error based on the expansion
c       coefficients in a given basis
c       
c       input
c        nd - integer
c          number of functions
c        n -  integer
c           number of points at which function is tabulated
c        fcoefs: double precision(nd,n) 
c          tensor product legendre coeffs of func
c        rmask: double precision(n)
c           coefficients which will be accounted for in the error
c        iptype: integer
c           type of error to be computed
c           iptype = 0 - linf error
c           iptype = 1 - l1 error
c           iptype = 2 - l2 error
c        rscale: double precision
c          scaling factor
c
c       output
c         err: double precision
c           max scaled error in the functions
c
        implicit none
        integer n,i,iptype,idim,nd
        real *8 rscale,erra
        real *8 fcoefs(nd,n),rmask(n)
        real *8, allocatable :: errtmp(:),ftmp(:,:)
        real *8 alpha, beta

        allocate(errtmp(nd),ftmp(nd,n))

        alpha = 1.0d0
        beta = 0.0d0
   
        erra = 0
        do idim=1,nd
          errtmp(idim) = 0
        enddo
        if(iptype.eq.0) then
          do i=1,n
            if(rmask(i).gt.0.5d0) then
              do idim = 1,nd
                if(errtmp(idim).lt.abs(fcoefs(idim,i))) 
     1             errtmp(idim)=abs(fcoefs(idim,i))
              enddo
            endif
          enddo
        endif

        if(iptype.eq.1) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = abs(fcoefs(idim,i))
            enddo
          enddo
          call dgemv_f77('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)
        endif


        if(iptype.eq.2) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = fcoefs(idim,i)**2
            enddo
          enddo
          call dgemv_f77('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)

          do idim=1,nd
            errtmp(idim) = sqrt(errtmp(idim))
          enddo
        endif

        erra = 0
        do idim=1,nd
          if(errtmp(idim).gt.erra) erra = errtmp(idim)
        enddo

        erra = erra*rscale

        return
        end

c
c
c
c       
c
      subroutine update_rints(ndim,nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts,rsc,rintbs,rint)
c
c------------------------
c  This subroutine updates the integrals of the function to
c  be resolved on the computational domain. It subtracts the
c  integral of the boxes which have been refined and adds
c  in the integrals corresponding to the function values 
c  tabulated at the children
c
c  Input arguments:
c  
c    - nd: integer
c        number of functions
c    - npbox: integer
c        number of points per box where the function is tabulated
c    - nbmax: integer
c        max number of boxes
c    - fvals: real *8 (nd,npbox,nbmax)
c        tabulated function values
c    - ifirstbox: integer
c        first box in the list of boxes to be processed
c    - nbloc: integer
c        number of boxes to be processed
c    - iptype: integer
c        Lp version of the scheme
c        * iptype = 0, linf
c        * iptype = 1, l1
c        * iptype = 2, l2
c    - nchild: integer(nbmax)
c        number of children 
c    - ichild: integer(8,nbmax)
c        list of children
c    - wts: real *8 (npbox)
c        quadrature weights for intgegrating functions 
c    - rsc: real *8
c        scaling parameter for computing integrals
c  
c  Inout arguemnts:
c
c     - rintbs: real *8(nbmax)
c         the integral for the new boxes cretated will be updated
c     - rint: real *8
c         the total integral will be updated
c    
c  
c      
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndim,nd,npbox,nbmax
      real *8, intent(in) :: fvals(nd,npbox,nbmax)
      integer, intent(in) :: ifirstbox,nbloc,iptype
      integer, intent(in) :: nchild(nbmax),ichild(2**ndim,nbmax)
      real *8, intent(in) :: wts(npbox),rsc
      real *8, intent(inout) :: rintbs(nbmax),rint

      integer mc

      mc = 2**ndim
c
c
c      compute the integrals for the newly formed boxes
c   and update the overall integral
c
      if(iptype.eq.0) then
        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,mc
              jbox = ichild(j,ibox)
              rintbs(jbox) = maxval(fvals(1:nd,1:npbox,jbox))
              if(rintbs(jbox).gt.rint) rint = rintbs(jbox)
            enddo
          endif
        enddo
      endif

      if(iptype.eq.1) then
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c     subtract contribution of ibox from rint
            rint = rint - rintbs(ibox) 
          endif
        enddo

c
c     add back contribution of children
c
        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,mc
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               abs(fvals(idim,l,jbox))*wts(l)*rsc
                enddo
              enddo
              rint = rint + rintbs(jbox)
            enddo
          endif
        enddo
      endif

      if(iptype.eq.2) then
        rintsq = rint**2
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c
c    note that if iptype = 2, then rintbs stores squares
c    of the integral on the box
c
             rintsq = rintsq - rintbs(ibox)
          endif
        enddo

        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,mc
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               fvals(idim,l,jbox)**2*wts(l)*rsc
                enddo
              enddo
              rintsq = rintsq + rintbs(jbox)
            enddo
          endif
        enddo
        rint = sqrt(rintsq)
      endif
          

      return
      end
c
c
c
c
      subroutine get_children_qwts(ndim,norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wts(norder),qwts(npc)

      if (ndim.eq.1) call get_children_qwts_1d(norder,npc,wts,qwts)
      if (ndim.eq.2) call get_children_qwts_2d(norder,npc,wts,qwts)
      if (ndim.eq.3) call get_children_qwts_3d(norder,npc,wts,qwts)

      return
      end
c
c
c
c
c
c
c
      subroutine get_children_qwts_1d(norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wts(norder),qwts(npc)
      
      do i=1,norder
         qwts(i) = wts(i)/2
      enddo

      do j=1,norder
         qwts(norder+j) = qwts(j)
      enddo


      return
      end
c
c
c
c
c
c
c
      subroutine get_children_qwts_2d(norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wts(norder),qwts(npc)
      
      ipt = 1
      do i=1,norder
        do j=1,norder
           qwts(ipt) = wts(i)*wts(j)/4
           ipt = ipt+1
        enddo
      enddo

      npbox = norder**2
      do i=1,3
        do j=1,npbox
          qwts(i*npbox+j) = qwts(j)
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
c
c
c
c
      subroutine get_children_qwts_3d(norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wts(norder),qwts(npc)
      
      ipt = 1
      do i=1,norder
        do j=1,norder
          do k=1,norder
            qwts(ipt) = wts(i)*wts(j)*wts(k)/8
            ipt = ipt+1
          enddo
        enddo
      enddo

      npbox = norder**3
      do i=1,7
        do j=1,npbox
          qwts(i*npbox+j) = qwts(j)
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
c     
       subroutine vol_tree_copy(ndim,nd,nb,npb,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)

       implicit none
       integer ndim,nd,nb,npb
       real *8 centers(ndim,nb),centers2(ndim,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(2**ndim,nb),ichild2(2**ndim,nb)
       real *8 fvals(nd,npb,nb),fvals2(nd,npb,nb)

       integer i,j,nel,incx,incy,mc,k

       mc=2**ndim
       
       nel = nd*npb*nb
       incx = 1
       incy = 1
       call dcopy_f77(nel,fvals,incx,fvals2,incy)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         do j=1,ndim
            centers2(j,i) = centers(j,i)
         enddo
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         do j=1,mc
           ichild2(j,i) = ichild(j,i)
         enddo
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c
c-------------------------------------------------------------      
      subroutine vol_tree_fix_lr(ndim,fun,nd,dpars,zpars,ipars,norder,
     1       npbox,fvals,grid,centers,nlevels,nboxes,boxsize,
     2       nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,
     3       nbors)
c
c
c       convert an adaptive tree into a level restricted tree
c
      implicit none
      integer nd,ipars(*),norder,npbox,nlevels,nboxes,nlmax
      integer nbmax,ndim
      real *8 dpars(*),fvals(nd,npbox,nbmax),grid(ndim,npbox)
      real *8 centers(ndim,nbmax),boxsize(0:nlmax)
      complex *16 zpars(*)
      integer laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax)
      integer nchild(nbmax),ichild(2**ndim,nbmax),nnbors(nbmax)
      integer nbors(3**ndim,nbmax)
      integer laddrtail(2,0:nlmax),isum
      integer, allocatable :: iflag(:)

      integer i,j,k,l,ibox,jbox,kbox,ilev,idad,igranddad
      integer nbloc,ict,iper,mc,mnbors,ifnbor
      real *8 dis,distest

      external fun

      mc=2**ndim
      mnbors=3**ndim
      
      allocate(iflag(nbmax))

c     Initialize flag array
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         iflag(i) = 0
      enddo
C$OMP END PARALLEL DO     



c     Flag boxes that violate level restriction by "1"
c     Violatioin refers to any box that is directly touching
c     a box that is more than one level finer
c
c     Method:
c     1) Carry out upward pass. For each box B, look at
c     the colleagues of B's grandparent
c     2) See if any of those colleagues are childless and in
c     contact with B.
c
c     Note that we only need to get up to level two, as
c     we will not find a violation at level 0 and level 1
c
c     For such boxes, we set iflag(i) = 1
c
      do ilev=nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by two levels are touching
         distest = 1.05d0*(boxsize(ilev-1) + boxsize(ilev-2))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,igranddad,i,jbox)         
C$OMP$ PRIVATE(ict,dis,k)
         do ibox = laddr(1,ilev),laddr(2,ilev) 
            idad = iparent(ibox)
            igranddad = iparent(idad)
            
c           Loop over colleagues of granddad            
            do i=1,nnbors(igranddad)
               jbox = nbors(i,igranddad)
c              Check if the colleague of grandad
c              is a leaf node. This automatically
c              eliminates the granddad
               if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                  ict = 0
                  do k=1,ndim
                     dis = centers(k,jbox) - centers(k,idad)
                     if(abs(dis).le.distest) ict = ict + 1
                  enddo
                  if(ict.eq.ndim) then
                     iflag(jbox) = 1
                  endif
               endif
c              End of checking criteria for the colleague of
c              granddad
            enddo
c           End of looping over colleagues of
c           granddad
         enddo
c        End of looping over boxes at ilev         
C$OMP END PARALLEL DO
      enddo
c     End of looping over levels and flagging boxes


c     Find all boxes that need to be given a flag+
c     A flag+ box will be denoted by setting iflag(box) = 2
c     This refers to any box that is not already flagged and
c     is bigger than and is contacting a flagged box
c     or another box that has already been given a flag +.
c     It is found by performing an upward pass and looking
c     at the flagged box's parents colleagues and a flag+
c     box's parents colleagues and seeing if they are
c     childless and present the case where a bigger box 
c     is contacting a flagged or flag+ box.

      do ilev = nlevels,1,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,i,jbox,dis)
C$OMP$PRIVATE(ict,k)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if(iflag(ibox).eq.1.or.iflag(ibox).eq.2) then
               idad = iparent(ibox)
c              Loop over dad's colleagues               
               do i=1,nnbors(idad)
                  jbox = nbors(i,idad)
c                 Check if the colleague of dad
c                 is a leaf node. This automatically
c                 eliminates the dad
                  if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                     ict = 0
                     do k=1,ndim
                        dis = centers(k,jbox) - centers(k,ibox)
                        if(abs(dis).le.distest) ict = ict + 1
                     enddo
                     if(ict.eq.ndim) then
                        iflag(jbox) = 2
                     endif
                  endif
c                 End of checking criteria for the colleague of
c                dad
               enddo
c              End of looping over dad's colleagues               
            endif
c           End of checking if current box is relevant for
c           flagging flag+ boxes
         enddo
c        End of looping over boxes at ilev        
C$OMP END PARALLEL DO 
      enddo
c     End of looping over levels

c     Subdivide all flag and flag+ boxes. Flag all the children
c     of flagged boxes as flag++. Flag++ boxes are denoted
c     by setting iflag(box) = 3. The flag++ boxes need 
c     to be checked later to see which of them need further
c     refinement. While creating new boxes, we will
c     need to update all the tree structures as well.
c     Note that all the flagged boxes live between
c     levels 1 and nlevels - 2. We process the boxes via a
c     downward pass. We first determine the number of boxes
c     that are going to be subdivided at each level and 
c     everything else accordingly
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo
 
      do ilev = 1,nlevels-2
c        First subdivide all the flag and flag+
c        boxes with boxno nboxes+1, nboxes+ 2
c        and so on. In the second step, we reorganize
c        all the structures again to bring it back
c        in the standard format

         laddrtail(1,ilev+1) = nboxes+1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call vol_tree_refine_boxes_flag(ndim,iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call vol_tree_reorg(ndim,nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO     
      iper = 0
      call computecoll(ndim,nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)

c     Processing of flag and flag+ boxes is done
c     Start processing flag++ boxes. We will use a similar
c     strategy as before. We keep checking the flag++
c     boxes that require subdivision if they still
c     violate the level restriction criterion, create
c     the new boxes, append them to the end of the list to begin
c     with and in the end reorganize the tree structure.
c     We shall accomplish this via a downward pass
c     as new boxes that get added in the downward pass
c     will also be processed simultaneously.
c     We shall additionally also need to keep on updating
c     the colleague information as we proceed in the 
c     downward pass

c     Reset the flags array to remove all the flag and flag+
c     cases. This is to ensure reusability of the subdivide
c     _flag routine to handle the flag++ case

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox)
      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
C$OMP END PARALLEL DO      
 
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo

      do ilev = 2,nlevels-2

c     Step 1: Determine which of the flag++ boxes need
c     further division. In the even a flag++ box needs
c     further subdivision then flag the box with iflag(box) = 1
c     This will again ensure that the subdivide_flag routine
c     will take care of handling the flag++ case
         call vol_updateflags(ndim,ilev,nboxes,nlevels,laddr,
     1       nchild,ichild,nnbors,nbors,centers,boxsize,iflag)

         call vol_updateflags(ndim,ilev,nboxes,nlevels,laddrtail,
     1       nchild,ichild,nnbors,nbors,centers,boxsize,iflag)
         
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call vol_tree_refine_boxes_flag(ndim,iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         nbloc = laddrtail(2,ilev)-laddrtail(1,ilev)+1
         call vol_tree_refine_boxes_flag(ndim,iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddrtail(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

          laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,idad,jbox,j,kbox)
C$OMP$PRIVATE(k,ifnbor)
          do ibox = laddrtail(1,ilev+1),laddrtail(2,ilev+1)
            nnbors(ibox) = 0
c           Find the parent of the current box         
            idad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out colleagues
            do i=1,nnbors(idad)
                jbox = nbors(i,idad)
                do j=1,mc
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c     Check if kbox is a nearest neighbor or in list 2
                      ifnbor=1
                      do k=1,ndim
                         if((abs(centers(k,kbox)-centers(k,ibox)).gt.
     1                   1.05*boxsize(ilev+1))) then
                            ifnbor=0
                            exit
                         endif
                      enddo
                      if (ifnbor.eq.1) then
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

c     Reorganize tree once again and we are all done      
      call vol_tree_reorg(ndim,nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO    

      call computecoll(ndim,nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)
      

      return
      end
      

c-------------------------------------------------------------      
      subroutine vol_tree_reorg(ndim,nboxes,nd,npbox,centers,
     1    nlevels,laddr,laddrtail,ilevel,iparent,nchild,ichild,
     2    fvals,iflag)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    nd             in: integer
c                   number of real value functions
c
c    npbox          in: integer
c                   number of grid points per function
c
c    centers        in/out: double precision(3,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    iflag        in/out: integer(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer ndim,nboxes,nlevels,npbox,nd
      double precision centers(ndim,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer laddrtail(2,0:nlevels)
      integer ilevel(nboxes)
      integer iparent(nboxes)
      integer nchild(nboxes)
      integer ichild(2**ndim,nboxes)
      integer iflag(nboxes)
      double precision fvals(nd,npbox,nboxes)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:),tiflag(:)
      integer, allocatable :: iboxtocurbox(:),ilevptr(:),ilevptr2(:)

      double precision, allocatable :: tfvals(:,:,:),tcenters(:,:)



c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox,idim,nblev,mc

      mc=2**ndim
      
      allocate(tilevel(nboxes),tiparent(nboxes),tnchild(nboxes))
      allocate(tichild(mc,nboxes),tiflag(nboxes),iboxtocurbox(nboxes))
      allocate(tfvals(nd,npbox,nboxes),tcenters(ndim,nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo
      call vol_tree_copy(ndim,nd,nboxes,npbox,centers,ilevel,
     1    iparent,nchild,ichild,fvals,tcenters,tilevel,
     2    tiparent,tnchild,tichild,tfvals)

      do ibox=1,nboxes
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
         enddo
      enddo

      allocate(ilevptr(nlevels+1),ilevptr2(nlevels))

      ilevptr(2) = laddr(1,2)


      do ilev=2,nlevels
        nblev = laddr(2,ilev)-laddr(1,ilev)+1
        ilevptr2(ilev) = ilevptr(ilev) + nblev
        nblev = laddrtail(2,ilev)-laddrtail(1,ilev)+1
        ilevptr(ilev+1) = ilevptr2(ilev) + nblev
      enddo

      curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,curbox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo
            nchild(curbox) = tnchild(ibox)
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,curbox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         do i=1,mc
            if(tichild(i,ibox).eq.-1) ichild(i,iboxtocurbox(ibox)) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo

      return
      end
c--------------------------------------------------------------------      
      subroutine vol_updateflags(ndim,curlev,nboxes,nlevels,laddr,
     1    nchild,ichild,nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(3,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer ndim,curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(2**ndim,nboxes)
      integer nnbors(nboxes), nbors(3**ndim,nboxes)
      integer iflag(nboxes)
      double precision centers(ndim,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict,mc
      double precision distest,dis

      mc=2**ndim
      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,dis,k)
C$OMP$PRIVATE(ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,mc
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        ict = 0
                        do k=1,ndim
                           dis = centers(k,kbox) - centers(k,ibox) 
                           if(abs(dis).le.distest) ict = ict + 1
                        enddo
                        if(ict.eq.ndim) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c


      subroutine vol_tree_refine_boxes_flag(ndim,iflag,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer ndim,nd,npbox,nboxes
      real *8 fvals(nd,npbox,nboxes)
      integer nbloc,nbctr,nlctr
      real *8 centers(ndim,nboxes),bs,grid(ndim,npbox),xyz(ndim)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(2**ndim,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox,ilastbox
      integer, allocatable :: isum(:)
      real *8 dpars(*)
      integer ipars(*)
      complex *16 zpars(*)

      integer i,ibox,nel,j,l,jbox,nbl,ii,k,mc

      real *8 bsh
      integer isgn(ndim,2**ndim)

      external fun

      call get_child_box_sign(ndim,isgn)

      mc=2**ndim

      ilastbox = ifirstbox+nbloc-1


      bsh = bs/2

      allocate(isum(nbloc))

      call cumsum_nz(nbloc,iflag(ifirstbox),isum)

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j,jbox,nbl,l)
C$OMP$PRIVATE(xyz,k)
      do ibox = ifirstbox,ilastbox
        if(iflag(ibox).gt.0) then
          nchild(ibox) = mc
          nbl = nbctr + (isum(ibox-ifirstbox+1)-1)*mc
          do j=1,mc
            jbox = nbl+j
            do k=1,ndim
               centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo

            do l=1,npbox
              do k=1,ndim
                 xyz(k) = centers(k,jbox) + grid(k,l)*bs
              enddo
              call fun(nd,xyz,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,mc
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO
      
      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*mc


      return
      end
c
c
c
c
c
      subroutine cumsum(n,a,b)
c
c        this subroutine computes the cumulative sum of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer(n)
c              input array
c
c       output:
c         b - integer(n)
c            b(i) = sum_{j=1}^{i} a(j)
c
      implicit none
      integer a(n),b(n),n,i,isum
      isum = 0


      do i=1,n
        isum = isum + a(i)
        b(i) = isum
      enddo
      
      return
      end
c
c
c
c
c
      subroutine cumsum_nz(n,a,b)
c
c        this subroutine computes the cumulative sum of postive
c        entries of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer(n)
c              input array
c
c       output:
c         b - integer(n)
c            b(i) = sum_{j=1}^{i} I_{a(j)>0}
c
      implicit none
      integer a(n),b(n),n,i,isum

      isum = 0
      do i=1,n
        if(a(i).gt.0) isum = isum+1
        b(i) = isum
      enddo
      
      return
      end
c
c
c
c
c
c
      subroutine fgt_vol_tree_reorg(ndim,nboxes,nd,npbox,nblock,
     1    nboxid,nnewboxes,nboxes0,
     1    centers,nlevels,laddr,ilevel,iparent,nchild,ichild,
     2    fvals)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    nd             in: integer
c                   number of real value functions
c
c    npbox          in: integer
c                   number of grid points per function
c
c    nblock         in: integer(0:nlevels)
c                   number of new boxes on each level
c
c
c    centers        in/out: double precision(3,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
      implicit real *8 (a-h,o-z)
c     Calling sequence variables and temporary variables
      integer ndim,nboxes,nlevels,npbox,nd
      double precision centers(ndim,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer ilevel(nboxes)
      integer iparent(nboxes)
      integer nchild(nboxes)
      integer ichild(2**ndim,nboxes)
      integer nblock(0:nlevels),nboxid(nnewboxes)
      double precision fvals(nd,npbox,nboxes)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:)
      integer, allocatable :: iboxtocurbox(:)

      double precision, allocatable :: tfvals(:,:,:),tcenters(:,:)

c     Temporary variables
      integer i,j,k,l,ilevstart(0:nlevels+1)
      integer ibox,ilev, curbox,idim,nblev,mc

      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
      mc=2**ndim
      
      allocate(tilevel(nboxes),tiparent(nboxes),tnchild(nboxes))
      allocate(tichild(mc,nboxes),iboxtocurbox(nboxes))
      allocate(tfvals(nd,npbox,nboxes),tcenters(ndim,nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo
      call vol_tree_copy(ndim,nd,nboxes,npbox,centers,ilevel,
     1    iparent,nchild,ichild,fvals,tcenters,tilevel,
     2    tiparent,tnchild,tichild,tfvals)

     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
         enddo
      enddo

      curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,curbox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
c        now add new boxes to level ilev
         do i=1,nblock(ilev)
            ibox=nboxid(ilevstart(ilev)+i)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo
            do j=1,npbox
              do idim=1,nd
                fvals(idim,j,curbox) = tfvals(idim,j,ibox)
              enddo
            enddo
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1       iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         do i=1,mc
            if(tichild(i,ibox).eq.-1) ichild(i,iboxtocurbox(ibox)) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo

      return
      end
c
c
c
c
      subroutine vol_tree_coarsen(nd,ndim,eps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    fvals,fcoefs,ifdelete,idelboxid,ndelete)

c    This subroutine marks all boxes that can be coarsened
c    and computes the expansion coefficients for parent boxes
c    when its children are to be deleted in a subsequent call
c    to tree_reorg.
c
c    input:
c         nd - integer
c            number of right hand sides
c         ndim - integer
c            dimension of the underlying space
c         eps - double precision
c            precision requested
c         ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c         norder - integer
c           order of expansions for input function value array
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c         nboxes - integer
c            number of boxes
c       input and output
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
c     input/output:
c         fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c         fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c
c     output:
c         ifdelete - integer nboxes
c                    1 - box can be deleted
c                    0 - box should be kept
c         idelboxid - indices of the boxes to be removed
c
c         ndelete - total of number of removable boxes
c
c
c
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer nboxes,nlevels
      integer iptr(8),ltree,ifdelete(nboxes),idelboxid(nboxes)
      integer itree(ltree),norder,npbox
      real *8 fvals(nd,npbox,nboxes),fcoefs(nd,npbox,nboxes)

      real *8 centers(ndim,nboxes),boxsize(0:nlevels)

      real *8 fvals2(nd,npbox/2**ndim)
      real *8, allocatable :: polyv(:,:,:,:) 

      real *8 umat_nd(norder,norder,ndim)

      real *8 rmask(npbox)
      
      integer isgn(ndim,2**ndim)

      mc=2**ndim
      n2=norder/2
      
      call get_child_box_sign(ndim,isgn)
      allocate(polyv(norder,norder/2,ndim,mc))
      call get_c2p_interp_matrices(ndim,ipoly,norder,isgn,polyv)
      call get_val2coefs_matrices(ndim,ipoly,norder,umat_nd)

      do i=1,nboxes
         ifdelete(i)=0
      enddo
      
      iptype=0
      eta=1.0d0

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0

c     compute the L_{iptype} norm of the pot by going throught each leaf box
      call treedata_lpnorm(ndim,iptype,ipoly,nd,nlevels,itree,
     1    iptr,boxsize,norder,npbox,fvals,rnorm,nleaf)
      call prin2('l2 norm of the function=*',rnorm,1)
      
      rsc = rsc*rnorm

      ndelete=0
      
      itype=0
      do ilev = nlevels-1,0,-1
         sc = boxsize(ilev)/2
         rscale=sc**eta
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4) + ibox-1)
            if(nchild.eq.mc) then
               do j=1,mc
                  jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                  call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,jbox),
     1                n2,fvals2,polyv(1,1,1,j))
                  call tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,
     1                fvals(1,1,ibox),j)
               enddo
               call ortho_trans_nd(ndim,nd,itype,norder,fvals(1,1,ibox),
     1             fcoefs(1,1,ibox),umat_nd)
               call fun_err(nd,npbox,fcoefs(1,1,ibox),
     1            rmask,iptype,rscale,erra)
cccc               print *, 'erra=',erra, eps*rsc
               
               erra = erra/rsum
               if(erra.le.eps*rsc) then
                  do j=1,mc
                     ndelete=ndelete+1

                     jbox=itree(iptr(5)+mc*(ibox-1)+j-1)

                     idelboxid(ndelete)=jbox
                     
c                    parent
                     itree(iptr(3)+jbox-1)=-1
c                    nchild of jbox
                     itree(iptr(4)+jbox-1)=0
c                    jbox's children
                     do l=1,mc
                        itree(iptr(5)+mc*(jbox-1)+l-1) = -1
                     enddo
                     
                     ifdelete(jbox)=1
                  enddo
c                 nchild of ibox
                  itree(iptr(4)+ibox-1)=0
                  do j=1,mc
                     itree(iptr(5)+mc*(ibox-1)+j-1)=-1
                  enddo
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
      subroutine get_child_box_sign(ndim,isgn)
c     This subroutine computes the signs of all child boxes
c     The convention is as follows.
c      
c     1d - 1:m 2:p
c      
c     2d - 1:mm 
c          2:pm
c          3:mp
c          4:pp
c      
c     3d - 1:mmm
c          2:pmm
c          3:mpm
c          4:ppm
c          5:mmp
c          6:pmp
c          7:mpp
c          8:ppp
c
c     input:
c     ndim - dimension of the underlying space
c
c     output:
c     isgn - integer(ndim,2**ndim)
c            the signs of the center coordinates of each child box,
c            when the center of the parent box is at the origin
c
      implicit real *8 (a-h,o-z)
      integer isgn(ndim,2**ndim)

      mc = 2**ndim
      do j=1,ndim
         isgn(j,1)=-1
      enddo

      do j=1,ndim
         do i=1,mc,2**(j-1)
            if (i.gt.1) isgn(j,i)=-isgn(j,i-2**(j-1))
            do k=1,2**(j-1)-1
               isgn(j,i+k)=isgn(j,i)
            enddo
         enddo
      enddo

      return
      end
c     
c      
c      
c      
c
      subroutine tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,fvals,k)
c     This subroutine copies the function values on the kth child box
c     to the proper locations of arrays for the parent box
c      
c     input:
c     ndim - dimension of the underlying space
c     nd - number of vectors
c     norder - number of point along each dimension for the data on the 
c              parent box
c     k - child box index
c     fvals2 - function values on a subset of the tensor grid of the parent box
c              that are in the kth child box
c
c     output:
c     fvals - function values on the tensor product grid of the parent box
c
      implicit real *8 (a-h,o-z)
      real *8 fvals2(nd,(norder/2)**ndim)
      real *8 fvals(nd,norder**ndim)

      if (ndim.eq.1) then
         call tree_data_c2p_copy_1d(nd,norder,fvals2,fvals,k)
      elseif (ndim.eq.2) then
         call tree_data_c2p_copy_2d(nd,norder,fvals2,fvals,k)
      elseif (ndim.eq.3) then
         call tree_data_c2p_copy_3d(nd,norder,fvals2,fvals,k)
      endif

      return
      end
c     
c      
c      
c      
c
      subroutine tree_data_c2p_copy_1d(nd,norder,fvals2,fvals,k)
      implicit real *8 (a-h,o-z)
      real *8 fvals2(nd,norder/2)
      real *8 fvals(nd,norder)

      n=norder/2

      if (k.eq.1) then
         do i=1,n
         do ind=1,nd
            fvals(ind,i)=fvals2(ind,i)
         enddo
         enddo
      elseif (k.eq.2) then
         do i=1,n
         do ind=1,nd
            fvals(ind,i+n)=fvals2(ind,i)
         enddo
         enddo
      endif

      return
      end
c     
c      
c      
c      
c
      subroutine tree_data_c2p_copy_2d(nd,norder,fvals2,fvals,k)
      implicit real *8 (a-h,o-z)
      real *8 fvals2(nd,norder/2,norder/2)
      real *8 fvals(nd,norder,norder)

      n=norder/2

      if (k.eq.1) then
         do i=1,n
         do j=1,n
         do ind=1,nd
            fvals(ind,j,i)=fvals2(ind,j,i)
         enddo
         enddo
         enddo
      elseif (k.eq.2) then
         do i=1,n
         do j=1,n
         do ind=1,nd
            fvals(ind,j+n,i)=fvals2(ind,j,i)
         enddo
         enddo
         enddo
      elseif (k.eq.3) then
         do i=1,n
         do j=1,n
         do ind=1,nd
            fvals(ind,j,i+n)=fvals2(ind,j,i)
         enddo
         enddo
         enddo
      elseif (k.eq.4) then
         do i=1,n
         do j=1,n
         do ind=1,nd
            fvals(ind,j+n,i+n)=fvals2(ind,j,i)
         enddo
         enddo
         enddo
      endif

      return
      end
c     
c      
c      
c      
c      
c      
c
      subroutine tree_data_c2p_copy_3d(nd,norder,fvals2,fvals,k)
      implicit real *8 (a-h,o-z)
      real *8 fvals2(nd,norder/2,norder/2,norder/2)
      real *8 fvals(nd,norder,norder,norder)

      n=norder/2

      if (k.eq.1) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.2) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.3) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j+n,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.4) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j+n,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.5) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.6) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.7) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j+n,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (k.eq.8) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j+n,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end
c     
c      
c      
c      
c
      subroutine vol_tree_p2c_interp(nd,ndim,ipoly,norder,
     1    npbox,itree,iptr,centers,boxsize,isgn,polyv,umat_nd,
     2    fvals,fcoefs,iparentbox,listcid)

c     This subroutine computes the data on each child box using the data on the parent box
c     by interpolation and update parent,child,level,centers info of the tree.
c
c     input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c
c     centers - double precision (ndim,nboxes)
c           xyz coordintes of boxes in the tree structure
c     boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     isgn - integer (ndim,2**ndim)
c           contains the signs of center coordinates of child boxes
c           when the parent box is centered at the origin
c     polyv - precomputed matrices for carrying out interpolation from the parent
c                 box to its child boxes
c     umat_nd - precomputed matrices for converting function values on a tensor
c                 grid to its polynomial expansion coefficients
c
c     iparentbox - the index of the parent box to be refined
c     listcid - the indices of the child boxes
c
c     input and output:
c     itree - integer(ltree)
c            array containing the tree structure
c
c     input/output:
c         fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c         fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer iptr(8)
      integer itree(*),norder,npbox
      integer listcid(2**ndim)
      
      real *8 fvals(nd,npbox,nboxes),fcoefs(nd,npbox,nboxes)

      real *8 centers(ndim,nboxes),boxsize(0:nlevels)

      integer isgn(ndim,2**ndim)
      real *8 polyv(norder,norder/2,ndim,2**ndim)
      real *8 umat_nd(norder,norder,ndim)
      
      mc=2**ndim

      itype=0
      ilev = itree(iptr(2)+iparentbox-1)
c     nchild of ibox
      itree(iptr(4)+iparentbox-1)=mc
      jlev = ilev+1
      bsh = boxsize(jlev)/2
      
      do j=1,mc
         jbox=listcid(j)
         itree(iptr(2)+jbox-1) = jlev
         do k=1,ndim
            centers(k,jbox) = centers(k,iparentbox)+isgn(k,i)*bsh
         enddo
         
         call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,iparentbox),
     1       norder,fvals(1,1,jbox),polyv(1,1,1,j))

         call ortho_trans_nd(ndim,nd,itype,norder,fvals(1,1,jbox),
     1       fcoefs(1,1,jbox),umat_nd)
                     
c        parent
         itree(iptr(3)+jbox-1)=ibox
c        nchild of jbox
         itree(iptr(4)+jbox-1)=0
c        jbox's children
         do l=1,mc
            itree(iptr(5)+mc*(jbox-1)+l-1) = -1
         enddo
c        ichild of ibox
         itree(iptr(5)+mc*(iparentbox-1)+j-1)=jbox
      enddo
      
      return
      end
c
c
c
c
      subroutine get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyv)
c     This subroutine returns to the user a list of matrices needed to
c     carry out interpolation from the parent box to its 2**ndim children.
c
c     input:
c     ndim - integer
c            dimension of the underlying space
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     isgn - integer (ndim,2**ndim)
c           contains the signs of center coordinates of child boxes
c           when the parent box is centered at the origin
c
c     output:
c     polyv - double precision (norder,norder,ndim,2**ndim)
c             polyv(1,1,1,j) contains ndim matrices where each of them contains
c             polynomial values in a tensor grid of the jth child box along
c             each dimension
c
      implicit real *8 (a-h,o-z)
      real *8 polyv(norder,norder,ndim,2**ndim) 

      real *8 xs(norder)
      real *8 xp(norder),xm(norder)
      real *8 polyvp(norder,norder),polyvm(norder,norder)
      
      integer isgn(ndim,2**ndim)

      mc=2**ndim
      norder2=norder*norder
      
      itype=0
      if (ipoly.eq.0) then
         call legeexps(itype,norder,xs,u,v,ws)
      elseif (ipoly.eq.1) then
         call chebexps(itype,norder,xs,u,v,ws)
      endif

      do i=1,norder
         xm(i) = xs(i)/2-0.5d0
         xp(i) = xs(i)/2+0.5d0
      enddo

      do i=1,norder
         if (ipoly.eq.0) then
            call legepols(xp(i),norder-1,polyvp(1,i))
            call legepols(xm(i),norder-1,polyvm(1,i))
         elseif (ipoly.eq.1) then
            call chebpols(xp(i),norder-1,polyvp(1,i))
            call chebpols(xm(i),norder-1,polyvm(1,i))
         endif
      enddo

      do j=1,mc
         do k=1,ndim
            if (isgn(k,j).lt.0) then
               call dcopy_f77(norder2,polyvm,1,polyv(1,1,k,j),1)
            else
               call dcopy_f77(norder2,polyvp,1,polyv(1,1,k,j),1)
            endif
         enddo
      enddo

      return
      end
c      
c            
c      
      subroutine get_c2p_interp_matrices(ndim,ipoly,norder,isgn,polyv)
c     This subroutine returns to the user a list of matrices needed to
c     carry out interpolation from child boxes to the parent box.
c
c     input:
c     ndim - integer
c            dimension of the underlying space
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     isgn - integer (ndim,2**ndim)
c           contains the signs of center coordinates of child boxes
c           when the parent box is centered at the origin
c
c     output:
c     polyv - double precision (norder,norder/2,ndim,2**ndim)
c             the parent box had a norder**ndim tensor grid. When norder is even,
c             this tensor grid has (norder/2)**ndim tensor grid points in each of
c             its child box.
c             polyv(1,1,1,j) contains ndim matrices where each of them contains
c             polynomial values in the small tensor grid of the jth child box along
c             each dimension
c
c
      implicit real *8 (a-h,o-z)
      real *8 polyv(norder,norder/2,ndim,2**ndim) 

      real *8 xs(norder)
      real *8 xp(norder/2),xm(norder/2)
      real *8 polyvp(norder,norder/2),polyvm(norder,norder/2)
      
      integer isgn(ndim,2**ndim)

      mc=2**ndim
      n2=norder/2
      
      itype=0
      if (ipoly.eq.0) then
         call legeexps(itype,norder,xs,u,v,ws)
      elseif (ipoly.eq.1) then
         call chebexps(itype,norder,xs,u,v,ws)
      endif

      do i=1,n2
         xm(i) = xs(i)*2+1.0d0
         xp(i) = xs(n2+i)*2-1.0d0
      enddo

      do i=1,n2
         if (ipoly.eq.0) then
            call legepols(xp(i),norder-1,polyvp(1,i))
            call legepols(xm(i),norder-1,polyvm(1,i))
         elseif (ipoly.eq.1) then
            call chebpols(xp(i),norder-1,polyvp(1,i))
            call chebpols(xm(i),norder-1,polyvm(1,i))
         endif
      enddo

      do j=1,mc
         do k=1,ndim
            if (isgn(k,j).lt.0) then
               call dcopy_f77(norder*n2,polyvm,1,polyv(1,1,k,j),1)
            else
               call dcopy_f77(norder*n2,polyvp,1,polyv(1,1,k,j),1)
            endif
         enddo
      enddo

      return
      end
c      
c            
c
      subroutine get_val2coefs_matrices(ndim,ipoly,norder,umat_nd)
c     This subroutine returns to the user transformation matrices
c     converting function values to polynomial expansion coefficients.
c
c     input:
c     ndim - integer
c            dimension of the underlying space
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c
c     output:
c     umat_nd - double precision (norder,norder,ndim)
c             umat_nd(1,1,k) is the transformation matrix converting function
c             values on a tensor grid to its polynomial expansion coefficients
c             along the kth dimension.
c
c
      implicit real *8 (a-h,o-z)
      real *8 umat_nd(norder,norder,ndim)

      real *8 xs(norder),ws(norder)
      real *8 umat(norder,norder),vmat(norder,norder)
      
      itype=2
      if (ipoly.eq.0) then
         call legeexps(itype,norder,xs,umat,vmat,ws)
      elseif (ipoly.eq.1) then
         call chebexps(itype,norder,xs,umat,vmat,ws)
      endif

      n2=norder*norder
      do k=1,ndim
         call dcopy_f77(n2,umat,1,umat_nd(1,1,k),1)
      enddo

      return
      end
c      
c            
c            
