c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c     this is the end of the debugging code and the beginning
c     of the legendre expansion routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c     This file contains a set of subroutines for the handling
c     of a function defined on an adaptive tree. 
c     Following is a brief description of these
c     subroutines.
c
      subroutine treedata_trans2d(nd,itype,nlevels,itree,iptr,boxsize,
     1    norder,fin,fout,uxmat,uymat)
c
c     This code converts an input tree data given on a tensor product grid on each leaf box
c     to the output tree data on each leaf box. 
c     Depending on the 1d transformation matrix, it could be:  
c     (1) converting function values to the coefficients of orthogonal polynomial expansion
c     coefficients;
c     (2) if umax converts coefficients to function values, then it's the inverse map of (1).
c     (3) converting fcoefs to the coefficients of its derivatives.
c 
c     It's user's responsibility to ensure that uxmat,uymat are the correct 1D
c     transformation matrices.
c     
c 
c 
c       input
c
c         nd - integer,   number of functions
c         itype - itype = 1 - transformation is only along x-direction, for example, calculating
c                         the expansion coefficients of the x-derivative
c                 itype = 2 - transformation is only along y-direction, for example, calculating
c                         the expansion coefficients of the y-derivative
c                 itype = 0 - transformation is done on both directions, for example, 
c                         val2coefs or coef2vals
c         nlevels - integer
c            number of levels
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
c         fin - double (nd,norder**2,nboxes)
c            input data on the tree
c         uxmat - 1D transformation matrix along the x-direction
c         uymat - 1D transformation matrix along the y-direction
c     output:
c         fout - double precision (nd,norder**2,nboxes)
c            putput data given on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fin(nd,norder*norder,*)
      real *8 fout(nd,norder*norder,*)

      real *8 uxmat(norder,norder)
      real *8 uymat(norder,norder)


      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_trans2d(nd,itype,norder,fin(1,1,ibox),
     1           fout(1,1,ibox),uxmat,uymat)
             if (itype.gt.0) then
                do j=1,norder*norder
                   do ind=1,nd
                      fout(ind,j,ibox)=fout(ind,j,ibox)*sc
                   enddo
                enddo
             endif
          endif
        enddo
      enddo


      end
c
c
c
c
      subroutine treedata_evalg2d(nd,ipoly,nlevels,itree,iptr,boxsize,
     1    norder,coefs,grad)
c
c     This code evaluates the gradient at the tensor grid on an adaptive tree
c     given the orthogonal polynomial expansion coefficients at each leaf box.
c 
c       input
c
c         nd - integer,   number of functions
c         nlevels - integer
c            number of levels
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
c         coefs - double (nd,norder**2,nboxes)
c           expansion coefficients on quad tree
c
c     output:
c         grad - double precision (nd,2,norder**2,nboxes)
c            gradient values on tensor grid on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefs(nd,norder*norder,*)
      real *8 grad(nd,2,norder*norder,*)

      real *8, allocatable :: vmat(:,:),vpmat(:,:),vppmat(:,:)
      real *8, allocatable :: xs(:),ws(:),umat(:,:)

      allocate(xs(norder),ws(norder),umat(norder,norder))
      allocate(vmat(norder,norder))
      allocate(vpmat(norder,norder))
      allocate(vppmat(norder,norder))

      itype=2
      if (ipoly.eq.0) then
         call legeexps2(itype,norder,xs,umat,vmat,ws,vpmat,vppmat)
      elseif (ipoly.eq.1) then
         call chebexps2(itype,norder,xs,umat,vmat,ws,vpmat,vppmat)
      endif

      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_evalg2d(nd,norder,coefs(1,1,ibox),sc,
     1           grad(1,1,1,ibox),vmat,vpmat)             
          endif
        enddo
      enddo


      end
c
c
c
c
      subroutine treedata_evalh2d(nd,ipoly,nlevels,itree,iptr,boxsize,
     1    norder,coefs,hess)
c
c     This code evaluates and hessian at the tensor 
c     grid on an adaptive tree given the orthogonal polynomial expansion coefficients 
c     at each leaf box.
c
c       input
c
c         nd - integer,   number of functions
c         ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c         nlevels - integer
c            number of levels
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
c         coefs - double (nd,norder**2,nboxes)
c           expansion coefficients on quad tree
c
c     output:
c         hess - double precision (nd,3,norder**2,nboxes)
c            hessian values on tensor grid on each leaf box
c            in the order of uxx, uxy, uyy
c      
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefs(nd,norder*norder,*)
      real *8 hess(nd,3,norder*norder,*)

      real *8, allocatable :: vmat(:,:),vpmat(:,:),vppmat(:,:)
      real *8, allocatable :: xs(:),ws(:),umat(:,:)

      allocate(xs(norder),ws(norder),umat(norder,norder))
      allocate(vmat(norder,norder))
      allocate(vpmat(norder,norder))
      allocate(vppmat(norder,norder))

      itype=2
      if (ipoly.eq.0) then
         call legeexps2(itype,norder,xs,umat,vmat,ws,vpmat,vppmat)
      elseif (ipoly.eq.1) then
         call chebexps2(itype,norder,xs,umat,vmat,ws,vpmat,vppmat)
      endif
      
      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_evalh2d(nd,norder,coefs(1,1,ibox),sc,
     1           hess(1,1,1,ibox),
     2           vmat,vpmat,vppmat)             
          endif
        enddo
      enddo


      end
c
c
c
c
      subroutine treedata_coefs_p_to_g2d(nd,nlevels,itree,iptr,
     1    boxsize,norder,coefsp,coefsg,umat)
c
c     This code converts expansion coefficients of a tree data given on a tensor product
c     grid on each leaf box to the expansion coefficients of its gradient
c
c       input
c
c         nd - integer,   number of functions
c         nlevels - integer
c            number of levels
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
c         coefsp - double (nd,norder**2,nboxes)
c            expansion coefficients of the potential
c     output:
c         coefsg - double precision (nd,2,norder**2,nboxes)
c            expansion coefficients of the gradient
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefsp(nd,norder*norder,*)
      real *8 coefsg(nd,2,norder*norder,*)

      real *8 umat(norder,norder)

      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_coefsg2d(nd,norder,coefsp(1,1,ibox),
     1           coefsg(1,1,1,ibox),umat)
             do j=1,norder*norder
                do i=1,2
                   do ind=1,nd
                      coefsg(ind,i,j,ibox)=coefsg(ind,i,j,ibox)*sc
                   enddo
                enddo
             enddo
          endif
        enddo
      enddo


      end
c
c
c
c

      subroutine treedata_coefs_p_to_gh2d(nd,nlevels,itree,iptr,
     1    boxsize,norder,coefsp,coefsg,coefsh,umat)
c
c     This code converts expansion coefficients of a tree data given on a tensor product
c     grid on each leaf box to the expansion coefficients of its gradient
c
c       input
c
c         nd - integer,   number of functions
c         nlevels - integer
c            number of levels
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
c         coefsp - double (nd,norder**2,nboxes)
c            expansion coefficients of the potential
c     output:
c         coefsg - double precision (nd,2,norder**2,nboxes)
c            expansion coefficients of the gradient
c         coefsh - double precision (nd,3,norder**2,nboxes)
c            expansion coefficients of the hessian
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefsp(nd,norder*norder,*)
      real *8 coefsg(nd,2,norder*norder,*)
      real *8 coefsh(nd,3,norder*norder,*)

      real *8 umat(norder,norder)

      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
        sc2 = sc*sc 
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_coefsgh2d(nd,norder,coefsp(1,1,ibox),
     1           coefsg(1,1,1,ibox),coefsh(1,1,1,ibox),umat)
             do j=1,norder*norder
                do ind=1,nd
                   coefsg(ind,1,j,ibox)=coefsg(ind,1,j,ibox)*sc
                   coefsg(ind,2,j,ibox)=coefsg(ind,2,j,ibox)*sc

                   coefsh(ind,1,j,ibox)=coefsh(ind,1,j,ibox)*sc2
                   coefsh(ind,2,j,ibox)=coefsh(ind,2,j,ibox)*sc2
                   coefsh(ind,3,j,ibox)=coefsh(ind,3,j,ibox)*sc2
                enddo
             enddo
          endif
        enddo
      enddo


      end
c
c
c
c
      subroutine treedata_evalt2d(nd,ipoly,itree,ltree,nboxes,nlevels,
     1    iptr,tcenters,boxsize,norder,fcoefs,nt,targ,fvals)
c
c       This code evaluates function values at nt target points where the function is given 
c       as the orthogonal polynomial expansion coefficients on each leaf box
c 
c       input
c
c         nd - integer,   number of functions
c         ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c         itree - integer(ltree)
c            array containing the tree structure
c         ltree - integer
c            length of array containing the tree structure
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
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
c         tcenters - double (2,nboxes) coordinates of box centers
c         boxsize - double(nboxes) box size of each box
c      
c      
c         norder - integer
c           order of expansions for input coefficients array
c         fcoefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients on each leaf box
c         nt - number of target points
c         targ - (2,nt) coordinates of targets
c      
c      
c     output:
c         fvals - double precision (nd,nt) function values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 tcenters(2,nboxes),boxsize(0:nlevels)
      
      real *8 fcoefs(nd,norder*norder,nboxes)
      real *8 targ(2,nt),fvals(nd,nt)
      real *8 xy(2)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      ndim=2
      
      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(2,nt),potsort(nd,nt))
c      
c     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)
c
c     
c     
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1

               if (npts.gt.0) then
                  cx = tcenters(1,ibox)
                  cy = tcenters(2,ibox)

                  do i=istart,iend
                     xy(1) = (targsort(1,i)-cx)*sc
                     xy(2) = (targsort(2,i)-cy)*sc
                     call orth_evalt2d(nd,ipoly,norder,fcoefs(1,1,ibox),
     1                   xy,potsort(1,i))
                  enddo
               endif
            endif
        enddo
      enddo
c
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,fvals,itarg)

      end
c
c
c
c
      subroutine treedata_evalpgt2d(nd,ipoly,itree,ltree,nboxes,nlevels,
     1    iptr,tcenters,boxsize,norder,coefs,
     2    nt,targ,pot,grad)
c
c     This code evaluates potential and gradient at nt target points
c     where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c       input
c
c         nd - integer,   number of functions
c         ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c         itree - integer(ltree)
c            array containing the tree structure
c         ltree - integer
c            length of array containing the tree structure
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
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
c         tcenters - double (2,nboxes) coordinates of box centers
c         boxsize - double(nboxes) box size of each box
c      
c      
c         norder - integer
c           order of expansions for input coefficients array
c         coefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients of the potential on each leaf box
c         nt - number of target points
c         targ - (2,nt) coordinates of targets
c      
c      
c     output:
c         pot - double precision (nd,nt) potential values
c         grad - double precision (nd,2,nt) gradient values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      
      integer itree(ltree)
      real *8 tcenters(2,nboxes),boxsize(0:nlevels)
      
      real *8 coefs(nd,norder*norder,nboxes)
      real *8 targ(2,nt),pot(nd,nt),grad(nd,2,nt)
      real *8 xy(2)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)
      real *8, allocatable :: gradsort(:,:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      ndim=2
      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(2,nt),potsort(nd,nt),gradsort(nd,2,nt))
c      
c     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)
c
c     
c     
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1

               if (npts.gt.0) then
                  cx = tcenters(1,ibox)
                  cy = tcenters(2,ibox)

                  do i=istart,iend
                     xy(1) = (targsort(1,i)-cx)*sc
                     xy(2) = (targsort(2,i)-cy)*sc
                     call orth_evalpgt2d(nd,ipoly,norder,
     1                   coefs(1,1,ibox),sc,xy,
     2                   potsort(1,i),gradsort(1,1,i))
                  enddo
               endif
            endif
        enddo
      enddo
c
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,pot,itarg)
      call dreorderi(nd*2,nt,gradsort,grad,itarg)

      end
c
c
c
c
      subroutine treedata_evalpght2d(nd,ipoly,itree,ltree,nboxes,
     1    nlevels,iptr,tcenters,boxsize,norder,coefs,
     2    nt,targ,pot,grad,hess)
c
c     This code evaluates potential, gradient, and hessian at nt target points
c     where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c       input
c
c         nd - integer,   number of functions
c         ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c         itree - integer(ltree)
c            array containing the tree structure
c         ltree - integer
c            length of array containing the tree structure
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
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
c         tcenters - double (2,nboxes) coordinates of box centers
c         boxsize - double(nboxes) box size of each box
c      
c      
c         norder - integer
c           order of expansions for input coefficients array
c         coefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients of the potential on each leaf box
c         nt - number of target points
c         targ - (2,nt) coordinates of targets
c      
c      
c     output:
c         pot - double precision (nd,nt) potential values
c         grad - double precision (nd,2,nt) gradient values
c         hess - double precision (nd,3,nt) gradient values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 tcenters(2,nboxes),boxsize(0:nlevels)
      
      real *8 coefs(nd,norder*norder,nboxes)
      real *8 targ(2,nt),pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)
      real *8 xy(2)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)
      real *8, allocatable :: gradsort(:,:,:)
      real *8, allocatable :: hesssort(:,:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      ndim=2
      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(2,nt),potsort(nd,nt),gradsort(nd,2,nt))
      allocate(hesssort(nd,3,nt))
c      
c     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)
c
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1

               if (npts.gt.0) then
                  cx = tcenters(1,ibox)
                  cy = tcenters(2,ibox)

                  do i=istart,iend
                     xy(1) = (targsort(1,i)-cx)*sc
                     xy(2) = (targsort(2,i)-cy)*sc
                     call orth_evalpght2d(nd,ipoly,norder,
     1                   coefs(1,1,ibox),sc,xy,
     2                   potsort(1,i),gradsort(1,1,i),hesssort(1,1,i))
                  enddo
               endif
            endif
        enddo
      enddo
c     
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,pot,itarg)
      call dreorderi(nd*2,nt,gradsort,grad,itarg)
      call dreorderi(nd*3,nt,hesssort,hess,itarg)

      end
c
c
c
c
      subroutine treedata_derror(nd,nlevels,itree,iptr,
     1    npbox,fex,fcomp,abserr,rnorm,nleaf)
c
c     computes the absolute l2 error of two tree data fcomp
c     fex, the l2 norm of fex, and the total number of leaf boxes
c 
c       input
c
c         nd - integer,   number of functions
c         nlevels - integer
c            number of levels
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
c         npbox - integer
c           number of points in each leaf box
c         fex - double (nd,npbox,nboxes)
c           values of the reference tree data on each leaf box
c         fcomp - double (nd,npbox,nboxes)
c           values of the numerical tree data on each leaf box
c     output:
c         abserr - absolute l2 error
c         rnorm - l2 norm of fex
c         nleaf  - total number of leaf boxes
c      
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 fex(nd,npbox,*)
      real *8 fcomp(nd,npbox,*)

      abserr=0
      rnorm=0
      nleaf=0

      do ilev = 0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             nleaf=nleaf+1
             do j=1,npbox
                do ind=1,nd
                   rnorm=rnorm+fex(ind,j,ibox)**2
                   abserr=abserr+(fex(ind,j,ibox)-fcomp(ind,j,ibox))**2
                enddo
             enddo
c             call prin2('fex=*',fex(1,1,ibox),nd*npbox)
c             call prin2('fcomp=*',fcomp(1,1,ibox),nd*npbox)
          endif
        enddo
      enddo

      abserr=sqrt(abserr)
      rnorm=sqrt(rnorm)
      
      end
c
c
c
c
      
