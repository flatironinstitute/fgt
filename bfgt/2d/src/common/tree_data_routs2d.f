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
     1    norder,fin,fout,umat)
c
c     This code converts an input tree data given on a tensor product grid on each leaf box
c     to the output tree data on each leaf box. 
c     Depending on the 1d transformation matrix, it could be:  
c     (1) converting function values to the coefficients of orthogonal polynomial expansion
c     coefficients;
c     (2) if umax converts coefficients to function values, then it's the inverse map of (1).
c     (3) converting fcoefs to the coefficients of its derivatives.
c 
c     It's user's responsibility to ensure that umat is the correct 1D transformation matrix.
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

      real *8 umat(norder,norder)
      real *8, allocatable :: work(:,:,:)

      allocate(work(nd,norder,norder))

      do ilev = 0,nlevels
        sc = 2/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_trans2d(nd,itype,norder,fin(1,1,ibox),
     1           fout(1,1,ibox),umat,work)
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
      subroutine treedata_coefs_p_to_g2d(nd,nlevels,itree,iptr,
     1    boxsize,norder,coefsp,coefsg,umat)
c
c     This code converts expansion coefficients of a tree data given on a tensor product
c     grid on each leaf box to the expansion coefficients of its gradient
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
        sc = 2/boxsize(ilev)
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
        sc = 2/boxsize(ilev)
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
      subroutine treedata_eval2d(nd,itype,itree,ltree,nboxes,nlevels,
     1    iptr,tcenters,boxsize,norder,fcoefs,nt,targ,fvals)
c
c       This code evaluates function values at nt target points where the function is given 
c       as the orthogonal polynomial expansion coefficients on each leaf box
c 
c       input
c
c         nd - integer,   number of functions
c         itype - 1: Legendre polynomials
c                 2: Chebyshev polynomials
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
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8)
      real *8 tcenters(2,nboxes),boxsize(0:nlevels)
      
      real *8 fcoefs(nd,norder*norder,nboxes)
      real *8 targ(2,nt),fvals(nd,nt)
      real *8 xy(2)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
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
         bs = boxsize(ilev)/2
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
                     xy(1) = (targsort(1,i)-cx)/bs
                     xy(2) = (targsort(2,i)-cy)/bs
                     call orth_eval2d(nd,itype,norder,fcoefs(1,1,ibox),
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
      subroutine treedata_evalpg2d(nd,itype,itree,ltree,nboxes,nlevels,
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
c         itype - 1: Legendre polynomials
c                 2: Chebyshev polynomials
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
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8)
      real *8 tcenters(2,nboxes),boxsize(0:nlevels)
      
      real *8 coefs(nd,norder*norder,nboxes)
      real *8 targ(2,nt),pot(nd,nt),grad(nd,2,nt)
      real *8 xy(2)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)
      real *8, allocatable :: gradsort(:,:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
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
         bs = boxsize(ilev)/2
         sc = 2/boxsize(ilev)
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
                     xy(1) = (targsort(1,i)-cx)/bs
                     xy(2) = (targsort(2,i)-cy)/bs
                     call orth_evalpg2d(nd,itype,norder,
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
      subroutine treedata_evalpgh2d(nd,itype,itree,ltree,nboxes,nlevels,
     1    iptr,tcenters,boxsize,norder,coefs,
     2    nt,targ,pot,grad,hess)
c
c     This code evaluates potential, gradient, and hessian at nt target points
c     where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c       input
c
c         nd - integer,   number of functions
c         itype - 1: Legendre polynomials
c                 2: Chebyshev polynomials
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
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8)
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
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(2,nt),potsort(nd,nt),gradsort(nd,2,nt))
      allocate(hesssort(nd,3,nt))
c      
c     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)
c
c     
c
      call cpu_time(t1)
      do ilev = 0,nlevels
         bs = boxsize(ilev)/2
         sc = 1/bs
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
                     xy(1) = (targsort(1,i)-cx)/bs
                     xy(2) = (targsort(2,i)-cy)/bs
                     call orth_evalpgh2d(nd,itype,norder,
     1                   coefs(1,1,ibox),sc,xy,
     2                   potsort(1,i),gradsort(1,1,i),hesssort(1,1,i))
                  enddo
               endif
            endif
        enddo
      enddo
      call cpu_time(t2)
      call prin2('time on extra targ pot eval=*',t2-t1,1)
      call prin2('speed in pps=*',(nt+0.0d0)/(t2-t1),1)

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
      subroutine orth_trans2d(nd,itype,norder,fin,fout,umat,work)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases. The polynomial could be any orthogonal
c     polyomial.
c
c     1. if umat converts function values to orthogonal polynomial expansion coefficients,  
c        then it converts function values on the tensor product grid to orthogonal
c        polynomial expansion coefficients in 2d.
c     2. if umax converts coefficients to function values, then it's the inverse map of (1).
c     3. if umax is the differentiation matrix in the coefficient space, then it converts
c        fcoefs to the coefficients of its derivative. 
c
c     input:
c     nd - number of input data
c     itype - itype = 1 - transformation is only along x-direction, for example, calculating
c                         the expansion coefficients of the x-derivative
c             itype = 2 - transformation is only along y-direction, for example, calculating
c                         the expansion coefficients of the y-derivative
c             itype = 0 - transformation is done on both directions, for example, 
c                         val2coefs or coef2vals
c      
c     norder - order of polynomial expansion
c     fin - input data
c     umat - 1D transformation matrix
c     work - work array
c      
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,norder,norder)
      real *8 fout(nd,norder,norder),umat(norder,norder)
      real *8 work(nd,norder,norder)

      if (itype.eq.0) then
c        transform in x
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fin(ind,k1,j)
               enddo
               work(ind,k,j)=dd
            enddo
         enddo
         enddo
c        transform in y
         do i=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do i1=1,norder
                  dd=dd+umat(i,i1)*work(ind,k,i1)
               enddo
               fout(ind,k,i)=dd
            enddo
         enddo
         enddo

      elseif (itype.eq. 1) then
c        transform in x
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fin(ind,k1,j)
               enddo
               fout(ind,k,j)=dd
            enddo
         enddo
         enddo
      
      elseif (itype.eq. 2) then
c        transform in y
         do k=1,norder
         do j=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fin(ind,j,k1)
               enddo
               fout(ind,j,k)=dd
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
      subroutine orth_coefsg2d(nd,norder,coefsp,coefsg,umat)
c
c     this subroutine computes the expansion coefficients for the gradient
c     given the expansion coefficients of the potential
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefsp - expansion coefficients of the potential
c     umat - 1D transformation matrix
c     work - work array
c      
c     output:
c     coefsg - expansion coefficients of the gradient
c
      implicit real *8 (a-h,o-z)
      real *8 coefsp(nd,norder,norder)
      real *8 coefsg(nd,2,norder,norder),umat(norder,norder)

c     transform in x
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1,j)
            enddo
            coefsg(ind,1,k,j)=dd
         enddo
      enddo
      enddo
      
c     transform in y
      do k=1,norder
      do j=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,j,k1)
            enddo
            coefsg(ind,2,j,k)=dd
         enddo
      enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine orth_coefsgh2d(nd,norder,coefsp,coefsg,coefsh,umat)
c
c     this subroutine computes the expansion coefficients for the gradient
c     and the hessian given the expansion coefficients of the potential
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefsp - expansion coefficients of the potential
c     umat - 1D transformation matrix
c     work - work array
c      
c     output:
c     coefsg - expansion coefficients of the gradient
c     coefsh - expansion coefficients of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefsp(nd,norder,norder)
      real *8 coefsg(nd,2,norder,norder),umat(norder,norder)
      real *8 coefsh(nd,3,norder,norder)

c     compute the expansion coefficients of the gradient
c     differentiation in x
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1,j)
            enddo
            coefsg(ind,1,k,j)=dd
         enddo
      enddo
      enddo
      
c     differentiation in y
      do k=1,norder
      do j=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,j,k1)
            enddo
            coefsg(ind,2,j,k)=dd
         enddo
      enddo
      enddo

c     compute the expansion coefficients of the hessian
c     u_{xx}
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsg(ind,1,k1,j)
            enddo
            coefsh(ind,1,k,j)=dd
         enddo
      enddo
      enddo

c     u_{yx}
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsg(ind,2,k1,j)
            enddo
            coefsh(ind,2,k,j)=dd
         enddo
      enddo
      enddo
      
c     u_{yy}
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*coefsg(ind,2,k,j1)
            enddo
            coefsh(ind,3,k,j)=dd
         enddo
      enddo
      enddo

      
      return
      end
c
c
c
c
      subroutine orth_eval2d(nd,itype,norder,fcoefs,targ,fval)
c
c     this subroutine evaluates the orthogonal polynomial expansion at the given target point
c
c     input:
c     nd - number of input data
c     itype - itype=1 Legendre polynomials
c                   2 Chebyshev polynomials
c     norder - order of polynomial expansion
c     fcoefs - orthogonal polynomial expansion coefficients
c     targ - xy coordinates of the target point
c      
c     output:
c     fval - value of the function
c
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,norder,norder)
      real *8 targ(2),fval(nd)

      real *8 px(norder),py(norder),tmp(nd,norder)

      x=targ(1)
      y=targ(2)
      
      if (itype .eq. 1) then
         call legepols(x,norder-1,px)
         call legepols(y,norder-1,py)
      elseif (itype .eq. 2) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
      endif
c     
      do j=1,norder
         do ind=1,nd
            dd=0
            do k=1,norder
               dd=dd+fcoefs(ind,k,j)*px(k)
            enddo
            tmp(ind,j)=dd
         enddo
      enddo
c
      do ind=1,nd
         dd=0
         do j=1,norder
            dd=dd+tmp(ind,j)*py(j)
         enddo
         fval(ind)=dd
      enddo

      return
      end
c
c
c
c
      subroutine orth_evalpg2d(nd,itype,norder,coefs,sc,
     1    targ,pot,grad)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     itype - itype=1 Legendre polynomials
c                   2 Chebyshev polynomials
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     targ - xy coordinates of the target point
c      
c     output:
c     pot - value of the potential
c     grad - value of the potential
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder)
      real *8 targ(2),pot(nd),grad(nd,2)

      real *8 px(norder),py(norder)
      real *8 pxp(norder),pyp(norder)
      
      real *8, allocatable :: tmp(:,:,:)

      allocate(tmp(nd,2,norder))
      
      x=targ(1)
      y=targ(2)
      
      if (itype .eq. 1) then
         call legepolders(x,px,pxp,norder-1)
         call legepolders(y,py,pyp,norder-1)
      elseif (itype .eq. 2) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
      endif
c     
      do j=1,norder
         do ind=1,nd
            pp=0
            d1=0
            do k=1,norder
               pp=pp+coefs(ind,k,j)*px(k)
               d1=d1+coefs(ind,k,j)*pxp(k)
            enddo
            tmp(ind,1,j)=pp
            tmp(ind,2,j)=d1
         enddo
      enddo
c
      do ind=1,nd
         pp=0
         d1=0
         d2=0
         do j=1,norder
            pp=pp+tmp(ind,1,j)*py(j)
            d1=d1+tmp(ind,2,j)*py(j)
            d2=d2+tmp(ind,1,j)*pyp(j)
         enddo
         pot(ind)=pp
         grad(ind,1)=d1*sc
         grad(ind,2)=d2*sc
      enddo

      return
      end
c
c
c
c
      subroutine orth_evalpgh2d(nd,itype,norder,coefs,sc,
     1    targ,pot,grad,hess)
c
c     this subroutine evaluates the potential, gradient and hessian at the given target point
c
c     input:
c     nd - number of input data
c     itype - itype=1 Legendre polynomials
c                   2 Chebyshev polynomials
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     targ - xy coordinates of the target point
c      
c     output:
c     pot - value of the potential
c     grad - value of the gradient
c     hess - value of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder)
      real *8 targ(2),pot(nd),grad(nd,2),hess(nd,3)

      real *8 px(norder),py(norder)
      real *8 pxp(norder),pyp(norder)
      real *8 pxpp(norder),pypp(norder)
      
      real *8, allocatable :: tmp(:,:,:)

      allocate(tmp(nd,3,norder))

      sc2=sc*sc
      
      x=targ(1)
      y=targ(2)
      
      if (itype .eq. 1) then
         call legepolders2(x,px,pxp,pxpp,norder-1)
         call legepolders2(y,py,pyp,pypp,norder-1)
         do i=1,norder
            pxp(i)=pxp(i)*sc
            pyp(i)=pyp(i)*sc
            pxpp(i)=pxpp(i)*sc2
            pypp(i)=pypp(i)*sc2
         enddo
      elseif (itype .eq. 2) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
      endif
c     
      do j=1,norder
         do ind=1,nd
            pp=0
            d1=0
            d2=0
            do k=1,norder
               pp=pp+coefs(ind,k,j)*px(k)
               d1=d1+coefs(ind,k,j)*pxp(k)
               d2=d2+coefs(ind,k,j)*pxpp(k)
            enddo
            tmp(ind,1,j)=pp
            tmp(ind,2,j)=d1
            tmp(ind,3,j)=d2
         enddo
      enddo
c
      do ind=1,nd
         pp=0
         d1=0
         d2=0
         h1=0
         h2=0
         h3=0
         do j=1,norder
            pp=pp+tmp(ind,1,j)*py(j)
            d2=d2+tmp(ind,1,j)*pyp(j)
            h3=h3+tmp(ind,1,j)*pypp(j)

            d1=d1+tmp(ind,2,j)*py(j)
            h2=h2+tmp(ind,2,j)*pyp(j)

            h1=h1+tmp(ind,3,j)*py(j)
         enddo
         pot(ind)=pp
         grad(ind,1)=d1
         grad(ind,2)=d2
         
         hess(ind,1)=h1
         hess(ind,2)=h2
         hess(ind,3)=h3
      enddo

      return
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
      
