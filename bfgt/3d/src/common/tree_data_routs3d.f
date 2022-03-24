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
      subroutine treedata_trans3d(nd,itype,nlevels,itree,iptr,boxsize,
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
c                 itype = 3 - transformation is only along z-direction, for example, calculating
c                         the expansion coefficients of the z-derivative
c                 itype = 0 - transformation is done on all directions, for example, 
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
c         fin - double (nd,norder**3,nboxes)
c            input data on the tree
c     output:
c         fout - double precision (nd,norder**3,nboxes)
c            putput data given on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fin(nd,norder*norder*norder,nboxes)
      real *8 fout(nd,norder*norder*norder,nboxes)

      real *8 umat(norder,norder)
      real *8, allocatable :: work(:,:,:)

      allocate(work(nd,norder,norder,norder))

      do ilev = 0,nlevels
        sc = 2/boxsize(ilev)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call orth_trans3d(nd,itype,norder,fin(1,1,ibox),
     1           fout(1,1,ibox),umat,work)
             if (itype.gt.0) then
                do j=1,norder*norder*norder
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
      subroutine treedata_eval3d(nd,itype,itree,ltree,nboxes,nlevels,
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
c         tcenters - double (3,nboxes) coordinates of box centers
c         boxsize - double(nboxes) box size of each box
c      
c      
c         norder - integer
c           order of expansions for input coefficients array
c         fcoefs - double (nd,norder*norder*norder,nboxes)
c                  expansion coefficients on each leaf box
c         nt - number of target points
c         targ - (3,nt) coordinates of targets
c      
c      
c     output:
c         fvals - double precision (nd,nt) function values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8)
      real *8 tcenters(3,nboxes),boxsize(0:nlevels)
      
      real *8 fcoefs(nd,norder*norder*norder,nboxes)
      real *8 targ(3,nt),fvals(nd,nt)
      real *8 xy(3)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)

      allocate(itarg(nt),itargse(3,nboxes))
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(3,nt),potsort(nd,nt))
c      
c     reorder targets
c
      call dreorderf(3,nt,targ,targsort,itarg)
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
                  cz = tcenters(3,ibox)

                  do i=istart,iend
                     xy(1) = (targsort(1,i)-cx)/bs
                     xy(2) = (targsort(2,i)-cy)/bs
                     xy(3) = (targsort(3,i)-cz)/bs
                     call orth_eval3d(nd,itype,norder,fcoefs(1,1,ibox),
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
      subroutine orth_trans3d(nd,itype,norder,fin,fout,umat,work)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases. The polynomial could be any orthogonal
c     polyomial.
c
c     1. if umat converts function values to orthogonal polynomial expansion coefficients,  
c        then it converts function values on the tensor product grid to orthogonal
c        polynomial expansion coefficients in 3d.
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
c             itype = 3 - transformation is only along z-direction, for example, calculating
c                         the expansion coefficients of the z-derivative
c             itype = 0 - transformation is done on all directions, for example, 
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
      real *8 fin(nd,norder,norder,norder)
      real *8 fout(nd,norder,norder,norder),umat(norder,norder)
      real *8, allocatable :: work(:,:,:,:),work2(:,:,:,:)

      allocate(work(nd,norder,norder,norder))
      allocate(work2(nd,norder,norder,norder))

      if (itype.eq.0) then
c        transform in x
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fin(ind,k1,j,i)
               enddo
               work(ind,k,j,i)=dd
            enddo
         enddo
         enddo
         enddo
c        transform in y
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do j1=1,norder
                  dd=dd+umat(j,j1)*work(ind,k,j1,i)
               enddo
               work2(ind,k,j,i)=dd
            enddo
         enddo
         enddo
         enddo
c        transform in z
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do i1=1,norder
                  dd=dd+umat(i,i1)*work2(ind,k,j,i1)
               enddo
               fout(ind,k,j,i)=dd
            enddo
         enddo
         enddo
         enddo
      elseif (itype.eq. 1) then
c        transform in x
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fin(ind,k1,j,i)
               enddo
               fout(ind,k,j,i)=dd
            enddo
         enddo
         enddo
         enddo
      elseif (itype.eq. 2) then
c        transform in y
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do j1=1,norder
                  dd=dd+umat(j,j1)*fin(ind,k,j1,i)
               enddo
               fout(ind,k,j,i)=dd
            enddo
         enddo
         enddo
         enddo
      elseif (itype.eq. 3) then
c        transform in z
         do i=1,norder
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do i1=1,norder
                  dd=dd+umat(i,i1)*fin(ind,k,j,i1)
               enddo
               fout(ind,k,j,i)=dd
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
      subroutine orth_eval3d(nd,itype,norder,fcoefs,targ,fval)
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
      real *8 fcoefs(nd,norder,norder,norder)
      real *8 targ(3),fval(nd)

      real *8 px(norder),py(norder),pz(norder)
      real *8, allocatable :: tmp(:,:,:), tmp2(:,:)

      allocate(tmp(nd,norder,norder))
      allocate(tmp2(nd,norder))

      x=targ(1)
      y=targ(2)
      z=targ(3)
      
      if (itype .eq. 1) then
         call legepols(x,norder-1,px)
         call legepols(y,norder-1,py)
         call legepols(z,norder-1,pz)
      elseif (itype .eq. 2) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
         call chebpols(z,norder-1,pz)
      endif
c     tmp_{ij}=sum_k f_{ijk} p_k(x)
      do i=1,norder
      do j=1,norder
         do ind=1,nd
            dd=0
            do k=1,norder
               dd=dd+fcoefs(ind,k,j,i)*px(k)
            enddo
            tmp(ind,j,i)=dd
         enddo
      enddo
      enddo
c     tmp2_{i} = sum_j tmp_{ij} p_j(y)
      do i=1,norder
      do ind=1,nd
         dd=0
         do j=1,norder
            dd=dd+tmp(ind,j,i)*py(j)
         enddo
         tmp2(ind,i)=dd
      enddo
      enddo
c     fval = sum_i tmp2_{i} p_i(z)
      do ind=1,nd
         dd=0
         do i=1,norder
            dd=dd+tmp2(ind,i)*pz(j)
         enddo
         fval(ind)=dd
      enddo
      

      return
      end
c
c
c
c
      
