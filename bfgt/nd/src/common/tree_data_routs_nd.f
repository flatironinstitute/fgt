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
      subroutine treedata_trans_nd(ndim,nd,itype,nlevels,itree,iptr,
     1    boxsize,norder,fin,fout,umat)
c
c     This code converts an input tree data given on a tensor product grid on each leaf box
c     to the output tree data on each leaf box. 
c     Depending on the 1d transformation matrices, it could be:  
c     (1) converting function values to the coefficients of orthogonal polynomial expansion
c     coefficients;
c     (2) if umax converts coefficients to function values, then it's the inverse map of (1).
c     (3) converting fcoefs to the coefficients of its derivatives.
c 
c     It's user's responsibility to ensure that umat are the correct 1D
c     transformation matrices.
c 
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     itype - itype = i>0 - transformation is only along ith coordinates, 
c                         for example, calculating the expansion coefficients of the x-derivative
c             itype = 0 - transformation is done on both directions, for example, 
c                         val2coefs or coef2vals
c     nlevels - integer
c              number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     fin - double (nd,norder**2,nboxes)
c            input data on the tree
c     umat - 1D transformation matrices along each direction
c      
c     output:
c     fout - double precision (nd,norder**ndim,nboxes)
c            output data given on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fin(nd,norder**ndim,*)
      real *8 fout(nd,norder**ndim,*)

      real *8 umat(norder,norder,ndim)


      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call ortho_trans_nd(ndim,nd,itype,norder,
     1           fin(1,1,ibox),fout(1,1,ibox),umat)
             if (itype.gt.0) then
                do j=1,norder**ndim
                   do ind=1,nd
                      fout(ind,j,ibox)=fout(ind,j,ibox)*sc
                   enddo
                enddo
             endif
          endif
        enddo
C$OMP END PARALLEL DO
      enddo


      end
c
c
c
c
      subroutine treedata_evalg_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1    boxsize,norder,coefs,grad)
c
c     This code evaluates the gradient at the tensor grid on an adaptive tree
c     given the orthogonal polynomial expansion coefficients at each leaf box.
c 
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     coefs - double (nd,norder**2,nboxes)
c           expansion coefficients on quad tree
c
c     output:
c     grad - double precision (nd,2,norder**2,nboxes)
c            gradient values on tensor grid on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefs(nd,norder**ndim,*)
      real *8 grad(nd,ndim,norder**ndim,*)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)

      call ortho_eval_tables(ipoly,norder,umat,vmat,vpmat,vppmat)

      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call ortho_evalg_nd(ndim,nd,norder,coefs(1,1,ibox),sc,
     1           grad(1,1,1,ibox),vmat,vpmat)             
          endif
        enddo
C$OMP END PARALLEL DO
      enddo


      end
c
c
c
c
      subroutine treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1    boxsize,norder,coefs,grad,hess)
c
c     This code evaluates and hessian at the tensor 
c     grid on an adaptive tree given the orthogonal polynomial expansion coefficients 
c     at each leaf box.
c
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     coefs - double (nd,norder**2,nboxes)
c           expansion coefficients on quad tree
c
c     output:
c     grad - double precision (nd,ndim,norder**ndim,nboxes)
c            gradient values on tensor grid on each leaf box
c     hess - double precision (nd,nhess,norder**ndim,nboxes)
c            hessian values on tensor grid on each leaf box
c            in the order of uxx, uxy, uyy in 2d
c            and uxx,uyy,uzz,uxy,uxz,uyz in 3d. nhess=ndim*(ndim+1)/2
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefs(nd,norder**ndim,*)
      real *8 grad(nd,ndim,norder**ndim,*)
      real *8 hess(nd,ndim*(ndim+1)/2,norder**ndim,*)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)

      call ortho_eval_tables(ipoly,norder,umat,vmat,vpmat,vppmat)
      
      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call ortho_evalgh_nd(ndim,nd,norder,coefs(1,1,ibox),sc,
     1           grad(1,1,1,ibox),hess(1,1,1,ibox),
     2           vmat,vpmat,vppmat)
          endif
        enddo
C$OMP END PARALLEL DO
      enddo

      return
      end
c
c
c
c
      subroutine treedata_eval_laplacian_nd(ndim,nd,ipoly,nlevels,
     1    itree,iptr,boxsize,norder,coefs,rlap)
c
c     This code evaluates the laplacian at the tensor grid on an adaptive tree
c     given the orthogonal polynomial expansion coefficients at each leaf box.
c 
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     coefs - double (nd,norder**2,nboxes)
c           expansion coefficients on quad tree
c
c     output:
c     laplacian - double precision (nd,norder**ndim,nboxes)
c            gradient values on tensor grid on each leaf box
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefs(nd,norder**ndim,*)
      real *8 rlap(nd,norder**ndim,*)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)

      call ortho_eval_tables(ipoly,norder,umat,vmat,vpmat,vppmat)
      
      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call ortho_eval_laplacian_nd(ndim,nd,norder,
     1           coefs(1,1,ibox),sc,rlap(1,1,ibox),vmat,vppmat)             
          endif
        enddo
C$OMP END PARALLEL DO
      enddo

      return
      end
c
c
c
c
      subroutine treedata_eval_pot_nd_asym(ndim,nd,delta,ipoly,nasym,
     1    nlevels,itree,iptr,boxsize,norder,fvals,pot)
c
c     This code evaluates the potential value at the tensor grid on an adaptive tree
c     using asympotitic expansions 
c 
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     nasym - order of asympotitic expansion
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     fvals - double (nd,norder**ndim,nboxes)
c           function values on the tree
c
c     output:
c     pot - double precision (nd,norder**ndim,nboxes)
c           potential values on the tree
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fvals(nd,norder**ndim,*)
      real *8 pot(nd,norder**ndim,*)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)
      real *8 umat_nd(norder,norder,ndim)

      real *8, allocatable :: fcoefs(:,:)
      real *8, allocatable :: flvals(:,:)
      real *8, allocatable :: flcoefs(:,:)
      real *8, allocatable :: fl2vals(:,:)

      call ortho_eval_tables(ipoly,norder,umat,vmat,vpmat,vppmat)
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
      
      npbox=norder**ndim
      allocate(fcoefs(nd,npbox))
      allocate(flvals(nd,npbox))
      allocate(flcoefs(nd,npbox))
      allocate(fl2vals(nd,npbox))
      
      sqrtpi = sqrt(4*atan(1.0d0))
      d0 = sqrtpi*sqrt(delta)
      d2 = delta/4
      d4 = delta*delta/32

      c0 = d0**ndim
      c2 = c0*d2
      c4 = c0*d4

      itype=0
      do ilev = 0,nlevels
        sc=2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             if (nasym.eq.2) then
                call ortho_trans_nd(ndim,nd,itype,norder,
     1              fvals(1,1,ibox),fcoefs,umat_nd)
                call ortho_eval_laplacian_nd(ndim,nd,norder,fcoefs,
     1              sc,flvals,vmat,vppmat)
                do j=1,npbox
                do ind=1,nd
                   pot(ind,j,ibox)=c0*fvals(ind,j,ibox)
     1                 +c2*flvals(ind,j)
                enddo
                enddo
             elseif (nasym.eq.3) then
                call ortho_trans_nd(ndim,nd,itype,norder,
     1              fvals(1,1,ibox),fcoefs,umat_nd)
                call ortho_eval_laplacian_nd(ndim,nd,norder,fcoefs,
     1              sc,flvals,vmat,vppmat)
                call ortho_trans_nd(ndim,nd,itype,norder,flvals,
     1              flcoefs,umat_nd)
                call ortho_eval_laplacian_nd(ndim,nd,norder,flcoefs,
     1              sc,fl2vals,vmat,vppmat)
                do j=1,npbox
                do ind=1,nd
                   pot(ind,j,ibox)=c0*fvals(ind,j,ibox)
     1                 +c2*flvals(ind,j)+c4*fl2vals(ind,j)
                enddo
                enddo
             endif
          endif
        enddo
C$OMP END PARALLEL DO
      enddo

      return
      end
c
c
c
c
      subroutine treedata_coefs_p_to_g2d(ndim,nd,nlevels,itree,iptr,
     1    boxsize,norder,coefsp,coefsg,umat)
c
c     This code converts expansion coefficients of a tree data given on a tensor product
c     grid on each leaf box to the expansion coefficients of its gradient
c
c     input:
c     nd - integer,   number of functions
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     coefsp - double (nd,norder**2,nboxes)
c            expansion coefficients of the potential
c
c     output:
c     coefsg - double precision (nd,2,norder**2,nboxes)
c            expansion coefficients of the gradient
c
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 coefsp(nd,norder**ndim,*)
      real *8 coefsg(nd,ndim,norder**ndim,*)

      real *8 umat(norder,norder)

      do ilev = 0,nlevels
        sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,i,ind)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) then
             call ortho_coefsg_nd(ndim,nd,norder,coefsp(1,1,ibox),
     1           coefsg(1,1,1,ibox),umat)
             do j=1,norder**ndim
                do i=1,ndim
                   do ind=1,nd
                      coefsg(ind,i,j,ibox)=coefsg(ind,i,j,ibox)*sc
                   enddo
                enddo
             enddo
          endif
        enddo
C$OMP END PARALLEL DO
      enddo

      return
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
c     input:
c
c     nd - integer,   number of functions
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     norder - integer
c           order of expansions for input coefficients array
c     coefsp - double (nd,norder**2,nboxes)
c            expansion coefficients of the potential
c     output:
c     coefsg - double precision (nd,2,norder**2,nboxes)
c            expansion coefficients of the gradient
c     coefsh - double precision (nd,3,norder**2,nboxes)
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
             call ortho_coefsgh_2d(nd,norder,coefsp(1,1,ibox),
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

      return
      end
c
c
c
c
      subroutine treedata_evalt_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1    ltree,itree,iptr,tcenters,boxsize,fcoefs,nt,targ,fvals)
c
c     This code evaluates function values at nt target points where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c     input:
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input coefficients array
c     nboxes - integer
c            number of boxes
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
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
c     tcenters - double (2,nboxes) coordinates of box centers
c     boxsize - double(nboxes) box size of each box
c      
c      
c     fcoefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients on each leaf box
c     nt - number of target points
c     targ - (2,nt) coordinates of targets
c      
c      
c     output:
c     fvals - double precision (nd,nt) function values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 tcenters(ndim,nboxes),boxsize(0:nlevels)
      
      real *8 fcoefs(nd,norder**ndim,nboxes)
      real *8 targ(ndim,nt),fvals(nd,nt)
      real *8 xyz(ndim),cen(ndim)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      allocate(targsort(ndim,nt),potsort(nd,nt))
c      
c     reorder targets
c
      call dreorderf(ndim,nt,targ,targsort,itarg)
c
c     
c     
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts,j,i,cen,xyz)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
cccc            call prinf('ibox=*',ibox,1)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1
cccc               print *, ibox, istart,iend,npts
               if (npts.gt.0) then
                  do j=1,ndim
                     cen(j)=tcenters(j,ibox)
                  enddo
cccc                  print *, cen(1),cen(2)
                  do i=istart,iend
                     do j=1,ndim
                        xyz(j) = (targsort(j,i)-cen(j))*sc
                     enddo
                     call ortho_evalt_nd(ndim,nd,ipoly,norder,
     1                   fcoefs(1,1,ibox),xyz,potsort(1,i))
                  enddo
               endif
            endif
        enddo
C$OMP END PARALLEL DO
      enddo
c
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,fvals,itarg)

      return
      end
c
c
c
c
      subroutine treedata_evalpgt_nd(ndim,nd,ipoly,norder,nboxes,
     1    nlevels,ltree,itree,iptr,tcenters,boxsize,coefs,
     2    nt,targ,pot,grad)
c
c     This code evaluates potential and gradient at nt target points
c     where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c     input:
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input coefficients array
c     nboxes - integer
c            number of boxes
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
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
c     tcenters - double (2,nboxes) coordinates of box centers
c     boxsize - double(nboxes) box size of each box
c      
c      
c     coefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients of the potential on each leaf box
c     nt - number of target points
c     targ - (2,nt) coordinates of targets
c      
c      
c     output:
c     pot - double precision (nd,nt) potential values
c     grad - double precision (nd,2,nt) gradient values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      
      integer itree(ltree)
      real *8 tcenters(ndim,nboxes),boxsize(0:nlevels)
      
      real *8 coefs(nd,norder**ndim,nboxes)
      real *8 targ(ndim,nt),pot(nd,nt),grad(nd,ndim,nt)
      real *8 xyz(ndim),cen(ndim)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)
      real *8, allocatable :: gradsort(:,:,:)

      allocate(itarg(nt),itargse(2,nboxes))
      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(ndim,nt),potsort(nd,nt),gradsort(nd,ndim,nt))
c      
c     reorder targets
c
      call dreorderf(ndim,nt,targ,targsort,itarg)
c
c     
c     
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts,j,i,cen,xyz)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1

               if (npts.gt.0) then
                  do j=1,ndim
                     cen(j)=tcenters(j,ibox)
                  enddo
                  do i=istart,iend
                     do j=1,ndim
                        xyz(j) = (targsort(j,i)-cen(j))*sc
                     enddo
                     call ortho_evalpgt_nd(ndim,nd,ipoly,norder,
     1                   coefs(1,1,ibox),sc,xyz,
     2                   potsort(1,i),gradsort(1,1,i))
                  enddo
               endif
            endif
        enddo
C$OMP END PARALLEL DO
      enddo
c
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,pot,itarg)
      call dreorderi(nd*ndim,nt,gradsort,grad,itarg)

      return
      end
c
c
c
c
      subroutine treedata_evalpght_nd(ndim,nd,ipoly,norder,
     1    nboxes,nlevels,ltree,itree,iptr,tcenters,boxsize,coefs,
     2    nt,targ,pot,grad,hess)
c
c     This code evaluates potential, gradient, and hessian at nt target points
c     where the function is given 
c     as the orthogonal polynomial expansion coefficients on each leaf box
c 
c     input:
c
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input coefficients array
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
c     nboxes - integer
c            number of boxes
c     nlevels - integer
c            number of levels
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
c     tcenters - double (2,nboxes) coordinates of box centers
c     boxsize - double(nboxes) box size of each box
c      
c      
c     coefs - double (nd,norder*norder,nboxes)
c                  expansion coefficients of the potential on each leaf box
c     nt - number of target points
c     targ - (2,nt) coordinates of targets
c      
c      
c     output:
c     pot - double precision (nd,nt) potential values
c     grad - double precision (nd,2,nt) gradient values
c     hess - double precision (nd,3,nt) gradient values
c            
      implicit real *8 (a-h,o-z)
      integer nd
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 tcenters(ndim,nboxes),boxsize(0:nlevels)
      
      real *8 coefs(nd,norder**ndim,nboxes)
      real *8 targ(ndim,nt),pot(nd,nt),grad(nd,ndim,nt)
      real *8 hess(nd,ndim*(ndim+1)/2,nt)
      real *8 xyz(ndim),cen(ndim)
      
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: potsort(:,:)
      real *8, allocatable :: gradsort(:,:,:)
      real *8, allocatable :: hesssort(:,:,:)

      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(ndim,nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,itarg,itargse)
      
      allocate(targsort(ndim,nt),potsort(nd,nt),gradsort(nd,ndim,nt))
      allocate(hesssort(nd,ndim*(ndim+1)/2,nt))
c      
c     reorder targets
c
      call dreorderf(ndim,nt,targ,targsort,itarg)
c
      do ilev = 0,nlevels
         sc = 2.0d0/boxsize(ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,cen,xyz,j,i)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if (nchild.eq.0) then
               istart = itargse(1,ibox) 
               iend = itargse(2,ibox) 
               npts = iend-istart+1

               if (npts.gt.0) then
                  do j=1,ndim
                     cen(j)=tcenters(j,ibox)
                  enddo
                  do i=istart,iend
                     do j=1,ndim
                        xyz(j) = (targsort(j,i)-cen(j))*sc
                     enddo
                     call ortho_evalpght_nd(ndim,nd,ipoly,norder,
     1                   coefs(1,1,ibox),sc,xyz,
     2                   potsort(1,i),gradsort(1,1,i),hesssort(1,1,i))
                  enddo
               endif
            endif
        enddo
C$OMP END PARALLEL DO
      enddo
c     
c     resort the output arrays in input order
c
      call dreorderi(nd,nt,potsort,pot,itarg)
      call dreorderi(nd*ndim,nt,gradsort,grad,itarg)
      call dreorderi(nd*ndim*(ndim+1)/2,nt,hesssort,hess,itarg)

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
c     input:
c
c     nd - integer,   number of functions
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     npbox - integer
c           number of points in each leaf box
c     fex - double (nd,npbox,nboxes)
c           values of the reference tree data on each leaf box
c     fcomp - double (nd,npbox,nboxes)
c           values of the numerical tree data on each leaf box
c     output:
c     abserr - absolute l2 error
c     rnorm - l2 norm of fex
c     nleaf  - total number of leaf boxes
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
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
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
C$OMP END PARALLEL DO
      enddo

      abserr=sqrt(abserr)
      rnorm=sqrt(rnorm)

      return
      end
c
c
c
c
      
      subroutine treedata_lpnorm(ndim,iptype,ipoly,nd,nlevels,itree,
     1    iptr,boxsize,norder,npbox,fvals,rnorm,nleaf)
c
c     computes the lp norm of the tree date fvals
c 
c     input:
c
c     nd - integer,   number of functions
c     nlevels - integer
c            number of levels
c     itree - integer(ltree)
c            array containing the tree structure
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
c     npbox - integer
c           number of points in each leaf box
c     fvals - double (nd,npbox,nboxes)
c           values of the tree data on each leaf box
c
c     output:
c     rnorm - lp norm of fex
c     nleaf  - total number of leaf boxes
c      
      implicit real *8 (a-h,o-z)
      integer nd
      integer nlevels
      integer itree(*),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fvals(nd,npbox,*)

      real *8 wts(npbox),xs(ndim,npbox)

      itype = 1
      call polytens_exps_nd(ndim,ipoly,itype,norder,'f',xs,
     1    utmp,1,vtmp,1,wts)
      
      rnorm=0
      nleaf=0

      if (iptype.eq.0) then
         do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox = itree(2*ilev+1),itree(2*ilev+2)
               nchild = itree(iptr(4) + ibox-1)
               if(nchild.eq.0) then
                  nleaf=nleaf+1
                  rtmp=maxval(fvals(1:nd,1:npbox,ibox))
                  if (rtmp.gt.rnorm) rnorm=rtmp
               endif
            enddo
C$OMP END PARALLEL DO
         enddo

      else
         do ilev = 0,nlevels
            sc = (boxsize(ilev)/2)**ndim
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,j,ind)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox = itree(2*ilev+1),itree(2*ilev+2)
               nchild = itree(iptr(4) + ibox-1)
               if(nchild.eq.0) then
                  nleaf=nleaf+1
                  do i=1,npbox
                     do ind=1,nd
                        rnorm=rnorm
     1                      +abs(fvals(ind,i,ibox))**iptype*wts(i)*sc
                     enddo
                  enddo
               endif
            enddo
C$OMP END PARALLEL DO
         enddo
         rnorm=rnorm**(1.0d0/iptype)
      endif
      
      return
      end
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
      subroutine tree_data_c2p_copy_3d(nd,norder,fvals2,fvals,kk)
      implicit real *8 (a-h,o-z)
      real *8 fvals2(nd,norder/2,norder/2,norder/2)
      real *8 fvals(nd,norder,norder,norder)

      n=norder/2

      if (kk.eq.1) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.2) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.3) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j+n,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.4) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j+n,i)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.5) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.6) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k+n,j,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.7) then
         do i=1,n
         do j=1,n
         do k=1,n
         do ind=1,nd
            fvals(ind,k,j+n,i+n)=fvals2(ind,k,j,i)
         enddo
         enddo
         enddo
         enddo
      elseif (kk.eq.8) then
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
c     fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c     fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer iptr(8)
      integer itree(*),norder,npbox
      integer listcid(2**ndim)
      
      real *8 fvals(nd,npbox,*),fcoefs(nd,npbox,*)

      real *8 centers(ndim,*),boxsize(0:*)

      integer isgn(ndim,2**ndim)
      real *8 polyv(norder,norder,ndim,2**ndim)
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
      subroutine vol_tree_c2p_interp(nd,ndim,norder,npbox,polyv,umat_nd,
     1    fvals,fcoefs,iparentbox,listcid)

c     This subroutine computes the data on the parent box using the data from its
c     child boxes by interpolation.
c
c     input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     norder - integer
c           order of expansions for input function value array
c     polyv - precomputed matrices for carrying out interpolation from the parent
c                 box to its child boxes
c     umat_nd - precomputed matrices for converting function values on a tensor
c                 grid to its polynomial expansion coefficients
c
c     iparentbox - the index of the parent box to be refined
c     listcid - the indices of child boxes
c
c     input/output:
c     fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c     fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c
c
      implicit real *8 (a-h,o-z)
      integer listcid(2**ndim)
      
      real *8 fvals(nd,npbox,*),fcoefs(nd,npbox,*)

      real *8 polyv(norder,norder/2,ndim,2**ndim)
      real *8 umat_nd(norder,norder,ndim)
      real *8 fvals2(nd,npbox/2**ndim)
      
      mc=2**ndim
      n2=norder/2
      
      do j=1,mc
         jbox=listcid(j)
         call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,jbox),
     1       n2,fvals2,polyv(1,1,1,j))
         call tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,
     1       fvals(1,1,iparentbox),j)
      enddo
      call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,iparentbox),
     1    fcoefs(1,1,iparentbox),umat_nd)
      
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
      
