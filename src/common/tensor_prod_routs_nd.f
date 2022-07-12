c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c     this is the end of the debugging code and the beginning
c     of the tensor product routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c     This file contains a set of subroutines for the handling
c     of tensor product in n dimensions.
c     Following is a brief description of these
c     subroutines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     general transform (val2coef,coef2val,ceof2der, etc.)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_trans_nd(ndim,nd,itype,norder,fin,fout,umat)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases. 
c
c     input:
c     nd - number of input data
c     itype - itype = 1 - transformation is only along x-direction, for example, calculating
c                         the expansion coefficients of the x-derivative
c             itype = 2 - transformation is only along y-direction, for example, calculating
c                         the expansion coefficients of the y-derivative
c             itype = 3 - transformation is only along z-direction, for example, calculating
c                         the expansion coefficients of the z-derivative
c             itype = 0 - transformation is done on both directions, for example, 
c                         val2coefs or coef2vals
c      
c     norder - order of the expansion in each direction
c     fin - input data
c     umat - 1D transformation matrix along each direction
c      
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,norder**ndim)
      real *8 fout(nd,norder**ndim)
      real *8 umat(norder,norder,ndim)

      if (ndim.eq.1) then
         call ortho_trans_1d(nd,itype,norder,fin,fout,umat)
      elseif (ndim.eq.2) then
         call ortho_trans_2d(nd,itype,norder,fin,fout,umat)
      elseif (ndim.eq.3) then
         call ortho_trans_3d(nd,itype,norder,fin,fout,umat)
      endif
      
      return
      end
c
c
c
c
      subroutine ortho_trans_1d(nd,itype,norder,fin,fout,umat)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases - value to coefficients, coefficients
c     to value, coefficients of the function to coefficients of the derivative, etc.
c
c     input:
c     nd - number of input data
c     itype - not used in 1d case
c      
c     norder - order of the expansion in each direction
c     fin - input data
c     umat - 1D transformation matrix 
c      
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,norder)
      real *8 fout(nd,norder)
      real *8 umat(norder,norder)

      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*fin(ind,k1)
            enddo
            fout(ind,k)=dd
         enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine ortho_trans_2d(nd,itype,norder,fin,fout,umat)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases. 
c
c     1. if uxmat,uymat convert function values to expansion coefficients,  
c        then it converts function values on the tensor product grid to 
c        expansion coefficients in 2d.
c     2. if uxmax,uymat convert coefficients to function values, then it's the inverse map of (1).
c     3. if uxmax,uymat are differentiation matrices in the coefficient space, then it converts
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
c     norder - order of the expansion in each direction
c     fin - input data
c     uxmat - 1D transformation matrix along the x-direction
c     uymat - 1D transformation matrix along the y-direction
c      
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,norder,norder)
      real *8 fout(nd,norder,norder)
      real *8 umat(norder,norder,2)
      real *8, allocatable :: work(:,:,:)
      
      allocate(work(nd,norder,norder))

      if (itype.eq.0) then
c        transform in x
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1,1)*fin(ind,k1,j)
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
                  dd=dd+umat(i,i1,2)*work(ind,k,i1)
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
                  dd=dd+umat(k,k1,1)*fin(ind,k1,j)
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
                  dd=dd+umat(k,k1,2)*fin(ind,j,k1)
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
      subroutine ortho_trans_3d(nd,itype,norder,fin,fout,umat)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. Depending on the 1d transformation matrix umat, it can
c     be used in the following several cases. 
c
c     input:
c     nd - number of input data
c     itype - itype = 1 - transformation is only along x-direction, for example, calculating
c                         the expansion coefficients of the x-derivative
c             itype = 2 - transformation is only along y-direction, for example, calculating
c                         the expansion coefficients of the y-derivative
c             itype = 3 - transformation is only along z-direction, for example, calculating
c                         the expansion coefficients of the z-derivative
c             itype = 0 - transformation is done on both directions, for example, 
c                         val2coefs or coef2vals
c      
c     norder - order of the expansion in each direction
c     fin - input data
c     umat - 1D transformation matrix along each direction
c      
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,norder,norder,norder)
      real *8 fout(nd,norder,norder,norder)
      real *8 umat(norder,norder,3)
      real *8, allocatable :: work(:,:,:,:)
      real *8, allocatable :: work2(:,:,:,:)
      
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
                  dd=dd+umat(k,k1,1)*fin(ind,k1,j,i)
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
                  dd=dd+umat(j,j1,2)*work(ind,k,j1,i)
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
                  dd=dd+umat(i,i1,3)*work2(ind,k,j,i1)
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
                  dd=dd+umat(k,k1,1)*fin(ind,k1,j,i)
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
                  dd=dd+umat(j,j1,2)*fin(ind,k,j1,i)
               enddo
               fout(ind,j,k,i)=dd
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
                  dd=dd+umat(i,i1,3)*fin(ind,k,j,i1)
               enddo
               fout(ind,j,k,i)=dd
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate function on a general tensor grid given its coefficients 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_eval_nd(ndim,nd,n1,fin,n2,fout,umat)
c
c     this subroutine transform the input data given on a tensor product grid
c     to the output data. 
c
c     input:
c     nd - number of input data
c     n1 - number of points in each direction for input data
c     fin - input data
c     umat - 1D transformation matrices along each direction
c     n2 - number of points in each direction for output data
c     output:
c     fout - the transform data, see the above explanation.
c
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,n1**ndim)
      real *8 fout(nd,n2**ndim)
      real *8 umat(n1,n2,ndim)

      if (ndim.eq.1) then
         call ortho_eval_1d(nd,n1,fin,n2,fout,umat)
      elseif (ndim.eq.2) then
         call ortho_eval_2d(nd,n1,fin,n2,fout,umat)
      elseif (ndim.eq.3) then
         call ortho_eval_3d(nd,n1,fin,n2,fout,umat)
      endif
      
      return
      end
c
c
c
c
      subroutine ortho_eval_1d(nd,n1,fin,n2,fout,umat)
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,n1)
      real *8 fout(nd,n2)
      real *8 umat(n1,n2)

      do k=1,n2
         do ind=1,nd
            dd=0
            do k1=1,n1
               dd=dd+umat(k1,k)*fin(ind,k1)
            enddo
            fout(ind,k)=dd
         enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine ortho_eval_2d(nd,n1,fin,n2,fout,umat)
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,n1,n1)
      real *8 fout(nd,n2,n2)
      real *8 umat(n1,n2,2)
      real *8, allocatable :: work(:,:,:)
      
      allocate(work(nd,n2,n1))

c     transform in x
      do j=1,n1
      do k=1,n2
         do ind=1,nd
            dd=0
            do k1=1,n1
               dd=dd+umat(k1,k,1)*fin(ind,k1,j)
            enddo
            work(ind,k,j)=dd
         enddo
      enddo
      enddo
c     transform in y
      do i=1,n2
      do k=1,n2
         do ind=1,nd
            dd=0
            do i1=1,n1
               dd=dd+umat(i1,i,2)*work(ind,k,i1)
            enddo
            fout(ind,k,i)=dd
         enddo
      enddo
      enddo

      return
      end
c
c
c
c
      subroutine ortho_eval_3d(nd,n1,fin,n2,fout,umat)
      implicit real *8 (a-h,o-z)
      real *8 fin(nd,n1,n1,n1)
      real *8 fout(nd,n2,n2,n2)
      real *8 umat(n1,n2,3)
      real *8, allocatable :: work(:,:,:,:)
      real *8, allocatable :: work2(:,:,:,:)
      
      allocate(work(nd,n2,n1,n1))
      allocate(work2(nd,n2,n2,n1))

c     transform in x
      do i=1,n1
      do j=1,n1
      do k=1,n2
         do ind=1,nd
            dd=0
            do k1=1,n1
               dd=dd+umat(k1,k,1)*fin(ind,k1,j,i)
            enddo
            work(ind,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo
c     transform in y
      do i=1,n1
      do j=1,n2
      do k=1,n2
         do ind=1,nd
            dd=0
            do j1=1,n1
               dd=dd+umat(j1,j,2)*work(ind,k,j1,i)
            enddo
            work2(ind,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo
c     transform in z
      do i=1,n2
      do j=1,n2
      do k=1,n2
         do ind=1,nd
            dd=0
            do i1=1,n1
               dd=dd+umat(i1,i,3)*work2(ind,k,j,i1)
            enddo
            fout(ind,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo
      
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate gradient at tensor grid given potential coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ortho_evalg_nd(ndim,nd,norder,coefs,sc,
     1    grad,vmat,vpmat)
c
c     this subroutine evaluates the gradient at the tensor grid
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     vmat - matrix converting coefficients into function values
c     vpmat - matrix converting coefficients into values of the first derivative
c      
c     output:
c     grad - value of the gradient
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder**ndim)
      real *8 grad(nd,ndim,norder**ndim)

      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)

      if (ndim.eq.1) then
         call ortho_evalg_1d(nd,norder,coefs,sc,grad,vmat,vpmat)
      elseif (ndim.eq.2) then
         call ortho_evalg_2d(nd,norder,coefs,sc,grad,vmat,vpmat)
      elseif (ndim.eq.3) then
         call ortho_evalg_3d(nd,norder,coefs,sc,grad,vmat,vpmat)
      endif
      
      return
      end
c
c
c
c
C*********************************************************************C
      subroutine ortho_evalg_1d(nd,n,fcoefs,sc,grad,
     1    vmat,vpmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n),grad(nd,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)
c
      
      do ind = 1,nd
c        transform in x
         do k1=1,n
            cdx=0
            do j1=1,n
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1)
            enddo
            grad(ind,k1)=cdx*sc
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
C*********************************************************************C
      subroutine ortho_evalg_2d(nd,n,fcoefs,sc,grad,
     1    vmat,vpmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n),pot(nd,n,n),grad(nd,2,n,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)
      real *8 ff(n,n)
      real *8 ffx(n,n)
c
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=1,n
            cd=0
            cdx=0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2)
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
            ffx(k1,j2)=cdx
         enddo
         enddo
c        transfrom in y
         do k2=1,n
         do k1=1,n
c            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            do j2=1,n
               cdy=cdy+vpmat(k2,j2)*ff(k1,j2)
               cdx=cdx+vmat(k2,j2)*ffx(k1,j2)
            enddo
            grad(ind,1,k1,k2)=cdx*sc
            grad(ind,2,k1,k2)=cdy*sc
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c
C
c
C*********************************************************************C
      subroutine ortho_evalg_3d(nd,n,fcoefs,sc,grad,
     1    vmat,vpmat)
C*********************************************************************C
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)

      real *8 ff(n,n,n),ff2(n,n,n)
      real *8 ffx(n,n,n),ff2x(n,n,n)
      real *8 ff2y(n,n,n)
c
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=1,n
            cd=0
            cdx=0.0d0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2,j3)
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
            ffx(k1,j2,j3)=cdx
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=1,n
         do k1=1,n
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            do j2=1,n
               cd=cd+vmat(k2,j2)*ff(k1,j2,j3)
               cdy=cdy+vpmat(k2,j2)*ff(k1,j2,j3)
               
               cdx=cdx+vmat(k2,j2)*ffx(k1,j2,j3)
            enddo
            ff2(k1,k2,j3)=cd
            ff2x(k1,k2,j3)=cdx
            ff2y(k1,k2,j3)=cdy
         enddo
         enddo
         enddo

c        transform in z
         do k3=1,n
         do k2=1,n
         do k1=1,n
            cdx = 0.0d0
            cdy = 0.0d0
            cdz = 0.0d0
            do j3=1,n
               cdz=cdz+vpmat(k3,j3)*ff2(k1,k2,j3)
               
               cdx=cdx+vmat(k3,j3)*ff2x(k1,k2,j3)
               cdy=cdy+vmat(k3,j3)*ff2y(k1,k2,j3)
            enddo
            grad(ind,1,k1,k2,k3)=cdx*sc
            grad(ind,2,k1,k2,k3)=cdy*sc
            grad(ind,3,k1,k2,k3)=cdz*sc
         enddo
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate gradient and hessian at tensor grid given potential coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_evalgh_nd(ndim,nd,norder,coefs,sc,grad,
     1    hess,vmat,vpmat,vppmat)
c
c     this subroutine evaluates the hessian at the tensor grid
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     vmat - matrix converting coefficients into function values
c     vpmat - matrix converting coefficients into values of the first derivative
c     vppmat - matrix converting coefficients into values of the second derivative
c
c     output:
c     grad - value of the gradient
c     hess - value of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder**ndim)
      real *8 grad(nd,ndim,norder**ndim)
      real *8 hess(nd,ndim*(ndim+1)/2,norder**ndim)

      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)

      if (ndim.eq.1) then
         call ortho_evalgh_1d(nd,norder,coefs,sc,grad,hess,vmat,
     1       vpmat,vppmat)
      elseif (ndim.eq.2) then
         call ortho_evalgh_2d(nd,norder,coefs,sc,grad,hess,vmat,
     1       vpmat,vppmat)
      elseif (ndim.eq.3) then
         call ortho_evalgh_3d(nd,norder,coefs,sc,grad,hess,vmat,
     1       vpmat,vppmat)
      endif
      
      return
      end
c
c
c
c
C*********************************************************************C
      subroutine ortho_evalgh_1d(nd,n,fcoefs,sc,grad,hess,
     1    vmat,vpmat,vppmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n),grad(nd,n),hess(nd,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)
      real *8 vppmat(n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do k1=1,n
c            cd=0
            cdx=0
            cdxx=0
            do j1=1,n
c               cd=cd+vmat(k1,j1)*fcoefs(ind,j1)
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1)
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1)
            enddo
            grad(ind,k1)=cdx*sc
            hess(ind,k1)=cdxx*sc2
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
      subroutine ortho_evalgh_2d(nd,n,fcoefs,sc,grad,hess,
     1    vmat,vpmat,vppmat)
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n),grad(nd,2,n,n),hess(nd,3,n,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)
      real *8 vppmat(n,n)

      real *8 ff(n,n)
      real *8 ffx(n,n)
      real *8 ffxx(n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=1,n
            cd=0
            cdx=0
            cdxx=0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2)
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1,j2)
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
            ffx(k1,j2)=cdx
            ffxx(k1,j2)=cdxx
         enddo
         enddo
c        transfrom in y
         do k2=1,n
         do k1=1,n
c            cd = 0.0d0
            cdx = 0.0d0
            cdy = 0.0d0
            cdxx = 0.0d0
            cdxy = 0.0d0
            cdyy = 0.0d0
            do j2=1,n
c               cd=cd+vmat(k2,j2)*ff(k1,j2)
               cdy=cdy+vpmat(k2,j2)*ff(k1,j2)
               cdyy=cdyy+vppmat(k2,j2)*ff(k1,j2)

               cdx=cdx+vmat(k2,j2)*ffx(k1,j2)
               cdxy=cdxy+vpmat(k2,j2)*ffx(k1,j2)

               cdxx=cdxx+vmat(k2,j2)*ffxx(k1,j2)
            enddo
            grad(ind,1,k1,k2)=cdx*sc
            grad(ind,2,k1,k2)=cdy*sc

            hess(ind,1,k1,k2)=cdxx*sc2
            hess(ind,2,k1,k2)=cdxy*sc2
            hess(ind,3,k1,k2)=cdyy*sc2
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c
C
c
      subroutine ortho_evalgh_3d(nd,n,fcoefs,sc,grad,hess,
     1    vmat,vpmat,vppmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n,n)
c      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 hess(nd,6,n,n,n)

      real *8 vmat(n,n)
      real *8 vpmat(n,n)
      real *8 vppmat(n,n)

      real *8 ff(n,n,n)
      real *8 ffx(n,n,n)
      real *8 ffxx(n,n,n)
      
      real *8 ff2(n,n,n)
      real *8 ff2x(n,n,n)
      real *8 ff2xy(n,n,n)
      real *8 ff2xx(n,n,n)
      
      real *8 ff2y(n,n,n)
      real *8 ff2yy(n,n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=1,n
            cd=0
            cdx=0.0d0
            cdxx=0.0d0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2,j3)
               cdx=cdx+vpmat(k1,j1)*fcoefs(ind,j1,j2,j3)
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
            ffx(k1,j2,j3)=cdx
            ffxx(k1,j2,j3)=cdxx
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=1,n
         do k1=1,n
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            cdxx = 0.0d0
            cdxy = 0.0d0
            cdyy = 0.0d0
            do j2=1,n
               cd   = cd   +   vmat(k2,j2)*ff(k1,j2,j3)
               cdy  = cdy  +  vpmat(k2,j2)*ff(k1,j2,j3)
               cdyy = cdyy + vppmat(k2,j2)*ff(k1,j2,j3)

               cdx  = cdx  +   vmat(k2,j2)*ffx(k1,j2,j3)
               cdxy = cdxy +  vpmat(k2,j2)*ffx(k1,j2,j3)
               
               cdxx = cdxx +   vmat(k2,j2)*ffxx(k1,j2,j3)
            enddo
            ff2(k1,k2,j3)=cd
            ff2x(k1,k2,j3)=cdx
            ff2xx(k1,k2,j3)=cdxx
            ff2y(k1,k2,j3)=cdy
            ff2yy(k1,k2,j3)=cdyy
            ff2xy(k1,k2,j3)=cdxy
         enddo
         enddo
         enddo

c        transform in z
         do k3=1,n
         do k2=1,n
         do k1=1,n
c            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            cdz = 0.0d0
            cdxx = 0.0d0
            cdyy = 0.0d0
            cdzz = 0.0d0
            cdxy = 0.0d0
            cdxz = 0.0d0
            cdyz = 0.0d0
            do j3=1,n
c               cd   = cd   +   vmat(k3,j3)*ff2(k1,k2,j3)
               cdz  = cdz  +  vpmat(k3,j3)*ff2(k1,k2,j3)
               cdzz = cdzz + vppmat(k3,j3)*ff2(k1,k2,j3)

               cdx  = cdx  +  vmat(k3,j3)*ff2x(k1,k2,j3)
               cdxz = cdxz + vpmat(k3,j3)*ff2x(k1,k2,j3)

               cdy  = cdy  +  vmat(k3,j3)*ff2y(k1,k2,j3)
               cdyz = cdyz + vpmat(k3,j3)*ff2y(k1,k2,j3)
               
               cdxx = cdxx + vmat(k3,j3)*ff2xx(k1,k2,j3)
               cdxy = cdxy + vmat(k3,j3)*ff2xy(k1,k2,j3)
               cdyy = cdyy + vmat(k3,j3)*ff2yy(k1,k2,j3)
            enddo
c            pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
            grad(ind,1,k1,k2,k3)=cdx*sc
            grad(ind,2,k1,k2,k3)=cdy*sc
            grad(ind,3,k1,k2,k3)=cdz*sc

            hess(ind,1,k1,k2,k3)=cdxx*sc2
            hess(ind,2,k1,k2,k3)=cdyy*sc2
            hess(ind,3,k1,k2,k3)=cdzz*sc2

            hess(ind,4,k1,k2,k3)=cdxy*sc2
            hess(ind,5,k1,k2,k3)=cdxz*sc2
            hess(ind,6,k1,k2,k3)=cdyz*sc2
         enddo
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c      
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate the laplacian at tensor grid given potential coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_eval_laplacian_nd(ndim,nd,norder,coefs,sc,rlap,
     1    vmat,vppmat)
c
c     this subroutine evaluates the laplacian at the tensor grid
c
c     input:
c     ndim - dimension of the underlying space
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     vmat - matrix converting coefficients into function values
c     vppmat - matrix converting coefficients into values of the second derivative
c
c     output:
c     rlap - value of the laplacian on the tensor grid
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder**ndim)
      real *8 rlap(nd,norder**ndim)

      real *8 vmat(norder,norder)
      real *8 vppmat(norder,norder)

      if (ndim.eq.1) then
         call ortho_eval_laplacian_1d(nd,norder,coefs,sc,rlap,
     1    vmat,vppmat)
      elseif (ndim.eq.2) then
         call ortho_eval_laplacian_2d(nd,norder,coefs,sc,rlap,
     1    vmat,vppmat)
      elseif (ndim.eq.3) then
         call ortho_eval_laplacian_3d(nd,norder,coefs,sc,rlap,
     1    vmat,vppmat)
      endif

      return
      end
c
c
c
c
C*********************************************************************C
      subroutine ortho_eval_laplacian_1d(nd,n,fcoefs,sc,rlap,
     1    vmat,vppmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n),rlap(nd,n)
      real *8 vmat(n,n)
      real *8 vppmat(n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do k1=1,n
            cdxx=0
            do j1=1,n
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1)
            enddo
            rlap(ind,k1)=cdxx*sc2
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
      subroutine ortho_eval_laplacian_2d(nd,n,fcoefs,sc,rlap,
     1    vmat,vppmat)
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n),rlap(nd,n,n)
      real *8 vmat(n,n)
      real *8 vppmat(n,n)

      real *8 ff(n,n)
      real *8 ffxx(n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=1,n
            cd=0
            cdxx=0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2)
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
            ffxx(k1,j2)=cdxx
         enddo
         enddo
c        transfrom in y
         do k2=1,n
         do k1=1,n
            cdxx = 0.0d0
            cdyy = 0.0d0
            do j2=1,n
               cdxx=cdxx+vmat(k2,j2)*ffxx(k1,j2)
               cdyy=cdyy+vppmat(k2,j2)*ff(k1,j2)
            enddo
            rlap(ind,k1,k2)=(cdxx+cdyy)*sc2
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c
C
c
      subroutine ortho_eval_laplacian_3d(nd,n,fcoefs,sc,rlap,
     1    vmat,vppmat)
C*********************************************************************C
      implicit real *8 (a-h,o-z)
      real *8 fcoefs(nd,n,n,n)
      real *8 rlap(nd,n,n,n)

      real *8 vmat(n,n)
      real *8 vppmat(n,n)

      real *8 ff(n,n,n)
      real *8 ffxx(n,n,n)
      
      real *8 ff2(n,n,n)
      real *8 ff2xx(n,n,n)
      real *8 ff2yy(n,n,n)
c
      sc2=sc*sc
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=1,n
            cd=0
            cdxx=0.0d0
            do j1=1,n
               cd=cd+vmat(k1,j1)*fcoefs(ind,j1,j2,j3)
               cdxx=cdxx+vppmat(k1,j1)*fcoefs(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
            ffxx(k1,j2,j3)=cdxx
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=1,n
         do k1=1,n
            cd=0
            cdxx = 0.0d0
            cdyy = 0.0d0
            do j2=1,n
               cd   = cd   +   vmat(k2,j2)*ff(k1,j2,j3)
               cdyy = cdyy + vppmat(k2,j2)*ff(k1,j2,j3)
               cdxx = cdxx +   vmat(k2,j2)*ffxx(k1,j2,j3)
            enddo
            ff2(k1,k2,j3)=cd
            ff2xx(k1,k2,j3)=cdxx
            ff2yy(k1,k2,j3)=cdyy
         enddo
         enddo
         enddo

c        transform in z
         do k3=1,n
         do k2=1,n
         do k1=1,n
            cdxx = 0.0d0
            cdyy = 0.0d0
            cdzz = 0.0d0
            do j3=1,n
               cdxx = cdxx + vmat(k3,j3)*ff2xx(k1,k2,j3)
               cdyy = cdyy + vmat(k3,j3)*ff2yy(k1,k2,j3)
               cdzz = cdzz + vppmat(k3,j3)*ff2(k1,k2,j3)
            enddo
            rlap(ind,k1,k2,k3)=(cdxx+cdyy+cdzz)*sc2
         enddo
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c      
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     potential coefficients to gradient coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_coefsg_nd(ndim,nd,norder,coefsp,coefsg,umat)
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
      real *8 coefsp(nd,norder**ndim)
      real *8 coefsg(nd,ndim,norder**ndim),umat(norder,norder)

      if (ndim.eq.1) then
         call ortho_coefsg_1d(nd,norder,coefsp,coefsg,umat)         
      elseif (ndim.eq.2) then
         call ortho_coefsg_2d(nd,norder,coefsp,coefsg,umat)         
      elseif (ndim.eq.3) then
         call ortho_coefsg_3d(nd,norder,coefsp,coefsg,umat)         
      endif
      
      return
      end
c
c
c
c
      subroutine ortho_coefsg_1d(nd,norder,coefsp,coefsg,umat)
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
      real *8 coefsp(nd,norder)
      real *8 coefsg(nd,norder),umat(norder,norder)

c     transform in x
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1)
            enddo
            coefsg(ind,k)=dd
         enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine ortho_coefsg_2d(nd,norder,coefsp,coefsg,umat)
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
      subroutine ortho_coefsg_3d(nd,norder,coefsp,coefsg,umat)
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
      real *8 coefsp(nd,norder,norder,norder)
      real *8 coefsg(nd,3,norder,norder,norder),umat(norder,norder)

c     transform in x
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1,j,i)
            enddo
            coefsg(ind,1,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo
      
c     transform in y
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*coefsp(ind,k,j1,i)
            enddo
            coefsg(ind,2,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo

c     transform in z
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do i1=1,norder
               dd=dd+umat(i,i1)*coefsp(ind,k,j,i1)
            enddo
            coefsg(ind,3,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo

      
      return
      end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     potential coefficients to gradient and hessian coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_coefsgh_nd(ndim,nd,norder,coefsp,
     1    coefsg,coefsh,umat)
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
c     coefsh - expansion coefficients of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefsp(nd,norder**ndim)
      real *8 coefsg(nd,ndim,norder**ndim),umat(norder,norder)
      real *8 coefsh(nd,ndim*(ndim+1)/2,norder**ndim)

      if (ndim.eq.1) then
         call ortho_coefsgh_1d(nd,norder,coefsp,coefsg,coefsh,umat)         
      elseif (ndim.eq.2) then
         call ortho_coefsgh_2d(nd,norder,coefsp,coefsg,coefsh,umat)         
      elseif (ndim.eq.3) then
         call ortho_coefsgh_3d(nd,norder,coefsp,coefsg,coefsh,umat)         
      endif
      
      return
      end
c      
c      
c      
c      
c      
      subroutine ortho_coefsgh_1d(nd,norder,coefsp,coefsg,coefsh,umat)
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
      real *8 coefsp(nd,norder)
      real *8 coefsg(nd,norder),umat(norder,norder)
      real *8 coefsh(nd,norder)

c     compute the expansion coefficients of the gradient
c     differentiation in x
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1)
            enddo
            coefsg(ind,k)=dd
         enddo
      enddo

c     compute the expansion coefficients of the hessian
c     u_{xx}
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsg(ind,k1)
            enddo
            coefsh(ind,k)=dd
         enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine ortho_coefsgh_2d(nd,norder,coefsp,coefsg,coefsh,umat)
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
      subroutine ortho_coefsgh_3d(nd,norder,coefsp,coefsg,coefsh,umat)
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
      real *8 coefsp(nd,norder,norder,norder)
      real *8 coefsg(nd,3,norder,norder,norder),umat(norder,norder)
      real *8 coefsh(nd,6,norder,norder,norder)

c     compute the expansion coefficients of the gradient
c     differentiation in x
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsp(ind,k1,j,i)
            enddo
            coefsg(ind,1,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo
      
c     differentiation in y
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*coefsp(ind,k,j1,i)
            enddo
            coefsg(ind,2,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     differentiation in z
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do i1=1,norder
               dd=dd+umat(i,i1)*coefsp(ind,k,j,i1)
            enddo
            coefsg(ind,3,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     compute the expansion coefficients of the hessian
c     u_{xx}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*coefsg(ind,1,k1,j,i)
            enddo
            coefsh(ind,1,k,j,i)=dd
         enddo
      enddo
      enddo
      enddo

c     u_{yy}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*coefsg(ind,2,k,j1,i)
            enddo
            coefsh(ind,2,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     u_{zz}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do i1=1,norder
               dd=dd+umat(i,i1)*coefsg(ind,3,k,j,i1)
            enddo
            coefsh(ind,3,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     u_{xy}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*coefsg(ind,1,k,j1,i)
            enddo
            coefsh(ind,4,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     u_{xz}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do i1=1,norder
               dd=dd+umat(i,i1)*coefsg(ind,1,k,j,i1)
            enddo
            coefsh(ind,5,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

c     u_{yz}
      do i=1,norder
      do j=1,norder
      do k=1,norder
         do ind=1,nd
            dd=0
            do i1=1,norder
               dd=dd+umat(i,i1)*coefsg(ind,2,k,j,i1)
            enddo
            coefsh(ind,6,j,k,i)=dd
         enddo
      enddo
      enddo
      enddo

      
      
      return
      end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate potential at arbitrary target 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_evalt_nd(ndim,nd,ipoly,norder,coefs,targ,pot)
c
c     this subroutine evaluates the orthogonal polynomial expansion at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - itype=0 Legendre polynomials
c                   1 Chebyshev polynomials
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients
c     targ - xy coordinates of the target point
c      
c     output:
c     pot - value of the function at the target point
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder**ndim)
      real *8 targ(ndim),pot(nd)

      if (ndim.eq.1) then
         call ortho_evalt_1d(nd,ipoly,norder,coefs,targ,pot)
      elseif (ndim.eq.2) then
         call ortho_evalt_2d(nd,ipoly,norder,coefs,targ,pot)
      elseif (ndim.eq.3) then
         call ortho_evalt_3d(nd,ipoly,norder,coefs,targ,pot)
      endif
      
      return
      end
c
c
c
      subroutine ortho_evalt_1d(nd,ipoly,norder,coefs,
     1    targ,pot)
c
c     this subroutine evaluates the potential at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     targ - xy coordinates of the target point
c      
c     output:
c     pot - value of the potential
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder)
      real *8 targ,pot(nd)

      real *8 px(norder)
      
      x=targ
      
      if (ipoly .eq. 0) then
         call legepols(x,norder-1,px)
      elseif (ipoly .eq. 1) then
         call chebpols(x,norder-1,px)
      endif
c     
      do ind=1,nd
         pp=0
         do k=1,norder
            pp=pp+coefs(ind,k)*px(k)
         enddo
         pot(ind)=pp
      enddo

      return
      end
c
c
c
c
      subroutine ortho_evalt_2d(nd,ipoly,norder,coefs,
     1    targ,pot)
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder)
      real *8 targ(2),pot(nd)

      real *8 px(norder),py(norder)
      
      real *8, allocatable :: tmp(:,:)

      allocate(tmp(nd,norder))
      
      x=targ(1)
      y=targ(2)
      
      if (ipoly .eq. 0) then
         call legepols(x,norder-1,px)
         call legepols(y,norder-1,py)
      elseif (ipoly .eq. 1) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
      endif
c     
      do j=1,norder
         do ind=1,nd
            pp=0
            do k=1,norder
               pp=pp+coefs(ind,k,j)*px(k)
            enddo
            tmp(ind,j)=pp
         enddo
      enddo
c
      do ind=1,nd
         pp=0
         do j=1,norder
            pp=pp+tmp(ind,j)*py(j)
         enddo
         pot(ind)=pp
      enddo

      return
      end
c
c
c
c
      subroutine ortho_evalt_3d(nd,ipoly,norder,coefs,
     1    targ,pot)
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder,norder)
      real *8 targ(3),pot(nd)

      real *8 px(norder),py(norder),pz(norder)
      
      real *8, allocatable :: tmp(:,:,:),tmp2(:,:)

      allocate(tmp(nd,norder,norder))
      allocate(tmp2(nd,norder))
      
      x=targ(1)
      y=targ(2)
      z=targ(3)
      
      if (ipoly .eq. 0) then
         call legepols(x,norder-1,px)
         call legepols(y,norder-1,py)
         call legepols(z,norder-1,pz)
      elseif (ipoly .eq. 1) then
         call chebpols(x,norder-1,px)
         call chebpols(y,norder-1,py)
         call chebpols(z,norder-1,pz)
      endif
c     
      do i=1,norder
      do j=1,norder
         do ind=1,nd
            pp=0
            do k=1,norder
               pp=pp+coefs(ind,k,j,i)*px(k)
            enddo
            tmp(ind,j,i)=pp
         enddo
      enddo
      enddo
c
      do i=1,norder
      do ind=1,nd
         pp=0
         do j=1,norder
            pp=pp+tmp(ind,j,i)*py(j)
         enddo
         tmp2(ind,i)=pp
      enddo
      enddo
c
      do ind=1,nd
         pp=0
         do j=1,norder
            pp=pp+tmp2(ind,j)*pz(j)
         enddo
         pot(ind)=pp
      enddo
      
      return
      end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate potential and gradient at arbitrary target 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_evalpgt_nd(ndim,nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      real *8 coefs(nd,norder**ndim)
      real *8 targ(ndim),pot(nd),grad(nd,ndim)

      if (ndim.eq.1) then
         call ortho_evalpgt_1d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
      elseif (ndim.eq.2) then
         call ortho_evalpgt_2d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
      elseif (ndim.eq.3) then
         call ortho_evalpgt_3d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
      endif
      
      return
      end
c
c
c
c
      subroutine ortho_evalpgt_1d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      real *8 coefs(nd,norder)
      real *8 targ,pot(nd),grad(nd)

      real *8 px(norder)
      real *8 pxp(norder)
      
      x=targ
      
      if (ipoly .eq. 0) then
         call legepolders(x,px,pxp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders(x,px,pxp,norder-1)
      endif
      do i=1,norder
         pxp(i)=pxp(i)*sc
      enddo
c     
      do ind=1,nd
         pp=0
         d1=0
         do k=1,norder
            pp=pp+coefs(ind,k)*px(k)
            d1=d1+coefs(ind,k)*pxp(k)
         enddo
         pot(ind)=pp
         grad(ind)=d1
      enddo

      return
      end
c
c
c
c
      subroutine ortho_evalpgt_2d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      
      if (ipoly .eq. 0) then
         call legepolders(x,px,pxp,norder-1)
         call legepolders(y,py,pyp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders(x,px,pxp,norder-1)
         call chebpolders(y,py,pyp,norder-1)
      endif
      do i=1,norder
         pxp(i)=pxp(i)*sc
         pyp(i)=pyp(i)*sc
      enddo
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
         grad(ind,1)=d1
         grad(ind,2)=d2
      enddo

      return
      end
c
c
c
c
      subroutine ortho_evalpgt_3d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      real *8 coefs(nd,norder,norder,norder)
      real *8 targ(3),pot(nd),grad(nd,3)

      real *8 px(norder),py(norder),pz(norder)
      real *8 pxp(norder),pyp(norder),pzp(norder)
      
      real *8, allocatable :: tmp(:,:,:,:),tmp2(:,:,:)

      allocate(tmp(nd,2,norder,norder))
      allocate(tmp2(nd,3,norder))
      
      x=targ(1)
      y=targ(2)
      z=targ(3)
      
      if (ipoly .eq. 0) then
         call legepolders(x,px,pxp,norder-1)
         call legepolders(y,py,pyp,norder-1)
         call legepolders(z,pz,pzp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders(x,px,pxp,norder-1)
         call chebpolders(y,py,pyp,norder-1)
         call chebpolders(z,pz,pzp,norder-1)
      endif

      do i=1,norder
         pxp(i)=pxp(i)*sc
         pyp(i)=pyp(i)*sc
         pzp(i)=pzp(i)*sc
      enddo
c     
      do i=1,norder
      do j=1,norder
         do ind=1,nd
            pp=0
            d1=0
            do k=1,norder
               pp=pp+coefs(ind,k,j,i)*px(k)
               d1=d1+coefs(ind,k,j,i)*pxp(k)
            enddo
            tmp(ind,1,j,i)=pp
            tmp(ind,2,j,i)=d1
         enddo
      enddo
      enddo
c
      do i=1,norder
      do ind=1,nd
         pp=0
         d1=0
         d2=0
         do j=1,norder
            pp=pp+tmp(ind,1,j,i)*py(j)
            d1=d1+tmp(ind,2,j,i)*py(j)
            d2=d2+tmp(ind,1,j,i)*pyp(j)
         enddo
         tmp2(ind,1,i)=pp
         tmp2(ind,2,i)=d1
         tmp2(ind,3,i)=d2
      enddo
      enddo
c
      do ind=1,nd
         pp=0
         d1=0
         d2=0
         d3=0
         do j=1,norder
            pp=pp+tmp2(ind,1,j)*pz(j)
            d1=d1+tmp2(ind,2,j)*pz(j)
            d2=d2+tmp2(ind,3,j)*pz(j)
            d3=d3+tmp2(ind,1,j)*pzp(j)
         enddo
         pot(ind)=pp
         grad(ind,1)=d1
         grad(ind,2)=d2
         grad(ind,3)=d3
      enddo
      
      return
      end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     evaluate potential, gradient, and hessian at arbitrary target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_evalpght_nd(ndim,nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
c
c     this subroutine evaluates the potential and gradient at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     targ - xy coordinates of the target point
c      
c     output:
c     pot - value of the potential
c     grad - value of the potential
c     hess - value of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder**ndim)
      real *8 targ(ndim),pot(nd),grad(nd,ndim),hess(nd,ndim*(ndim+1)/2)

      if (ndim.eq.1) then
         call ortho_evalpght_1d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
      elseif (ndim.eq.2) then
         call ortho_evalpght_2d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
      elseif (ndim.eq.3) then
         call ortho_evalpght_3d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
      endif
      
      return
      end
c
c
c
c
      subroutine ortho_evalpght_1d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
c
c     this subroutine evaluates the potential, gradient and hessian at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      real *8 coefs(nd,norder)
      real *8 targ,pot(nd),grad(nd),hess(nd)

      real *8 px(norder)
      real *8 pxp(norder)
      real *8 pxpp(norder)
      
      sc2=sc*sc
      
      x=targ
      
      if (ipoly .eq. 0) then
         call legepolders2(x,px,pxp,pxpp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders2(x,px,pxp,pxpp,norder-1)
      endif
      do i=1,norder
         pxp(i)=pxp(i)*sc
         pxpp(i)=pxpp(i)*sc2
      enddo
c     
      do ind=1,nd
         pp=0
         dd=0
         hh=0
         do k=1,norder
            pp=pp+coefs(ind,k)*px(k)
            dd=dd+coefs(ind,k)*pxp(k)
            hh=hh+coefs(ind,k)*pxpp(k)
         enddo
         pot(ind)=pp
         grad(ind)=dd
         hess(ind)=hh
      enddo

      return
      end
c
c
c
c
      subroutine ortho_evalpght_2d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
c
c     this subroutine evaluates the potential, gradient and hessian at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      
      if (ipoly .eq. 0) then
         call legepolders2(x,px,pxp,pxpp,norder-1)
         call legepolders2(y,py,pyp,pypp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders2(x,px,pxp,pxpp,norder-1)
         call chebpolders2(y,py,pyp,pypp,norder-1)
      endif
      do i=1,norder
         pxp(i)=pxp(i)*sc
         pyp(i)=pyp(i)*sc
         pxpp(i)=pxpp(i)*sc2
         pypp(i)=pypp(i)*sc2
      enddo
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
      subroutine ortho_evalpght_3d(nd,ipoly,norder,coefs,sc,
     1    targ,pot,grad,hess)
c
c     this subroutine evaluates the potential, gradient and hessian at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      real *8 coefs(nd,norder,norder,norder)
      real *8 targ(3),pot(nd),grad(nd,3),hess(nd,6)

      real *8 px(norder),py(norder),pz(norder)
      real *8 pxp(norder),pyp(norder),pzp(norder)
      real *8 pxpp(norder),pypp(norder),pzpp(norder)
      
      real *8, allocatable :: tmp(:,:,:,:),tmp2(:,:,:)

      allocate(tmp(nd,3,norder,norder))
      allocate(tmp2(nd,6,norder))

      sc2=sc*sc
      
      x=targ(1)
      y=targ(2)
      z=targ(3)
      
      if (ipoly .eq. 0) then
         call legepolders2(x,px,pxp,pxpp,norder-1)
         call legepolders2(y,py,pyp,pypp,norder-1)
         call legepolders2(z,pz,pzp,pzpp,norder-1)
      elseif (ipoly .eq. 1) then
         call chebpolders2(x,px,pxp,pxpp,norder-1)
         call chebpolders2(y,py,pyp,pypp,norder-1)
         call chebpolders2(z,pz,pzp,pzpp,norder-1)
      endif
      do i=1,norder
         pxp(i)=pxp(i)*sc
         pyp(i)=pyp(i)*sc
         pzp(i)=pzp(i)*sc
         pxpp(i)=pxpp(i)*sc2
         pypp(i)=pypp(i)*sc2
         pzpp(i)=pzpp(i)*sc2
      enddo
c     
      do i=1,norder
      do j=1,norder
         do ind=1,nd
            pp=0
            dx=0
            dxx=0
            do k=1,norder
               pp=pp  + coefs(ind,k,j,i)*px(k)
               dx=dx  + coefs(ind,k,j,i)*pxp(k)
               dxx=dxx+ coefs(ind,k,j,i)*pxpp(k)
            enddo
            tmp(ind,1,j,i)=pp
            tmp(ind,2,j,i)=dx
            tmp(ind,3,j,i)=dxx
         enddo
      enddo
      enddo
c
      do i=1,norder
      do ind=1,nd
         pp=0
         dy=0
         dyy=0
         dx=0
         dxy=0
         dxx=0
         do j=1,norder
            pp=pp+tmp(ind,1,j,i)*py(j)
            dy=dy+tmp(ind,1,j,i)*pyp(j)
            dyy=dyy+tmp(ind,1,j,i)*pypp(j)

            dx=dx+tmp(ind,2,j,i)*py(j)
            dxy=dxy+tmp(ind,2,j,i)*pyp(j)

            dxx=dxx+tmp(ind,3,j,i)*py(j)
         enddo
         tmp2(ind,1,i)=pp
         tmp2(ind,2,i)=dx
         tmp2(ind,3,i)=dy
         tmp2(ind,4,i)=dxx
         tmp2(ind,5,i)=dxy
         tmp2(ind,6,i)=dyy
      enddo
      enddo

c
      do ind=1,nd
         pp=0
         dx=0
         dy=0
         dz=0
         dxx=0
         dyy=0
         dzz=0
         dxy=0
         dxz=0
         dyz=0
         do j=1,norder
            pp  = pp  + tmp2(ind,1,j)*pz(j)
            dx  = dx  + tmp2(ind,2,j)*pz(j)
            dy  = dy  + tmp2(ind,3,j)*pz(j)

            dxx = dxx + tmp2(ind,4,j)*pz(j)
            dxy = dxy + tmp2(ind,5,j)*pz(j)
            dyy = dyy + tmp2(ind,6,j)*pz(j)

            dz  = dz  + tmp2(ind,1,j)*pzp(j)
            dxz = dxz + tmp2(ind,2,j)*pzp(j)
            dyz = dyz + tmp2(ind,3,j)*pzp(j)

            dzz = dzz + tmp2(ind,1,j)*pzpp(j)
         enddo
         pot(ind)=pp
         grad(ind,1)=dx
         grad(ind,2)=dy
         grad(ind,3)=dz
         
         hess(ind,1)=dxx
         hess(ind,2)=dyy
         hess(ind,3)=dzz
         
         hess(ind,4)=dxy
         hess(ind,5)=dxz
         hess(ind,6)=dyz
      enddo
      
      return
      end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     make 1D tables for computing function values and its derivatives
c     from its expansion coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ortho_eval_tables(ipoly,n,umat,vmat,vpmat,vppmat)
c
c     returns 1d transformation tables 
c
c     input:
c     ipoly - ipoly=0 Legendre polynomials
c                   1 Chebyshev polynomials
c     n - order of polynomial expansion
c      
c     output:
c     umat - nxn matrix converting function values to expansion coefficients 
c     vmat - nxn matrix converting expansion coefficients to function values
c     vpmat - nxn matrix converting expansion coefficients to first derivatives
c     vppmat - nxn matrix converting expansion coefficients to second derivatives
c
      implicit real *8 (a-h,o-z)
      real *8 umat(n,n)
      real *8 vmat(n,n)
      real *8 vpmat(n,n)
      real *8 vppmat(n,n)
      real *8 xs(n),ws(n)

      itype=2
      if (ipoly.eq.0) then
         call legeexps2(itype,n,xs,umat,vmat,ws,vpmat,vppmat)
      elseif (ipoly.eq.1) then
         call chebexps2(itype,n,xs,umat,vmat,ws,vpmat,vppmat)
      endif

      return
      end
c
c
c
c
      
