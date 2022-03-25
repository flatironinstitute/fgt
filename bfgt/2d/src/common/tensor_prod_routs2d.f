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
c     of tensor product expansions in two dimensions.
c     Following is a brief description of these
c     subroutines.
c
      subroutine orth_trans2d(nd,itype,norder,fin,fout,uxmat,uymat)
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
      real *8 uxmat(norder,norder),uymat(norder,norder)
      real *8, allocatable :: work(:,:,:)
      
      allocate(work(nd,norder,norder))

      if (itype.eq.0) then
c        transform in x
         do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+uxmat(k,k1)*fin(ind,k1,j)
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
                  dd=dd+uymat(i,i1)*work(ind,k,i1)
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
                  dd=dd+uxmat(k,k1)*fin(ind,k1,j)
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
                  dd=dd+uymat(k,k1)*fin(ind,j,k1)
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
      subroutine orth_evalg2d(nd,norder,coefs,sc,
     1    grad,vmat,vpmat)
c
c     this subroutine evaluates the gradient at the tensor grid
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     vmat - matrix converting coefficient into function values
c     vpmat - matrix converting coefficient into values of the first derivative
c      
c     output:
c     grad - value of the gradient
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder)
      real *8 grad(nd,2,norder,norder)

      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)

      real *8, allocatable :: gradx(:,:,:),grady(:,:,:)

      allocate(gradx(nd,norder,norder))
      allocate(grady(nd,norder,norder))

      itype=0
      call orth_trans2d(nd,itype,norder,coefs,gradx,vpmat,vmat)
      call orth_trans2d(nd,itype,norder,coefs,grady,vmat,vpmat)

      do i=1,norder
         do j=1,norder
            do ind=1,nd
               grad(ind,1,j,i)=gradx(ind,j,i)*sc
               grad(ind,2,j,i)=grady(ind,j,i)*sc
            enddo
         enddo
      enddo
      
      return
      end
c
c
c
c
      subroutine orth_evalh2d(nd,norder,coefs,sc,
     1    hess,vmat,vpmat,vppmat)
c
c     this subroutine evaluates the hessian at the tensor grid
c
c     input:
c     nd - number of input data
c     norder - order of polynomial expansion
c     coefs - orthogonal polynomial expansion coefficients for the potential
c     sc - scaling factor for the derivative
c     vmat - matrix converting coefficient into function values
c     vpmat - matrix converting coefficient into values of the first derivative
c     vppmat - matrix converting coefficient into values of the second derivative
c
c     output:
c     hess - value of the hessian
c
      implicit real *8 (a-h,o-z)
      real *8 coefs(nd,norder,norder)
      real *8 hess(nd,3,norder,norder)

      real *8 vmat(norder,norder)
      real *8 vpmat(norder,norder)
      real *8 vppmat(norder,norder)

      real *8, allocatable :: hess1(:,:,:),hess2(:,:,:)
      real *8, allocatable :: hess3(:,:,:)

      allocate(hess1(nd,norder,norder))
      allocate(hess2(nd,norder,norder))
      allocate(hess3(nd,norder,norder))

      itype=0
      call orth_trans2d(nd,itype,norder,coefs,hess1,vppmat,vmat)
      call orth_trans2d(nd,itype,norder,coefs,hess2,vpmat,vpmat)
      call orth_trans2d(nd,itype,norder,coefs,hess3,vmat,vppmat)

      sc2=sc*sc
      do i=1,norder
         do j=1,norder
            do ind=1,nd
               hess(ind,1,j,i)=hess1(ind,j,i)*sc2
               hess(ind,2,j,i)=hess2(ind,j,i)*sc2
               hess(ind,3,j,i)=hess3(ind,j,i)*sc2
            enddo
         enddo
      enddo
      
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
      subroutine orth_evalt2d(nd,ipoly,norder,fcoefs,targ,fval)
c
c     this subroutine evaluates the orthogonal polynomial expansion at the given target point
c
c     input:
c     nd - number of input data
c     ipoly - itype=0 Legendre polynomials
c                   1 Chebyshev polynomials
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
      subroutine orth_evalpgt2d(nd,ipoly,norder,coefs,sc,
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
      subroutine orth_evalpght2d(nd,ipoly,norder,coefs,sc,
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
