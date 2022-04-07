c     This file contains a set of subroutines that carry out transforms
c     for 2D box FGT.
c
c     legtrans2d: converts function values on a Legendre tensor product grid
c                 to Legendre expansion coefficients.
c
c     leg2d_to_pw: converts 2D Legendre expansion coefficients to 2D planewave
c                 expansion coefficients.
c
c     g2d_pw2pot: converts 2D planewave expansion coefficients to potential values
c                 on a Legendre tensor product grid in the same box.
c
c     leg2d_to_potloc: converts 2D Legendre expansion coefficients in the source box 
c                 to potential values on a Legendre tensor product grid in 
c                 a target box that is in the list1 of the source box.
c
c     pw_translation_matrices: precomputes all translation matrices for 2D PW
C                 expansions for mp to loc at the cutoff level.
c
c     merge_split_pw_matrices: precomputes all translation matrices for 
c                 2D PW translations from child to parent or vice versa at all 
c                 needed levels.
c
c     g2dshiftpw_vec: translates a 2D PW expansion (pwexp1) about
C                 the center (CENT1) into another PW expansion (pwexp2) about 
C                 (CENT2) using precomputed translation matrix wshift.
c
c     g2dcopypwexp_vec: adds one 2D PW expansion (pwexp1) 
C                 to another 2D PW expansion (pwexp2).
c
c     g2dpwzero_vec: sets a vector 2D planewave expansions to zero.
c
c
C*********************************************************************C
c
c     convert function values on a tensor product grid in each leaf
c     box to the plane wave expansion
c
C*********************************************************************C
      subroutine gnd_tens_prod_to_pw(ndim,nd,n,fvals,npw,
     1    tab_poly2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     tenosr product function values in general n dimensions
c
c     INPUT:
c     ndim         dimension of underlying space
c     nd           vector length (for multiple RHS)
c     n            number of points on tensor grid along each dimension
c     fvals        function values at a tensor product grid
c     npw          number of plane waves along each dimension
c                  NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_poly2pw  precomputed table of 1D conversion factors
c                  (n,j) entry is: ws(j)*(D/2)* 
c                  int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                  where D is box dimension at current level in
c                  tree hierarchy.
c     OUTPUT:
c     pwexp        plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n**ndim)
      complex *16 tab_poly2pw(n,npw)
      complex *16 pwexp(npw**ndim/2,nd)

      if (ndim.eq.1) then
         call g1d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
      elseif (ndim.eq.2) then
         call g2d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
      elseif (ndim.eq.3) then
         call g3d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
      endif
      return
      end subroutine
c
c
C*********************************************************************C
      subroutine g1d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     function values on a tensor product grid (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     fvals    function values on tensor product grid
c     npw      number of plane waves
c                 NOTE 1D convention is pwexp(npw/2)
c     tab_poly2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n)
      complex *16 tab_poly2pw(n,npw)
      complex *16 pwexp(npw/2,nd),cd
c
      do ind = 1,nd
         do k1 = 1,npw/2
            cd = 0.0d0
            do m1 = 1,n
               cd = cd+tab_poly2pw(m1,k1)*fvals(ind,m1)
            enddo
            pwexp(k1,ind) = cd
         enddo
      enddo
      
      return
      end subroutine
c
c
c
c
c
c
c
C*********************************************************************C
      subroutine g2d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     function values on a tensor product grid (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     fvals    function values on tensor grid
c     npw      number of plane waves
c                 NOTE 2D convention is pwexp(npw,npw,npw/2)
c     tab_poly2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n)
      complex *16 tab_poly2pw(n,npw)
      complex *16 pwexp(npw,npw/2,nd)
      complex *16 ff(n,npw/2),cd
c
      do ind = 1,nd
         do m1 = 1,n
            do k2 = 1,npw/2
               cd = 0.0d0
               do m2 = 1,n
                  cd = cd+tab_poly2pw(m2,k2)*fvals(ind,m1,m2)
               enddo
               ff(m1,k2) = cd
            enddo
         enddo
c
         do k2 = 1,npw/2
            do k1 = 1,npw
               cd = 0.0d0
               do m1 = 1,n
                  cd = cd+tab_poly2pw(m1,k1)*ff(m1,k2)
               enddo
               pwexp(k1,k2,ind) = cd
            enddo
         enddo
c
      enddo
      
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine g3d_tens_prod_to_pw(nd,n,fvals,npw,tab_poly2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     tenosr product function values.
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     fvals    function values at Legendre tensor product grid
c
c     npw      number of plane waves
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_poly2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n)
      complex *16 tab_poly2pw(n,npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      complex *16 ff(n,n,npw/2)
      complex *16 ff2(n,npw,npw/2)
      complex *16 cd
c
      do ind = 1,nd
         do k3 = 1,npw/2
         do m2 = 1,n
         do m1 = 1,n
            cd = 0.0d0
            do m3 = 1,n
               cd = cd+tab_poly2pw(m3,k3)*fvals(ind,m1,m2,m3)
            enddo
            ff(m1,m2,k3) = cd
         enddo
         enddo
         enddo
c
         do k3 = 1,npw/2
         do k2 = 1,npw
         do m1 = 1,n
            cd = 0.0d0
            do m2 = 1,n
               cd = cd+tab_poly2pw(m2,k2)*ff(m1,m2,k3)
            enddo
            ff2(m1,k2,k3) = cd
         enddo
         enddo
         enddo
c
         do k3 = 1,npw/2
         do k2 = 1,npw
         do k1 = 1,npw
            cd = 0.0d0
            do m1 = 1,n
               cd = cd+tab_poly2pw(m1,k1)*ff2(m1,k2,k3)
            enddo
            pwexp(k1,k2,k3,ind) = cd
         enddo
         enddo
         enddo
c
      enddo
      
      return
      end subroutine
c
c
c
c
c
c*******************************************************************************
c
c     Plane-wave expansion to potential subroutines
c     
C*********************************************************************C
      subroutine gnd_pw2pgh(ndim,nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,tab_pw2dxx,ifpgh,pot,grad,hess)
C*********************************************************************C
c     This routine computes the potential, gradient, and hessian on a tensor 
c     grid from the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     ndim     dimension of underlying space
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane wave expansion along each dimension
c     pwexp    plane wave expansion
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw**ndim/2,nd)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 tab_pw2dxx(npw,n)
      real *8 pot(nd,n**ndim)
      real *8 grad(nd,ndim,n**ndim)
      real *8 hess(nd,ndim*(ndim+1)/2,n**ndim)
c
      if (ifpgh.eq.1) then
        if (ndim.eq.1) then
          call g1d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
        elseif (ndim.eq.2) then
          call g2d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
        elseif (ndim.eq.3) then
          call g3d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
        endif
      elseif (ifpgh.eq.2) then
        if (ndim.eq.1) then
          call g1d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,pot,grad)
        elseif (ndim.eq.2) then
          call g2d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,pot,grad)
        elseif (ndim.eq.3) then
          call g3d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,pot,grad)
        endif
      elseif (ifpgh.eq.3) then
        if (ndim.eq.1) then
          call g1d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,
     1         tab_pw2dxx,pot,grad,hess)
        elseif (ndim.eq.2) then
          call g2d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,
     1         tab_pw2dxx,pot,grad,hess)
        elseif (ndim.eq.3) then
          call g3d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,tab_pw2der,
     1         tab_pw2dxx,pot,grad,hess)
        endif
      endif
      
      return
      end subroutine
c
c
c
c************************************************************************
c
c     1d pw2pot,pg,pgh
c
c*************************************************************************
      subroutine g1d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
C*********************************************************************C
c     This routine computes the potential on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 1D convention is pwexp(npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 pwexp(npw/2,nd),cd
c
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            do m1 = 1,npw/2
               cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,ind)
            enddo
            pot(ind,k1) = pot(ind,k1)+dreal(cd)*2
         enddo
      enddo
c
      return
      end subroutine
c      
c      
c      
cC*********************************************************************C
      subroutine g1d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,pot,grad)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 2D convention is pwexp(npw,(npw+1)/2)
c     tab_pw2pot    precomputed table of 1D conversion 
c     tab_pw2der  precomputed table of deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n)
      real *8 grad(nd,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 pwexp(npw/2,nd),cd,cdx
c
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            do m1 = 1,npw/2
               cd  = cd  + tab_pw2pot(m1,k1)*pwexp(m1,ind)
               cdx = cdx + tab_pw2der(m1,k1)*pwexp(m1,ind)
            enddo
            pot(ind,k1)=pot(ind,k1)+dreal(cd)*2
            grad(ind,k1)=grad(ind,k1)+dreal(cdx)*2
         enddo
      enddo
c
      return
      end subroutine
c
c
c
c      
c      
c
c
C*********************************************************************C
      subroutine g1d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,tab_pw2dxx,pot,grad,hess)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c     tab_pw2pot    precomputed table of 1D conversion 
c     tab_pw2der  precomputed table of deriv of 1D conversion 
c     tab_pw2dxx  precomputed table of second deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c     hess     hess values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n)
      real *8 grad(nd,n)
      real *8 hess(nd,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 tab_pw2dxx(npw,n)
      complex *16 pwexp(npw/2,nd),cd,cdx,cdxx
c
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            cdxx=0.0d0
            do m1 = 1,npw/2
               cd   = cd   + tab_pw2pot(m1,k1)*pwexp(m1,ind)
               cdx  = cdx  + tab_pw2der(m1,k1)*pwexp(m1,ind)
               cdxx = cdxx + tab_pw2dxx(m1,k1)*pwexp(m1,ind)
            enddo
            pot(ind,k1)  = pot(ind,k1)  + dreal(cd)*2
            grad(ind,k1) = grad(ind,k1) + dreal(cdx)*2
            hess(ind,k1) = hess(ind,k1) + dreal(cdxx)*2
         enddo
      enddo
c
      return
      end subroutine
c
c
c
c
c
c************************************************************************
c
c     2d pw2pot,pg,pgh
c
C*********************************************************************C
      subroutine g2d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
C*********************************************************************C
c     This routine computes the potential on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 2D convention is pwexp(npw,npw/2)
c     ff       complex workspace (n,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 pwexp(npw,npw/2,nd),cd
      complex *16 ff(n,npw/2)
c
      do ind = 1,nd
         do m2 = 1,npw/2
         do k1 = 1,n
            cd=0.0d0
            do m1 = 1,npw
               cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,m2,ind)
            enddo
            ff(k1,m2) = cd
         enddo
         enddo
c
         do k2 = 1,n
         do k1 = 1,n
            cd = 0.0d0
            do m2 = 1,npw/2
               cd = cd+tab_pw2pot(m2,k2)*ff(k1,m2)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+dreal(cd)*2
         enddo
         enddo
c
      enddo

      return
      end subroutine
c      
c      
c      
c      
cC*********************************************************************C
      subroutine g2d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,pot,grad)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 2D convention is pwexp(npw,(npw+1)/2)
c     tab_pw2pot    precomputed table of 1D conversion 
c     tab_pw2der  precomputed table of deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n)
      real *8 grad(nd,2,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 ff(npw/2,n)
      complex *16 ffx(npw/2,n)
      complex *16 pwexp(npw,npw/2,nd),cd,cdx,cdy
c
      do ind = 1,nd
         do m2 = 1,npw/2
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            do m1 = 1,npw
               cd  = cd  + tab_pw2pot(m1,k1)*pwexp(m1,m2,ind)
               cdx = cdx + tab_pw2der(m1,k1)*pwexp(m1,m2,ind)
            enddo
            ff(m2,k1) = cd
            ffx(m2,k1) = cdx
         enddo
         enddo
c
         do k2 = 1,n
         do k1 = 1,n
            cd = 0.0d0
            cdx = 0.0d0
            cdy = 0.0d0
            do m2 = 1,npw/2
               cd  = cd  + tab_pw2pot(m2,k2)*ff(m2,k1)
               cdy = cdy + tab_pw2der(m2,k2)*ff(m2,k1)
               
               cdx = cdx + tab_pw2pot(m2,k2)*ffx(m2,k1)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+dreal(cd)*2
            grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+dreal(cdx)*2
            grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+dreal(cdy)*2
         enddo
         enddo
c
      enddo
c
      return
      end subroutine
c
c
c
c      
c      
c
c
C*********************************************************************C
      subroutine g2d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,tab_pw2dxx,pot,grad,hess)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd           vector length (for multiple RHS)
c     n            number of nodes
c     npw          number of plane waves
c     pwexp        plane wave expansion
c                  NOTE 2D convention is pwexp(npw,(npw+1)/2)
c     tab_pw2pot   precomputed table of 1D conversion 
c     tab_pw2der   precomputed table of deriv of 1D conversion 
c     tab_pw2dxx   precomputed table of second deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c     hess     hess values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n)
      real *8 grad(nd,2,n,n)
      real *8 hess(nd,3,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 tab_pw2dxx(npw,n)
      complex *16 ff(npw/2,n)
      complex *16 ffx(npw/2,n)
      complex *16 ffxx(npw/2,n)
      complex *16 pwexp(npw,npw/2,nd),cd,cdx,cdy
      complex *16 cdxx,cdxy,cdyy
c
      do ind = 1,nd
         do m2 = 1,npw/2
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            cdxx=0.0d0
            do m1 = 1,npw
               cd   = cd   + tab_pw2pot(m1,k1)*pwexp(m1,m2,ind)
               cdx  = cdx  + tab_pw2der(m1,k1)*pwexp(m1,m2,ind)
               cdxx = cdxx + tab_pw2dxx(m1,k1)*pwexp(m1,m2,ind)
            enddo
            ff(m2,k1) = cd
            ffx(m2,k1) = cdx
            ffxx(m2,k1) = cdxx
         enddo
         enddo
c
         do k2 = 1,n
         do k1 = 1,n
            cd = 0.0d0
            cdx = 0.0d0
            cdy = 0.0d0
            cdxx = 0.0d0
            cdxy = 0.0d0
            cdyy = 0.0d0
            do m2 = 1,npw/2
               cd   = cd   + tab_pw2pot(m2,k2)*ff(m2,k1)
               cdy  = cdy  + tab_pw2der(m2,k2)*ff(m2,k1)
               cdyy = cdyy + tab_pw2dxx(m2,k2)*ff(m2,k1)

               cdx  = cdx  + tab_pw2pot(m2,k2)*ffx(m2,k1)
               cdxy = cdxy + tab_pw2der(m2,k2)*ffx(m2,k1)

               cdxx = cdxx + tab_pw2pot(m2,k2)*ffxx(m2,k1)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+dreal(cd)*2

            grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+dreal(cdx)*2
            grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+dreal(cdy)*2

            hess(ind,1,k1,k2)=hess(ind,1,k1,k2)+dreal(cdxx)*2
            hess(ind,2,k1,k2)=hess(ind,2,k1,k2)+dreal(cdxy)*2
            hess(ind,3,k1,k2)=hess(ind,3,k1,k2)+dreal(cdyy)*2
         enddo
         enddo
c
      enddo
c
      return
      end subroutine
c
c
c
c
c
c************************************************************************
c
c     3d pw2pot,pg,pgh
c
C*********************************************************************C
      subroutine g3d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
C*********************************************************************C
c     This routine computes the potential on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw,npw,npw/2,nd)
      complex *16 tab_pw2pot(npw,n)
      real *8 pot(nd,n,n,n)

      complex *16 ff(n,npw,npw/2)
      complex *16 ff2(n,n,npw/2)
      complex *16 cd
c
      do ind = 1,nd
         do m3 = 1,npw/2
         do m2 = 1,npw
         do k1 = 1,n
            cd=0.0d0
            do m1 = 1,npw
               cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,m2,m3,ind)
            enddo
            ff(k1,m2,m3) = cd
         enddo
         enddo
         enddo
c
         do m3 = 1,npw/2
         do k2 = 1,n
         do k1 = 1,n
            cd = 0.0d0
            do m2 = 1,npw
               cd = cd+tab_pw2pot(m2,k2)*ff(k1,m2,m3)
            enddo
            ff2(k1,k2,m3) = cd
         enddo
         enddo
         enddo
c
         do k3 = 1,n
         do k2 = 1,n
         do k1 = 1,n
            cd = 0.0d0
            do m3 = 1,npw/2
               cd = cd+tab_pw2pot(m3,k3)*ff2(k1,k2,m3)
            enddo
            pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+dreal(cd)*2
         enddo
         enddo
         enddo
c
      enddo
      
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine g3d_pw2pg(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,pot,grad)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c     tab_pw2deriv  precomputed table of deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw,npw,npw/2,nd)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)

      complex *16 ff(n,npw,npw/2)
      complex *16 ffx(n,npw,npw/2)
      
      complex *16 ff2(n,n,npw/2)
      complex *16 ff2x(n,n,npw/2)
      complex *16 ff2y(n,n,npw/2)
      complex *16 cd,cdx,cdy,cdz
c
      do ind = 1,nd
        do m3 = 1,npw/2
        do m2 = 1,npw
        do k1 = 1,n
           cd=0.0d0
           cdx=0.0d0
           do m1 = 1,npw
              cd  = cd  + tab_pw2pot(m1,k1)*pwexp(m1,m2,m3,ind)
              cdx = cdx + tab_pw2der(m1,k1)*pwexp(m1,m2,m3,ind)
           enddo
           ff(k1,m2,m3) = cd
           ffx(k1,m2,m3) = cdx
        enddo
        enddo
        enddo
c
        do m3 = 1,npw/2
        do k2 = 1,n
        do k1 = 1,n
           cd = 0.0d0
           cdx = 0.0d0
           cdy = 0.0d0
           do m2 = 1,npw
              cd = cd   + tab_pw2pot(m2,k2)*ff(k1,m2,m3)
              cdy = cdy + tab_pw2der(m2,k2)*ff(k1,m2,m3)

              cdx = cdx + tab_pw2pot(m2,k2)*ffx(k1,m2,m3)
           enddo
           ff2(k1,k2,m3) = cd
           ff2x(k1,k2,m3) = cdx
           ff2y(k1,k2,m3) = cdy
        enddo
        enddo
        enddo
c
        do k3 = 1,n
        do k2 = 1,n
        do k1 = 1,n
           cd = 0.0d0
           cdx = 0.0d0
           cdy = 0.0d0
           cdz = 0.0d0
           do m3 = 1,npw/2
              cd  = cd  + tab_pw2pot(m3,k3)*ff2(k1,k2,m3)
              cdz = cdz + tab_pw2der(m3,k3)*ff2(k1,k2,m3)
              
              cdx = cdx + tab_pw2pot(m3,k3)*ff2x(k1,k2,m3)
              cdy = cdy + tab_pw2pot(m3,k3)*ff2y(k1,k2,m3)
           enddo
           pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+dreal(cd)*2
           grad(ind,1,k1,k2,k3)=grad(ind,1,k1,k2,k3)+dreal(cdx)*2
           grad(ind,2,k1,k2,k3)=grad(ind,2,k1,k2,k3)+dreal(cdy)*2
           grad(ind,3,k1,k2,k3)=grad(ind,3,k1,k2,k3)+dreal(cdz)*2
        enddo
        enddo
        enddo
c
      enddo
      
      return
      end subroutine
c
c
c
c
C*********************************************************************C
      subroutine g3d_pw2pgh(nd,n,npw,pwexp,tab_pw2pot,
     1           tab_pw2der,tab_pw2dxx,pot,grad,hess)
C*********************************************************************C
c     This routine computes the potential, gradient and hessian on a 
c     tensor grid from the plane wave expansion coefficients
c     (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c     tab_pw2der  precomputed table of first deriv of 1D conversion 
c     tab_pw2dxx  precomputed table of second deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c     hess     hess values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw,npw,npw/2,nd)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 tab_pw2dxx(npw,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 hess(nd,6,n,n,n)

      complex *16 ff(n,npw,npw/2)
      complex *16 ffx(n,npw,npw/2)
      complex *16 ffxx(n,npw,npw/2)
      
      complex *16 ff2(n,n,npw/2)
      complex *16 ff2x(n,n,npw/2)
      complex *16 ff2xy(n,n,npw/2)
      complex *16 ff2xx(n,n,npw/2)
      complex *16 ff2y(n,n,npw/2)
      complex *16 ff2yy(n,n,npw/2)
      
      complex *16 cd,cdx,cdy,cdz,cdxx,cdyy,cdzz,cdxy,cdxz,cdyz
c
      do ind = 1,nd
c       transformation in x
        do m3 = 1,npw/2
        do m2 = 1,npw
        do k1 = 1,n
           cd=0.0d0
           cdx=0.0d0
           cdxx=0.0d0
           do m1 = 1,npw
              cd =   cd   + tab_pw2pot(m1,k1)*pwexp(m1,m2,m3,ind)
              cdx  = cdx  + tab_pw2der(m1,k1)*pwexp(m1,m2,m3,ind)
              cdxx = cdxx + tab_pw2dxx(m1,k1)*pwexp(m1,m2,m3,ind)
           enddo
           ff(k1,m2,m3) = cd
           ffx(k1,m2,m3) = cdx
           ffxx(k1,m2,m3) = cdxx
        enddo
        enddo
        enddo
c       transformation in y
        do m3 = 1,npw/2
        do k2 = 1,n
        do k1 = 1,n
           cd = 0.0d0
           cdx = 0.0d0
           cdy = 0.0d0
           cdxx = 0.0d0
           cdxy = 0.0d0
           cdyy = 0.0d0
           do m2 = 1,npw
              cd   = cd   + tab_pw2pot(m2,k2)*ff(k1,m2,m3)
              cdy  = cdy  + tab_pw2der(m2,k2)*ff(k1,m2,m3)
              cdyy = cdyy + tab_pw2dxx(m2,k2)*ff(k1,m2,m3)

              cdx  = cdx  + tab_pw2pot(m2,k2)*ffx(k1,m2,m3)
              cdxy = cdxy + tab_pw2der(m2,k2)*ffx(k1,m2,m3)

              cdxx = cdxx  + tab_pw2pot(m2,k2)*ffxx(k1,m2,m3)
           enddo
           ff2(k1,k2,m3) = cd
           ff2x(k1,k2,m3) = cdx
           ff2xx(k1,k2,m3) = cdxx
           ff2y(k1,k2,m3) = cdy
           ff2yy(k1,k2,m3) = cdyy
           ff2xy(k1,k2,m3) = cdxy
        enddo
        enddo
        enddo
c       transformation in z
        do k3 = 1,n
        do k2 = 1,n
        do k1 = 1,n
           cd = 0.0d0
           cdx = 0.0d0
           cdy = 0.0d0
           cdz = 0.0d0
           cdxx = 0.0d0
           cdyy = 0.0d0
           cdzz = 0.0d0
           cdxy = 0.0d0
           cdxz = 0.0d0
           cdyz = 0.0d0
           do m3 = 1,npw/2
              cd   = cd  +  tab_pw2pot(m3,k3)*ff2(k1,k2,m3)
              cdz  = cdz +  tab_pw2der(m3,k3)*ff2(k1,k2,m3)
              cdzz = cdzz + tab_pw2dxx(m3,k3)*ff2(k1,k2,m3)

              cdx  = cdx  + tab_pw2pot(m3,k3)*ff2x(k1,k2,m3)
              cdxz = cdxz + tab_pw2der(m3,k3)*ff2x(k1,k2,m3)

              cdy  = cdy  + tab_pw2pot(m3,k3)*ff2y(k1,k2,m3)
              cdyz = cdyz + tab_pw2der(m3,k3)*ff2y(k1,k2,m3)

              cdxx = cdxx + tab_pw2pot(m3,k3)*ff2xx(k1,k2,m3)
              cdxy = cdxy + tab_pw2pot(m3,k3)*ff2xy(k1,k2,m3)
              cdyy = cdyy + tab_pw2pot(m3,k3)*ff2yy(k1,k2,m3)
           enddo
           pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+dreal(cd)*2
           
           grad(ind,1,k1,k2,k3)=grad(ind,1,k1,k2,k3)+dreal(cdx)*2
           grad(ind,2,k1,k2,k3)=grad(ind,2,k1,k2,k3)+dreal(cdy)*2
           grad(ind,3,k1,k2,k3)=grad(ind,3,k1,k2,k3)+dreal(cdz)*2

           hess(ind,1,k1,k2,k3)=hess(ind,1,k1,k2,k3)+dreal(cdxx)*2
           hess(ind,2,k1,k2,k3)=hess(ind,2,k1,k2,k3)+dreal(cdyy)*2
           hess(ind,3,k1,k2,k3)=hess(ind,3,k1,k2,k3)+dreal(cdzz)*2

           hess(ind,4,k1,k2,k3)=hess(ind,4,k1,k2,k3)+dreal(cdxy)*2
           hess(ind,5,k1,k2,k3)=hess(ind,5,k1,k2,k3)+dreal(cdxz)*2
           hess(ind,6,k1,k2,k3)=hess(ind,6,k1,k2,k3)+dreal(cdyz)*2
        enddo
        enddo
        enddo
c
      enddo
      
      return
      end subroutine
c
c
c
c
c*******************************************************************************
c
c     local direct interactions by precomputed 1D tables
c     (1) without a uniform nd user interface
c     (2) without speedup for small delta
c
C*********************************************************************C
      subroutine g1d_tens_prod_to_potloc_old(nd,n,fvals,pot,tabx)
C*********************************************************************C
c     This routine computes the volume Gauss transform over a 
c     single box source distribution given as a Legendre series.
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays (tabx, taby). 
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in 1D, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x or y in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, or y.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, or y.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       Legendre coefficients
c                 f = sum coeff(n) P_n(x)
c     tabx    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in x.
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n),pot(nd,n)
      real *8 tabx(n,n)
c
      do ind = 1,nd
c        transform in x
         do k1=1,n
            cd=0
            do j1=1,n
               cd=cd+tabx(j1,k1)*fvals(ind,j1)
            enddo
            pot(ind,k1)=pot(ind,k1)+cd
         enddo
c     end of the ind loop
      enddo

      return
      end subroutine
c
c
C
c
C
C*********************************************************************C
      subroutine g2d_tens_prod_to_potloc_old(nd,n,fvals,pot,tabx,taby)
C*********************************************************************C
c     This routine computes the volume Gauss transform over a 
c     single box source distribution given as a Legendre series.
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays (tabx, taby). 
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in 2D, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x or y in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, or y.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, or y.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     fvals       function values at tensor grids
c     tabx    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in x.
c     taby    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in y.
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n),pot(nd,n,n)
      real *8 ff(n,n),tabx(n,n),taby(n,n)
c
      do ind = 1,nd
c        transform in x
         do j2=1,n
            do k1=1,n
               cd=0
               do j1=1,n
                  cd=cd+tabx(j1,k1)*fvals(ind,j1,j2)
               enddo
               ff(k1,j2)=cd
            enddo
         enddo
c        transfrom in y
         do k2=1,n
            do k1=1,n
               cd=0
               do j2=1,n
                  cd=cd+taby(j2,k2)*ff(k1,j2)
               enddo
               pot(ind,k1,k2)=pot(ind,k1,k2)+cd
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
      subroutine g3d_tens_prod_to_potloc_old(nd,n,fvals,pot,
     1    tabx,taby,tabz)
C*********************************************************************C
c     This routine computes the volume Gauss transform over a 
c     single box source distribution given as a Legendre series.
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays (tabx, taby,tabz). 
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in 3D, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x, y, or z in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, y, or z.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, y, or z.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     fvals       functions values at Legendre tensor product grid
c                 
c     ff          workspace
c     ff2          workspace
c     tabx    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in x.
c     taby    precomputed table of 1D integrals
c                 int_{Source box} P_n(x) exp( (\xi_j -x)^2/delta)
c                 for targets at current level in tree hierarchy with
c                 desired offset in y.
c     tabz    precomputed table of 1D integrals
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n),pot(nd,n,n,n)
      real *8 ff(n,n,n),ff2(n,n,n),tabx(n,n),taby(n,n),tabz(n,n)
c
      do ind = 1,nd
c        transform in x
         do j3=1,n
            do j2=1,n
               do k1=1,n
                  cd=0
                  do j1=1,n
                     cd=cd+tabx(j1,k1)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(k1,j2,j3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do j3=1,n
            do k2=1,n            
               do k1=1,n
                  cd=0
                  do j2=1,n
                     cd=cd+taby(j2,k2)*ff(k1,j2,j3)
                  enddo
                  ff2(k1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transfrom in z
         do k3=1,n
            do k2=1,n
               do k1=1,n
                  cd=0
                  do j3=1,n
                     cd=cd+tabz(j3,k3)*ff2(k1,k2,j3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo

      return
      end subroutine
c
c
C**************************************************************************
c
C     local direction interactions by precomputed 1D tables
C
C*********************************************************************C
      subroutine gnd_tens_prod_to_pghloc(ndim,nd,n,fvals,
     1    ifpgh,pot,grad,hess,
     2    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes the volume Gauss transform in general n dimensions
c     over a single box source distribution given as function values on a
c     tensor product grid.
c      
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays tab_loc.
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in higher dimensions, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x or y in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, or y.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, or y.
c
c     INPUT:
c     ndim          dimension of the underlying space
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     tabxx_loc     precomputed tables for second derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot           potential values on tensor product grid
c     grad          gradient values on tensor product grid
c     hess          hessian values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n**ndim)
      real *8 pot(nd,n**ndim)
      real *8 grad(nd,ndim,n**ndim)
      real *8 hess(nd,ndim*(ndim+1)/2,n**ndim)

      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      real *8 tabxx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(ndim)

      if (ifpgh.eq.1) goto 1100
      if (ifpgh.eq.2) goto 2200
      if (ifpgh.eq.3) goto 3300

 1100 continue
      if (ndim.eq.1) then
         call g1d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      elseif (ndim.eq.2) then
         call g2d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      elseif (ndim.eq.3) then
         call g3d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      endif
      return
      
 2200 continue
      if (ndim.eq.1) then
         call g1d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ixyz)
      elseif (ndim.eq.2) then
         call g2d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ixyz)
      elseif (ndim.eq.3) then
         call g3d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ixyz)
      endif
      return

 3300 continue
      if (ndim.eq.1) then
         call g1d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixyz)
      elseif (ndim.eq.2) then
         call g2d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixyz)
      elseif (ndim.eq.3) then
         call g3d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixyz)
      endif
         
      return
      end subroutine
c
c
C
c
C
C
c************************************************************************
c
c     1d tens_prod_to_potloc,pgloc,pghloc
c
C*********************************************************************C
      subroutine g1d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ix)
C*********************************************************************C
c     This routine computes the 1D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot           potential values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n),pot(nd,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ix
c
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1

      if (nx.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1)
            enddo
            pot(ind,k1)=pot(ind,k1)+cd
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
C*********************************************************************C
      subroutine g1d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ix)
C*********************************************************************C
c     This routine computes the 1D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot          potential values on tensor product grid
c     grad         gradient values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n),pot(nd,n),grad(nd,n)
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ix
c
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1

      if (nx.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1)
            enddo
            pot(ind,k1)=pot(ind,k1)+cd
            grad(ind,k1)=grad(ind,k1)+cdx
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
C*********************************************************************C
      subroutine g1d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ix)
C*********************************************************************C
c     This routine computes the 1D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     tabxx_loc     precomputed tables for second derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot          potential values on tensor product grid
c     grad         gradient values on tensor product grid
c     hess         hessian values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n),pot(nd,n),grad(nd,n),hess(nd,n)
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      real *8 tabxx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ix
c
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1

      if (nx.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx=0
            cdxx=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1)
               cdxx=cdxx+tabxx_loc(j1,k1,ix)*fvals(ind,j1)
            enddo
            pot(ind,k1)=pot(ind,k1)+cd
            grad(ind,k1)=grad(ind,k1)+cdx
            hess(ind,k1)=hess(ind,k1)+cdxx
         enddo
      enddo
      
      return
      end subroutine
c
c
C
c
c************************************************************************
c
c     2d tens_prod_to_potloc,pgloc,pghloc
c
C*********************************************************************C
      subroutine g2d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixy)
C*********************************************************************C
c     This routine computes the 2D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot           potential values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n),pot(nd,n,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixy(2)
      real *8 ff(n,n)
c
      ix = ixy(1)
      iy = ixy(2)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1

      if (nx.eq.0 .or. ny.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
         enddo
         enddo
c        transfrom in y
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
               cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+cd
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
      subroutine g2d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ixy)
C*********************************************************************C
c     This routine computes the 2D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot          potential values on tensor product grid
c     grad         gradient values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n),pot(nd,n,n),grad(nd,2,n,n)
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixy(2)
      real *8 ff(n,n)
      real *8 ffx(n,n)
c
      ix = ixy(1)
      iy = ixy(2)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1

      if (nx.eq.0 .or. ny.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
            ffx(k1,j2)=cdx
         enddo
         enddo
c        transfrom in y
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
               cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2)
               cdy=cdy+tabx_loc(j2,k2,iy)*ff(k1,j2)
               cdx=cdx+tab_loc(j2,k2,iy)*ffx(k1,j2)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+cd
            grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+cdx
            grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+cdy
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
      subroutine g2d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixy)
C*********************************************************************C
c     This routine computes the 2D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     tabxx_loc     precomputed tables for second derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot           potential values on tensor product grid
c     grad          gradient values on tensor product grid
c     hess          hessian values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n),pot(nd,n,n),grad(nd,2,n,n),hess(nd,3,n,n)
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      real *8 tabxx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixy(2)
      real *8 ff(n,n)
      real *8 ffx(n,n)
      real *8 ffxx(n,n)
c
      ix = ixy(1)
      iy = ixy(2)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1

      if (nx.eq.0 .or. ny.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx=0
            cdxx=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1,j2)
               cdxx=cdxx+tabxx_loc(j1,k1,ix)*fvals(ind,j1,j2)
            enddo
            ff(k1,j2)=cd
            ffx(k1,j2)=cdx
            ffxx(k1,j2)=cdxx
         enddo
         enddo
c        transfrom in y
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd = 0.0d0
            cdx = 0.0d0
            cdy = 0.0d0
            cdxx = 0.0d0
            cdxy = 0.0d0
            cdyy = 0.0d0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
               cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2)
               cdy=cdy+tabx_loc(j2,k2,iy)*ff(k1,j2)
               cdyy=cdyy+tabxx_loc(j2,k2,iy)*ff(k1,j2)

               cdx=cdx+tab_loc(j2,k2,iy)*ffx(k1,j2)
               cdxy=cdxy+tabx_loc(j2,k2,iy)*ffx(k1,j2)

               cdxx=cdxx+tab_loc(j2,k2,iy)*ffxx(k1,j2)
            enddo
            pot(ind,k1,k2)=pot(ind,k1,k2)+cd
            grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+cdx
            grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+cdy

            hess(ind,1,k1,k2)=hess(ind,1,k1,k2)+cdxx
            hess(ind,2,k1,k2)=hess(ind,2,k1,k2)+cdxy
            hess(ind,3,k1,k2)=hess(ind,3,k1,k2)+cdyy
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
c************************************************************************
c
c     3d tens_prod_to_potloc,pgloc,pghloc
c
C*********************************************************************C
      subroutine g3d_tens_prod_to_potloc(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes 3D volume Gauss transform over a 
c     single box source distribution given as function values
c     on a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     ind_loc       nonzero patterns of tab_loc
c     ixyz          pointers to local table, specify which local table should
c                   be used
c
c     OUTPUT:
c     pot           potential values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n)
      real *8 pot(nd,n,n,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(3)

      real *8 ff(n,n,n),ff2(n,n,n)
c

      ix=ixyz(1)
      iy=ixyz(2)
      iz=ixyz(3)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1
      nz = ind_loc(2,n+1,iz)-ind_loc(1,n+1,iz)+1

      if (nx.eq.0 .or. ny.eq.0 .or. nz.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
               cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2,j3)
            enddo
            ff2(k1,k2,j3)=cd
         enddo
         enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
               cd=cd+tab_loc(j3,k3,iz)*ff2(k1,k2,j3)
            enddo
            pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
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
C
c
C*********************************************************************C
      subroutine g3d_tens_prod_to_pgloc(nd,n,fvals,pot,grad,
     1    tab_loc,tabx_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes 3D volume Gauss transform over a 
c     single box source distribution given as function values
c     on a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot          potential values on tensor product grid
c     grad         gradient values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(3)

      real *8 ff(n,n,n),ff2(n,n,n)
      real *8 ffx(n,n,n),ff2x(n,n,n)
      real *8 ff2y(n,n,n)
c

      ix=ixyz(1)
      iy=ixyz(2)
      iz=ixyz(3)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1
      nz = ind_loc(2,n+1,iz)-ind_loc(1,n+1,iz)+1

      if (nx.eq.0 .or. ny.eq.0 .or. nz.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx=0.0d0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
            ffx(k1,j2,j3)=cdx
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
               cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2,j3)
               cdy=cdy+tabx_loc(j2,k2,iy)*ff(k1,j2,j3)
               
               cdx=cdx+tab_loc(j2,k2,iy)*ffx(k1,j2,j3)
            enddo
            ff2(k1,k2,j3)=cd
            ff2x(k1,k2,j3)=cdx
            ff2y(k1,k2,j3)=cdy
         enddo
         enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            cdz = 0.0d0
            do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
               cd=cd+tab_loc(j3,k3,iz)*ff2(k1,k2,j3)
               cdz=cdz+tabx_loc(j3,k3,iz)*ff2(k1,k2,j3)
               
               cdx=cdx+tab_loc(j3,k3,iz)*ff2x(k1,k2,j3)
               cdy=cdy+tab_loc(j3,k3,iz)*ff2y(k1,k2,j3)
            enddo
            pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
            grad(ind,1,k1,k2,k3)=grad(ind,1,k1,k2,k3)+cdx
            grad(ind,2,k1,k2,k3)=grad(ind,2,k1,k2,k3)+cdy
            grad(ind,3,k1,k2,k3)=grad(ind,3,k1,k2,k3)+cdz
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
      subroutine g3d_tens_prod_to_pghloc(nd,n,fvals,pot,grad,hess,
     1    tab_loc,tabx_loc,tabxx_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes 3D volume Gauss transform over a 
c     single box source distribution given as function values
c     on a tensor product grid.
c
c     INPUT:
c     nd            vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     fvals         function values at tensor grid
c     tab_loc       precomputed tables of 1D integrals
c     tabx_loc      precomputed tables for first derivatives
c     tabxx_loc     precomputed tables for second derivatives
c     ind_loc       precomputed nonzero pattern of tables
c     ixy           pointers to local table, specify which local table should
c                   be used
c     OUTPUT:
c     pot           potential values on tensor product grid
c     grad          gradient values on tensor product grid
c     hess          hessian values on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 hess(nd,6,n,n,n)

      real *8 tab_loc(n,n,-6:6)
      real *8 tabx_loc(n,n,-6:6)
      real *8 tabxx_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(3)

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

      ix=ixyz(1)
      iy=ixyz(2)
      iz=ixyz(3)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1
      nz = ind_loc(2,n+1,iz)-ind_loc(1,n+1,iz)+1

      if (nx.eq.0 .or. ny.eq.0 .or. nz.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do j3=1,n
         do j2=1,n
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
c         do k1=1,n
            cd=0
            cdx=0.0d0
            cdxx=0.0d0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
c            do j1=1,n
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
               cdx=cdx+tabx_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
               cdxx=cdxx+tabxx_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
            enddo
            ff(k1,j2,j3)=cd
            ffx(k1,j2,j3)=cdx
            ffxx(k1,j2,j3)=cdxx
         enddo
         enddo
         enddo

c        transform in y
         do j3=1,n
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
c         do k2=1,n
c         do k1=1,n
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            cdxx = 0.0d0
            cdxy = 0.0d0
            cdyy = 0.0d0
            do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
c            do j2=1,n
               cd   = cd   +   tab_loc(j2,k2,iy)*ff(k1,j2,j3)
               cdy  = cdy  +  tabx_loc(j2,k2,iy)*ff(k1,j2,j3)
               cdyy = cdyy + tabxx_loc(j2,k2,iy)*ff(k1,j2,j3)

               cdx  = cdx  +   tab_loc(j2,k2,iy)*ffx(k1,j2,j3)
               cdxy = cdxy +  tabx_loc(j2,k2,iy)*ffx(k1,j2,j3)
               
               cdxx = cdxx +   tab_loc(j2,k2,iy)*ffxx(k1,j2,j3)
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
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
c         do k3=1,n
c         do k2=1,n
c         do k1=1,n
            cd=0
            cdx = 0.0d0
            cdy = 0.0d0
            cdz = 0.0d0
            cdxx = 0.0d0
            cdyy = 0.0d0
            cdzz = 0.0d0
            cdxy = 0.0d0
            cdxz = 0.0d0
            cdyz = 0.0d0
            do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
c            do j3=1,n
               cd   = cd   +   tab_loc(j3,k3,iz)*ff2(k1,k2,j3)
               cdz  = cdz  +  tabx_loc(j3,k3,iz)*ff2(k1,k2,j3)
               cdzz = cdzz + tabxx_loc(j3,k3,iz)*ff2(k1,k2,j3)

               cdx  = cdx  +  tab_loc(j3,k3,iz)*ff2x(k1,k2,j3)
               cdxz = cdxz + tabx_loc(j3,k3,iz)*ff2x(k1,k2,j3)

               cdy  = cdy  +  tab_loc(j3,k3,iz)*ff2y(k1,k2,j3)
               cdyz = cdyz + tabx_loc(j3,k3,iz)*ff2y(k1,k2,j3)
               
               cdxx = cdxx + tab_loc(j3,k3,iz)*ff2xx(k1,k2,j3)
               cdxy = cdxy + tab_loc(j3,k3,iz)*ff2xy(k1,k2,j3)
               cdyy = cdyy + tab_loc(j3,k3,iz)*ff2yy(k1,k2,j3)
            enddo
            pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
            grad(ind,1,k1,k2,k3)=grad(ind,1,k1,k2,k3)+cdx
            grad(ind,2,k1,k2,k3)=grad(ind,2,k1,k2,k3)+cdy
            grad(ind,3,k1,k2,k3)=grad(ind,3,k1,k2,k3)+cdz

            hess(ind,1,k1,k2,k3)=hess(ind,1,k1,k2,k3)+cdxx
            hess(ind,2,k1,k2,k3)=hess(ind,2,k1,k2,k3)+cdyy
            hess(ind,3,k1,k2,k3)=hess(ind,3,k1,k2,k3)+cdzz

            hess(ind,4,k1,k2,k3)=hess(ind,4,k1,k2,k3)+cdxy
            hess(ind,5,k1,k2,k3)=hess(ind,5,k1,k2,k3)+cdxz
            hess(ind,6,k1,k2,k3)=hess(ind,6,k1,k2,k3)+cdyz
         enddo
         enddo
         enddo
c     end of the ind loop
      enddo
      
      return
      end subroutine
c
c******************************************************************************
C
c     fast version - use sparse patterns of local tables
C
c******************************************************************************
      subroutine gnd_tens_prod_to_potloc_fast(ndim,nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes the volume Gauss transform in general n dimensions
c     over a single box source distribution given as function values on a
c     tensor product grid.
c      
c     The target points have a fixed location w.r.t. source box
c     and the integrals of Gaussians times Legendre polynomials at 
c     those points is assumed to have been precomputed and stored 
c     in arrays tab_loc.
c     Thus, the specific geometric relation of the source and target
c     boxes are IMPLICITLY contained in these arrays.
c     There are many such relations in higher dimensions, but only a few one-dimensional
c     tables are needed corresponding to the range of possible shifts
c     of the box center in any single dimension.
c
c
c     Case 1: same level
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |     |    | 
c         |     |  T  |    |    target points in T
c         |_____|_____|____|    source box has offset in x and y.  
c         |     |     |    |    Because of separation of variables,
c         |     |     |    |    we can use 1D tables for desired 
c         |_____|_____|____|    offsets in x or y in range (-1,0,1). 
c
c     Case 2: different levels
c          _____ _____ ____  
c         |     |     |    | 
c         |     |     |    | 
c         |_____|_____|____| 
c         |     |A |  |    | 
c         |     |--|--| B  |   for target points in small box A, of 
c         |_____|__|__|____|   dimension D, adjacent large boxes can be   
c         |     |     |    |   offset by one of -3D/2,-D/2,D/2,3D/2
c         |     |     |    |   in either x, or y.
c         |_____|_____|____|   
c                              For target points in large box B, of
c                              dimension D, adjacent small boxes can be
c                              offset by one of -3D/4,-D/4,D/4,3D/4
c                              in either x, or y.
c
c     INPUT:
c     ndim        dimension of the underlying space
c     nd          vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     n           dimension of coeff array
c     fvals       function values at tensor grid
c     tab_loc     precomputed tables of 1D integrals
c     ind_loc     precomputed nonzero pattern of tables
c                 This is used to speed up the matrix-vector multiplication,
c                 especially when the Gaussian is sharply peaked, i.e., small delta.
c
c     ixyz        pointers to local table, specify which local table should
c                 be used along each direction
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n**ndim),pot(nd,n*ndim)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(ndim)

      if (ndim.eq.1) then
         call g1d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      elseif (ndim.eq.2) then
         call g2d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      elseif (ndim.eq.3) then
         call g3d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
      endif
      
      return
      end subroutine
c
c
C
c
C
C
C*********************************************************************C
      subroutine g1d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ix)
C*********************************************************************C
c     This routine computes the 1D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     n           dimension of coeff array
c     fvals       function values at tensor grid
c     tab_loc     precomputed tables of 1D integrals
c     ind_loc     precomputed nonzero pattern of tables
c     ixy         pointers to local table, specify which local table should
c                 be used
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n),pot(nd,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ix
c
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1

      if (nx.eq.0) return
      
      do ind = 1,nd
c        transform in x
         do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
            cd=0
            do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
               cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1)
            enddo
            pot(ind,k1)=pot(ind,k1)+cd
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
      subroutine g2d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixy)
C*********************************************************************C
c     This routine computes the 2D volume Gauss transform over a 
c     single box source distribution given as function values on
c     a tensor product grid.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     n           dimension of coeff array
c     fvals       function values at tensor grid
c     tab_loc     precomputed tables of 1D integrals
c     ind_loc     precomputed nonzero pattern of tables
c     ixy         pointers to local table, specify which local table should
c                 be used
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n),pot(nd,n,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixy(2)
      real *8 ff(n,n)
c
      ix = ixy(1)
      iy = ixy(2)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1

      if (nx.eq.0 .or. ny.eq.0) return
      
      if (nx.le.ny) then
      
      do ind = 1,nd
c        transform in x
         do j2=1,n
            do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
               cd=0
               do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                  cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2)
               enddo
               ff(k1,j2)=cd
            enddo
         enddo
c        transfrom in y
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
            do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
               cd=0
               do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                  cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2)
               enddo
               pot(ind,k1,k2)=pot(ind,k1,k2)+cd
            enddo
         enddo
c     end of the ind loop
      enddo
      else
      do ind = 1,nd
c        transform in y
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
            do j1=1,n
               cd=0
               do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                  cd=cd+tab_loc(j2,k2,iy)*fvals(ind,j1,j2)
               enddo
               ff(j1,k2)=cd
            enddo
         enddo
c        transfrom in x
         do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
            do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
               cd=0
               do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                  cd=cd+tab_loc(j1,k1,ix)*ff(j1,k2)
               enddo
               pot(ind,k1,k2)=pot(ind,k1,k2)+cd
            enddo
         enddo
c     end of the ind loop
      enddo
      endif
      
      return
      end subroutine
c
c
C
c
C*********************************************************************C
      subroutine g3d_tens_prod_to_potloc_fast(nd,n,fvals,pot,
     1    tab_loc,ind_loc,ixyz)
C*********************************************************************C
c     This routine computes 3D volume Gauss transform over a 
c     single box source distribution given as function values
c     on a tensor product grid.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n             number of nodes along each dimension
c     n           dimension of coeff array
c     fvals       function values at tensor grid
c     tab_loc     precomputed tables of 1D integrals
c     ind_loc     nonzero patterns of tab_loc
c     ixyz        pointers to local table, specify which local table should
c                 be used
c
c     OUTPUT:
c     pot         output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,n,n,n)
      real *8 pot(nd,n,n,n)
      real *8 tab_loc(n,n,-6:6)
      integer ind_loc(2,n+1,-6:6)
      integer ixyz(3)

      real *8 ff(n,n,n),ff2(n,n,n)
c

      ix=ixyz(1)
      iy=ixyz(2)
      iz=ixyz(3)
      
      nx = ind_loc(2,n+1,ix)-ind_loc(1,n+1,ix)+1
      ny = ind_loc(2,n+1,iy)-ind_loc(1,n+1,iy)+1
      nz = ind_loc(2,n+1,iz)-ind_loc(1,n+1,iz)+1

      if (nx.eq.0 .or. ny.eq.0 .or. nz.eq.0) return
      
      if (nx.le.ny .and. ny.le.nz) then
      do ind = 1,nd
c        transform in x
         do j3=1,n
            do j2=1,n
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(k1,j2,j3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do j3=1,n
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*ff(k1,j2,j3)
                  enddo
                  ff2(k1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*ff2(k1,k2,j3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo
      elseif (nx .le. nz .and. nz.lt.ny) then
      do ind = 1,nd
c        transform in x
         do j3=1,n
            do j2=1,n
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(k1,j2,j3)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do j2=1,n
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*ff(k1,j2,j3)
                  enddo
                  ff2(k1,j2,k3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*ff2(k1,j2,k3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo
      elseif (ny .lt. nx .and. nx.le.nz) then
      do ind = 1,nd
c        transform in y
         do j3=1,n
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do j1=1,n
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(j1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transform in x
         do j3=1,n
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*ff(j1,k2,j3)
                  enddo
                  ff2(k1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*ff2(k1,k2,j3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo
      elseif (ny .le. nz .and. nz.lt.nx) then
      do ind = 1,nd
c        transform in y
         do j3=1,n
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do j1=1,n
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(j1,k2,j3)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do j1=1,n
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*ff(j1,k2,j3)
                  enddo
                  ff2(j1,k2,k3)=cd
               enddo
            enddo
         enddo

c        transform in x
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*ff2(j1,k2,k3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo      
      elseif (nz .lt. nx .and. nx.le.ny) then
      do ind = 1,nd
c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do j2=1,n
               do j1=1,n
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(j1,j2,k3)=cd
               enddo
            enddo
         enddo

c        transform in x
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do j2=1,n
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*ff(j1,j2,k3)
                  enddo
                  ff2(k1,j2,k3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*ff2(k1,j2,k3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo
      elseif (nz .lt. ny .and. ny.lt.nx) then
      do ind = 1,nd
c        transform in z
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do j2=1,n
               do j1=1,n
                  cd=0
                  do j3=ind_loc(1,k3,iz),ind_loc(2,k3,iz)
                     cd=cd+tab_loc(j3,k3,iz)*fvals(ind,j1,j2,j3)
                  enddo
                  ff(j1,j2,k3)=cd
               enddo
            enddo
         enddo

c        transform in y
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do j1=1,n
                  cd=0
                  do j2=ind_loc(1,k2,iy),ind_loc(2,k2,iy)
                     cd=cd+tab_loc(j2,k2,iy)*ff(j1,j2,k3)
                  enddo
                  ff2(j1,k2,k3)=cd
               enddo
            enddo
         enddo

c        transform in x
         do k3=ind_loc(1,n+1,iz),ind_loc(2,n+1,iz)
            do k2=ind_loc(1,n+1,iy),ind_loc(2,n+1,iy)
               do k1=ind_loc(1,n+1,ix),ind_loc(2,n+1,ix)
                  cd=0
                  do j1=ind_loc(1,k1,ix),ind_loc(2,k1,ix)
                     cd=cd+tab_loc(j1,k1,ix)*ff2(j1,k2,k3)
                  enddo
                  pot(ind,k1,k2,k3)=pot(ind,k1,k2,k3)+cd
               enddo
            enddo
         enddo
c     end of the ind loop
      enddo      
      endif

      
      return
      end subroutine
c
c
C
c
C
      subroutine gnd_find_loctab_ind(ndim,iperiod,tcenter,scenter,
     1    sboxsize,bs0,ixyz)
c     returns an index arrary used in gnd_tens_prod_to_potloc
c
c     input:
c     ndim - dimension of the underlying space
c     iperiod - 0: free space; 1: doubly periodic
c     tcenter - target box center
c     scenter - source box center
c     sboxsize - source box size
c     bs0 - root box size
c
c     output
c     ixyz - an index array determining which local table 
c            should be used along each dimension
c
      implicit real *8 (a-h,o-z)
      real *8 tcenter(ndim),scenter(ndim)
      integer ixyz(ndim),i

      bs=sboxsize/4

      do i=1,ndim
         dx = (tcenter(i)-scenter(i))/bs
         if (iperiod .eq. 1) then
            dxp1=dx-bs0/bs
            dxm1=dx+bs0/bs
            if (abs(dx).gt.abs(dxp1)) dx=dxp1
            if (abs(dx).gt.abs(dxm1)) dx=dxm1
         endif
         ixyz(i)=dx
      enddo
      
      
      return
      end subroutine
c
c
C
c
C
c************************************************************************
C
C     make plane-wave translation matrices at the cutoff level (mp to loc)
C
C************************************************************************
      subroutine gnd_mk_pw_translation_matrices(ndim,xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C     ndim    = dimension of underlying space
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw**ndim/2,(2*nmax+1)**ndim)

      if (ndim.eq.1) then
         call  g1d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1       wshift)
      elseif (ndim.eq.2) then
         call  g2d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1       wshift)
      elseif (ndim.eq.3) then
         call  g3d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1       wshift)
      endif
      
      return
      end
c
c
c     
c
      subroutine g1d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw/2,(2*nmax+1))
      
      complex *16,allocatable:: ww(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(ww(npw,-nmax:nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         ww(j1,0)=1
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
            ww(j1,-k1) = dconjg(ww(j1,k1))
         enddo
      enddo

      k=0
      do k1=-nmax,nmax
         k=k+1
         j=0   
         do j1=1,npw/2
            j=j+1
            wshift(j,k) = ww(j1,k1)
         enddo
      enddo
c
      return
      end
c
c
c     
c
      subroutine g2d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw*npw/2,(2*nmax+1)**2)
      
      complex *16,allocatable:: ww(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(ww(npw,-nmax:nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         ww(j1,0)=1
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
            ww(j1,-k1) = dconjg(ww(j1,k1))
         enddo
      enddo

      k=0
      do k1=-nmax,nmax
      do k2=-nmax,nmax
         k=k+1
         j=0   
         do j1=1,npw/2
         do j2=1,npw
            j=j+1
            wshift(j,k) = ww(j2,k2)*ww(j1,k1)
         enddo
         enddo
      enddo
      enddo
c
      return
      end
c
c
c     
c
      subroutine g3d_mk_pw_translation_matrices(xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw*npw*npw/2,(2*nmax+1)**3)
      
      complex *16,allocatable:: ww(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(ww(npw,-nmax:nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         ww(j1,0)=1
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
            ww(j1,-k1) = dconjg(ww(j1,k1))
         enddo
      enddo

      k=0
      do k1=-nmax,nmax
      do k2=-nmax,nmax
      do k3=-nmax,nmax
         k=k+1
         j=0   
         do j1=1,npw/2
         do j2=1,npw
         do j3=1,npw
            j=j+1
            wshift(j,k) = ww(j3,k3)*ww(j2,k2)*ww(j1,k1)
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
c
      return
      end
c
c
c     
c
c************************************************************************
C
C     make plane-wave merge-split matrices from children to parent
c     and vice versa. That is, mp to mp from children to parent
c     and local to local from parent to children.
C
C************************************************************************
      subroutine gnd_mk_merge_split_pw_matrices(ndim,xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C     ndim     = dimension of the underlying space
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ts    = real *8, 1d PW expansion nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      integer ndim
      real *8 xmin
      complex *16 wshift(npw**ndim/2,2**ndim,nmax)
      real *8 ts(npw)

      if (ndim.eq.1) then
         call  g1d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1       wshift)
      elseif (ndim.eq.2) then
         call  g2d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1       wshift)
      elseif (ndim.eq.3) then
         call  g3d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1       wshift)
      endif
      
      return
      end
c
c
c     
c
      subroutine g1d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ws,ts    = real *8, 1d PW expansion weights and nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw/2,2,nmax)
      real *8 ts(npw)
      complex *16 ztmp
      complex *16 eye
      complex *16, allocatable:: ww(:,:)
C
      eye =dcmplx(0,1)
      
      allocate(ww(npw,nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo
      
      do k1=1,nmax
         do j1=1,npw/2
c           p    
            wshift(j1,1,k1) = ww(j1,k1)
c           m
            wshift(j1,2,k1) = conjg(ww(j1,k1))
         enddo
      enddo

      return
      end
c
c
c
c      
      subroutine g2d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ws,ts    = real *8, 1d PW expansion weights and nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw*npw/2,4,nmax)
      real *8 ts(npw)
      complex *16 ztmp
      complex *16 eye
      complex *16, allocatable:: ww(:,:)
C
      eye =dcmplx(0,1)
      
      allocate(ww(npw,nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo
      
      do k1=1,nmax
         j=0
         do j1=1,npw/2
         do j2=1,npw
            j=j+1
c           pp            
            wshift(j,1,k1) = ww(j2,k1)*ww(j1,k1)
c           mp
            wshift(j,2,k1) = conjg(ww(j2,k1))*ww(j1,k1)
c           pm
            wshift(j,3,k1) = ww(j2,k1)*conjg(ww(j1,k1))
c           mm
            wshift(j,4,k1) = conjg(ww(j2,k1))*conjg(ww(j1,k1))
         enddo
         enddo
      enddo

      return
      end
c
c
c
c      
      subroutine g3d_mk_merge_split_pw_matrices(xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ws,ts    = real *8, 1d PW expansion weights and nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw*npw*npw/2,8,nmax)
      real *8 ts(npw)
      complex *16 ztmp
      complex *16 eye
      complex *16, allocatable:: ww(:,:)
C
      eye =dcmplx(0,1)
      
      allocate(ww(npw,nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo
      
      do k1=1,nmax
         j=0
         do j1=1,npw/2
         do j2=1,npw
         do j3=1,npw
            j=j+1
c           ppp              
            wshift(j,1,k1) = ww(j3,k1)*ww(j2,k1)
     1          *ww(j1,k1)
c           mpp
            wshift(j,2,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *ww(j1,k1)
c           pmp
            wshift(j,3,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *ww(j1,k1)
c           mmp
            wshift(j,4,k1) = conjg(ww(j3,k1))*conjg(ww(j2,k1))
     1          *ww(j1,k1)
            
c           ppm              
            wshift(j,5,k1) = ww(j3,k1)*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           mpm
            wshift(j,6,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           pmm
            wshift(j,7,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *conjg(ww(j1,k1))
c           mmm
            wshift(j,8,k1) = conjg(ww(j3,k1))*conjg(ww(j2,k1))
     1          *conjg(ww(j1,k1))
            
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
c
c************************************************************************
C
C     translate plane-wave expansion at the cutoff level (mp to loc)
C
C************************************************************************
      subroutine gnd_shiftpw(nd,nexp,pwexp1,
     1              pwexp2,wshift)
C
C     This subroutine converts the PW expansion (pwexp1) about
C     the center (CENT1) into an PW expansion (pwexp2) about 
C     (CENT2) using precomputed translation matrix wshift.
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     delta   = Gaussian variance
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C     wshift  = precomputed PW exp translation matrix 
C
C     OUTPUT:
C
C     pwexp2 = shifted expansion 
C
      implicit none
      integer nd,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)
      complex *16 wshift(nexp)

C
      do ind=1,nd
      do j=1,nexp
         pwexp2(j,ind) = pwexp2(j,ind) + pwexp1(j,ind)*wshift(j)
      enddo
      enddo
c
      return
      end
c
C
c
      subroutine gnd_find_pwshift_ind(ndim,iperiod,tcenter,scenter,
     1   bs0,xmin,nmax,ind)
c     returns an index used in gnd_shiftpw
c
c     input:
c     ndim - dimension of the underlying space
c     iperiod - 0: free space; 1: doubly periodic
c     tcenter - target box center
c     scenter - source box center
c     bs0 - root box size
c     xmin - box size at npwlevel
c     nmax - (2*nmax+1) is the number of possible translations along
c            each direction
c     output
c     ind - index to determine which precomputed pwshift matrix 
c           should be used in plane-wave translation at the cut-off
c           level
c
      implicit real *8 (a-h,o-z)
      real *8 tcenter(ndim),scenter(ndim)
      integer jxyz(ndim),i

      do i=1,ndim
         jx= nint((tcenter(i) - scenter(i))/xmin)
         
         if (iperiod .eq. 1) then
            jxp1=nint((tcenter(i) - scenter(i) - bs0)/xmin)
            jxm1=nint((tcenter(i) - scenter(i) + bs0)/xmin)
            if (abs(jx).gt.abs(jxp1)) jx=jxp1
            if (abs(jx).gt.abs(jxm1)) jx=jxm1
         endif

         jxyz(i)=jx+nmax
      enddo

      n=2*nmax+1
      ind=1
      inc=1
      do i=1,ndim
         ind=ind+jxyz(i)*inc
         inc=inc*n
      enddo
      
      return
      end subroutine
c
c
C
c
C
C
c************************************************************************
C
C     translate plane-wave expansion from children to parent
c     (merge mp to mp) or from parent to children (split loc to loc)
C
C************************************************************************
      subroutine gnd_shiftpw_loc(nd,nexp,pwexp1,
     1              pwexp2,wshift)
C
C     This subroutine converts the PW expansion (pwexp1) about
C     the center (CENT1) into an PW expansion (pwexp2) about 
C     (CENT2) using precomputed translation matrix wshift.
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     delta   = Gaussian variance
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C     wshift  = precomputed PW exp translation matrix 
C
C     OUTPUT:
C
C     pwexp2 = shifted expansion 
C
      implicit none
      integer nd,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)
      complex *16 wshift(nexp)

C
      do ind=1,nd
      do j=1,nexp
         pwexp2(j,ind) = pwexp1(j,ind)*wshift(j)
      enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine gnd_find_msshift_ind(ndim,tcen,scen,k)
C
C     This subroutine returns the correct index so that 
C     the correct precomputed plane-wave merge-split
C     matrix is used. To be used either in merge multipole
c     or split local, i.e., parent-children translations.
C
C     INPUT
C     ndim - dimension of the underlysing space
c     tcen - target box center
C     scen - original box center
C
C     OUTPUT:
C
C     k - index for plane-wave merge-split shift table 
C
      implicit real *8 (a-h,o-z)
      real *8 tcen(ndim),scen(ndim)
      integer k, i

C
      k=1
      inc=1
      do i=1,ndim
         dx = tcen(i)-scen(i)
         if (dx.lt.0) k=k+inc
         inc=inc*2
      enddo
c
      return
      end
c
C
c
c************************************************************************
C
C     translate the plane-wave expansion at the cutoff level (mp to loc)
c     from self to self
C
C************************************************************************
C
      subroutine gnd_copy_pwexp(nd,nexp,pwexp1,pwexp2)
C
C     This subroutine copy one PW expansion (pwexp1) 
C     to an PW expansion (pwexp2).
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C
C     OUTPUT:
C
C     pwexp2 = copied expansion 
C
      implicit none
      integer nd,nn,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)

C
      do ind=1,nd
      do j=1,nexp
         pwexp2(j,ind) = pwexp1(j,ind)
      enddo
      enddo
c
      return
      end
c
C
c
C
cC***********************************************************************
      subroutine gnd_pwzero(nd,pwexp,nexp)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector planewave expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     nexp   :   number of terms in the planewave expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     pwexp  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer nd,nexp,n,ind
      complex *16 pwexp(nexp,nd)
c
      do ind=1,nd
      do n=1,nexp
         pwexp(n,ind)=0.0d0
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
C*********************************************************************C
      subroutine ortho_to_g(nd,n,coeff,ff,ff2,grad,
     1           tabf,tabfx,boxsize)
C*********************************************************************C
c     This routine takes an orthogonal polynomial expansion for a 
c     potential on the standard box [-1,1]^2 and computes the 
c     correspomnding gradient and Hessian on the corresponding
c     tensor product grid.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       orthogonal polynomial coefficients
c                 f = sum coeff(n,m) F_n(x) F_m(y)
c     ff          workspace
c     ff2         workspace
c     ff3         workspace
c     tabf        precomputed table  tabf(n,i) =   F_n(x_i)
c     tabfx       precomputed table  tabfx(n,i) =  F'_n(x_i)
c     tabfxx      precomputed table  tabfxx(n,i) = F''_n(x_i)
c
c     OUTPUT:
c     grad        output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(nd,n,n)
      real *8 grad(nd,2,n,n)
      real *8 ff2(n,n)
      real *8 ff(n,n),tabf(n,n)
      real *8 tabfx(n,n)
c
      b2 = boxsize**2
      do ind = 1,nd
c        transform in x
         do j2=1,n
            do k1=1,n
               cd=0
               cdx=0
               do j1=1,n
                  cd=cd+tabf(j1,k1)*coeff(ind,j1,j2)
                  cdx=cdx+tabfx(j1,k1)*coeff(ind,j1,j2)
               enddo
               ff(j2,k1)=cd
               ff2(j2,k1)=cdx
            enddo
         enddo
c        transform in y
         do k2=1,n
            do k1=1,n
               cd=0
               cdx=0
               cdy=0
               do j2=1,n
                  cdx=cdx+tabf(j2,k2)*ff2(j2,k1)
                  cdy=cdy+tabfx(j2,k2)*ff(j2,k1)
               enddo
               grad(ind,1,k1,k2)=cdx*2/boxsize
               grad(ind,2,k1,k2)=cdy*2/boxsize
            enddo
         enddo
c     end of the ind loop
      enddo
c
c
      return
      end subroutine
c
c
c
c
c
c
C*********************************************************************C
      subroutine ortho_to_gh(nd,n,coeff,ff,ff2,ff3,grad,
     1           hess,tabf,tabfx,tabfxx,boxsize)
C*********************************************************************C
c     This routine takes an orthogonal polynomial expansion for a 
c     potential on the standard box [-1,1]^2 and computes the 
c     correspomnding gradient and Hessian on the corresponding
c     tensor product grid.
c
c     INPUT:
c     nd          vector length (for multiple RHS)
c     n           dimension of coeff array
c     coeff       orthogonal polynomial coefficients
c                 f = sum coeff(n,m) F_n(x) F_m(y)
c     ff          workspace
c     ff2         workspace
c     ff3         workspace
c     tabf        precomputed table  tabf(n,i) =   F_n(x_i)
c     tabfx       precomputed table  tabfx(n,i) =  F'_n(x_i)
c     tabfxx      precomputed table  tabfxx(n,i) = F''_n(x_i)
c
c     OUTPUT:
c     grad        output on tensor product grid
c     hess        output on tensor product grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(nd,n,n)
      real *8 grad(nd,2,n,n)
      real *8 hess(nd,3,n,n)
      real *8 ff3(n,n)
      real *8 ff2(n,n)
      real *8 ff(n,n),tabf(n,n)
      real *8 tabfx(n,n),tabfxx(n,n)
c
      b2 = boxsize**2
      do ind = 1,nd
c        transform in x
         do j2=1,n
            do k1=1,n
               cd=0
               cdx=0
               cdxx=0
               do j1=1,n
                  cd=cd+tabf(j1,k1)*coeff(ind,j1,j2)
                  cdx=cdx+tabfx(j1,k1)*coeff(ind,j1,j2)
                  cdxx=cdxx+tabfxx(j1,k1)*coeff(ind,j1,j2)
               enddo
               ff(j2,k1)=cd
               ff2(j2,k1)=cdx
               ff3(j2,k1)=cdxx
            enddo
         enddo
c        transform in y
         do k2=1,n
            do k1=1,n
               cd=0
               cdx=0
               cdy=0
               cdxx=0
               cdxy=0
               cdyy=0
               do j2=1,n
                  cdx=cdx+tabf(j2,k2)*ff2(j2,k1)
                  cdxx=cdxx+tabf(j2,k2)*ff3(j2,k1)
                  cdy=cdy+tabfx(j2,k2)*ff(j2,k1)
                  cdxy=cdxy+tabfx(j2,k2)*ff2(j2,k1)
                  cdyy=cdyy+tabfxx(j2,k2)*ff(j2,k1)
               enddo
               grad(ind,1,k1,k2)=cdx*2/boxsize
               grad(ind,2,k1,k2)=cdy*2/boxsize
               hess(ind,1,k1,k2)=cdxx*4/b2
               hess(ind,2,k1,k2)=cdxy*4/b2
               hess(ind,3,k1,k2)=cdyy*4/b2
            enddo
         enddo
c     end of the ind loop
      enddo
c
c
      return
      end subroutine
c
c
C
c 
c 
c 
      subroutine mk_poly_tables(norder,xnodes,tabf,tabfx,tabfxx)
      implicit real *8 (a-h,o-z)
      integer norder
      real *8 xnodes(norder)
      real *8 tabf(norder,norder)
      real *8 tabfx(norder,norder)
      real *8 tabfxx(norder,norder)
C 
C     This subroutine computes the values and the derivatives
c     of first norder Legendre polynomials at the Legendre nodes
C     in interval [-1,1] as well as their first and second derivatives.
c 
c     INPUT:
c 
C     norder = order of expansion
C     xnodes = Legendre nodes
c 
c     OUTPUT:
c 
C     tabf:    tabf(i,j) = P_i(x_j)
C     tabfx:   tabf(i,j) = P'_i(x_j)
C     tabfxx:  tabf(i,j) = P''_i(x_j)
C 
C 
      do k = 1,norder
         x = xnodes(k)
         done=1
         pjm2=1
         pjm1=x
         dxjm2=0
         dxjm1=1
         dxxjm2=0
         dxxjm1=0
c 
         tabf(1,k)=1
         tabfx(1,k)=0
         tabfxx(1,k)=0
c 
         tabf(2,k)=x
         tabfx(2,k)=1
         tabfxx(2,k)=0
c 
         do j = 2,norder-1
            pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
            dxj=(2*j-1)*(pjm1+x*dxjm1)-(j-1)*dxjm2
            dxxj=(2*j-1)*(2.0d0*dxjm1+x*dxxjm1)-(j-1)*dxxjm2
c 
            dxj=dxj/j
            dxxj=dxxj/j
  
            tabf(j+1,k)=pj
            tabfx(j+1,k)=dxj
            tabfxx(j+1,k)=dxxj
c 
            pjm2=pjm1
            pjm1=pj
            dxjm2=dxjm1
            dxjm1=dxj
            dxxjm2=dxxjm1
            dxxjm1=dxxj
         enddo
      enddo
      return
      end
c 
c 
c 
