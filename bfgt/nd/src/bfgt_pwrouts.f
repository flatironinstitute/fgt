C     This file contains routines for handling plane wave expansions
c     used in the box fgt in general n dimensions, including tensor
c     product function values to plane wave expansions, plane wave
C     expansions to potential/gradient/hessian values at tensor grid
c     or arbritrary targets, and plane wave merge (mp to mp), split
c     (loc to loc), and translation (mp to loc).
c
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
c                  NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
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
      complex *16 pwexp(((npw+1)/2)*npw**(ndim-1),nd)

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
c                 NOTE 1D convention is pwexp(((npw+1)/2))
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
      complex *16 pwexp(((npw+1)/2),nd),cd
c
      do ind = 1,nd
         do k1 = 1,((npw+1)/2)
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
c                 NOTE 2D convention is pwexp(npw,npw,((npw+1)/2))
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
      complex *16 pwexp(npw,((npw+1)/2),nd)
      complex *16 ff(n,((npw+1)/2)),cd
c
      do ind = 1,nd
         do m1 = 1,n
            do k2 = 1,((npw+1)/2)
               cd = 0.0d0
               do m2 = 1,n
                  cd = cd+tab_poly2pw(m2,k2)*fvals(ind,m1,m2)
               enddo
               ff(m1,k2) = cd
            enddo
         enddo
c
         do k2 = 1,((npw+1)/2)
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
c                 NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
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
      complex *16 pwexp(npw,npw,((npw+1)/2),nd)

      complex *16 ff(n,n,((npw+1)/2))
      complex *16 ff2(n,npw,((npw+1)/2))
      complex *16 cd
c
      do ind = 1,nd
         do k3 = 1,((npw+1)/2)
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
         do k3 = 1,((npw+1)/2)
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
         do k3 = 1,((npw+1)/2)
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
c                 NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(((npw+1)/2)*npw**(ndim-1),nd)
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
c                 NOTE 1D convention is pwexp(((npw+1)/2))
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 pwexp(((npw+1)/2),nd),cd
c
      npw2=npw/2
      
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            do m1 = 1,npw2
               cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,ind)
            enddo
c     when npw is an odd number, zero frequency needs special treatment
            m1=((npw+1)/2)
            if (m1.gt.npw2) then
               cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,ind)/2
            endif
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
      complex *16 pwexp(((npw+1)/2),nd),cd,cdx
c
      npw2=npw/2
      
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            do m1 = 1,npw2
               cd  = cd  + tab_pw2pot(m1,k1)*pwexp(m1,ind)
               cdx = cdx + tab_pw2der(m1,k1)*pwexp(m1,ind)
            enddo
            
            m1=((npw+1)/2)
            if (m1.gt.npw2) then
               cd  = cd  + tab_pw2pot(m1,k1)*pwexp(m1,ind)/2
               cdx = cdx + tab_pw2der(m1,k1)*pwexp(m1,ind)/2
            endif
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
      complex *16 pwexp(((npw+1)/2),nd),cd,cdx,cdxx
c
      npw2=npw/2
      do ind = 1,nd
         do k1 = 1,n
            cd=0.0d0
            cdx=0.0d0
            cdxx=0.0d0
            do m1 = 1,npw2
               cd   = cd   + tab_pw2pot(m1,k1)*pwexp(m1,ind)
               cdx  = cdx  + tab_pw2der(m1,k1)*pwexp(m1,ind)
               cdxx = cdxx + tab_pw2dxx(m1,k1)*pwexp(m1,ind)
            enddo
            m1=((npw+1)/2)
            if (m1.gt.npw2) then
               cd   = cd   + tab_pw2pot(m1,k1)*pwexp(m1,ind)/2
               cdx  = cdx  + tab_pw2der(m1,k1)*pwexp(m1,ind)/2
               cdxx = cdxx + tab_pw2dxx(m1,k1)*pwexp(m1,ind)/2
            endif
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
c                 NOTE 2D convention is pwexp(npw,((npw+1)/2))
c     ff       complex workspace (n,((npw+1)/2))
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 pwexp(npw,((npw+1)/2),nd),cd
      complex *16 ff(n,((npw+1)/2))
c
      npw2=npw/2
      
      do ind = 1,nd
         do m2 = 1,((npw+1)/2)
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
            do m2 = 1,npw2
               cd = cd+tab_pw2pot(m2,k2)*ff(k1,m2)
            enddo
            m2=((npw+1)/2)
            if (m2.gt.npw2) then
               cd = cd+tab_pw2pot(m2,k2)*ff(k1,m2)/2
            endif
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
      complex *16 ff(((npw+1)/2),n)
      complex *16 ffx(((npw+1)/2),n)
      complex *16 pwexp(npw,((npw+1)/2),nd),cd,cdx,cdy
c
      npw2=npw/2
      
      do ind = 1,nd
         do m2 = 1,((npw+1)/2)
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
            do m2 = 1,npw2
               cd  = cd  + tab_pw2pot(m2,k2)*ff(m2,k1)
               cdy = cdy + tab_pw2der(m2,k2)*ff(m2,k1)
               
               cdx = cdx + tab_pw2pot(m2,k2)*ffx(m2,k1)
            enddo
            m2=((npw+1)/2)
            if (m2.gt.npw2) then
               cd  = cd  + tab_pw2pot(m2,k2)*ff(m2,k1)/2
               cdy = cdy + tab_pw2der(m2,k2)*ff(m2,k1)/2
               
               cdx = cdx + tab_pw2pot(m2,k2)*ffx(m2,k1)/2
            endif
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
      complex *16 ff(((npw+1)/2),n)
      complex *16 ffx(((npw+1)/2),n)
      complex *16 ffxx(((npw+1)/2),n)
      complex *16 pwexp(npw,((npw+1)/2),nd),cd,cdx,cdy
      complex *16 cdxx,cdxy,cdyy
c
      npw2=npw/2
      
      do ind = 1,nd
         do m2 = 1,((npw+1)/2)
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
            do m2 = 1,npw2
               cd   = cd   + tab_pw2pot(m2,k2)*ff(m2,k1)
               cdy  = cdy  + tab_pw2der(m2,k2)*ff(m2,k1)
               cdyy = cdyy + tab_pw2dxx(m2,k2)*ff(m2,k1)

               cdx  = cdx  + tab_pw2pot(m2,k2)*ffx(m2,k1)
               cdxy = cdxy + tab_pw2der(m2,k2)*ffx(m2,k1)

               cdxx = cdxx + tab_pw2pot(m2,k2)*ffxx(m2,k1)
            enddo
            m2=((npw+1)/2)
            if (m2.gt.npw2) then
               cd   = cd   + tab_pw2pot(m2,k2)*ff(m2,k1)/2
               cdy  = cdy  + tab_pw2der(m2,k2)*ff(m2,k1)/2
               cdyy = cdyy + tab_pw2dxx(m2,k2)*ff(m2,k1)/2

               cdx  = cdx  + tab_pw2pot(m2,k2)*ffx(m2,k1)/2
               cdxy = cdxy + tab_pw2der(m2,k2)*ffx(m2,k1)/2

               cdxx = cdxx + tab_pw2pot(m2,k2)*ffxx(m2,k1)/2
            endif
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
c                 NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw,npw,((npw+1)/2),nd)
      complex *16 tab_pw2pot(npw,n)
      real *8 pot(nd,n,n,n)

      complex *16 ff(n,npw,((npw+1)/2))
      complex *16 ff2(n,n,((npw+1)/2))
      complex *16 cd
c
      npw2=npw/2
      
      do ind = 1,nd
         do m3 = 1,((npw+1)/2)
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
         do m3 = 1,((npw+1)/2)
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
            do m3 = 1,npw2
               cd = cd+tab_pw2pot(m3,k3)*ff2(k1,k2,m3)
            enddo
            m3=((npw+1)/2)
            if (m3.gt.npw2) then
               cd = cd+tab_pw2pot(m3,k3)*ff2(k1,k2,m3)/2
            endif
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
c                 NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
c     tab_pw2pot  precomputed table of 1D conversion 
c     tab_pw2deriv  precomputed table of deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      complex *16 pwexp(npw,npw,((npw+1)/2),nd)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)

      complex *16 ff(n,npw,((npw+1)/2))
      complex *16 ffx(n,npw,((npw+1)/2))
      
      complex *16 ff2(n,n,((npw+1)/2))
      complex *16 ff2x(n,n,((npw+1)/2))
      complex *16 ff2y(n,n,((npw+1)/2))
      complex *16 cd,cdx,cdy,cdz
c
      npw2=npw/2
      
      do ind = 1,nd
        do m3 = 1,((npw+1)/2)
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
        do m3 = 1,((npw+1)/2)
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
           do m3 = 1,npw2
              cd  = cd  + tab_pw2pot(m3,k3)*ff2(k1,k2,m3)
              cdz = cdz + tab_pw2der(m3,k3)*ff2(k1,k2,m3)
              
              cdx = cdx + tab_pw2pot(m3,k3)*ff2x(k1,k2,m3)
              cdy = cdy + tab_pw2pot(m3,k3)*ff2y(k1,k2,m3)
           enddo
           m3=((npw+1)/2)
           if (m3.gt.npw2) then
              cd  = cd  + tab_pw2pot(m3,k3)*ff2(k1,k2,m3)/2
              cdz = cdz + tab_pw2der(m3,k3)*ff2(k1,k2,m3)/2
              
              cdx = cdx + tab_pw2pot(m3,k3)*ff2x(k1,k2,m3)/2
              cdy = cdy + tab_pw2pot(m3,k3)*ff2y(k1,k2,m3)/2
           endif
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
c                 NOTE 3D convention is pwexp(npw,npw,((npw+1)/2))
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
      complex *16 pwexp(npw,npw,((npw+1)/2),nd)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2der(npw,n)
      complex *16 tab_pw2dxx(npw,n)
      real *8 pot(nd,n,n,n)
      real *8 grad(nd,3,n,n,n)
      real *8 hess(nd,6,n,n,n)

      complex *16 ff(n,npw,((npw+1)/2))
      complex *16 ffx(n,npw,((npw+1)/2))
      complex *16 ffxx(n,npw,((npw+1)/2))
      
      complex *16 ff2(n,n,((npw+1)/2))
      complex *16 ff2x(n,n,((npw+1)/2))
      complex *16 ff2xy(n,n,((npw+1)/2))
      complex *16 ff2xx(n,n,((npw+1)/2))
      complex *16 ff2y(n,n,((npw+1)/2))
      complex *16 ff2yy(n,n,((npw+1)/2))
      
      complex *16 cd,cdx,cdy,cdz,cdxx,cdyy,cdzz,cdxy,cdxz,cdyz
c
      npw2=npw/2
      
      do ind = 1,nd
c       transformation in x
        do m3 = 1,((npw+1)/2)
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
        do m3 = 1,((npw+1)/2)
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
           do m3 = 1,npw2
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
           m3=((npw+1)/2)
           if (m3.gt.npw2) then
              cd   = cd  +  tab_pw2pot(m3,k3)*ff2(k1,k2,m3)/2
              cdz  = cdz +  tab_pw2der(m3,k3)*ff2(k1,k2,m3)/2
              cdzz = cdzz + tab_pw2dxx(m3,k3)*ff2(k1,k2,m3)/2

              cdx  = cdx  + tab_pw2pot(m3,k3)*ff2x(k1,k2,m3)/2
              cdxz = cdxz + tab_pw2der(m3,k3)*ff2x(k1,k2,m3)/2

              cdy  = cdy  + tab_pw2pot(m3,k3)*ff2y(k1,k2,m3)/2
              cdyz = cdyz + tab_pw2der(m3,k3)*ff2y(k1,k2,m3)/2

              cdxx = cdxx + tab_pw2pot(m3,k3)*ff2xx(k1,k2,m3)/2
              cdxy = cdxy + tab_pw2pot(m3,k3)*ff2xy(k1,k2,m3)/2
              cdyy = cdyy + tab_pw2pot(m3,k3)*ff2yy(k1,k2,m3)/2
           endif
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
      
      complex *16 wshift(((npw+1)/2)*npw**(ndim-1),(2*nmax+1)**ndim)

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
      
      complex *16 wshift(((npw+1)/2),(2*nmax+1))
      
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
         do j1=1,((npw+1)/2)
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
      
      complex *16 wshift(npw*((npw+1)/2),(2*nmax+1)**2)
      
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
         do j1=1,((npw+1)/2)
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
      
      complex *16 wshift(npw*npw*((npw+1)/2),(2*nmax+1)**3)
      
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
         do j1=1,((npw+1)/2)
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
      complex *16 wshift(((npw+1)/2)*npw**(ndim-1),2**ndim,nmax)
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
      complex *16 wshift(((npw+1)/2),2,nmax)
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
         do j1=1,((npw+1)/2)
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
      complex *16 wshift(npw*((npw+1)/2),4,nmax)
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
         do j1=1,((npw+1)/2)
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
      complex *16 wshift(npw*npw*((npw+1)/2),8,nmax)
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
         do j1=1,((npw+1)/2)
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
