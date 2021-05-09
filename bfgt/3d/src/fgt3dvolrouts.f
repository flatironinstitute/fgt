C*********************************************************************C
      subroutine legtrans3d(nd,n,fdat,f,texp,umat,coeff)
C*********************************************************************C
c     3d Legendre transformation of data on tensor product
c     Legendre grid
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        grid size (and expansion order)
c     fdat     tensor product function values
c              fdat(i,j,k) = f(x_i,y_j,z_k)
c     f, texp  work arrays of length n
c     umat     1D transform matrix from legeexps
c
c     OUTPUT:
c     coeff    Legendre coefficients
c              f = sum coeff(n,m,l) P_n(x) P_m(y) P_l(z)
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 fdat(n,n,n,nd),coeff(n,n,n,nd),umat(n,n)
      real *8 f(n),texp(n)
c
c     transform rows
c
      do ind = 1,nd
         do j3 = 1,n
         do j2 = 1,n
            do j1 = 1,n
               f(j1) = fdat(j1,j2,j3,ind)
               texp(j1) = 0.0d0
            enddo
            do k = 1,n
            do i = 1,n
               texp(i) = texp(i) + umat(i,k)*f(k)
            enddo
            enddo
            do j1 = 1,n
               coeff(j1,j2,j3,ind) = texp(j1)
            enddo
         enddo
         enddo
c
c     transform columns
c
         do j3 = 1,n
         do j1 = 1,n
            do j2 = 1,n
               f(j2) = coeff(j1,j2,j3,ind)
               texp(j2) = 0.0d0
            enddo
            do k = 1,n
            do i = 1,n
               texp(i) = texp(i) + umat(i,k)*f(k)
            enddo
            enddo
            do j2 = 1,n
               coeff(j1,j2,j3,ind) = texp(j2)
            enddo
         enddo
         enddo
c
c     transform planes
c
         do j2 = 1,n
         do j1 = 1,n
            do j3 = 1,n
               f(j3) = coeff(j1,j2,j3,ind)
               texp(j3) = 0.0d0
            enddo
            do k = 1,n
            do i = 1,n
               texp(i) = texp(i) + umat(i,k)*f(k)
            enddo
            enddo
            do j3 = 1,n
               coeff(j1,j2,j3,ind) = texp(j3)
            enddo
         enddo
         enddo
      enddo
      return
      end subroutine
c
c
c
C*********************************************************************C
      subroutine leg3deval(n,coeff,xx,yy,zz,val,w1,w2,w3)
C*********************************************************************C
c
c     3d Legendre series evaluation at single point.
c
c     INPUT:
c     n          grid size (and expansion order)
c     coeff      Legendre coefficients
c                f = sum coeff(n,m,l) P_n(x) P_m(y) P_l(z)
c     xx,yy,zz   target coordinates
c     w1,w2,w3   work arrays of length n
c
c     OUTPUT:
c     val      series value at (xx,yy,zz)
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,n)
      real *8 f,w1(n),w2(n),w3(n)
c
      do k = 1,n
         do j = 1,n
            do i = 1,n
               w1(i) = coeff(i,j,k)
            enddo
            call legeexev(xx,val,w1,n)
            w2(j) = val
         enddo
         call legeexev(yy,val,w2,n)
         w3(k) = val
      enddo
      call legeexev(zz,val,w3,n)
c
      return
      end subroutine
c
c
c
c
c
C*********************************************************************C
      subroutine leg3d_to_pw(nd,n,coeff,npw,ff,ff2,tab_leg2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     Legendre series coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     coeff    Legendre coefficients
c              f = sum coeff(n,m,l) P_n(x) P_m(y) P_l(z)
c     npw      number of plane waves
c                 NOTE 3D convention is pwexp(npw,npw,npw/2)
c     ff       complex workspace (npw,n,n)
c     ff2      complex workspace (npw,npw,n)
c     tab_leg2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,n,nd)
      complex *16 ff(npw,n,n),tab_leg2pw(n,npw)
      complex *16 ff2(npw,npw,n)
      complex *16 pwexp(npw,npw,npw/2,nd),cd
c
      do ind = 1,nd
         do m3 = 1,n
         do m2 = 1,n
            do k1 = 1,npw
               cd = 0.0d0
               do m1 = 1,n
                  cd = cd+tab_leg2pw(m1,k1)*coeff(m1,m2,m3,ind)
               enddo
               ff(k1,m2,m3) = cd
            enddo
         enddo
         enddo
c
         do m3 = 1,n
         do k2 = 1,npw
            do k1 = 1,npw
               cd = 0.0d0
               do m2 = 1,n
                  cd = cd+tab_leg2pw(m2,k2)*ff(k1,m2,m3)
               enddo
               ff2(k1,k2,m3) = cd
            enddo
         enddo
         enddo
c
         do k3 = 1,npw/2
         do k2 = 1,npw
         do k1 = 1,npw
               cd = 0.0d0
               do m3 = 1,n
                  cd = cd+tab_leg2pw(m3,k3)*ff2(k1,k2,m3)
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
c
c
c
c
C*********************************************************************C
      subroutine g3d_pw2pot(nd,n,npw,pwexp,ff,ff2,tab_pw2pot,pot)
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
c     ff       complex workspace (n,npw,npw/2)
c     ff2      complex workspace (n,n,npw/2)
c     tab_pw2pot  precomputed table of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 ff(n,npw,npw/2)
      complex *16 ff2(n,n,npw/2)
      complex *16 pwexp(npw,npw,npw/2,nd),cd
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
