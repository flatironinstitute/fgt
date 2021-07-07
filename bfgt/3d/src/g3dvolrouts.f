c     This file contains a set of subroutines that carry out transforms
c     for 3D box FGT.
c
c     legtrans3d: converts function values on a Legendre tensor product grid
c                 to Legendre expansion coefficients.
c
c     leg3d_to_pw: converts 3D Legendre expansion coefficients to 3D planewave
c                 expansion coefficients.
c
c     g3d_pw2pot: converts 3D planewave expansion coefficients to potential values
c                 on a Legendre tensor product grid in the same box.
c
c     leg3d_to_potloc: converts 3D Legendre expansion coefficients in the source box 
c                 to potential values on a Legendre tensor product grid in 
c                 a target box that is in the list1 of the source box.
c
c     pw_translation_matrices: precomputes all translation matrices for 3D PW
C                 expansions for mp to loc at the cutoff level.
c
c     merge_split_pw_matrices: precomputes all translation matrices for 
c                 3D PW translations from child to parent or vice versa at all 
c                 needed levels.
c
c     g3dshiftpw_vec: translates a 3D PW expansion (pwexp1) about
C                 the center (CENT1) into another PW expansion (pwexp2) about 
C                 (CENT2) using precomputed translation matrix wshift.
c
c     g3dcopypwexp_vec: adds one 3D PW expansion (pwexp1) 
C                 to another 3D PW expansion (pwexp2).
c
c     g3dpwzero_vec: sets a vector 3D planewave expansions to zero.
c
c
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
      complex *16 ff(n,n,npw/2),tab_leg2pw(n,npw)
      complex *16 ff2(n,npw,npw/2)
      complex *16 pwexp(npw,npw,npw/2,nd),cd
c
      do ind = 1,nd
         do k3 = 1,npw/2
         do m2 = 1,n
         do m1 = 1,n
            cd = 0.0d0
            do m3 = 1,n
               cd = cd+tab_leg2pw(m3,k3)*coeff(m1,m2,m3,ind)
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
               cd = cd+tab_leg2pw(m2,k2)*ff(m1,m2,k3)
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
               cd = cd+tab_leg2pw(m1,k1)*ff2(m1,k2,k3)
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
C*********************************************************************C
      subroutine leg3d_to_potloc(nd,n,coeff,ff,ff2,pot,tabx,taby,tabz)
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
c     coeff       Legendre coefficients
c                 f = sum coeff(n,m,k) P_n(x) P_m(y) P_k(z)
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
      real *8 coeff(n,n,n,nd),pot(nd,n,n,n)
      real *8 ff(n,n,n),ff2(n,n,n),tabx(n,n),taby(n,n),tabz(n,n)
c
      do ind = 1,nd
c        transform in x
         do j3=1,n
            do j2=1,n
               do k1=1,n
                  cd=0
                  do j1=1,n
                     cd=cd+tabx(j1,k1)*coeff(j1,j2,j3,ind)
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
C
c
C
c*********************************************************************
C
C shift PW expansions (mp to loc, mp to mp, loc to loc)
C
C*********************************************************************
      subroutine pw_translation_matrices(xmin,npw,ts,nmax,
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
      
      complex *16 wshift(npw*npw*npw/2,-nmax:nmax,-nmax:nmax,-nmax:nmax)
      
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
      
      do k1=-nmax,nmax
      do k2=-nmax,nmax
      do k3=-nmax,nmax
         j=0   
         do j1=1,npw/2
         do j2=1,npw
         do j3=1,npw
            j=j+1
            wshift(j,k3,k2,k1) = ww(j3,k3)*ww(j2,k2)*ww(j1,k1)
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
      subroutine merge_split_pw_matrices(xmin,npw,ts,nmax,
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
c           pp p              
            wshift(j,1,k1) = ww(j3,k1)*ww(j2,k1)
     1          *ww(j1,k1)
c           pm p
            wshift(j,2,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *ww(j1,k1)
c           mp p
            wshift(j,3,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *ww(j1,k1)
c           mm p
            wshift(j,4,k1) = conjg(ww(j3,k1))*conjg(ww(j2,k1))
     1          *ww(j1,k1)
            
c           pp m              
            wshift(j,5,k1) = ww(j3,k1)*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           pm m
            wshift(j,6,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *conjg(ww(j1,k1))
c           mp m
            wshift(j,7,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           mm m
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
      subroutine g3dshiftpw_vec(nd,nexp,pwexp1,
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
            pwexp2(j,ind) = pwexp2(j,ind)
     1          +pwexp1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g3dshiftpw_loc_vec(nd,nexp,pwexp1,
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
C     Note: there is no incrementation here!
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
      subroutine g3dcopypwexp_vec(nd,nexp,pwexp1,
     1              pwexp2)
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
      subroutine g3dpwzero_vec(nd,pwexp,npw)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector planewave expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     npw    :   number of terms in 1d planewave expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     pwexp  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,npw,nd,ii
      complex *16 pwexp(npw*npw*npw/2,nd)
c
      do ii=1,nd
         do n=1,npw*npw*npw/2
            pwexp(n,ii)=0.0d0
         enddo
      enddo
      
      return
      end
c      
c      
c      
c      
c
