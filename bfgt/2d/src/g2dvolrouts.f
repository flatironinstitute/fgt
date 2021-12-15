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
      subroutine legtrans2d(nd,norder,fvals,fcoefs,umat)
c
c     converts function values on a 2d Legendre tensor product grid
c     to legendre expansion coefficients.
c
c
c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,norder,norder)
      real *8 fcoefs(norder,norder,nd),umat(norder,norder)
      real *8, allocatable:: fcv(:,:,:)

      allocate(fcv(nd,norder,norder))
c     transform in x
      do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fvals(ind,k1,j)
               enddo
               fcv(ind,k,j)=dd
            enddo
         enddo
      enddo
c     transform in y
      do i=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do i1=1,norder
                  dd=dd+umat(i,i1)*fcv(ind,k,i1)
               enddo
               fcoefs(k,i,ind)=dd
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
C*********************************************************************C
      subroutine leg2d_to_pw(nd,n,coeff,npw,ff,tab_leg2pw,pwexp)
cccc      subroutine leg2d_to_pw(nd,n,coeff,npw,tab_leg2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     Legendre series coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     coeff    Legendre coefficients
c              f = sum coeff(n,m,l) P_n(x) P_m(y) 
c     npw      number of plane waves
c                 NOTE 2D convention is pwexp(npw,npw,npw/2)
c     ff       complex workspace (npw,n)
c     tab_leg2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,n,nd)
      complex *16 ff(n,npw/2),tab_leg2pw(n,npw)
      complex *16 pwexp(npw,npw/2,nd),cd
c
      do ind = 1,nd
         do m1 = 1,n
            do k2 = 1,npw/2
               cd = 0.0d0
               do m2 = 1,n
                  cd = cd+tab_leg2pw(m2,k2)*coeff(m1,m2,ind)
               enddo
               ff(m1,k2) = cd
            enddo
         enddo
c
         do k2 = 1,npw/2
            do k1 = 1,npw
               cd = 0.0d0
               do m1 = 1,n
                  cd = cd+tab_leg2pw(m1,k1)*ff(m1,k2)
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
c
c
c
c
c
C*********************************************************************C
      subroutine g2d_pw2pot(nd,n,npw,pwexp,ff,tab_pw2pot,pot)
C*********************************************************************C
c     This routine computes the potential on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
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
      complex *16 ff(n,npw/2)
      complex *16 pwexp(npw,npw/2,nd),cd
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
C*********************************************************************C
      subroutine leg2d_to_potloc(nd,n,coeff,ff,pot,tabx,taby)
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
c     coeff       Legendre coefficients
c                 f = sum coeff(n,m) P_n(x) P_m(y)
c     ff          workspace
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
      real *8 coeff(n,n,nd),pot(nd,n,n)
      real *8 ff(n,n),tabx(n,n),taby(n,n)
c
      do ind = 1,nd
c        transform in x
         do j2=1,n
            do k1=1,n
               cd=0
               do j1=1,n
                  cd=cd+tabx(j1,k1)*coeff(j1,j2,ind)
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
      
      complex *16 wshift(npw*npw/2,-nmax:nmax,-nmax:nmax)
      
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
         j=0   
         do j1=1,npw/2
         do j2=1,npw
            j=j+1
            wshift(j,k2,k1) = ww(j2,k2)*ww(j1,k1)
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
c           pm
            wshift(j,2,k1) = ww(j2,k1)*conjg(ww(j1,k1))
c           mp
            wshift(j,3,k1) = conjg(ww(j2,k1))*ww(j1,k1)
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
      subroutine g2dshiftpw_vec(nd,nexp,pwexp1,
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
      subroutine g2dshiftpw_loc_vec(nd,nexp,pwexp1,
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
            pwexp2(j,ind) = 
     1          pwexp1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g2dcopypwexp_vec(nd,nexp,pwexp1,
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
      subroutine g2dpwzero_vec(nd,pwexp,npw)
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
      complex *16 pwexp(npw*((npw+1)/2),nd)
c
      do ii=1,nd
         do n=1,npw*((npw+1)/2)
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
c
C*********************************************************************C
      subroutine g2d_pw2pg(nd,n,npw,pwexp,ff,ff2,tab_pw2pot,
     1           tab_pw2deriv,pot,grad)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 2D convention is pwexp(npw,(npw+1)/2)
c     ff       complex workspace (n,(npw+1)/2)
c     ff2       complex workspace (n,(npw+1)/2)
c     tab_pw2pot    precomputed table of 1D conversion 
c     tab_pw2deriv  precomputed table of deriv of 1D conversion 
c
c     OUTPUT:
c     pot      potential values on the tensor grid
c     grad     grad values on the tensor grid
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,n,n)
      real *8 grad(nd,2,n,n)
      complex *16 tab_pw2pot(npw,n)
      complex *16 tab_pw2deriv(npw,n)
ccc      complex *16 ff(n,(npw+1)/2)
ccc      complex *16 ff2(n,(npw+1)/2)
      complex *16 ff(npw/2,n)
      complex *16 ff2(npw/2,n)
      complex *16 pwexp(npw,npw/2,nd),cd,cdx,cdy
c
ccc      call prin2(' tab_pwpot *',tab_pw2pot,2*npw*n)
ccc      call prin2(' tab_pwderiv *',tab_pw2deriv,2*npw*n)
      do ind = 1,nd
         do m2 = 1,npw/2
            do k1 = 1,n
               cd=0.0d0
               cdx=0.0d0
               do m1 = 1,npw
                  cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,m2,ind)
                  cdx = cdx+tab_pw2deriv(m1,k1)*pwexp(m1,m2,ind)
               enddo
               ff(m2,k1) = cd
               ff2(m2,k1) = cdx
            enddo
         enddo
c
         do k2 = 1,n
            do k1 = 1,n
               cd = 0.0d0
               cdx = 0.0d0
               cdy = 0.0d0
               do m2 = 1,npw/2
                  cd = cd+tab_pw2pot(m2,k2)*ff(m2,k1)
                  cdx = cdx+tab_pw2pot(m2,k2)*ff2(m2,k1)
                  cdy = cdy+tab_pw2deriv(m2,k2)*ff(m2,k1)
               enddo
               pot(ind,k1,k2)=pot(ind,k1,k2)+dreal(cd)*2
               grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+dreal(cdx)
               grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+dreal(cdy)
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
      subroutine g2d_pw2pgh(nd,n,npw,pwexp,ff,ff2,ff3,tab_pw2pot,
     1           tab_pw2deriv,tab_pw2dxx,pot,grad,hess)
C*********************************************************************C
c     This routine computes the potential and gradient on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
c     npw      number of plane waves
c     pwexp    plane wave expansion
c                 NOTE 2D convention is pwexp(npw,(npw+1)/2)
c     ff       complex workspace (n,(npw+1)/2)
c     ff2       complex workspace (n,(npw+1)/2)
c     ff3       complex workspace (n,(npw+1)/2)
c     tab_pw2pot    precomputed table of 1D conversion 
c     tab_pw2deriv  precomputed table of deriv of 1D conversion 
c     tab_pw2dxx  precomputed table of second deriv of 1D conversion 
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
      complex *16 tab_pw2deriv(npw,n)
      complex *16 tab_pw2dxx(npw,n)
ccc      complex *16 ff(n,(npw+1)/2)
ccc      complex *16 ff2(n,(npw+1)/2)
      complex *16 ff(npw/2,n)
      complex *16 ff2(npw/2,n)
      complex *16 ff3(npw/2,n)
      complex *16 pwexp(npw,npw/2,nd),cd,cdx,cdy
      complex *16 cdxx,cdxy,cdyy
c
ccc      call prin2(' tab_pwpot *',tab_pw2pot,2*npw*n)
ccc      call prin2(' tab_pwderiv *',tab_pw2deriv,2*npw*n)
      do ind = 1,nd
         do m2 = 1,npw/2
            do k1 = 1,n
               cd=0.0d0
               cdx=0.0d0
               cdxx=0.0d0
               do m1 = 1,npw
                  cd = cd+tab_pw2pot(m1,k1)*pwexp(m1,m2,ind)
                  cdx = cdx+tab_pw2deriv(m1,k1)*pwexp(m1,m2,ind)
                  cdxx = cdxx+tab_pw2dxx(m1,k1)*pwexp(m1,m2,ind)
               enddo
               ff(m2,k1) = cd
               ff2(m2,k1) = cdx
               ff3(m2,k1) = cdxx
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
                  cd = cd+tab_pw2pot(m2,k2)*ff(m2,k1)
                  cdx = cdx+tab_pw2pot(m2,k2)*ff2(m2,k1)
                  cdxx = cdxx+tab_pw2pot(m2,k2)*ff3(m2,k1)
                  cdy = cdy+tab_pw2deriv(m2,k2)*ff(m2,k1)
                  cdyy = cdyy+tab_pw2dxx(m2,k2)*ff(m2,k1)
                  cdxy = cdxy+tab_pw2deriv(m2,k2)*ff2(m2,k1)
               enddo
               pot(ind,k1,k2)=pot(ind,k1,k2)+dreal(cd)*2
               grad(ind,1,k1,k2)=grad(ind,1,k1,k2)+dreal(cdx)
               grad(ind,2,k1,k2)=grad(ind,2,k1,k2)+dreal(cdy)
               hess(ind,1,k1,k2)=hess(ind,1,k1,k2)+dreal(cdxx)
               hess(ind,2,k1,k2)=hess(ind,2,k1,k2)+dreal(cdxy)
               hess(ind,3,k1,k2)=hess(ind,3,k1,k2)+dreal(cdyy)
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
      real *8 coeff(n,n,nd)
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
                  cd=cd+tabf(j1,k1)*coeff(j1,j2,ind)
                  cdx=cdx+tabfx(j1,k1)*coeff(j1,j2,ind)
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
      real *8 coeff(n,n,nd)
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
                  cd=cd+tabf(j1,k1)*coeff(j1,j2,ind)
                  cdx=cdx+tabfx(j1,k1)*coeff(j1,j2,ind)
                  cdxx=cdxx+tabfxx(j1,k1)*coeff(j1,j2,ind)
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
      subroutine mklegtables(norder,xnodes,tabf,tabfx,tabfxx)
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

