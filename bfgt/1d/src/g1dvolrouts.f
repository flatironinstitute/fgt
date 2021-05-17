c     This file contains a set of subroutines that carry out transforms
c     for 1D box FGT.
c
c     legtrans1d: converts function values on a Legendre tensor product grid
c                 to Legendre expansion coefficients.
c
c     leg1d_to_pw: converts 1D Legendre expansion coefficients to 1D planewave
c                 expansion coefficients.
c
c     g1d_pw2pot: converts 1D planewave expansion coefficients to potential values
c                 on a Legendre tensor product grid in the same box.
c
c     leg1d_to_potloc: converts 1D Legendre expansion coefficients in the source box 
c                 to potential values on a Legendre tensor product grid in 
c                 a target box that is in the list1 of the source box.
c
c     pw_translation_matrices: precomputes all translation matrices for 1D PW
C                 expansions for mp to loc at the cutoff level.
c
c     merge_split_pw_matrices: precomputes all translation matrices for 
c                 1D PW translations from child to parent or vice versa at all 
c                 needed levels.
c
c     g1dshiftpw_vec: translates a 1D PW expansion (pwexp1) about
C                 the center (CENT1) into another PW expansion (pwexp2) about 
C                 (CENT2) using precomputed translation matrix wshift.
c
c     g1dcopypwexp_vec: adds one 1D PW expansion (pwexp1) 
C                 to another 1D PW expansion (pwexp2).
c
c     g1dpwzero_vec: sets a vector 1D planewave expansions to zero.
c
c
      subroutine legtrans1d(nd,norder,fvals,fcoefs,umat)
c
c     converts function values on a 1d Legendre tensor product grid
c     to legendre expansion coefficients.
c
      implicit real *8 (a-h,o-z)
      real *8 fvals(nd,norder)
      real *8 fcoefs(norder,nd),umat(norder,norder)

c     transform in x
      do ind=1,nd
         do k=1,norder
            dd=0
            do k1=1,norder
               dd=dd+umat(k,k1)*fvals(ind,k1)
            enddo
            fcoefs(k,ind)=dd
         enddo
      enddo

      return
      end
c
c
c
c
C*********************************************************************C
      subroutine leg1d_to_pw(nd,n,coeff,npw,tab_leg2pw,pwexp)
C*********************************************************************C
c     This routine computes the plane wave expansion from
c     Legendre series coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        dimension of coeff array
c     coeff    Legendre coefficients
c              f = sum coeff(n) P_n(x)
c     npw      number of plane waves
c                 NOTE 1D convention is pwexp(npw/2)
c     tab_leg2pw  precomputed table of 1D conversion factors
c                 (n,j) entry is: ws(j)*(D/2)* 
c                 int_{-1}^1 P_n(x)exp(- i ts(j)Dx/(2 \sqrt{delta}))dx
c                 where D is box dimension at current level in
c                 tree hierarchy.
c     OUTPUT:
c     pwexp    plane wave expansion
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
      real *8 coeff(n,nd)
      complex *16 tab_leg2pw(n,npw)
      complex *16 pwexp(npw/2,nd),cd
c
      do ind = 1,nd
         do k1 = 1,npw/2
            cd = 0.0d0
            do m1 = 1,n
               cd = cd+tab_leg2pw(m1,k1)*coeff(m1,ind)
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
c
c
C*********************************************************************C
      subroutine g1d_pw2pot(nd,n,npw,pwexp,tab_pw2pot,pot)
C*********************************************************************C
c     This routine computes the potential on a tensor grid from 
c     the plane wave expansion coefficients (implicitly about the box center).
c
c     INPUT:
c     nd       vector length (for multiple RHS)
c     n        number of Legendre nodes
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
c
C*********************************************************************C
      subroutine leg1d_to_potloc(nd,n,coeff,pot,tabx)
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
      real *8 coeff(n,nd),pot(nd,n)
      real *8 tabx(n,n)
c
      do ind = 1,nd
c        transform in x
         do k1=1,n
            cd=0
            do j1=1,n
               cd=cd+tabx(j1,k1)*coeff(j1,ind)
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
      
      complex *16 wshift(npw/2,-nmax:nmax)
      
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
         do j1=1,npw/2
            wshift(j1,k1) = ww(j1,k1)
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
      subroutine g1dshiftpw_vec(nd,nexp,pwexp1,
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
      subroutine g1dcopypwexp_vec(nd,nexp,pwexp1,
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
            pwexp2(j,ind) = pwexp2(j,ind)+pwexp1(j,ind)
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
      subroutine g1dpwzero_vec(nd,pwexp,npw)
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
      complex *16 pwexp(npw/2,nd)
c
      do ii=1,nd
         do n=1,npw/2
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
