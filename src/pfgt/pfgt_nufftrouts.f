c     This file contains the following subroutines
c
c      
c     gnd_formpw : computes the PW expansion due to charges and/or dipoles
c                      type 1 NUFFTs are called 
c      
c     gnd_pweval : evaluates the PW expansion
c                      potential + gradient + hessian
c                      type 2 NUFFTs are called
c      
c     gnd_mk_full_translation_matrices : returns precomputed translation matrices
c                   for PW mp to loc translations on the cutoff level
c      
c     nufft_weights : weights for type 1 and type 2 NUFFT calls
c      
c*********************************************************************
C
C     form PW expansions (charge, dipole, charge & dipole) using NUFFT
C
C*********************************************************************
      subroutine gnd_formpw(nd,dim,delta,eps,sources,ns,
     1    ifcharge,charge,ifdipole,rnormal,dipstr,cent,
     2    hpw,nexp,wnufft,ffexp,fftplan)
C
C     This subroutine computes the PW expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
c     dim           = dimension of the underlying space
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     sources       = source locations
C     ns            = number of sources
c     ifcharge      = whether charge is present
C     charge        = strengths of sources
c     ifdipole      = whether dipole is present
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
C     cent          = center of the expansion
C     hpw           = step size in the Fourier space
C     nexp          = total number of terms in the plane-wave expansion
C     wnufft        = real *8 weights for nufft, tensor product of 1d weights
C     fftplan       = integer *8 pointer to the fftw plan for finufft calls
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,dim
      real *8 cent(dim),sources(dim,ns),dipstr(nd,ns)
      real *8 rnormal(dim,ns),charge(nd,ns)
      real *8 wnufft(nexp,dim+1)
      complex *16 ffexp(nexp,nd)
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)

c     to pass null pointers to unused arguments...
      real *8, pointer :: dummy => null()

c     this is what you use as the "opaque" ptr to ptr to finufft_plan...
      integer *8 fftplan
      integer *8 ns8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      if (ifcharge.eq.1 .and. ifdipole.eq.0) then
         kstart=1
         kend=1
      elseif (ifcharge.eq.0 .and. ifdipole.eq.1) then
         kstart=2
         kend=dim+1
      elseif (ifcharge.eq.1 .and. ifdipole.eq.1) then
         kstart=1
         kend=dim+1
      else
         return
      endif

      allocate(xj(ns))
      allocate(cj(ns,kstart:kend,nd))
      allocate(fk(nexp,kstart:kend,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = hpw/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
      enddo

      if (dim.ge.2) then
         allocate(yj(ns))
         do j=1,ns
            yj(j) = (sources(2,j) - cent(2))*dsq
         enddo
      endif

      if (dim.ge.3) then
         allocate(zj(ns))
         do j=1,ns
            zj(j) = (sources(3,j) - cent(3))*dsq
         enddo
      endif

      if (ifcharge.eq.1 .and. ifdipole.eq.0) then
         do ind = 1,nd
         do j=1,ns
            cj(j,kstart,ind) = charge(ind,j)
         enddo
         enddo
      elseif (ifcharge.eq.0 .and. ifdipole.eq.1) then
         do ind = 1,nd
         do j=1,ns
            ztmp = -eye*dipstr(ind,j)*dsq0
            do k=1,dim
               cj(j,k+1,ind) = ztmp*rnormal(k,j)
            enddo
         enddo
         enddo
      elseif (ifcharge.eq.1 .and. ifdipole.eq.1) then
         do ind = 1,nd
         do j=1,ns
            cj(j,kstart,ind) = charge(ind,j)
            ztmp = -eye*dipstr(ind,j)*dsq0
            do k=1,dim
               cj(j,k+1,ind) = ztmp*rnormal(k,j)
            enddo
         enddo
         enddo
      endif
      
      ns8=ns
      if (dim.eq.1) then
         call finufft_setpts(fftplan,ns8,xj,dummy,dummy,dummy,dummy,
     1       dummy,dummy,ier)
      elseif (dim.eq.2) then
         call finufft_setpts(fftplan,ns8,xj,yj,dummy,dummy,dummy,
     1       dummy,dummy,ier)
      elseif (dim.eq.3) then
         call finufft_setpts(fftplan,ns8,xj,yj,zj,dummy,dummy,
     1       dummy,dummy,ier)
      endif

      call finufft_execute(fftplan,cj,fk,ier)

      do ind=1,nd
         do j=1,nexp
            ztmp=0
            do k=kstart,kend
               ztmp=ztmp+wnufft(j,k)*fk(j,k,ind)
            enddo
            ffexp(j,ind)=ztmp
         enddo
      enddo

      return
      end
C*********************************************************************
c
c     evaluate plane-wave expansions
C
C*********************************************************************
      subroutine gnd_pweval(nd,dim,delta,eps,center,hpw,
     1    nexp,wnufftgh,pwexp,targ,nt,ifpgh,pot,grad,hess,fftplan)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient + hess
C
C     INPUT:
c     nd            = vector length (for multiple charges at same locations)
c     dim           = dimension of the underlying space
C     delta         = Gaussian variance
c     eps           = prescribed precision
C     center        = center of the expansion
C     hpw           = step size in the Fourier space
C     nexp          = number of terms in the plane-wave expansion
c     wnufftgh      = weights for computing gradient and hessian
C     pwexp         = pw expansions 
C     targ          = target locations
C     nt            = number of targets
c     ifpgh         = output flag: 1-> pot; 2-> pot+grad; 3-> pot+grad+hess
C     fftplan       = integer *8 pointer to the fftw plan for finufft calls
c      
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = hessian (or vectorized hessians) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,dim
      real *8 delta,center(dim),targ(dim,nt)
      real *8 pot(nd,nt),grad(nd,dim,*),hess(nd,dim*(dim+1)/2,*)
      real *8 wnufftgh(nexp,dim+dim*(dim+1)/2)
      complex *16 pwexp(nexp,nd)

      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:)

c     to pass null pointers to unused arguments...
      real *8, pointer :: dummy => null()

c     this is what you use as the "opaque" ptr to ptr to finufft_plan...
      integer *8 fftplan 
      integer*8 nt8
      complex *16 z,z2
      complex *16 eye
C
      eye = dcmplx(0,1)
      
      dsq = hpw/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      allocate(xj(nt))
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
      enddo

      if (dim.ge.2) then
         allocate(yj(nt))
         do j=1,nt
            yj(j) = (targ(2,j) - center(2))*dsq
         enddo
      endif

      if (dim.ge.3) then
         allocate(zj(nt))
         do j=1,nt
            zj(j) = (targ(3,j) - center(3))*dsq
         enddo
      endif

      z = eye/dsqrt(delta)
      z2= z*z

      kstart=0
      kend=0

      if (ifpgh.eq.1) kend=0
      if (ifpgh.eq.2) kend=dim
      if (ifpgh.eq.3) kend=dim+dim*(dim+1)/2
      allocate(cj(nt,nd,kstart:kend))

      if (ifpgh.eq.1) goto 1200
      allocate(fk(nexp,nd,kstart:kend))
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,0,ind)=pwexp(j,ind)
         enddo
      enddo

      if (ifpgh.ge.2) then
         do k=1,dim
         do ind = 1,nd
         do j=1,nexp
            fk(j,ind,k)=pwexp(j,ind)*wnufftgh(j,k)*z
         enddo
         enddo
         enddo
      endif
      
      if (ifpgh.eq.3) then
         do k=dim+1,dim+dim*(dim+1)/2
         do ind = 1,nd
         do j=1,nexp
            fk(j,ind,k)=pwexp(j,ind)*wnufftgh(j,k)*z2
         enddo
         enddo
         enddo
      endif

 1200 continue
      
      nt8=nt
      if (dim.eq.1) then
         call finufft_setpts(fftplan,nt8,xj,dummy,dummy,dummy,dummy,
     1       dummy,dummy,ier)
      elseif (dim.eq.2) then
         call finufft_setpts(fftplan,nt8,xj,yj,dummy,dummy,dummy,
     1       dummy,dummy,ier)
      elseif (dim.eq.3) then
         call finufft_setpts(fftplan,nt8,xj,yj,zj,dummy,dummy,
     1       dummy,dummy,ier)
      endif
      
      if (ifpgh.eq.1) then
         call finufft_execute(fftplan,cj,pwexp,ier)
      else
         call finufft_execute(fftplan,cj,fk,ier)
      endif
      
      do ind=1,nd
      do j=1,nt
         pot(ind,j)=pot(ind,j)+dreal(cj(j,ind,0))
      enddo
      enddo

      if (ifpgh.ge.2) then
         do ind=1,nd
         do j=1,nt
         do k=1,dim      
            grad(ind,k,j)=grad(ind,k,j)+dreal(cj(j,ind,k))
         enddo
         enddo
         enddo
      endif

      if (ifpgh.eq.3) then
         do ind=1,nd
         do j=1,nt
         do k=1,dim*(dim+1)/2   
            hess(ind,k,j)=hess(ind,k,j)+dreal(cj(j,ind,k+dim))
         enddo
         enddo
         enddo
      endif
c     
      return
      end
C
C
c
c************************************************************************
C
C     make plane-wave translation matrices at the cutoff level (mp to loc)
C
C************************************************************************
      subroutine gnd_mk_full_translation_matrices(dim,xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level for the point FGT
c      
C     INPUT
C     dim     = dimension of the underlying space
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c               always equal to 1 in the current implementation
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      integer dim
      real *8 ts(npw)
      
      complex *16 wshift(npw**dim,(2*nmax+1)**dim)
      
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

      if (dim.eq.1) then
         k=0
         do k1=-nmax,nmax
            k=k+1
            j=0   
            do j1=1,npw
               j=j+1
               wshift(j,k) = ww(j1,k1)
            enddo
         enddo
      elseif (dim.eq.2) then
         k=0
         do k1=-nmax,nmax
         do k2=-nmax,nmax
            k=k+1
            j=0   
            do j1=1,npw
            do j2=1,npw
               j=j+1
               wshift(j,k) = ww(j2,k2)*ww(j1,k1)
            enddo
            enddo
         enddo
         enddo
      elseif (dim.eq.3) then
         k=0
         do k1=-nmax,nmax
         do k2=-nmax,nmax
         do k3=-nmax,nmax
            k=k+1
            j=0   
            do j1=1,npw
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
      endif
c
      return
      end
c
c
c     
c
c*********************************************************************
C
C     weights for nufft calls
C
C*********************************************************************
      subroutine nufft_weights(dim,npw,ws,ts,
     1    nexp,wnufftcd,wnufftgh)
C
C     This subroutine precomputes nufft weights for "form mp" and 
C     "eval loc" stages
C
C     INPUT
C
c     dim     = dimension of the underlying space
C     npw     = number of terms in plane wave exp
c     ws,ts   = weights and nodes of the 1d plane-wave expansion
C     nexp    = number of terms in the full plane-wave expansion
c      
C     OUTPUT:
C
C     wnufftcd - weights for the "form mp" stage
C     wnufftgh - weights for computing grad and hessian
c
      implicit real *8 (a-h,o-z)
      integer dim
      real *8 ws(npw),ts(npw)
      
      real *8 wnufftcd(nexp,dim+1),wnufftgh(nexp,dim+dim*(dim+1)/2)
      
      j=0
      if (dim.eq.1) then
         do j1=1,npw
            j=j+1
            wnufftcd(j,1) = ws(j1)
            wnufftcd(j,2) = ws(j1)*ts(j1)
         enddo
      elseif (dim.eq.2) then
         do j2=1,npw
         do j1=1,npw
            j=j+1
            wnufftcd(j,1) = ws(j1)*ws(j2)
            wnufftcd(j,2) = wnufftcd(j,1)*ts(j1)
            wnufftcd(j,3) = wnufftcd(j,1)*ts(j2)
         enddo
         enddo
      elseif (dim.eq.3) then
         do j3=1,npw
         do j2=1,npw
         do j1=1,npw
            j=j+1
            wnufftcd(j,1) = ws(j1)*ws(j2)*ws(j3)
            wnufftcd(j,2) = wnufftcd(j,1)*ts(j1)
            wnufftcd(j,3) = wnufftcd(j,1)*ts(j2)
            wnufftcd(j,4) = wnufftcd(j,1)*ts(j3)
         enddo
         enddo
         enddo
      endif
c
      
      j=0
      if (dim.eq.1) then
         do j1=1,npw
            j=j+1
            wnufftgh(j,1) = ts(j1)
            wnufftgh(j,2) = ts(j1)*ts(j1)
         enddo
      elseif (dim.eq.2) then
         do j2=1,npw
         do j1=1,npw
            j=j+1
            wnufftgh(j,1) = ts(j1)
            wnufftgh(j,2) = ts(j2)
            
            wnufftgh(j,3) = ts(j1)*ts(j1)
            wnufftgh(j,4) = ts(j1)*ts(j2)
            wnufftgh(j,5) = ts(j2)*ts(j2)
         enddo
         enddo
      elseif (dim.eq.3) then
         do j3=1,npw
         do j2=1,npw
         do j1=1,npw
            j=j+1
            wnufftgh(j,1) = ts(j1)
            wnufftgh(j,2) = ts(j2)
            wnufftgh(j,3) = ts(j3)

            wnufftgh(j,4) = ts(j1)*ts(j1)
            wnufftgh(j,5) = ts(j2)*ts(j2)
            wnufftgh(j,6) = ts(j3)*ts(j3)

            wnufftgh(j,7) = ts(j1)*ts(j2)
            wnufftgh(j,8) = ts(j1)*ts(j3)
            wnufftgh(j,9) = ts(j2)*ts(j3)
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
