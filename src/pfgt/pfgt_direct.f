c     This file contains nine direction evaluation subroutines
c     for point Gauss transform in 1,2,3 dimensions.
c
c     Naming convention for subroutines:
c     1. all subroutines start with gnd_direct, indicating that 
c        they are for direct evaluation of Gauss transform in general dimensions;
c     2. they are followed by c (charge input only), d (dipole input only),
c        cd (charge+dipole input), p (potential output only), g (potential+gradient output),
c        h (potential+gradient+hessian output).
c
c     pot(ii,i)  = \sum_j charge(ii,j)*exp(-(t_i-s_j)^2/delta)
c                + \sum_j dipstr(ii,j)*(grad_s(exp(-(t_i-s_j)^2/delta) \cdot rnormal(s_j))
c
c     grad  = d(pot)
c     hess  = dd(pot)
c     Note: in 2d, hessian has the order dxx, dxy, dyy
c           in 3d, hessian has the order dxx, dyy, dzz, dxy, dxz, dyz
c     
c     Note: all output variables are incremented. So proper initialization are needed, but
c     the subroutines can be called many times to calculate the output variable due to all
c     sources.
c
c     The union of input and output arguments are as follows.
c
c     Input parameters:
c     nd: number of input vectors (charge, rnormal, dipstr) and output vectors (pot, grad, hess)
c     dim: dimension of the underlying space
c     delta: the Gaussian variance
c     dmax: the squared distance outside which the Gaussian kernel is regarded as 0
c     ns: number of sources
c     sources: (dim,ns) source coordinates
c     ntarg: number of targets
c     targ: (dim,ntarg) target coordinates
c     charge: (nd,ns) charge strengths
c     rnormal: (dim,ns) dipole orientation vectors
c     dipstr: (nd,ns) dipole strengths
c
c     Output parameters:
c     pot: (nd,ntarg) incremented potential at targets
c     grad: (nd,dim,ntarg) incremented gradient at targets
c     hess: (nd,dim*(dim+1)/2,ntarg) incremented hessian at targets
c***********************************************************************
c
c     charge to potential
c
c**********************************************************************
      subroutine gnd_directcp(nd,dim,delta,dmax,sources,ns,charge,
     $           targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),targ(dim,*)
      real *8 dr(dim)
      real *8 rtmp,delta
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns)
c
      do itarg=1,ntarg
      do i = 1,ns
         rr=0
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c***********************************************************************
c
c     charge to potential+gradient
c
c**********************************************************************
      subroutine gnd_directcg(nd,dim,delta,dmax,sources,ns,charge,
     1             targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),targ(dim,ntarg)
      real *8 dr(dim)
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 charge(nd,ns)
c
      rfac=-2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + dr(k)*rinc
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
c***********************************************************************
c
c     charge to potential+gradient+hessian
c
c**********************************************************************
      subroutine gnd_directch(nd,dim,delta,dmax,sources,ns,charge,
     1           targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ifexpon,ii,nd
      integer dim
      real *8 sources(dim,ns),targ(dim,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,delta
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 hess(nd,dim*(dim+1)/2,ntarg)
      real *8 dr(dim)
      real *8 dg(dim),dh(dim*(dim+1)/2)
      real *8 charge(nd,ns)
c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do k=1,dim
               dg(k)=rfac*dr(k)
            enddo

            if (dim.ne.2) then
               do k=1,dim
                  dh(k)=rfac+dg(k)*dg(k)
               enddo

               m=dim
               do k=1,dim-1
                  do j=k+1,dim
                     m=m+1
                     dh(m)=dg(k)*dg(j)
                  enddo
               enddo
            else
               dh(1)=rfac+dg(1)*dg(1)
               dh(2)=dg(1)*dg(2)
               dh(3)=rfac+dg(2)*dg(2)
            endif
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + dg(k)*rinc
               enddo
               do k=1,dim*(dim+1)/2
                  hess(ii,k,itarg) = hess(ii,k,itarg) + dh(k)*rinc
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
c***********************************************************************
c
c     dipole to potential
c      
c**********************************************************************
      subroutine gnd_directdp(nd,dim,delta,dmax,sources,ns,
     $           rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),rnormal(dim,ns),targ(dim,ntarg)
      real *8 dr(dim)
      real *8 pot(nd,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo
            rinc=0
            do k=1,dim
               rinc=rinc+dr(k)*rnormal(k,i)
            enddo
            rinc=rinc*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c***********************************************************************
c
c     dipole to potential+gradient
c      
c**********************************************************************
      subroutine gnd_directdg(nd,dim,delta,dmax,sources,ns,
     1             rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),rnormal(dim,ns),targ(dim,ntarg)
      real *8 dr(dim)
      real *8 rincg(dim),delta
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo

            rinc=0
            do k=1,dim
               rinc=rinc+dr(k)*rnormal(k,i)
            enddo
            do k=1,dim
               rincg(k)=rnormal(k,i)*rfac-rinc*dr(k)
            enddo
c     
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + rincg(k)*rtmp
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
c***********************************************************************
c
c     dipole to potential+gradient+hessian
c      
c**********************************************************************
      subroutine gnd_directdh(nd,dim,delta,dmax,sources,ns,
     1           rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ifexpon,ii,nd
      integer dim
      real *8 sources(dim,ns),rnormal(dim,ns),targ(dim,ntarg)
      real *8 dr(dim),rincg(dim),rinch(dim*(dim+1)/2)
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 hess(nd,dim*(dim+1)/2,ntarg)
      real *8 dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo

            rinc=0
            do k=1,dim
               rinc=rinc+dr(k)*rnormal(k,i)
            enddo
            
            do k=1,dim
               rincg(k)=rnormal(k,i)*rfac-rinc*dr(k)
            enddo

            if (dim.ne.2) then
               do k=1,dim
                  rinch(k)=-rfac*(dr(k)*rnormal(k,i)+rinc)
     1                - rincg(k)*dr(k)
               enddo

               m=dim
               do k=1,dim-1
                  do j=k+1,dim
                     m=m+1
                     rinch(m)=-rfac*dr(k)*rnormal(j,i) - rincg(k)*dr(j)
                  enddo
               enddo
            else
               rinch(1)=-rfac*(dr(1)*rnormal(1,i)+rinc)
     1             - rincg(1)*dr(1)
               rinch(2)=-rfac*dr(1)*rnormal(2,i) - rincg(1)*dr(2)
               rinch(3)=-rfac*(dr(2)*rnormal(2,i)+rinc)
     1             - rincg(2)*dr(2)
            endif
c     
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + rincg(k)*rtmp
               enddo
               do k=1,dim*(dim+1)/2
                  hess(ii,k,itarg) = hess(ii,k,itarg) + rinch(k)*rtmp
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
c***********************************************************************
c
c     charge+dipole to potential
c      
c**********************************************************************
      subroutine gnd_directcdp(nd,dim,delta,dmax,sources,ns,
     $           charge,rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),rnormal(dim,ns),targ(dim,ntarg)
      real *8 dr(dim)
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then         
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
c
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo
            rinc=0
            do k=1,dim
               rinc=rinc+dr(k)*rnormal(k,i)
            enddo
            rinc=rinc*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c***********************************************************************
c
c     charge+dipole to potential+gradient
c      
c**********************************************************************
      subroutine gnd_directcdg(nd,dim,delta,dmax,sources,ns,
     1             charge,rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ii,nd
      integer dim
      real *8 sources(dim,ns),rnormal(dim,ns),targ(dim,ntarg),delta
      real *8 dr(dim),dg(dim),rincg(dim)
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            do k=1,dim
               dr(k)=rfac*dr(k)
            enddo
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + dr(k)*rinc
               enddo
            enddo
c
            rinc=0
            do k=1,dim
               rinc=rinc-dr(k)*rnormal(k,i)
            enddo
            do k=1,dim
               rincg(k)=-rnormal(k,i)*rfac+rinc*dr(k)
            enddo
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + rincg(k)*rtmp
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
c***********************************************************************
c
c     charge+dipole to potential+gradient+hessian
c      
c**********************************************************************
      subroutine gnd_directcdh(nd,dim,delta,dmax,sources,ns,
     1           charge,rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
      integer i,ns,ifexpon,ii,nd
      integer dim
      real *8 sources(dim,ns),targ(dim,ntarg),rnormal(dim,ns)
      real *8 dr(dim),dg(dim),dh(dim*(dim+1)/2)
      real *8 rincg(dim),rinch(dim*(dim+1)/2)
      
      real *8 pot(nd,ntarg),grad(nd,dim,ntarg)
      real *8 hess(nd,dim*(dim+1)/2,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)


c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         do k=1,dim
            dr(k)=targ(k,itarg)-sources(k,i)
         enddo
         rr=0
         do k=1,dim
            rr=rr+dr(k)*dr(k)
         enddo
         if (rr .lt. dmax) then
            dtmp = dexp(-rr/delta)
            do k=1,dim
               dg(k)=rfac*dr(k)
            enddo
            
            if (dim.ne.2) then
               do k=1,dim
                  dh(k)=rfac+dg(k)*dg(k)
               enddo

               m=dim
               do k=1,dim-1
                  do j=k+1,dim
                     m=m+1
                     dh(m)=dg(k)*dg(j)
                  enddo
               enddo
            else
               dh(1)=rfac+dg(1)*dg(1)
               dh(2)=dg(1)*dg(2)
               dh(3)=rfac+dg(2)*dg(2)
            endif
c
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + dg(k)*rinc
               enddo
               do k=1,dim*(dim+1)/2
                  hess(ii,k,itarg) = hess(ii,k,itarg) + dh(k)*rinc
               enddo
            enddo

            rinc=0
            do k=1,dim
               rinc=rinc-dg(k)*rnormal(k,i)
            enddo
            
            do k=1,dim
               rincg(k)=-rnormal(k,i)*rfac+rinc*dg(k)
            enddo

            if (dim.ne.2) then
               do k=1,dim
                  rinch(k)=-rfac*(dg(k)*rnormal(k,i)+rinc)
     1                - rincg(k)*dg(k)
               enddo

               m=dim
               do k=1,dim-1
                  do j=k+1,dim
                     m=m+1
                     rinch(m)=-rfac*dg(k)*rnormal(j,i) - rincg(k)*dg(j)
                  enddo
               enddo
            else
               rinch(1)=-rfac*(dg(1)*rnormal(1,i)+rinc)
     1             - rincg(1)*dg(1)
               rinch(2)=-rfac*dg(1)*rnormal(2,i) - rincg(1)*dg(2)
               rinch(3)=-rfac*(dg(2)*rnormal(2,i)+rinc)
     1             - rincg(2)*dg(2)
            endif
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               do k=1,dim
                  grad(ii,k,itarg) = grad(ii,k,itarg) + rincg(k)*rtmp
               enddo
               do k=1,dim*(dim+1)/2
                  hess(ii,k,itarg) = hess(ii,k,itarg) + rinch(k)*rtmp
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
C
      subroutine pfgt_find_local_shift(ndim,tcenter,scenter,
     1    bs0,shift)
c     returns the center shift of the image cell for the periodic point FGT
c
c     input:
c     ndim - dimension of the underlying space
c     tcenter - target box center
c     scenter - source box center
c     bs0 - root box size
c
c     output
c     shift - source center shift
c
      implicit real *8 (a-h,o-z)
      real *8 tcenter(ndim),scenter(ndim),shift(ndim)


      do i=1,ndim
         dx = tcenter(i)-scenter(i)
         shift(i)=0
         dxp1=dx-bs0
         dxm1=dx+bs0
         if (abs(dx).gt.abs(dxp1)) then
            dx=dxp1
            shift(i)=bs0
         endif
         if (abs(dx).gt.abs(dxm1)) shift(i)=-bs0
      enddo
      
      
      return
      end subroutine
