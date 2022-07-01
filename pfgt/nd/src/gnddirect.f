c
c**********************************************************************
      subroutine gnd_directcp(nd,dim,delta,dmax,sources,ns,charge,
     $           targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(3,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
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
            rr=rr+ (targ(k,itarg)-sources(k,i))**2
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
c
c
c**********************************************************************
      subroutine gnd_directcg(nd,dim,delta,dmax,sources,ns,charge,
     1             targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target

c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)      : gradient is incremented
c---------------------------------------------------------------------
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
c
c
c
c
c
c**********************************************************************
      subroutine gnd_directch(nd,dim,delta,dmax,sources,ns,charge,
     1           targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c     grad = grad
c     hess = (dxx,dxy,dyy)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c     hess(nd,3)    : Hessian is incremented
c                     using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
      
c---------------------------------------------------------------------
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
c
c
c**********************************************************************
      subroutine gnd_directdp(nd,dim,delta,dmax,sources,ns,
     $           rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
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
c
c
c**********************************************************************
      subroutine gnd_directdg(nd,dim,delta,dmax,sources,ns,
     1             rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c---------------------------------------------------------------------
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
c
c
c
c
c
c**********************************************************************
      subroutine gnd_directdh(nd,dim,delta,dmax,sources,ns,
     1           rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c     hess(nd)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
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

            do k=1,dim
               rinch(k)=-rfac*(dr(k)*rnormal(k,i)+rinc)
     1             - rincg(k)*dr(k)
            enddo

            m=dim
            do k=1,dim-1
            do j=k+1,dim
               m=m+1
               rinch(m)=-rfac*dr(k)*rnormal(j,i) - rincg(k)*dr(j)
            enddo
            enddo
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
c
c
c**********************************************************************
      subroutine gnd_directcdp(nd,dim,delta,dmax,sources,ns,
     $           charge,rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges and dipoles at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
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
c
c
c**********************************************************************
      subroutine gnd_directcdg(nd,dim,delta,dmax,sources,ns,
     1             charge,rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
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
c
c
c
c
c
c**********************************************************************
      subroutine gnd_directcdh(nd,dim,delta,dmax,sources,ns,
     1           charge,rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd)      : gradient is incremented
c     hess(nd)      : Hessian is incremented
c---------------------------------------------------------------------
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

            do k=1,dim
               rinch(k)=rfac*(-dg(k)*rnormal(k,i)+rinc)
     1             + rincg(k)*dg(k)
            enddo

            m=dim
            do k=1,dim-1
            do j=k+1,dim
               m=m+1
               rinch(m)=-rfac*dg(k)*rnormal(j,i) + rincg(k)*dg(j)
            enddo
            enddo
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
