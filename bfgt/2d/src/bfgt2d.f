cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

      subroutine bfgt2d(nd,delta,eps,nboxes,nlevels,ltree,
     1   itree,iptr,norder,ncbox,ttype,fvals,centers,boxsize,npbox,
     2   pot,timeinfo,tprecomp)

c
c       This code applies the heat volume potential
c       to a collection of right hand sides
c 
c       input
c         nd - integer
c            number of right hand sides
c         delta - double precision
c            Gaussian variance 
c         eps - double precision
c            tolerance requested
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         norder - integer
c           order of expansions for input coefficients array
c         ncbox - integer
c           number of coefficients of expansions of functions
c           in each of the boxes
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fvals - double precision (npbox,nboxes)
c           function tabulated on a grid
c         centers - double precision (3,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c
c     output:
c         pot - double precision (npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),norder,ncbox,npbox
      character *1 ttype
      real *8 fvals(nd,npbox,nboxes)
      real *8 pot(nd,npbox,nboxes)
      real *8 centers(2,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(6),tprecomp(3)
c
cc     additional fgt variables
c
      integer *8 lmptot
      integer, allocatable :: iaddr(:,:)
      real *8, allocatable :: rmlexp(:)

c
cc      temporary variables
c
      integer npwlevel
      integer i,ilev,lmptmp,idim,ndim
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      real *8 dcutoff
      real *8 omp_get_wtime
      real *8 time1,time2,pi,done,pmax,bs0,cen0(2),bsize,pweps

c     cutoff length      
      dcutoff = sqrt(delta*log(1.0d0/eps))
cccc      call prin2(' dcutoff=*',dcutoff,1)
      
c     find the cutoff level
      npwlevel = nlevels+1
      do i=nlevels,0,-1
         if (boxsize(i).ge. dcutoff) then
            npwlevel=i
            exit
         endif
      enddo
c
      if (boxsize(0) .le. dcutoff) npwlevel=0
      call prinf(' npwlevel =*',npwlevel,1)
cccc      call prin2(' boxsize(npwlevel)=*',boxsize(npwlevel),1)
c
      
      done = 1
      pi = atan(done)*4.0d0
      
      ifprint = 1
c
c     compute the length of plane wave expansion
      npw=0
      if (npwlevel.ge.0.and.npwlevel.le.nlevels) then
         bsize=2*boxsize(npwlevel)
         call fgtpwterms(bsize,delta,eps,pmax,npw)
      endif
      call prinf(' npw =*',npw,1)
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c
c       ... allocate iaddr and temporary arrays
c
      allocate(iaddr(2,nboxes))
c
c
c     allocate memory need by multipole, local expansions at all levels
c     
c     irmlexp is pointer for workspace need by various expansions.
c
      call g2dmpalloc(nd,itree,iaddr,
     1    nlevels,npwlevel,lmptot,npw)
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      allocate(rmlexp(lmptot),stat=ier)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in allocating rmlexp=*',
     1   time2-time1,1)
      
cccc      do i=1,lmptot
cccc         rmlexp(i)=0
cccc      enddo
cccc      call cpu_time(time2)
ccccC$        time2=omp_get_wtime()
cccc      if( ifprint .eq. 1 ) call prin2('time in allocating rmlexp=*',
cccc     1   time2-time1,1)
      
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call bfgt2dmain(nd,delta,eps,nboxes,nlevels,ltree,
     1    itree,iptr,norder,ncbox,ttype,fvals,centers,boxsize,
     2    iaddr,rmlexp,npwlevel,pmax,npw,
     3    npbox,pot,timeinfo)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in box fgt main=*',
     1   time2-time1,1)

      return
      end
c
c
c
c
      subroutine bfgt2dmain(nd,delta,eps,nboxes,nlevels,ltree,
     1    itree,iptr,norder,ncbox,ttype,fvals,centers,boxsize,
     2    iaddr,rmlexp,npwlevel,pmax,npw,
     3    npbox,pot,timeinfo)

c
c       This code applies the heat volume potential
c       to a collection of right hand sides
c 
c       input
c
c         nd - integer,   number of box FGTs with the same tree
c         delta - double precision, Gaussian variance 
c         eps - double precision
c            tolerance requested
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         norder - integer
c           order of expansions for input coefficients array
c         ncbox - integer
c           number of coefficients of expansions of functions
c           in each of the boxes
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fvals - double complex (npbox,nboxes)
c            function values tabulated on the tree
c         centers - double precision (2,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         iaddr - (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c         rmlexp - double precision, stores multipole and local PW expansions
c                  for each box below the cutoff level npwlevel
c
c
c         npwlevel - integer
c             cutoff level at which the PW expansion is valid
c         pmax - double precision,  cutoff limit in the planewave expansion
c         npw  - integer,  length of planewave expansions
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c
c     output:
c         pot - double precision (nd,npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      integer nd
      real *8 delta,eps
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),ncbox,npbox
      real *8 fvals(nd,npbox,nboxes)
      real *8 pot(nd,npbox,nboxes)
      real *8 boxsize(0:nlevels),centers(2,nboxes)
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(*)

      character *1 ttype
c
cc        temporary list info
c
      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlistpw(:), listpw(:,:)
      integer, allocatable :: ifhung(:),iflocal(:)

      integer ndirect
      
      real *8, allocatable :: fcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)

      real *8 timelev(0:200)
      real *8 ws(100),ts(100)
      
      complex *16, allocatable :: wpwshift(:,:,:)
      complex *16, allocatable :: wpwmsshift(:,:,:)
      
      complex *16, allocatable :: tab_leg2pw(:,:,:),tab_pw2pot(:,:,:)
      complex *16, allocatable :: tab_leg2pw2(:,:,:)
      complex *16, allocatable :: ff(:,:)
      complex *16, allocatable :: gg(:,:)
      real *8, allocatable :: hh(:,:)

      real *8, allocatable :: tab_coll(:,:,:,:)
      real *8, allocatable :: tab_stob(:,:,:,:)
      real *8, allocatable :: tab_btos(:,:,:,:)

      ifprint = 1

      done = 1
      pi = atan(done)*4
c
c       compute coefs
c
      allocate(fcoefs(ncbox,nd,nboxes),xq(norder),umat(norder,norder),
     1   vmat(norder,norder),wts(norder))
     
      itype = 2
      call legeexps(itype,norder,xq,umat,vmat,wts)

      call cpu_time(ttt1)
      do ilev = 0,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0)
     1        call legtrans2d(nd,norder,fvals(1,1,ibox),
     2        fcoefs(1,1,ibox),umat)
        enddo
      enddo
      call cpu_time(ttt2)
      call prin2('val to coeffs time=*', ttt2-ttt1,1)
c
c       initialize potential
c 
      do i=1,nboxes
         do j=1,npbox
            do ind=1,nd
               pot(ind,j,i) = 0
            enddo
        enddo
      enddo
c      
      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      ndirect=0
      do ilev = 0,nlevend
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
c           Check if the current box is a nonempty leaf box            
            if(nchild.eq.0) then
               ifhung(ibox) = 1
               ndirect = ndirect+1
            endif
         enddo
      enddo
      call prinf('number of direct evaluation source boxes=*',
     1    ndirect,1)

      do i=1,6
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo

c
c
c       compute list info
c
      isep = 1
      iper = 0
      mnlist1 = 13
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))

c     modified list1 for direct evaluation
      call compute_modified_list1(nlevels,npwlevel,
     1  nboxes,itree,ltree,iptr,
     2  centers,boxsize,
     3  iper,mnlist1,nlist1,list1)

c     compute the tables converting Legendre polynomial expansion to potential
c     values, used in direct evaluation
      allocate(tab_coll(norder,norder,-1:1,0:nlevels))
      allocate(tab_stob(norder,norder,4,0:nlevels))
      allocate(tab_btos(norder,norder,4,0:nlevels))
      allocate(hh(norder,norder))
      call cpu_time(t1)

      nnodes=50
      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_coll(norder,nnodes,delta,boxsize(ilev),
     1       tab_coll(1,1,-1,ilev))
      enddo

      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_stob(norder,nnodes,delta,boxsize(ilev),
     1       tab_stob(1,1,1,ilev))
      enddo

      do ilev = 0,min(npwlevel,nlevels)
         call mk_loctab_btos(norder,nnodes,delta,boxsize(ilev),
     1       tab_btos(1,1,1,ilev))
      enddo

c
c     check whether we need to create and evaluate planewave expansions 
c     for boxes at the cutoff level
      allocate(iflocal(nboxes))
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
         do i=1,nboxes
            iflocal(i)=0
         enddo

         do ilev=npwlevel,npwlevel
            do ibox=itree(2*ilev+1),itree(2*ilev+2)
               nchild = itree(iptr(4)+ibox-1)
               if (nchild .eq. 0) then
                  ncoll = itree(iptr(6)+ibox-1)
                  do j=1,ncoll
                     jbox = itree(iptr(7) + (ibox-1)*9+j-1)
                     nchild = itree(iptr(4)+jbox-1)
                     jlev = itree(iptr(2)+jbox-1)
                     if (nchild .gt. 0 .and. jlev.eq.ilev) then
                        iflocal(ibox)=1
                        exit
                     endif
                  enddo
               endif
            enddo
         enddo
      endif
      
c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800

c      
c     get planewave nodes and weights
      call get_pwnodes(pmax,npw,ws,ts)

c     compute translation matrices for PW expansions
      nlevstart = max(npwlevel,0)

c     compute the tables converting Legendre polynomial expansion
c     to planewave expansion
      nnodes = 100
      allocate(tab_leg2pw(norder,npw,0:nlevels))
      allocate(tab_leg2pw2(norder,npw,0:nlevels))
      allocate(ff(norder,npw/2))

      do ilev=nlevstart,nlevels
cccc         call mk_leg2pw_old(norder,npw,nnodes,ws,ts,delta,boxsize(ilev),
cccc     1       tab_leg2pw2(1,1,ilev))
         call mk_leg2pw(norder,npw,nnodes,ws,ts,delta,boxsize(ilev),
     1       tab_leg2pw(1,1,ilev))
cccc         call prin2('leg2pw2=*',tab_leg2pw2(1,1,ilev),norder*npw*2)
cccc         call derr(tab_leg2pw(1,1,ilev),tab_leg2pw2(1,1,ilev),
cccc     1       norder*npw*2,rerr1)
cccc         print *, rerr1
      enddo
cccc      pause
c     compute the tables converting planewave expansions to potential values
      allocate(tab_pw2pot(npw,norder,0:nlevels))
      allocate(gg(norder,npw/2))
      
      do ilev=nlevstart,nlevels
         call mk_pw2pot(norder,npw,ts,xq,delta,boxsize(ilev),
     1       tab_pw2pot(1,1,ilev))      
      enddo
      
      nexp = npw*npw/2
      
      nmax = 1
      allocate(wpwshift(nexp,-nmax:nmax,-nmax:nmax))
      xmin  = boxsize(nlevstart)/sqrt(delta)
      call pw_translation_matrices(xmin,npw,ts,nmax,
     1    wpwshift)

      nmax = nlevels-max(npwlevel,0)      
      allocate(wpwmsshift(nexp,4,nmax))
      xmin2 = boxsize(nlevels)/sqrt(delta)/2
      call merge_split_pw_matrices(xmin2,npw,ts,nmax,wpwmsshift)
c     xmin is used in shiftpw subroutines to
c     determine the right translation matrices
c
      xmin  = boxsize(nlevstart)
      xmin2 = boxsize(nlevels)/2
      call cpu_time(t2)
      call prin2('precomputation in fgt =*', t2-t1,1)
      
c
c     compute list info
c
      call gt2d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1    iptr,centers,
     2    boxsize,iper,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call gt2d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,
     2    boxsize,itree(iptr(1)),
     3    mnlistpw,nlistpw,listpw)
      
c
c     ... set all multipole and local expansions to zero
c
cccc      do ilev = nlevstart,nlevels
ccccC$OMP PARALLEL DO DEFAULT (SHARED)
ccccC$OMP$PRIVATE(ibox)
cccc         do ibox = itree(2*ilev+1),itree(2*ilev+2)
cccc           call g3dpwzero_vec(nd,rmlexp(iaddr(1,ibox)),npw)
cccc           call g3dpwzero_vec(nd,rmlexp(iaddr(2,ibox)),npw)
cccc         enddo
ccccC$OMP END PARALLEL DO         
cccc       enddo
c
cccc      call prinf('finished initializing rmlexp*',i,0)
c

      call prinf('laddr=*',itree(iptr(1)),2*(nlevels+1))

      

c
c
c        step 1: convert coeffs to planewave expansions
c
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0)
      
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do 1100 ilev = nlevels,nlevstart,-1
cccc         nb=0
cccc         dt=0
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
cccc            if (ilev .eq. npwlevel .and. iflocal(ibox).eq.0) then
c              do nothing here
cccc            else    
            if (ilev .ne. npwlevel .or. iflocal(ibox).eq.1) then
               nchild = itree(iptr(4)+ibox-1)
c              Check if current box is a leaf box            
               if(nchild.eq.0) then
cccc                  nb=nb+1
cccc                  call cpu_time(t1)
c                 form PW expansion directly
                  call leg2d_to_pw(nd,norder,fcoefs(1,1,ibox),npw,
     1                ff,tab_leg2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
cccc     1                tab_leg2pw(1,1,ilev),rmlexp(iaddr(1,ibox)))
cccc                  call cpu_time(t2)
cccc                  dt=dt+t2-t1
               endif
            endif
         enddo
C$OMP END PARALLEL DO
ccc 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f5.2)         
cccc         write(6,111) ilev,nb,dt
c     end of ilev do loop
 1100 continue
      
      call cpu_time(time2)
C$       time2 = omp_get_wtime()
      timeinfo(1) = time2-time1

      
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do 1200 ilev=nlevels-1,nlevstart,-1
         klev = nlevels - ilev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,dx,dy,dz,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild .gt. 0) then
              call g2dpwzero_vec(nd,rmlexp(iaddr(1,ibox)),npw)
           endif
           do i=1,nchild
              jbox = itree(iptr(5)+4*(ibox-1)+i-1)

              dx= (centers(1,ibox) - centers(1,jbox))
              dy= (centers(2,ibox) - centers(2,jbox))
              if (dx.gt.0 .and. dy.gt.0) then
                 k=1
              elseif (dx.gt.0 .and. dy.lt.0) then
                 k=2
              elseif (dx.lt.0 .and. dy.gt.0) then
                 k=3
              elseif (dx.lt.0 .and. dy.lt.0) then
                 k=4
              endif

              call g2dshiftpw_vec(nd,nexp,rmlexp(iaddr(1,jbox)),
     1            rmlexp(iaddr(1,ibox)),wpwmsshift(1,k,klev))
           enddo
        enddo
C$OMP END PARALLEL DO    
 1200 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1
c
      
      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = nlevstart,nlevstart
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,j)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           call g2dcopypwexp_vec(nd,nexp,rmlexp(iaddr(1,ibox)),
     1         rmlexp(iaddr(2,ibox)))
c          shift PW expansions
           do j=1,nlistpw(ibox)
              jbox=listpw(j,ibox)
              jx= nint((centers(1,ibox) - centers(1,jbox))/xmin)
              jy= nint((centers(2,ibox) - centers(2,jbox))/xmin)

              call g2dshiftpw_vec(nd,nexp,rmlexp(iaddr(1,jbox)),
     1            rmlexp(iaddr(2,ibox)),wpwshift(1,jx,jy))
           enddo
        enddo
C$OMP END PARALLEL DO        
 1300 continue
c

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3) = time2-time1

cccc      call prin2('timeinfo4=*',time2-time1,1)

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      do 1400 ilev = nlevstart,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,dx,dy,dz,k)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if (nchild.gt.0) then
          do i=1,nchild
             jbox = itree(iptr(5)+4*(ibox-1)+i-1)
             dx= centers(1,jbox) - centers(1,ibox)
             dy= centers(2,jbox) - centers(2,ibox)
             if (dx.gt.0 .and. dy.gt.0) then
                k=1
             elseif (dx.gt.0 .and. dy.lt.0) then
                k=2
             elseif (dx.lt.0 .and. dy.gt.0) then
                k=3
             elseif (dx.lt.0 .and. dy.lt.0) then
                k=4
             endif
             call g2dshiftpw_loc_vec(nd,nexp,rmlexp(iaddr(2,ibox)),
     1           rmlexp(iaddr(2,jbox)),wpwmsshift(1,k,nlevels-ilev))
          enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1400 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1
      
      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()

      do 1500 ilev = nlevstart,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild)
C$OMP$SCHEDULE(DYNAMIC)
cccc         call cpu_time(t1)
cccc         nb=0
cccc         dt=0
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           if (ilev .eq. npwlevel .and. iflocal(ibox).eq.0) then
c            do nothing here
           else    
             nchild = itree(iptr(4)+ibox-1)
             if(nchild.eq.0) then
cccc                nb=nb+1
                call cpu_time(t1)
                call g2d_pw2pot(nd,norder,npw,rmlexp(iaddr(2,ibox)),
     1              gg,tab_pw2pot(1,1,ilev),pot(1,1,ibox))
cccc                call cpu_time(t2)
cccc                dt=dt+t2-t1
             endif
           endif
ccc    end of ibox loop        
        enddo
 222    format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f5.2)         
cccc        write(6,222) ilev,nb,dt
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(5) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nl1,bs,ix,iy,iz,jlev)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
c        ibox is the source box            
            if (ifhung(ibox) .eq. 1) then
               
               nl1 = nlist1(ibox)
               do j=1,nl1
cccc              jbox is the target box
                  jbox = list1(j,ibox)
                  jlev = itree(iptr(2)+jbox-1)
                  bs = boxsize(jlev)
c                 colleague                  
                  if (ilev .eq. jlev) then
                     if (iflocal(ibox).eq. 1 .and. jbox.eq.ibox) then
                     else
                     ix = (centers(1,jbox)-centers(1,ibox))/bs
                     iy = (centers(2,jbox)-centers(2,ibox))/bs

                     call leg2d_to_potloc(nd,norder,fcoefs(1,1,ibox),
     1                   hh,pot(1,1,jbox),
     2                   tab_coll(1,1,ix,jlev),
     3                   tab_coll(1,1,iy,jlev))
                     endif
c                 big source box to small target box                     
                  elseif (ilev .eq. jlev-1) then
                     ix = (centers(1,jbox)-centers(1,ibox))/bs+2.55d0
                     iy = (centers(2,jbox)-centers(2,ibox))/bs+2.55d0

                     call leg2d_to_potloc(nd,norder,fcoefs(1,1,ibox),
     1                   hh,pot(1,1,jbox),
     2                   tab_btos(1,1,ix,jlev),
     3                   tab_btos(1,1,iy,jlev))
c                 small source box to big target box 
                  elseif (ilev .eq. jlev+1) then
                     ix = (centers(1,jbox)-centers(1,ibox))/bs*2+2.55d0
                     iy = (centers(2,jbox)-centers(2,ibox))/bs*2+2.55d0
cccc                     print *, ilev, jlev, ix,iy,iz                     
                     call leg2d_to_potloc(nd,norder,fcoefs(1,1,ibox),
     1                   hh,pot(1,1,jbox),
     2                   tab_stob(1,1,ix,jlev),
     3                   tab_stob(1,1,iy,jlev))
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
cccc      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
c
c
c
c------------------------------------------------------------------    
      subroutine g2dmpalloc(nd,laddr,iaddr,
     1    nlevels,npwlevel,lmptot,npw)
c     This subroutine determines the size of the array
c     to be allocated for multipole/local expansions
c
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     npwlevel    in: Integer
c                 cutoff level where the plane wave expansion is
c                 valid at or below ilev = npwlevel
c
c     npw         in: Integer
c                 Number of terms in the plane wave expansion
c
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,npwlevel,npw,nd
      integer iaddr(2,*), laddr(2,0:nlevels)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart
c
      istart = 1
      if (npwlevel .eq. nlevels) then
         lmptot=0
         return
      endif
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      do i = nlevstart,nlevels

         nn = (npw*npw/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole PW expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i = nlevstart,nlevels

         nn = (npw*npw/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local PW expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(2,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
            
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
c
c
c
