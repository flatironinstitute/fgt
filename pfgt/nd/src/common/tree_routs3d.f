c
c
c    common routines for generating and processing
c     a level restricted quad tree in 2D
c   
c
c
c

c---------------------------------------------------------------      
      subroutine compute_modified_list1(nlevels,npwlevel,
     1                   nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,nlist1,
     3                   mnlist1,list1)
c     Compute max nuber of boxes in list1
      implicit none
      integer nlevels,npwlevel,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,isep
      integer nlist1(nboxes)
      integer list1(mnlist1,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO      
      if(nchild(1).eq.0) then
         nlist1(1) = 1
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif

      do ilev = 1,nlevels
        firstbox = laddr(1,ilev)
        lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,xdis,ydis,zdis,distest)
        do ibox = firstbox,lastbox
          dad = iparent(ibox)
c         Compute list1 of ibox if it is childless
          if(nchild(ibox).eq.0) then
            do i=1,nnbors(ibox)
              jbox = nbors(i,ibox)
c
cc            boxes in list 1 at the same level
c
cccc              if (ilev .gt. npwlevel) then
                 
cccc              elseif (ilev .eq. npwlevel) then
              if (ilev .eq. npwlevel) then
c     at the cutoff level, colleagues are directly added to list1
cccc                 if (nchild(jbox).eq. 0) then
cccc                    nlist1(ibox) = nlist1(ibox)+1
cccc                    list1(nlist1(ibox),ibox) = jbox
cccc                 endif

                nlist1(ibox) = nlist1(ibox)+1
                list1(nlist1(ibox),ibox) = jbox
              else
                if(nchild(jbox).eq.0) then
                   nlist1(ibox) = nlist1(ibox) + 1
                   list1(nlist1(ibox),ibox) = jbox
                else
c
c                  boxes in list1 at level ilev+1
c
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                 2.0d0*isep
                   do j=1,8
                      kbox = ichild(j,jbox)
                      if(kbox.gt.0) then
                         xdis = dabs(centers(1,kbox)-centers(1,ibox))
                         ydis = dabs(centers(2,kbox)-centers(2,ibox))
                         zdis = dabs(centers(3,kbox)-centers(3,ibox))

                         if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                       zdis.lt.distest) then
                            nlist1(ibox) = nlist1(ibox)+1
                            list1(nlist1(ibox),ibox) = kbox
                         endif
                      endif
                   enddo

                endif
              endif
            enddo
c     
cc          compute list1 at level ilev-1 
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                  distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                2.0d0*isep
                  xdis = dabs(centers(1,jbox)-centers(1,ibox))
                  ydis = dabs(centers(2,jbox)-centers(2,ibox))
                  zdis = dabs(centers(3,jbox)-centers(3,ibox))
                  if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                zdis.lt.distest) then
                     nlist1(ibox) = nlist1(ibox)+1
                     list1(nlist1(ibox),ibox) = jbox
                  endif
               endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c----------------------------------------------------------------
c      
c      
c      
c      
      subroutine gt3d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1   iptr,centers,
     1   boxsize,iper,mnlistpw)
c
c     determine maximum number of elements in listpw
c
c     NOTE: we use max values
c
      implicit real *8 (a-h,o-z)
      integer ltree, iptr(8)
      integer nlevels,nboxes,itree(ltree)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer mnlistpw

      mnlistpw = 26

      return
      end
c
c
c
c
c
      subroutine gt3d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,nchild,ichild,
     2    boxsize,laddr,mnlistpw,nlistpw,listpw)
c     this subroutine returns lists of pw interaction boxes
c
c     input parameters:
c
c     nlevels = number of levels, level 0 contains boxes of size 6\sqrt{delta}
c     npwlevel   = recursive SOE picture stops at the level npwlevel
c     nboxes = total number of boxes in the tree
c     itree - laddr: tree stuff
c     nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c     ichild         in: integer(8,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c
c     mnlistpw = maximum number of a particular type of PW expansions
c      
c     output parameters:
c
c     nlistpw (nboxes) = contains the number of boxes for the PW interaction
c                            for each box
c     listpw (mnlistpw,nboxes) = contains the box ID for the PW interaction
c                            for each box
c      
      implicit real *8 (a-h,o-z)
      integer iptr(8), ltree
      integer nlevels,nboxes
      integer itree(ltree),laddr(2,0:nlevels)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer nchild(nboxes),ichild(8,nboxes)
      integer mnlistpw
      integer nlistpw(nboxes), listpw(mnlistpw,nboxes)

      integer ibox

      do ibox=1,nboxes
         nlistpw(ibox)=0
      enddo

      nlevstart=max(npwlevel,0)
      do ilev=npwlevel,npwlevel
         do ibox = laddr(1,ilev),laddr(2,ilev)
            ncoll = itree(iptr(6)+ibox-1)
            do i=1,ncoll
               jbox = itree(iptr(7) + (ibox-1)*27+i-1)
cccc              jlev = itree(iptr(2)+jbox-1)
cccc              if (jbox.ne.ibox .and. ilev.eq.jlev .and. 
cccc     1            (nchild(jbox).gt.0 .or. nchild(ibox).gt.0)) then
               if (jbox.ne.ibox .and. nchild(jbox).gt.0) then                 
                 nlistpw(ibox)=nlistpw(ibox)+1
                 listpw(nlistpw(ibox),ibox) = jbox
               endif
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
      
