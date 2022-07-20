c
c
c     common routines for generating and processing
c     a level restricted quad tree in ndim=1,2,3 dimensions
c   
c     created by Shidong Jiang on 03/31/2022
c
c
      subroutine tree_refine_boxes(ndim,irefinebox,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer ndim,nboxes,nbloc,nbctr,nlctr
      real *8 centers(ndim,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(2**ndim,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)
      integer ii

      integer i,ibox,nel0,j,l,jbox,nel1,nbl,k,mc
      real *8 bsh
      integer isgn(ndim,2**ndim)

      call get_child_box_sign(ndim,isgn)

      mc=2**ndim
      bsh = bs/2

      allocate(isum(nbloc))
      if(nbloc.gt.0) call cumsum(nbloc,irefinebox,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*mc
          
          nchild(ibox) = mc
          do j=1,mc
            jbox = nbl+j
            do k=1,ndim
               centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,mc
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*mc


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(ndim,nb,centers,ilevel,iparent,nchild,
     1              ichild,centers2,ilevel2,iparent2,nchild2,ichild2)

       implicit none
       integer ndim,nd,nb,npb
       real *8 centers(ndim,nb),centers2(ndim,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(2**ndim,nb),ichild2(2**ndim,nb)

       integer i,j,nel,mc

       mc=2**ndim

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         do j=1,ndim
            centers2(j,i) = centers(j,i)
         enddo
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         do j=1,mc
            ichild2(j,i) = ichild(j,i)
         enddo
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c
      subroutine computecoll(ndim,nlevels,nboxes,laddr,boxsize,
     1    centers,iparent,nchild,ichild,iperiod,
     2    nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     ndim        in: integer
c                 dimension of the space
c
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iperiod     in: integer
c                 0: free space; 1: doubly periodic
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(3**ndim,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer ndim,nlevels,nboxes
      integer iperiod
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(ndim,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(2**ndim,nboxes)
      integer nnbors(nboxes)
      integer nbors(3**ndim,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox,mc,mnbors,k,ifnbor
      real *8 bs0,dis,dp1

      bs0=boxsize(0)
      
      mc= 2**ndim
      mnbors=3**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO     
      
c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,ifnbor,k,dp1,dis)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,mc
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c     Check if kbox is a nearest neighbor or in list 2
                      ifnbor=1
                      do k=1,ndim
                         dis=abs(centers(k,kbox)-centers(k,ibox))
                         if (iperiod.eq.1) then
                            dp1=bs0-dis
                            if (dp1.lt.dis) dis=dp1
                         endif
                         if (dis.gt.1.05*boxsize(ilev)) then
                            ifnbor=0
                            exit
                         endif
                      enddo
                         
                      if(ifnbor.eq.1) then
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
c--------------------------------------------------------------------      
      subroutine updateflags(ndim,curlev,nboxes,nlevels,laddr,
     1    nchild,ichild,nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer ndim,curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(2**ndim,nboxes)
      integer nnbors(nboxes), nbors(3**ndim,nboxes)
      integer iflag(nboxes)
      double precision centers(ndim,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict, mc
      double precision distest,dis

      mc=2**ndim
      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,dis,ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,mc
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        ict = 0
                        do k=1,ndim
                           dis = centers(k,kbox) - centers(k,ibox) 
                           if(abs(dis).le.distest) ict = ict + 1
                        enddo
                        if(ict.eq.ndim) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c

      subroutine tree_refine_boxes_flag(ndim,iflag,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer ndim,nboxes,nbloc,nbctr,nlctr
      real *8 centers(ndim,nboxes),bs,bsh
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(2**ndim,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox
      integer, allocatable :: isum(:),itmp(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl,k
      integer ii,mc
      integer isgn(ndim,2**ndim)

      call get_child_box_sign(ndim,isgn)
      
      bsh = bs/2
      mc = 2**ndim

      allocate(isum(nbloc),itmp(nbloc))
      do i=1,nbloc
        ibox = ifirstbox+i-1
        itmp(i) = 0
        if(iflag(ibox).gt.0) itmp(i) = 1
      enddo
      if(nbloc.gt.0) call cumsum(nbloc,itmp,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(iflag(ibox).gt.0) then
          nbl = nbctr + (isum(i)-1)*mc
          
          nchild(ibox) = mc
          do j=1,mc
            jbox = nbl+j
            do k=1,ndim
               centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,mc
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*mc


      return
      end
c
c
c
c
c
c
      subroutine print_tree2d(itree,ltree,nboxes,centers,boxsize,
     1    nlevels,iptr,fname)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(8)
c            pointer to various arrays inside itree
c          fname - character *
c            file name to which tree info is to be written
c 
c          output
c            file with name fname, which contains the tree info
c            file can be plotted using the python script
c              tree_plot.py containted in src/common
c
      implicit real *8 (a-h,o-z)
      integer ipts(8),ltree
      integer itree(ltree),nboxes,nlevels
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      character (len=*) fname

      open(unit=33,file=trim(fname))
      

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(4)+ibox-1).eq.0) then
           ilev = itree(iptr(2)+ibox-1)
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
         endif
      enddo

      close(33)

      return
      end
c
c
c
c
c
      subroutine compute_mnlist1(ndim,nboxes,nlevels,laddr,
     1    centers,boxsize,iparent,nchild,
     2    ichild,isep,nnbors,nbors,iper,mnlist1)
c     Compute max number of boxes in list1
      implicit none
      integer ndim,nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(ndim,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(2**ndim,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(3**ndim,nboxes)
      integer mnlist1
      integer nlist1(nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad,mc,ifnbor,iflist1
      double precision dis,distest,bs0,dp1

      bs0=boxsize(0)

      mnbors=3**ndim
      mc=2**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO
      
      nlist1(1) = 1

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,distest,dis,iflist1,ifnbor)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                     check for list1 at the same level
c
                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                  endif
c
cc                     check for list1 and list3 at one ilev+1
                  if(nchild(jbox).gt.0) then
                    distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                   2.0d0*isep
                    do j=1,mc
                      kbox = ichild(j,jbox)
                      if(kbox.gt.0) then
                        ifnbor=1
                        do k=1,ndim
                          dis = dabs(centers(k,kbox)-centers(k,ibox))
                          if (iper .eq. 1) then
                             dp1 = bs0-dis
                             if (dp1.lt.dis) dis=dp1
                          endif
                           
                          if (dis.gt.distest) then
                            ifnbor=0
                            exit
                          endif
                        enddo
                        if(ifnbor.eq.1) then
                           nlist1(ibox) = nlist1(ibox)+1
                        endif
                      endif
                    enddo
                  endif
               enddo
c
cc               compute list1 and list4 for boxes at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                    2.0d0*isep
                      iflist1=1
                      do k=1,ndim
                         dis = dabs(centers(k,jbox)-centers(k,ibox))
                         if (iper .eq. 1) then
                            dp1 = bs0-dis
                            if (dp1.lt.dis) dis=dp1
                         endif
                         
                         if (dis.ge.distest) then
                            iflist1=0
                            exit
                         endif
                      enddo
                      if(iflist1.eq.1) then
                         nlist1(ibox) = nlist1(ibox)+1
                      endif
                   endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo

      mnlist1 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
C$OMP$REDUCTION(max:mnlist1)      
      do i=1,nboxes
         if(nlist1(i).gt.mnlist1) mnlist1 = nlist1(i)
      enddo
C$OMP END PARALLEL DO      

      return
      end
c      
c      
c      
c---------------------------------------------------------------      
      subroutine computelists(ndim,nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)
c     Compute max number of boxes in list1,list2,list3,list4
      implicit none
      integer ndim,nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(ndim,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(2**ndim,nboxes)
      integer mnbors
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4,isep
      integer nlist1(nboxes),nlist2(nboxes),nlist3(nboxes)
      integer nlist4(nboxes)
      integer list1(mnlist1,nboxes),list2(mnlist2,nboxes)
      integer list3(mnlist3,nboxes),list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l,mc
      integer firstbox,lastbox,dad,iflist1,iflist2,iflist4
      double precision dis,distest

      mc=2**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
      enddo
C$OMP END PARALLEL DO      
      if(nchild(1).eq.0) then
         nlist1(1) = 1
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,dis,distest,k)
C$OMP$PRIVATE(iflist1,iflist2,iflist4)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               distest = 1.05d0*isep*boxsize(ilev)
               do j=1,mc
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     iflist2=0
                     do k=1,ndim
                        dis=centers(k,kbox)-centers(k,ibox)
                        if(abs(dis).ge.distest) then
                           iflist2=1
                           exit
                        endif
                     enddo
                     if (iflist2.eq.1) then
                        nlist2(ibox) = nlist2(ibox) + 1
                        list2(nlist2(ibox),ibox) = kbox
                     endif
                  endif
               enddo
            enddo
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)
c
cc                boxes in list 1 at the same level
c

                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                     list1(nlist1(ibox),ibox) = jbox
                  else
c
cc                     boxes in list1 and list3 at level ilev+1
c
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,mc
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                          iflist1=1
                          do k=1,ndim
                            dis = dabs(centers(k,kbox)-centers(k,ibox))
                            if (dis.ge.distest) then
                               iflist1=0
                               exit
                            endif
                          enddo
                          if(iflist1.eq.1) then
                              nlist1(ibox) = nlist1(ibox)+1
                              list1(nlist1(ibox),ibox) = kbox
                           else
                              nlist3(ibox) = nlist3(ibox)+1
                              list3(nlist3(ibox),ibox) = kbox
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                    2.0d0*isep
                      iflist1=1
                      do k=1,ndim
                         dis = dabs(centers(k,jbox)-centers(k,ibox))
                         if (dis.ge.distest) then
                            iflist1=0
                            exit
                         endif
                      enddo
                      if(iflist1.eq.1) then
                         nlist1(ibox) = nlist1(ibox)+1
                         list1(nlist1(ibox),ibox) = jbox
                      endif
                   endif
                enddo
            endif
c
cc           compute list 4 at level ilev-1
c
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                2.0d0*isep
                   iflist4=0
                   do k=1,ndim
                      dis = dabs(centers(k,jbox)-centers(k,ibox))
                      if (dis.gt.distest) then
                         iflist4=1
                         exit
                      endif
                   enddo
                   if (iflist4.eq.1)then
                      nlist4(ibox) = nlist4(ibox)+1
                      list4(nlist4(ibox),ibox)=jbox
                   endif
                endif
             enddo
          enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c      
c      
c      
      subroutine flag_box_refine(ndim,nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,
     3                   ifrefine,nrefine)
c     flag boxes that may require refinement for the box FGT
      implicit none
      integer ndim,nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(ndim,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(2**ndim,nboxes)
      integer mnbors
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer isep
      integer ifrefine(nboxes),nrefine

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l,mc
      integer firstbox,lastbox,iflist1
      double precision dis,distest

      mc=2**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         ifrefine(i)=0
      enddo
C$OMP END PARALLEL DO      

      nrefine=0
      do ilev = 1,nlevels-1
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,j,kbox,dis,distest,k,iflist1)
         do 1100 ibox = firstbox,lastbox
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)
                  if(nchild(jbox).ne.0) then
c     
cc                     boxes in list1 and list3 at level ilev+1
c
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                   2.0d0*isep
                     do j=1,mc
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                          iflist1=1
                          do k=1,ndim
                            dis = dabs(centers(k,kbox)-centers(k,ibox))
                            if (dis.ge.distest) then
                               iflist1=0
                               exit
                            endif
                          enddo
                          if(iflist1.eq.1) then
                             ifrefine(ibox)=1
                             nrefine=nrefine+1
                             goto 1100
                          endif
                       endif
                     enddo
                  endif
               enddo
            endif
 1100    continue
C$OMP END PARALLEL DO         
      enddo

      return
      end
c      
c      
c      
c      
c----------------------------------------------------------------
c      
      subroutine gnd_compute_modified_list1(ndim,npwlevel,ifpwexp,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,iperiod,
     2    mnlist1,nlist1,list1)
c
c
c     This subroutine computes the modified list1 of a given tree
c     structure
c
c     Note: new definition of list 1 - for ilev <= npwlevel,             
c     list1 of a leaf box contains all childless neighbors at or
c     above npwlevel. ifpwexp(ibox)=1, then ibox is excluded from 
c     list1 since then the self interaction is handled via
c     plane wave expansions.
c                  
c     Assume that the tree is level-restricted. Then the
c     new list1 of source ibox contains all target boxes 
c     that require the evaluation of direct interactions with ibox 
c     in the box FGT.  
c      
c      
c     haven't tested iper=1 case.
c      
c      
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     npwlevel    in: integer
c                 Cutoff level
c     ifpwexp     in: integer(nboxes)
c                 1: requires pwexp; 0: does not require pwexp
c      
c     nboxes      in: integer
c                 Total number of boxes
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     centers     in: real *8(2,nboxes)
c                 xy coordinates of centers of boxes
c   
c     boxsize     in: real *8(0:nlevels)
c                 Array of boxsizes
c   
c     iperiod     in: integer
c                 flag for periodic implementations. Currently not used.
c                 Feature under construction
c 
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i.
c                  
c      
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nlevels,npwlevel,nboxes,ndim
      integer iptr(8),ltree
      integer iperiod
      integer itree(ltree),ifpwexp(nboxes)
      real *8 boxsize(0:nlevels)
      real *8 centers(ndim,nboxes)
      integer mnlist1
      integer nlist1(nboxes), list1(mnlist1,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad,mnbors,mc
      integer i,j,ifirstbox,ilastbox,ii,iflist1,k
      real *8 distest,dis,bs0,dp1

      bs0=boxsize(0)

      mnbors=3**ndim
      mc=2**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
        nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO      


c     Setting parameters for level = 0
      if(itree(iptr(4)).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif

      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel
      do ilev = 1,nlevend
         ifirstbox = itree(2*ilev+1)
         ilastbox = itree(2*ilev+2)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,dis,distest,k)
C$OMP$PRIVATE(iflist1,dp1)
         do ibox = ifirstbox,ilastbox
            dad = itree(iptr(3)+ibox-1)
c           Compute list1 of ibox if it is childless
            if(itree(iptr(4)+ibox-1).eq.0) then
               do i=1,itree(iptr(6)+ibox-1)
                  jbox = itree(iptr(7)+mnbors*(ibox-1)+i-1)
c
cc                boxes in list 1 at the same level
c

                  if (itree(iptr(4)+jbox-1).eq.0) then
                     if ((ifpwexp(ibox).ne.1) .or. (jbox.ne.ibox)) then
                        nlist1(ibox) = nlist1(ibox) + 1
                        list1(nlist1(ibox),ibox) = jbox
                     endif
                  else
c
cc                     boxes in list1 at level ilev+1
c
                     if (ilev.lt.npwlevel) then
                        distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))
     1                      /2.0d0
                        do j=1,mc
                           kbox = itree(iptr(5)+mc*(jbox-1)+j-1)
                           if(itree(iptr(4)+kbox-1).eq.0) then
                              iflist1=1
                              do k=1,ndim
                                 dis = dabs(centers(k,kbox)
     1                               -centers(k,ibox))
                                 if (iperiod .eq. 1) then
                                    dp1 = bs0-dis
                                    if (dp1.lt.dis) dis=dp1
                                 endif
                                 
                                 if (dis.ge.distest) then
                                    iflist1=0
                                    exit
                                 endif
                              enddo
                              if(iflist1.eq.1) then
                                 nlist1(ibox) = nlist1(ibox)+1
                                 list1(nlist1(ibox),ibox) = kbox
                              endif
                           endif
                        enddo
                     endif
                  endif
               enddo
c
cc               compute list1 at level ilev-1 
               do i=1,itree(iptr(6)+dad-1)
                  jbox = itree(iptr(7)+mnbors*(dad-1)+i-1)
                  if(itree(iptr(4)+jbox-1).eq.0) then
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                   2.0d0
                     iflist1=1
                     do k=1,ndim
                        dis = dabs(centers(k,jbox)-centers(k,ibox))
                        if (iperiod .eq. 1) then
                           dp1 = bs0-dis
                           if (dp1.lt.dis) dis=dp1
                        endif
                        if (dis.ge.distest) then
                           iflist1=0
                           exit
                        endif
                     enddo
                     if(iflist1.eq.1) then
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
c      
c
c
c      
c      
c----------------------------------------------------------------
c      
      subroutine pfgt_compute_modified_list1(ndim,npwlevel,ifpwexp,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,iperiod,
     2    mnlist1,nlist1,list1)
c
c
c     This subroutine computes the modified list1 of a given tree
c     structure
c
c     Note: new definition of list 1 - for ilev <= npwlevel,             
c     list1 of a leaf box contains all childless neighbors at or
c     above npwlevel. ifpwexp(ibox)=1, then ibox is excluded from 
c     list1 since then the self interaction is handled via
c     plane wave expansions.
c                  
c     Assume that the tree is level-restricted. Then the
c     new list1 of source ibox contains all target boxes 
c     that require the evaluation of direct interactions with ibox 
c     in the box FGT.  
c      
c      
c     INPUT arguments
c     ndim        in: integer
c                 dimension of the space
c
c     nlevels     in: integer
c                 Number of levels
c
c     npwlevel    in: integer
c                 Cutoff level
c     ifpwexp     in: integer(nboxes)
c                 1: requires pwexp; 0: does not require pwexp
c      
c     nboxes      in: integer
c                 Total number of boxes
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     centers     in: real *8(ndim,nboxes)
c                 xy coordinates of centers of boxes
c   
c     boxsize     in: real *8(0:nlevels)
c                 Array of boxsizes
c   
c     iperiod     in: integer
c                 0: free space; 1: periodic
c 
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i.
c                  
c      
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nlevels,npwlevel,nboxes,ndim
      integer iptr(8),ltree
      integer iperiod
      integer itree(ltree),ifpwexp(nboxes)
      real *8 boxsize(0:nlevels)
      real *8 centers(ndim,nboxes)
      integer mnlist1
      integer nlist1(nboxes), list1(mnlist1,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad,mnbors,mc
      integer i,j,ifirstbox,ilastbox,ii,iflist1,k
      real *8 distest,dis,bs0,dp1

      bs0=boxsize(0)

      mnbors=3**ndim
      mc=2**ndim
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
        nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO      


c     Setting parameters for level = 0
      if(itree(iptr(4)).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif

      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel
      do ilev = 1,nlevend
         ifirstbox = itree(2*ilev+1)
         ilastbox = itree(2*ilev+2)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,dis,distest,k)
C$OMP$PRIVATE(iflist1,dp1)
         do ibox = ifirstbox,ilastbox
            dad = itree(iptr(3)+ibox-1)
c           Compute list1 of ibox if it is childless
            if(itree(iptr(4)+ibox-1).eq.0) then
               do i=1,itree(iptr(6)+ibox-1)
                  jbox = itree(iptr(7)+mnbors*(ibox-1)+i-1)
c
cc                boxes in list 1 at the same level
c
                  if (itree(iptr(4)+jbox-1).eq.0 .or. ilev.eq.npwlevel)
     1                then
c                     nlist1(ibox) = nlist1(ibox) + 1
c                     list1(nlist1(ibox),ibox) = jbox
                     nlist1(jbox) = nlist1(jbox) + 1
                     list1(nlist1(jbox),jbox) = ibox
                  else
c
cc                     boxes in list1 at level ilev+1
c
                     if (ilev.lt.npwlevel) then
                        distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))
     1                      /2.0d0
                        do j=1,mc
                           kbox = itree(iptr(5)+mc*(jbox-1)+j-1)
                           if(itree(iptr(4)+kbox-1).eq.0) then
                              iflist1=1
                              do k=1,ndim
                                 dis = dabs(centers(k,kbox)
     1                               -centers(k,ibox))
                                 if (iperiod .eq. 1) then
                                    dp1 = bs0-dis
                                    if (dp1.lt.dis) dis=dp1
                                 endif
                                 
                                 if (dis.ge.distest) then
                                    iflist1=0
                                    exit
                                 endif
                              enddo
                              if(iflist1.eq.1) then
c                                 nlist1(ibox) = nlist1(ibox)+1
c                                 list1(nlist1(ibox),ibox) = kbox
                                 nlist1(kbox) = nlist1(kbox)+1
                                 list1(nlist1(kbox),kbox) = ibox
                              endif
                           endif
                        enddo
                     endif
                  endif
               enddo
c
cc               compute list1 at level ilev-1 
               do i=1,itree(iptr(6)+dad-1)
                  jbox = itree(iptr(7)+mnbors*(dad-1)+i-1)
                  if(itree(iptr(4)+jbox-1).eq.0) then
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                   2.0d0
                     iflist1=1
                     do k=1,ndim
                        dis = dabs(centers(k,jbox)-centers(k,ibox))
                        if (iperiod .eq. 1) then
                           dp1 = bs0-dis
                           if (dp1.lt.dis) dis=dp1
                        endif
                        if (dis.ge.distest) then
                           iflist1=0
                           exit
                        endif
                     enddo
                     if(iflist1.eq.1) then
c                        nlist1(ibox) = nlist1(ibox)+1
c                        list1(nlist1(ibox),ibox) = jbox
                        nlist1(jbox) = nlist1(jbox)+1
                        list1(nlist1(jbox),jbox) = ibox
                     endif
                  endif
               enddo
            endif
          enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c      
c
c
cc      
c      
c      
c----------------------------------------------------------------
c      
      subroutine gnd_find_pwexp_boxes(ndim,npwlevel,nboxes,
     1    nlevels,ltree,itree,iptr,iper,ifpwexp)
c
c
c     Determine whether a box needs plane wave expansions
c     in the box fgt.
c
c     1. all boxes below npwlevel need plane wave expansions
c     2. at the cutoff level, i.e., npwlevel, a box needs
c     plane wave expansion if one of its colleagues is a nonleaf 
c     box. 
c                  
c     Thus, a leaf box at the cutoff level may or may not 
c     have plane wave expansion.
c     But if it has plane wave expansion, then the self interaction
c     is handled by the plane wave expansion instead of 
c     direct evaluation. 
c      
c      
c     INPUT arguments
c     npwlevel    in: integer
c                 Cutoff level
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     nlevels     in: integer
c                 Number of levels
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     iper        in: integer
c                 0: free space; 1: periodic
c 
c--------------------------------------------------------------
c     OUTPUT arguments:
c     ifpwexp     out: integer(nboxes)
c                 ifpwexp(ibox)=1, ibox needs plane wave expansion
c                 ifpwexp(ibox)=0, ibox does not need pwexp
c      
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nlevels,npwlevel,nboxes,ndim
      integer iptr(8),ltree
      integer iper
      integer itree(ltree)
      integer ifpwexp(nboxes)

      mnbors=3**ndim
      
      do i=1,nboxes
         ifpwexp(i)=0
      enddo

      if (npwlevel .le. 0) then
         ifpwexp(1)=1
         return
      endif
      
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
c     At the cutoff level, a box will have plane wave expansion if 1. it's a nonleaf
c     box; or 2. it's a leaf box but has a nonleaf colleague.
c     Note that we also use plane wave expansion to compute the self interaction of 
c     this box as well
         do ilev=npwlevel,npwlevel
            do 1000 ibox=itree(2*ilev+1),itree(2*ilev+2)
               nchild = itree(iptr(4)+ibox-1)
               if (nchild .eq. 0) then
                  ncoll = itree(iptr(6)+ibox-1)
                  do j=1,ncoll
                     jbox = itree(iptr(7) + (ibox-1)*mnbors+j-1)
                     nchild = itree(iptr(4)+jbox-1)
c                     jlev = itree(iptr(2)+jbox-1)
c                     if (nchild .gt. 0 .and. jlev.eq.ilev) then
                     if (nchild .gt. 0) then
                        ifpwexp(ibox)=1
                        exit
                     endif
                  enddo
               else
                  ifpwexp(ibox)=1
               endif
 1000       continue
         enddo
c     below the cutoff level, all boxes are handled via plane wave expansions
         do ilev=npwlevel+1,nlevels
            do ibox=itree(2*ilev+1),itree(2*ilev+2)
               ifpwexp(ibox)=1
            enddo
         enddo
      endif      

      return
      end
c      
c
c
c
c----------------------------------------------------------------
c      
      subroutine pfgt_find_pwexp_boxes(ndim,npwlevel,nboxes,
     1    nlevels,ltree,itree,iptr,ifpwexp)
c
c
c     Determine whether a box needs plane wave expansions
c     in the point fgt.
c
c     At the cutoff level, i.e., npwlevel, a box needs
c     plane wave expansion if one of its colleagues is a nonleaf 
c     box. 
c                  
c     Thus, a leaf box at the cutoff level may or may not 
c     have plane wave expansion.
c     But if it has plane wave expansion, then the self interaction
c     is handled by the plane wave expansion instead of 
c     direct evaluation. 
c      
c      
c     INPUT arguments
c     ndim        in: integer
c                 dimension of the space
c
c     npwlevel    in: integer
c                 Cutoff level
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     nlevels     in: integer
c                 Number of levels
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c--------------------------------------------------------------
c     OUTPUT arguments:
c     ifpwexp     out: integer(nboxes)
c                 ifpwexp(ibox)=1, ibox needs plane wave expansion
c                 ifpwexp(ibox)=0, ibox does not need pwexp
c      
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nlevels,npwlevel,nboxes,ndim
      integer iptr(8),ltree
      integer itree(ltree)
      integer ifpwexp(nboxes)

      mnbors=3**ndim
      
      do i=1,nboxes
         ifpwexp(i)=0
      enddo

      if (npwlevel .le. 0) then
         ifpwexp(1)=1
         return
      endif
      
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
c     At the cutoff level, a box will have plane wave expansion if 1. it's a nonleaf
c     box; or 2. it's a leaf box but has a nonleaf colleague.
c     Note that we also use plane wave expansion to compute the self interaction of 
c     this box as well
         do ilev=npwlevel,npwlevel
            do 1000 ibox=itree(2*ilev+1),itree(2*ilev+2)
               ncoll = itree(iptr(6)+ibox-1)
               do j=1,ncoll
                  jbox = itree(iptr(7) + (ibox-1)*mnbors+j-1)
                  nchild = itree(iptr(4)+jbox-1)
c                     jlev = itree(iptr(2)+jbox-1)
c                     if (nchild .gt. 0 .and. jlev.eq.ilev) then
                  if (nchild .gt. 0) then
                     ifpwexp(ibox)=1
                     exit
                  endif
               enddo
 1000       continue
         enddo
      endif      

      return
      end
c      
c
c
c
      subroutine gnd_compute_mnlistpw(ndim,nboxes,nlevels,ltree,itree,
     1   iptr,centers,boxsize,mnlistpw)
c
c     determine maximum number of elements in listpw
c
c     NOTE: we use max values
c
      implicit real *8 (a-h,o-z)
      integer iptr(8),ltree
      integer nlevels,nboxes,itree(ltree)
      real *8 centers(ndim,nboxes),boxsize(0:nlevels)
      integer mnlistpw

      mnlistpw = 3**ndim-1

      return
      end
c
c
c
c
c
      subroutine pfgt_compute_listpw(ndim,npwlevel,nboxes,nlevels,
     1    ltree,itree,iptr,centers,boxsize,laddr,
     2    mnlistpw,nlistpw,listpw)
c     this subroutine returns lists of pw interaction boxes for point FGT
c
c     input parameters:
c
c     nlevels = number of levels
c     npwlevel   = the cutoff level
c     nboxes = total number of boxes in the tree
c     itree - laddr: tree stuff
c
c     mnlistpw = maximum number of PW interaction boxes at the cutoff level
c      
c     output parameters:
c
c     nlistpw (nboxes) = contains the number of boxes for the PW interaction
c                            for each box
c     listpw (mnlistpw,nboxes) = contains the box ID for the PW interaction
c                            for each box
c      
      implicit real *8 (a-h,o-z)
      integer nlevels,nboxes
      integer iptr(8),ltree
      integer itree(ltree),laddr(2,0:nlevels)
      real *8 centers(ndim,nboxes),boxsize(0:nlevels)
      integer mnlistpw
      integer nlistpw(nboxes), listpw(mnlistpw,nboxes)

      integer ibox

      mnbors=3**ndim
      
      do ibox=1,nboxes
         nlistpw(ibox)=0
      enddo

      if (npwlevel.gt.nlevels) return
      
c     make sure ncutoff is >= 0
      ncutoff=max(npwlevel,0)
      do ilev=ncutoff,ncutoff
        do ibox = laddr(1,ilev),laddr(2,ilev)
           ncoll = itree(iptr(6)+ibox-1)
           do i=1,ncoll
              jbox = itree(iptr(7) + (ibox-1)*mnbors+i-1)
c              jlev = itree(iptr(2)+jbox-1)
c              if (jbox.ne.ibox .and. ilev.eq.jlev .and. 
              if (jbox.ne.ibox .and. 
     1            (itree(iptr(4)+ibox-1).gt.0)) then
c             The pw list of ibox at the cutoff level contains its
c             colleague jbox at the cutoff level if either of them
c             has children. If both of them are leaf boxes, their
c             interaction will be computed directly. Self interaction
c             is always directly copied from mp exp to the loc exp and
c             thus is not in the pw list.
c                 nlistpw(ibox)=nlistpw(ibox)+1
c                 listpw(nlistpw(ibox),ibox) = jbox
                 nlistpw(jbox)=nlistpw(jbox)+1
                 listpw(nlistpw(jbox),jbox) = ibox
              endif
           enddo
cccc           print *, ibox, ncoll, nlistpw(ibox)
        enddo
c     end of ilev do loop
      enddo

      return
      end
c
c
c
c
c
c
      subroutine gnd_compute_listpw(ndim,npwlevel,nboxes,nlevels,
     1    ltree,itree,iptr,centers,boxsize,laddr,
     2    mnlistpw,nlistpw,listpw)
c     this subroutine returns lists of pw interaction boxes for the box FGT
c
c     input parameters:
c
c     nlevels = number of levels
c     npwlevel   = the cutoff level
c     nboxes = total number of boxes in the tree
c     itree - laddr: tree stuff
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
      integer nlevels,nboxes
      integer iptr(8),ltree
      integer itree(ltree),laddr(2,0:nlevels)
      real *8 centers(ndim,nboxes),boxsize(0:nlevels)
      integer mnlistpw
      integer nlistpw(nboxes), listpw(mnlistpw,nboxes)

      integer ibox

      mnbors=3**ndim
      
      do ibox=1,nboxes
         nlistpw(ibox)=0
      enddo

      if (npwlevel.gt.nlevels) return
      
c     make sure ncutoff is >= 0
      ncutoff=max(npwlevel,0)
      do ilev=ncutoff,ncutoff
        do ibox = laddr(1,ilev),laddr(2,ilev)
           ncoll = itree(iptr(6)+ibox-1)
           do i=1,ncoll
              jbox = itree(iptr(7) + (ibox-1)*mnbors+i-1)
c              jlev = itree(iptr(2)+jbox-1)
c              if (jbox.ne.ibox .and. ilev.eq.jlev .and. 
              if (jbox.ne.ibox .and. 
     1            (itree(iptr(4)+jbox-1).gt.0
     2            .or. itree(iptr(4)+ibox-1).gt.0)) then
c             The pw list of ibox at the cutoff level contains its
c             colleague jbox at the cutoff level if either of them
c             has children. If both of them are leaf boxes, their
c             interaction will be computed directly. Self interaction
c             is always directly copied from mp exp to the loc exp and
c             thus is not in the pw list.
                 nlistpw(ibox)=nlistpw(ibox)+1
                 listpw(nlistpw(ibox),ibox) = jbox
              endif
           enddo
        enddo
c     end of ilev do loop
      enddo

      return
      end
c
c
c
c
c
      subroutine oldtree2newtree2d(nlev,levelbox,iparentbox,
     2    ichildbox,icolbox,irowbox,nboxes,nblevel,
     3    iboxlev,istartlev,cent0,xsize0,iperiod,
     4    ltree,nlevels,itree,iptr,centers,boxsize)
      implicit real *8 (a-h,o-z)
c
c     convert an old tree to the new tree
c
c     input parameters:
c     nlev - total number of  levels
c     levelbox - an array determining the level of each box
c     iparentbox - the parent of each box
c     ichildbox - the four children of each box
c     icolbox - the column of each box
c     irowbox - the row of each box
c     nboxes - integer
c          number of boxes
c     nblevel - the total number of boxes per level
c     iboxlev - the array in which the boxes are arranged
c     istartlev - the pointer to where each level
c               begins in the iboxlev array
c     cent0 - center of the root box
c     xsize0 - size of the root box
c     iperiod = 0 : free space
c               1: periodic
c     ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c
c     output:
c     nlevels - integer
c          number of levels
c     itree - integer(ltree)
c          tree info
c     iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c     centers - double precision (2,nboxes)
c          xy coordinates of box centers in the oct tree
c     boxsize - double precision (0:nlevels)
c          size of box at each of the levels
      integer *4  levelbox(1)
      integer *4  nlev, nboxes
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 cent0(2),xsize0
      real *8 centers(2,*),boxsize(0:1)
      integer iboxlevinv(nboxes)
      integer, allocatable :: icolleagbox(:,:)

      allocate(icolleagbox(9,nboxes))
      
      nlevels = nlev
      
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + 4*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + 9*nboxes



      boxsize(0) = xsize0

      centers(1,1) = cent0(1)
      centers(2,1) = cent0(2)
      
cccc      call prinf('iboxlev=*',iboxlev,nboxes)
cccc      call prinf('istartlev=*',istartlev,nlev+1)
      
      do i=1,nlevels
         boxsize(i)=boxsize(i-1)/2
      enddo

      do i=1,nboxes
         ibox=iboxlev(i)
         iboxlevinv(ibox)=i
      enddo
      
      do ilev=0,nlevels
         do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
            ibox=iboxlev(i)
            icol=icolbox(ibox)
            irow=irowbox(ibox)
         
            bs = boxsize(ilev)
         
            centers(1,i) = cent0(1)-xsize0/2 + (icol-0.5d0) * bs
            centers(2,i) = cent0(2)-xsize0/2 + (irow-0.5d0) * bs
         enddo
      enddo

cccc      call mkcolls(icolbox,
cccc     1      irowbox, icolleagbox, nboxes, nlev,
cccc     2      iparentbox, ichildbox, nblevel,
cccc     3      iboxlev, istartlev, iperiod)

      do ilev=0,nlevels
         itree(iptr(1)+2*ilev)=istartlev(ilev)
         itree(iptr(1)+2*ilev+1)=istartlev(ilev)+nblevel(ilev)-1
      enddo

      do ilev=0,nlevels
        do i=istartlev(ilev),istartlev(ilev)+nblevel(ilev)-1
          ibox=iboxlev(i)
          itree(iptr(2)+i-1) = levelbox(ibox)
cccc          print *, ilev, levelbox(ibox)
          
          jbox=iparentbox(ibox)
          if (jbox.gt.0) then
             itree(iptr(3)+i-1) = iboxlevinv(jbox)
          else
             itree(iptr(3)+i-1) = -1
          endif

          do j=1,4
             itree(iptr(5)+4*(i-1)+j-1)=-1
          enddo
          
          ichild=0
          do j=1,4
             jbox=ichildbox(j,ibox)
             if (jbox.gt.0) then
                itree(iptr(5)+4*(i-1)+ichild)=iboxlevinv(jbox)
                ichild=ichild+1
             endif 
          enddo
          itree(iptr(4)+i-1) = ichild

          do j=1,9
             itree(iptr(7)+9*(i-1)+j-1)=-1
          enddo
          
          icoll=0
          do j=1,9
             jbox=icolleagbox(j,ibox)
             if (jbox.gt.0) then
                itree(iptr(7)+9*(i-1)+icoll)=iboxlevinv(jbox)
                icoll=icoll+1
             endif
          enddo
          itree(iptr(6)+i-1)=icoll
        enddo
      enddo

cccc      call prinf('itree=*',itree,ltree)
      
      return
      end subroutine
      
