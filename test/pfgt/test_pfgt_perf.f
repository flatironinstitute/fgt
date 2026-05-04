c     test/pfgt/test_pfgt_perf.f
c
c     Performance harness for pfgt. Reads dimension, nsrc, eps from
c     environment variables FGT_DIM (2 or 3), FGT_NSRC, FGT_EPS, FGT_DELTA.
c     Prints one CSV line to stdout (last line of program output):
c       nthreads,nsrc,dim,eps,total_s
c     so the wrapper script can `tail -1` it cleanly.
c
c     Validation against direct sum is intentionally skipped here -
c     correctness is covered by test_pfgt_all.f. This driver only
c     exercises the dispatch + kernel path under varying OMP_NUM_THREADS.

      program test_pfgt_perf
      implicit real *8 (a-h,o-z)
      integer dim, nd, nsrc, ntarg, iperiod
      integer ifcharge, ifdipole, ifpgh, ifpghtarg
      integer nthreads, k, j
      real *8 eps, delta, t1, t2, dummy
      character*32 buf

      integer omp_get_max_threads
      real *8 omp_get_wtime
      external omp_get_max_threads, omp_get_wtime
      real *8 hkrand
      external hkrand

      real *8, allocatable :: sources(:,:), targ(:,:)
      real *8, allocatable :: charge(:,:), rnormal(:,:), dipstr(:,:)
      real *8, allocatable :: pot(:,:), grad(:,:,:), hess(:,:,:)
      real *8, allocatable :: pottarg(:,:), gradtarg(:,:,:),
     1                              hesstarg(:,:,:)

      real *8 bs0, cen0(3)

c     ---- read config from env ----
      call get_environment_variable('FGT_DIM', buf)
      if (len_trim(buf) .eq. 0) then
         dim = 3
      else
         read(buf, *) dim
      endif

      call get_environment_variable('FGT_NSRC', buf)
      if (len_trim(buf) .eq. 0) then
         nsrc = 1000000
      else
         read(buf, *) nsrc
      endif

      call get_environment_variable('FGT_EPS', buf)
      if (len_trim(buf) .eq. 0) then
         eps = 1d-6
      else
         read(buf, *) eps
      endif

      call get_environment_variable('FGT_DELTA', buf)
      if (len_trim(buf) .eq. 0) then
         delta = 1d-4
      else
         read(buf, *) delta
      endif

      nd = 1
      ntarg = 0
      iperiod = 0
      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 0

      bs0 = 1d0
      do k = 1, dim
         cen0(k) = 0.5d0
      enddo

c     ---- generate uniformly random sources in [0,1]^dim ----
      allocate(sources(dim, nsrc))
      allocate(charge(nd, nsrc))
      allocate(rnormal(dim, nsrc))
      allocate(dipstr(nd, nsrc))
      allocate(pot(nd, nsrc))
      allocate(grad(nd, dim, nsrc))
      allocate(hess(nd, dim*(dim+1)/2, nsrc))
      allocate(targ(dim, max(ntarg,1)))
      allocate(pottarg(nd, max(ntarg,1)))
      allocate(gradtarg(nd, dim, max(ntarg,1)))
      allocate(hesstarg(nd, dim*(dim+1)/2, max(ntarg,1)))

c     seed hkrand
      dummy = hkrand(7654321)
      do j = 1, nsrc
         do k = 1, dim
            sources(k, j) = hkrand(0)
         enddo
         charge(1, j) = hkrand(0) - 0.5d0
      enddo

c     ---- run pfgt and time it ----
      t1 = omp_get_wtime()

      call pfgt(nd, dim, delta, eps, iperiod, bs0, cen0,
     1     nsrc, sources,
     1     ifcharge, charge, ifdipole, rnormal, dipstr,
     2     ifpgh, pot, grad, hess, ntarg, targ,
     3     ifpghtarg, pottarg, gradtarg, hesstarg)

      t2 = omp_get_wtime()

c     ---- emit single CSV line at the very end ----
      nthreads = omp_get_max_threads()
      write(*, '(I0, A, I0, A, I0, A, ES9.2, A, F10.4)')
     $     nthreads, ',', nsrc, ',', dim, ',', eps, ',', t2 - t1

      stop
      end
