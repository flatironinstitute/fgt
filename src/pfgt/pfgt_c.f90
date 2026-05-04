! src/pfgt/pfgt_c.f90
!
! C-callable shim for the pfgt Fortran driver. Exposed symbol: fgt_pfgt_d
! (double precision, the only precision the underlying lib supports today;
! the _d suffix is forward-looking).
!
! Argument layout matches pfgt() exactly; the C-side allocator is responsible
! for sizing all output arrays. All real arrays are double precision; all
! integers are int32.

module fgt_c_pfgt
  use, intrinsic :: iso_c_binding
  implicit none

contains

  subroutine fgt_pfgt_d(                                              &
       nd, ndim, delta, eps, iperiod, bs0, cen0,                      &
       ns, sources,                                                   &
       ifcharge, charge, ifdipole, rnormal, dipstr,                   &
       ifpgh, pot, grad, hess,                                        &
       nt, targ, ifpghtarg, pottarg, gradtarg, hesstarg)              &
       bind(C, name="fgt_pfgt_d")
    integer(c_int), value :: nd, ndim, iperiod, ns, nt
    integer(c_int), value :: ifcharge, ifdipole, ifpgh, ifpghtarg
    real(c_double), value :: delta, eps, bs0
    real(c_double), intent(in)  :: cen0(ndim)
    real(c_double), intent(in)  :: sources(ndim, ns)
    real(c_double), intent(in)  :: charge(nd, max(ns, 1))
    real(c_double), intent(in)  :: rnormal(ndim, max(ns, 1))
    real(c_double), intent(in)  :: dipstr(nd, max(ns, 1))
    real(c_double), intent(out) :: pot(nd, max(ns, 1))
    real(c_double), intent(out) :: grad(nd, ndim, max(ns, 1))
    real(c_double), intent(out) :: hess(nd, ndim*(ndim+1)/2, max(ns, 1))
    real(c_double), intent(in)  :: targ(ndim, max(nt, 1))
    real(c_double), intent(out) :: pottarg(nd, max(nt, 1))
    real(c_double), intent(out) :: gradtarg(nd, ndim, max(nt, 1))
    real(c_double), intent(out) :: hesstarg(nd, ndim*(ndim+1)/2, max(nt, 1))

    call pfgt(nd, ndim, delta, eps, iperiod, bs0, cen0,                &
         ns, sources,                                                  &
         ifcharge, charge, ifdipole, rnormal, dipstr,                  &
         ifpgh, pot, grad, hess, nt, targ,                             &
         ifpghtarg, pottarg, gradtarg, hesstarg)
  end subroutine fgt_pfgt_d

end module fgt_c_pfgt
