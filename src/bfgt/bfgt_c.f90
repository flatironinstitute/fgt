! src/bfgt/bfgt_c.f90
!
! C-callable shim for the boxfgt Fortran driver. Exposed symbol: fgt_bfgt_d.
!
! For this initial wrapper we expose ONLY the post-tree-built signature.
! The Python wrapper is responsible for tree construction (or for now,
! exposes only pfgt; box-FGT Python support is a follow-up).

module fgt_c_bfgt
  use, intrinsic :: iso_c_binding
  implicit none

contains

  subroutine fgt_bfgt_d(                                              &
       nd, ndim, delta, eps, ipoly, iperiod, norder, npbox,           &
       nboxes, nlevels, ltree, itree, iptr,                           &
       centers, boxsize, fvals,                                       &
       ifpgh, pot, grad, hess,                                        &
       ifnewtree,                                                     &
       ntarg, targs, ifpghtarg, pote, grade, hesse,                   &
       timeinfo) bind(C, name="fgt_bfgt_d")
    integer(c_int), value :: nd, ndim, ipoly, iperiod, norder, npbox
    integer(c_int), value :: nboxes, nlevels, ltree
    integer(c_int), value :: ifpgh, ifnewtree, ntarg, ifpghtarg
    real(c_double), value :: delta, eps
    integer(c_int), intent(inout) :: itree(ltree)
    integer(c_int), intent(inout) :: iptr(8)
    real(c_double), intent(inout) :: centers(ndim, nboxes)
    real(c_double), intent(inout) :: boxsize(0:nlevels)
    real(c_double), intent(inout) :: fvals(nd, npbox, nboxes)
    real(c_double), intent(out)   :: pot(nd, npbox, nboxes)
    real(c_double), intent(out)   :: grad(nd, ndim, npbox, max(nboxes,1))
    real(c_double), intent(out)   :: hess(nd, ndim*(ndim+1)/2, npbox,    &
                                          max(nboxes,1))
    real(c_double), intent(in)    :: targs(ndim, max(ntarg,1))
    real(c_double), intent(out)   :: pote(nd, max(ntarg,1))
    real(c_double), intent(out)   :: grade(nd, ndim, max(ntarg,1))
    real(c_double), intent(out)   :: hesse(nd, ndim*(ndim+1)/2, max(ntarg,1))
    real(c_double), intent(out)   :: timeinfo(20)

    call boxfgt(nd, ndim, delta, eps, ipoly, iperiod, norder, npbox,   &
         nboxes, nlevels, ltree, itree, iptr,                          &
         centers, boxsize, fvals,                                      &
         ifpgh, pot, grad, hess, ifnewtree,                            &
         ntarg, targs, ifpghtarg, pote, grade, hesse, timeinfo)
  end subroutine fgt_bfgt_d

end module fgt_c_bfgt
