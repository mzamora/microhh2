! Wrapper inspired in ../src_rrtmg/rrtmg_lw_wrapper.f90
! which was imported from PyCLES
! Add approppriate copyrights

module fishpack_wrapper

use iso_c_binding, only: c_double, c_int
! USE fish
implicit none

contains

subroutine c_blktri &
            (iflg, np, n, an, bn, cn, mp, &
             m, am, bm, cm, idimy, y, &
             ierror,w,k) bind(c,name="c_blktri")
      integer(c_int), intent(in):: iflg       ! initialization flag
                                              ! 0: initialization only. work array w
                                              ! 1: initialized quantities used to solve x
      integer(c_int), intent(in):: np         ! 0: if an(1) and cn(n) are not zero: periodic BCs
                                              ! 1: if an(1) and cn(n) are zero
      integer(c_int), intent(in):: n          ! number of unknowns in the j direction
      real(c_double), intent(in):: an(n)      ! arrays of length n ~ coefficients
      real(c_double), intent(in):: bn(n)
      real(c_double), intent(in):: cn(n)
      real:: an_f(n)
      real:: bn_f(n)
      real:: cn_f(n)
      integer(c_int), intent(in):: mp         ! 0: if an(1) and cn(n) are not zero: periodic BCs
                                              ! 1: if an(1) and cn(n) are zero
      integer(c_int), intent(in):: m          ! number of unknowns in the i direction
      real(c_double), intent(in):: am(m)      ! arrays of length m ~ coefficients
      real(c_double), intent(in):: bm(m)
      real(c_double), intent(in):: cm(m)
      real:: am_f(m)
      real:: bm_f(m)
      real:: cm_f(m)
      integer(c_int), intent(in):: idimy      ! row dimension of the 2d array Y, at least m
      real(c_double), intent(inout):: y(m,n)  ! in: 2D array for the right hand side
                                              ! out: solution X
      real:: y_f(m,n)
      integer(c_int), intent(in) ::k
      real(c_double), intent(inout):: w(k)    ! IF NP=1 DEFINE K=INT(LOG2(N))+1 AND
                                              ! SET L=2**(K+1) THEN W MUST HAVE DIMENSION
                                              ! (K-2)*L+K+5+MAX(2N,6M)
      ! TYPE (fishworkspace), SAVE :: W
      real:: w_f(k)
      integer(c_int), intent(inout):: ierror     ! error flag
                                              ! 0  NO ERROR.
                                              ! 1  M IS LESS THAN 5
                                              ! 2  N IS LESS THAN 5
                                              ! 3  IDIMY IS LESS THAN M.
                                              ! 4  BLKTRI FAILED WHILE COMPUTING RESULTS
                                              !    THAT DEPEND ON THE COEFFICIENT ARRAYS
                                              !    AN, BN, CN.  CHECK THESE ARRAYS.
                                              ! 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.
    integer :: j, kk, ii

    an_f=real(an); bn_f=real(bn); cn_f=real(cn); y_f=real(y);
    am_f=real(am); bm_f=real(bm); cm_f=real(cm); w_f=real(w);
    ! PRINT *, "In Fortran, before calling blktri"
    ! PRINT *, "[fishpack_wrapper.f90] iflg=", iflg
    ! PRINT *, "[fishpack_wrapper.f90] np=", np
    ! PRINT *, "[fishpack_wrapper.f90] n=", n
    ! PRINT *, "[fishpack_wrapper.f90] an(1)=", an(1)
    ! PRINT *, "[fishpack_wrapper.f90] an(n)=", an(n)
    ! PRINT *, "[fishpack_wrapper.f90] bn(1)=", bn(1)
    ! PRINT *, "[fishpack_wrapper.f90] bn(n)=", bn(n)
    ! PRINT *, "[fishpack_wrapper.f90] cn(1)=", cn(1)
    ! PRINT *, "[fishpack_wrapper.f90] cn(n)=", cn(n)
!    PRINT *, "[fishpack_wrapper.f90] mp=", mp
!    PRINT *, "[fishpack_wrapper.f90] m=", m
!     PRINT *, "[fishpack_wrapper.f90] am(1)=", am(1)
!     PRINT *, "[fishpack_wrapper.f90] am(m)=", am(m)
!     PRINT *, "[fishpack_wrapper.f90] bm(1)=", bm(1)
!     PRINT *, "[fishpack_wrapper.f90] bm(m)=", bm(m)
!     PRINT *, "[fishpack_wrapper.f90] cm(1)=", cm(1)
!     PRINT *, "[fishpack_wrapper.f90] cm(m)=", cm(m)
! !    PRINT *, "[fishpack_wrapper.f90] idimy=", idimy
!     PRINT *, "[fishpack_wrapper.f90] y(1,1)=", y_f(1,1)
!     PRINT *, "[fishpack_wrapper.f90] y(m,1)=", y_f(m,1)
!     PRINT *, "[fishpack_wrapper.f90] y(1,n)=", y_f(1,n)
!     PRINT *, "[fishpack_wrapper.f90] y(m,n)=", y_f(m,n)
!    PRINT *, "[fishpack_wrapper.f90] ierror", ierror
!    PRINT *, "[fishpack_wrapper.f90] w(1)=", w_f(1)
!    PRINT *, "[fishpack_wrapper.f90] w(k)=", w_f(k)
!    PRINT *, "[fishpack_wrapper.f90] k", k

    !call blktri &
    !        (iflg, np, n, an, bn, cn, mp, &
    !         m, am, bm, cm, idimy, y, &
    !         ierror, w)
    call blktri &
            (iflg, np, n, an_f, bn_f, cn_f, mp, &
             m, am_f, bm_f, cm_f, idimy, y_f, &
             ierror, w_f)

     y=real(y_f,c_double);
     w=real(w_f,c_double);

    PRINT *, "In Fortran, after calling blktri"
!    PRINT *, "[fishpack_wrapper.f90] iflg=", iflg
!    PRINT *, "[fishpack_wrapper.f90] np=", np
!    PRINT *, "[fishpack_wrapper.f90] n=", n
!    PRINT *, "[fishpack_wrapper.f90] an(1)=", an(1)
!    PRINT *, "[fishpack_wrapper.f90] an(n)=", an(n)
!    PRINT *, "[fishpack_wrapper.f90] bn(1)=", bn(1)
!    PRINT *, "[fishpack_wrapper.f90] bn(n)=", bn(n)
!    PRINT *, "[fishpack_wrapper.f90] cn(1)=", cn(1)
!    PRINT *, "[fishpack_wrapper.f90] cn(n)=", cn(n)
!    PRINT *, "[fishpack_wrapper.f90] mp=", mp
!    PRINT *, "[fishpack_wrapper.f90] m=", m
!    PRINT *, "[fishpack_wrapper.f90] am(1)=", am(1)
!    PRINT *, "[fishpack_wrapper.f90] am(m)=", am(m)
!    PRINT *, "[fishpack_wrapper.f90] bm(1)=", bm(1)
!    PRINT *, "[fishpack_wrapper.f90] bm(m)=", bm(m)
!    PRINT *, "[fishpack_wrapper.f90] cm(1)=", cm(1)
!    PRINT *, "[fishpack_wrapper.f90] cm(m)=", cm(m)
!    PRINT *, "[fishpack_wrapper.f90] idimy=", idimy
    ! PRINT *, "[fishpack_wrapper.f90] y(1,1)=", y(1,1)
    ! PRINT *, "[fishpack_wrapper.f90] y(32,32)=", y(32,32)
    ! PRINT *, "[fishpack_wrapper.f90] y(1,n)=", y(1,n)
    ! PRINT *, "[fishpack_wrapper.f90] y(m,n)=", y(m,n)
    ! PRINT *, "[fishpack_wrapper.f90] ierror", ierror
    ! PRINT *, "[fishpack_wrapper.f90] w(1)=", w(1)
    ! PRINT *, "[fishpack_wrapper.f90] w(k)=", w(k)

end subroutine c_blktri

end module
