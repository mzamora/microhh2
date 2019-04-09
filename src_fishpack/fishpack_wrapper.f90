! Wrapper inspired in ../src_rrtmg/rrtmg_lw_wrapper.f90
! which was imported from PyCLES
! Add approppriate copyrights

module fishpack_wrapper

use iso_c_binding, only: c_double, c_int

implicit none

contains

subroutine c_blktri &
            (iflg, np, n, an, bn, cn, mp, &
             m, am, bm, cm, idimy, y, &
             ierror, w, k) bind(c,name="c_blktri")
      integer(c_int), intent(in):: iflg       ! initialization flag
                                              ! 0: initialization only. work array w 
                                              ! 1: initialized quantities used to solve x
      integer(c_int), intent(in):: np         ! 0: if an(1) and cn(n) are not zero: periodic BCs
                                              ! 1: if an(1) and cn(n) are zero
      integer(c_int), intent(in):: n          ! number of unknowns in the j direction
      real(c_double), intent(in):: an(n)      ! arrays of length n ~ coefficients
      real(c_double), intent(in):: bn(n) 
      real(c_double), intent(in):: cn(n) 
      integer(c_int), intent(in):: mp         ! 0: if an(1) and cn(n) are not zero: periodic BCs
                                              ! 1: if an(1) and cn(n) are zero
      integer(c_int), intent(in):: m          ! number of unknowns in the i direction
      real(c_double), intent(in):: am(m)      ! arrays of length m ~ coefficients
      real(c_double), intent(in):: bm(m)                                   
      real(c_double), intent(in):: cm(m) 
      integer(c_int), intent(in):: idimy      ! row dimension of the 2d array Y, at least m
      real(c_double), intent(inout):: y(m,n)  ! in: 2D array for the right hand side
                                              ! out: solution X
      integer(c_int), intent(in) ::k
      real(c_double), intent(inout):: w(k)    ! IF NP=1 DEFINE K=INT(LOG2(N))+1 AND          
                                              ! SET L=2**(K+1) THEN W MUST HAVE DIMENSION    
                                              ! (K-2)*L+K+5+MAX(2N,6M)       
      integer(c_int), intent(in):: ierror     ! error flag
                                              ! 0  NO ERROR.                                 
                                              ! 1  M IS LESS THAN 5                          
                                              ! 2  N IS LESS THAN 5                          
                                              ! 3  IDIMY IS LESS THAN M.                     
                                              ! 4  BLKTRI FAILED WHILE COMPUTING RESULTS     
                                              !    THAT DEPEND ON THE COEFFICIENT ARRAYS     
                                              !    AN, BN, CN.  CHECK THESE ARRAYS.          
                                              ! 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
      real:: an_f(n) !Extra arrays to avoid variable type issues for real and c_doubles
      real:: bn_f(n)
      real:: cn_f(n)
      real:: am_f(m)
      real:: bm_f(m)
      real:: cm_f(m)
      real:: y_f(m,n)
      real:: w_f(k)

    ! Transform c_double to real
    an_f=real(an); bn_f=real(bn); cn_f=real(cn); y_f=real(y);
    am_f=real(am); bm_f=real(bm); cm_f=real(cm); w_f=real(w);

    call blktri &
            (iflg, np, n, an_f, bn_f, cn_f, mp, &
             m, am_f, bm_f, cm_f, idimy, y_f, &
             ierror, w_f)

    ! Transform real to c_double (a,b,c arrays shouldn't matter)
    y=real(y_f,c_double);
    w=real(w_f,c_double);

end subroutine c_blktri

end module
