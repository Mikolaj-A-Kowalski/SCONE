module involuteUniverse_test
  use numPrecision
  use pFUnit_mod
  use involuteUniverse_class, only : phase_and_derivative, involute_newton

  implicit none


contains

! @Test
!   subroutine test_example()
!     real(defReal)     :: rb, a0, rl, theta_l
!     real(defReal)     :: d_min, d_max
!     real(defReal)     :: f, df, d
!     integer(shortInt) :: N, i

!     ! rb = 1.0_defReal
!     ! a0 = PI / 2.0
!     ! rl = 0.5_defReal
!     ! theta_l = PI / 4.0
!     ! N = 1000
!     ! d_max = 10.0_defReal
!     ! d_min = -10.0_defReal
!     ! print *, "HELLO"
!     ! do i = 1, N
!     !   d = d_min + (d_max - d_min) * (i - 1) / real(N - 1, defReal)
!     !   call phase_and_derivative(f, df, rb, a0, rl, theta_l, d)
!     !   print *, d, f, df
!     ! end do

!   end subroutine test_example

@Test
  subroutine test_newton_with_complex_gap()
    real(defReal), parameter :: rb = 1.5_defReal
    real(defReal), parameter :: a0 = PI / 2.0_defReal
    real(defReal), parameter :: rl = 0.5_defReal
    real(defReal), parameter :: theta_l = PI / 4.0_defReal
    real(defReal), parameter :: TOL = 1.0E-7_defReal
    real(defReal) :: d_guess, d_res, rhs

    !! Solve for the 0th Root
    rhs = ZERO
    d_guess = -1.0_defReal
    d_res = involute_newton(rb, a0, rl, theta_l, d_guess, rhs)
    @assertEqual(-0.5000_defReal, d_res,  TOL)

    !! Solve for the 2nd 0th root
    d_guess = -1.52_defReal
    d_res = involute_newton(rb, a0, rl, theta_l, d_guess, rhs)
    @assertEqual(-2.825952750312339_defReal, d_res,  TOL)

    !! Solve for the -2PI root on RHS of maximum
    rhs = -2.0_defReal * PI
    d_guess = 2.5_defReal
    d_res = involute_newton(rb, a0, rl, theta_l, d_guess, rhs)
    @assertEqual(8.18536592255579_defReal, d_res,  TOL)


    !! Solve for the -4PI root on LHS of maximum
    rhs = -4.0_defReal * PI
    d_guess = -2.5000_defReal
    d_res = involute_newton(rb, a0, rl, theta_l, d_guess, rhs)
    @assertEqual(-22.294137957839933_defReal, d_res,  TOL)

  end subroutine test_newton_with_complex_gap



end module involuteUniverse_test
