module involuteUniverse_test
  use numPrecision
  use universalVariables
  use pFUnit_mod
  use involuteUniverse_class, only : phase_and_derivative, involute_newton, involute_distance

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

@Test
  subroutine test_distance_planar()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: rb, a0, d, theta
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    !! Origin @ maximum
    rb = 1.0_defReal
    a0 = 0.0_defReal
    r = [-1.0_defReal, 0.0_defReal, 0.0_defReal]
    theta = 0.0_defReal
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(1.0_defReal, d, TOL)

    !! Tangent @ maximum
    rb = 2.0_defReal
    a0 = 0.0_defReal
    r = [-2.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = 0.0_defReal
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.141592653589793_defReal, d, TOL)

    !! Complex Gap Hit
    rb = 2.0_defReal
    a0 = 0.0_defReal
    r = [-2.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.15_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(4.405378529170532_defReal, d, TOL)

    !! Origin from the top
    rb = 2.0_defReal
    a0 = 0.0_defReal
    r = [0.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(2.0_defReal, d, TOL)

    !! Opposite involute
    rb = 1.0_defReal
    a0 = -PI
    r = [0.1_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.005326125763902_defReal, d, TOL)

    !! Opposite involute origin from the top
    !! (Could also be distance 2.0!)
    rb = 1.0_defReal
    a0 = -PI
    r = [0.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(4.971693870713802_defReal, d, TOL)

    !! Opposite involute +ve phase
    rb = 1.0_defReal
    a0 = PI
    r = [0.1_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.005326125763902_defReal, d, TOL)

    !! Hit in the complex gap
    rb = 1.0_defReal
    a0 = PI
    r = [4.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.87_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.035907396620697_defReal, d, TOL)

    !! Hit across the complex gap
    rb = 1.0_defReal
    a0 = PI
    r = [4.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.85_defReal * PI
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(6.242127085088714_defReal, d, TOL)

    !!<><><><><><><><><><><><><><><><><><><><>
    !! Edge Cases found during initial testing

    !! Unstable with Newton iteration only (close to maximum)
    !! Requires improved initial guess
    rb = 1.0720802059149177_defReal
    a0 = -1.4511564629610503_defReal
    r = [-0.3770951755705738_defReal, -1.2108722180281113_defReal, 0.0_defReal]
    theta = 1.8325433442424157_defReal
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(6.261406035145226_defReal, d, TOL)

    !! Hit in the complex gap but with start point on the RHS of maximum
    rb = 1.9454685012010147_defReal
    a0 = -1.0436493497633244_defReal
    r = [1.62573767830043_defReal, -1.330988362095857_defReal, 0.0_defReal]
    theta = 3.049211353314198_defReal
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(0.9033408415670273_defReal, d, TOL)

  end subroutine test_distance_planar

@Test
  subroutine test_distance_3D()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: rb, a0, d, theta
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    !! Catch a vertical particle
    rb = 1.0_defReal
    a0 = 0.0_defReal
    r = [2.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = 0.0_defReal
    u = [0.0_defReal, 0.0_defReal, 1.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(INF, d, TOL)

    !! Inclined hit in the complex gap
    rb = 2.0_defReal
    a0 = 0.0_defReal
    r = [-2.0_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.15_defReal * PI
    u = [cos(theta), sin(theta), 1.0_defReal]
    u = u / norm2(u)
    d = involute_distance(rb, a0, r, u)
    @assertEqual(4.405378529170532_defReal / (ONE - u(3)**2), d, TOL)

    !! Inclined hit outside the gap
    rb = 1.0_defReal
    a0 = PI
    r = [0.1_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 2.0_defReal]
    u = u / norm2(u)
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.005326125763902_defReal/  (ONE - u(3)**2), d, TOL)

  end subroutine test_distance_3D


end module involuteUniverse_test
