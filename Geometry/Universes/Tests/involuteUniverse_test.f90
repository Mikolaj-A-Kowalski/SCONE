module involuteUniverse_test
  use numPrecision
  use universalVariables
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use pFUnit_mod
  use involuteUniverse_class

  implicit none

  character(*), parameter :: DICT_INPUT = "&
  &  id 2;&
  &  type involuteUniverse;&
  &  baseRadius 1.0;&
  &  hubRadius  1.0;&
  &  numPlates 2;&
  &  plateThickness 0.1;&
  &  plateFill fuel;&
  &  hubFill Al;&
  &  channelFill water;"

  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats

contains

@Test
  subroutine distance_in_universe()
    type(involuteUniverse)      :: uni
    type(dictionary)            :: dict
    character(nameLen)          :: name
    integer(shortInt)           :: localId, cellIdx, surfIdx
    integer(shortInt), dimension(:), allocatable :: fill
    type(coord)                 :: coords
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: d
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Load Materials
    name = "water"
    call mats % add(name, 1)
    name = "fuel"
    call mats % add(name, 2)
    name = "Al"
    call mats % add(name, 3)

    ! Load dictionary
    call charToDict(dict, DICT_INPUT)

    ! Initialise universe
    call uni % init(fill, dict, cells, surfs, mats)

    !<><><><><><><><><>
    ! Proper Checks can begin

    ! Point very close to the surface and the hub
    r = [-1.009_defReal, -9.0E-4_defReal, -14.0_defReal]
    u = [-ONE, ZERO, ZERO]

    call uni % findCell(localId, cellIdx, r, u)

    ! Calculate the distance
    coords % r = r
    coords % dir = u
    coords % localId = localId
    call uni % distance(d, surfIdx, coords)

    @assertEqual(0.0006604491039443605_defReal / sqrt(ONE- u(3)**2), d, TOL)

    ! Cleanup
    call mats % kill()

  end subroutine distance_in_universe

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

    ! May not converge is tolerance is too small
    rb = 3.5000000000000000_defReal
    a0 = -2.7931476739450991_defReal
    r = [-3.7817149350289174_defReal, -1.6700800181359774_defReal, 0.0_defReal]
    theta = -4.296525884176383
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(1.0714829417698582e-06_defReal, d, TOL)

    rb = 3.5000000000000000_defReal
    a0 = 1.1995194940784073_defReal
    r = [-4.742309241034304_defReal, 6.686200553696211_defReal, 0.0_defReal]
    theta = 0.6169191679045625
    u = [cos(theta), sin(theta), 0.0_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(7.639077762594239e-08_defReal, d, TOL)

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
    @assertEqual(4.405378529170532_defReal / sqrt(ONE - u(3)**2), d, TOL)

    !! Inclined hit outside the gap
    rb = 1.0_defReal
    a0 = PI
    r = [0.1_defReal, 2.0_defReal, 0.0_defReal]
    theta = -0.5_defReal * PI
    u = [cos(theta), sin(theta), 2.0_defReal]
    u = u / norm2(u)
    d = involute_distance(rb, a0, r, u)
    @assertEqual(5.005326125763902_defReal / sqrt(ONE - u(3)**2), d, TOL)

    !! Regression Case
    !! Caused -ve distance if the guess was always pushed to the right of complex gap
    rb = 6.5_defReal
    a0 = -1.2788784253551366_defReal
    r = [2.8708518996940788_defReal, -7.3103261399936441_defReal, 10.468929715684157_defReal]
    u = [-0.61021412406725728_defReal, 0.77555085185490247_defReal, 0.16173929323440561_defReal]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(1.7637231237285702_defReal / sqrt(ONE - u(3)**2) , d, TOL)

    !! Regression case
    !! Caused -ve distance, point is extremely close to the surface
    rb = 0.999_defReal
    a0 = 3.1415926535897931_defReal
    r = [-1.0084672937865440_defReal, -8.7349079288844174E-004_defReal  ,-14.253067666806164_defReal]
    u = [[0.96015501731015551_defReal,      9.0670022704170905E-002_defReal ,-0.26435069456492838_defReal]]
    d = involute_distance(rb, a0, r, u)
    @assertEqual(0.009249182790172505_defReal / sqrt(ONE - u(3)**2) , d, TOL)

  end subroutine test_distance_3D


end module involuteUniverse_test
