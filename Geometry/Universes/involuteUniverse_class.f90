module involuteUniverse_class

    use numPrecision
    use universalVariables, only : OUTSIDE_MAT, INF, NUDGE
    use genericProcedures,  only : fatalError, numToChar
    use dictionary_class,   only : dictionary
    use coord_class,        only : coord
    use charMap_class,      only : charMap
    use surface_inter,      only : surface
    use surfaceShelf_class, only : surfaceShelf
    use cylinder_class,     only : cylinder
    use cell_inter,         only : cell
    use cellShelf_class,    only : cellShelf
    use universe_inter,     only : universe, kill_super => kill, charToFill
    use cylinder_class,     only : cylinder

    implicit none
    private

    !!
    !! Public functions for distance calculation
    !! Are made public mostly for testing
    !!
    public :: phase_and_derivative
    public :: involute_newton
    public :: involute_distance

    integer(shortInt), parameter :: SURF_CLOCKWISE = -1
    integer(shortInt), parameter :: SURF_ANTI_CLOCKWISE = -2
    integer(shortInt), parameter :: SURF_CYLINDER = -3

    !!
    !! A spiral of involute channels around a central cylindrical hub
    !!
    !! TODO:
    !!   Finish the documentation. Just wait until the input is finalised
    !!
    !!
    !! Sample input dictionary:
    !!
    !!  involute {
    !!    id 1;
    !!    type involuteUniverse;
    !!    #origin (0.0 0.0 0.0); #
    !!    #rotation (0.0 0.0 0.0); #
    !!    baseRadius 1.0;
    !!    hubRadius  1.0;
    !!    numPlates  10;
    !!    plateThickness 0.1;
    !!    plateFill fuel;
    !!    # cladThickness 0.05; #
    !!    # cladFill Al; #
    !!    channelFill water;
    !!    hubFill reflector;
    !!  };
    !!
    type, public, extends(universe) :: involuteUniverse
      real(defReal)     :: hubRadius
      real(defReal)     :: baseRadius
      integer(shortInt) :: numPlates
      real(defReal)     :: anglePitch
      real(defReal), dimension(:), allocatable :: angleThickness
      type(cylinder)    :: hubCylinder
    contains
      ! Superclass procedures
      procedure :: init
      procedure :: kill
      procedure :: findCell
      procedure :: distance
      procedure :: cross
      procedure :: cellOffset

    end type involuteUniverse


  contains

    !!
    !! Initialise Universe
    !!
    !! See universe_inter for details.
    !!
    !! Errors:
    !!   fatalError for invalid input
    !!
    subroutine init(self, fill, dict, cells, surfs, mats)
      class(involuteUniverse), intent(inout)                    :: self
      integer(shortInt), dimension(:), allocatable, intent(out) :: fill
      class(dictionary), intent(in)                             :: dict
      type(cellShelf), intent(inout)                            :: cells
      type(surfaceShelf), intent(inout)                         :: surfs
      type(charMap), intent(in)                                 :: mats
      integer(shortInt)                                         :: id, numRegions
      real(defReal), dimension(:), allocatable                  :: temp
      real(defReal)                                             :: circ, plateThickness, cladThickness
      character(nameLen)                                        :: temp_name
      logical(defBool)                                          :: hasCladding
      character(100), parameter :: Here = 'init (involuteUniverse_class.f90)'

      ! Load basic data
      call dict % get(id, 'id')
      if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: '//numToChar(id))
      call self % setId(id)

      ! Load origin
      if (dict % isPresent('origin')) then
        call dict % get(temp, 'origin')

        if (size(temp) /= 3) then
          call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(size(temp)))
        end if
        call self % setTransform(origin=temp)

      end if

      ! Load rotation
      if (dict % isPresent('rotation')) then
        call dict % get(temp, 'rotation')

        if (size(temp) /= 3) then
          call fatalError(Here, '3 rotation angles must be given. Has only: '//numToChar(size(temp)))
        end if
        call self % setTransform(rotation=temp)
      end if

      ! Data specific to this universe starts here
      ! Load radii
      call dict % get(self % baseRadius, 'baseRadius')
      if (self % baseRadius <= 0) then
        call fatalError(Here, 'Base radius must be +ve. Is: ' // numToChar(self % baseRadius))
      end if

      call dict % get(self % hubRadius, 'hubRadius')
      if (self % hubRadius < self % baseRadius) then
        call fatalError(Here, 'Hub radius must be >= base radius')
      end if

      ! Construct the central hub cylinder surface
      call self % hubCylinder % build(id=1, origin=[ZERO, ZERO, ZERO], type="zCylinder", radius=self % hubRadius)

      ! Load plate info
      call dict % get(self % numPlates, 'numPlates')
      if (self % numPlates <= 0) then
        call fatalError(Here, 'Number of plates must be +ve. Is: ' // numToChar(self % numPlates))
      end if

      call dict % get(plateThickness, 'plateThickness')
      if (plateThickness <= 0) then
        call fatalError(Here, 'Plate thickness must be +ve. Is: ' // numToChar(plateThickness))
      end if
      numRegions = 2

      ! Load cladding info of present
      cladThickness = ZERO
      hasCladding = dict % isPresent('cladThickness')
      if (hasCladding) then
        call dict % get(cladThickness, 'cladThickness')
        if (cladThickness <= 0) then
          call fatalError(Here, 'Cladding thickness must be +ve. Is: ' // numToChar(cladThickness))
        end if
        numRegions = 4
      end if

      ! Verify that there will be no overlap between plates
      circ = TWO * PI * self % baseRadius
      if ((plateThickness + TWO * cladThickness) > circ / self % numPlates) then
        call fatalError(Here, 'Plate thickness is too large. Plates will overlap.')
      end if

      !<><><><> Set the angular divisions
      self % anglePitch = TWO_PI / self % numPlates
      allocate(self % angleThickness(numRegions))

      if (hasCladding) then
        self % angleThickness(1) = cladThickness / circ * self % numPlates
        self % angleThickness(2) = plateThickness / circ * self % numPlates + self % angleThickness(1)
        self % angleThickness(3) = cladThickness / circ * self % numPlates + self % angleThickness(2)
        self % angleThickness(4) = ONE
      else
        self % angleThickness(1) = plateThickness / circ * self % numPlates
        self % angleThickness(2) = ONE
      end if

      ! Load fill info
      allocate( fill(1 + numRegions * self % numPlates))

      ! Load hub material
      call dict % get(temp_name, 'hubFill')
      fill(1) = charToFill(temp_name, mats, Here)

      if (hasCladding) then
        call dict % get(temp_name, 'cladFill')
        fill(2:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)
        fill(4:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)

        call dict % get(temp_name, 'plateFill')
        fill(3:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)

        call dict % get(temp_name, 'channelFill')
        fill(5:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)
      else
        call dict % get(temp_name, 'plateFill')
        fill(2:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)

        call dict % get(temp_name, 'channelFill')
        fill(3:1 + numRegions * self % numPlates: numRegions) = charToFill(temp_name, mats, Here)

      end if

    end subroutine init

    !!
    !! Find local cell ID given a point
    !!
    !! See universe_inter for details.
    !!
    subroutine findCell(self, localID, cellIdx, r, u)
      class(involuteUniverse), intent(inout)  :: self
      integer(shortInt), intent(out)          :: localID
      integer(shortInt), intent(out)          :: cellIdx
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal)                           :: radius, theta, ratio, theta_0, angle
      integer(shortInt)                       :: bin, plateIdx, i
      character(100), parameter               :: Here = 'findCell (involuteUniverse_class.f90)'

      ! There is no user defined cell
      cellIdx = 0

      ! If the particle is in -ve hub cylinder halfspace return
      if (.not. self % hubCylinder % halfspace(r, u)) then
        localID = 1
        return
      end if

      ! Place the particle
      radius = sqrt(r(1)**2 + r(2)**2)
      theta  = atan2(r(2), r(1))

      ratio =  self % baseRadius / radius
      theta_0 = sqrt(1 - ratio**2) / ratio - arcCos(ratio)

      ! Renormalise theta_0 to be in [-pi, pi) range
      theta_0 = theta_0 - int((theta_0 + PI) / TWO_PI) * TWO_PI

      ! We need to make sure that the distance is in [0, 2pi] range
      angle = theta - theta_0
      if (angle < 0) then
        angle = angle + TWO_PI
      end if
      bin = int(angle / self % anglePitch)

      ! Calculate normalised angle coordinate in the bin
      ! Takes values in [0; 1]
      angle = angle - bin * self % anglePitch
      angle = angle / self % anglePitch

      ! Find the bin index
      ! Do the counting search
      plateIdx = 1
      do i = 1, size(self % angleThickness)
        if (angle > self % angleThickness(i)) then
          plateIdx = plateIdx + 1
        end if
      end do

      localID = bin * size(self % angleThickness) + plateIdx + 1

    end subroutine findCell

    !!
    !! Return distance to the next boundary between local cells in the universe
    !!
    !! See universe_inter for details.
    !!
    subroutine distance(self, d, surfIdx, coords)
      class(involuteUniverse), intent(inout) :: self
      real(defReal), intent(out)             :: d
      integer(shortInt), intent(out)         :: surfIdx
      type(coord), intent(in)                :: coords
      real(defReal)                          :: d_cylinder, d_clockwise, d_aclockwise
      real(defReal)                          :: phase_clockwise, phase_aclockwise
      integer(shortInt)                      :: bin, plateIdx
      character(100), parameter :: Here = 'distance (involuteUniverse_class.f90)'

      d_cylinder = self % hubCylinder % distance(coords % r, coords % dir)

      ! Special case inside the hub region
      if (coords % localID == 1) then
        d = d_cylinder
        surfIdx = SURF_CYLINDER
        return
      end if

      ! Outer cell
      ! 1) Calculate the phases of the boundary surfaces (in [-pi;pi] range)
      ! 2) Calculate the distance to both and the cylinder
      ! 3) Select the smallest and set the right surdIdx

      ! Calculate the bin and plateIdx from localID
      bin = (coords % localID - 2) / size(self % angleThickness)
      plateIdx = (coords % localID - 1) - bin * size(self % angleThickness)

      ! Determine the phases of boundary involutes
      ! Phases are calculated in [0;2pi] range so we need to correctly map them to
      ! [-pi;pi] range
      phase_clockwise = bin * self % anglePitch
      if (plateIdx /= 1) then
        phase_clockwise = phase_clockwise + self % angleThickness(plateIdx - 1) * self % anglePitch
      end if
      phase_aclockwise = bin * self % anglePitch + self % angleThickness(plateIdx) * self % anglePitch

      if (phase_clockwise > PI) then
        phase_clockwise = phase_clockwise - TWO_PI
      end if
      if (phase_aclockwise > PI) then
        phase_aclockwise = phase_aclockwise - TWO_PI
      end if

      d_clockwise = involute_distance(self % baseRadius, phase_clockwise, coords % r, coords % dir)
      d_aclockwise = involute_distance(self % baseRadius, phase_aclockwise, coords % r, coords % dir)

      if (d_clockwise < ZERO) then
        print *, coords % r, coords % dir, phase_clockwise, self % baseRadius
        call fatalError(Here, 'Clockwise distance is negative: '//numToChar(d_clockwise))
      end if
      if (d_aclockwise < ZERO) then
        print *, coords % r, coords % dir, phase_aclockwise, self % baseRadius
        call fatalError(Here, 'Anti-clockwise distance is negative: '//numToChar(d_aclockwise))
      end if

      if (d_clockwise < d_aclockwise) then
        d = d_clockwise
        surfIdx = SURF_CLOCKWISE
      else
        d = d_aclockwise
        surfIdx = SURF_ANTI_CLOCKWISE
      end if

      ! Compare with the distance to the cylinder
      if (d_cylinder < d) then
        d = d_cylinder
        surfIdx = SURF_CYLINDER
      end if

    end subroutine distance

    !!
    !! Cross between local cells
    !!
    !! See universe_inter for details.
    !!
    !! Errors:
    !!   fatalError if surface from distance is not MOVING_IN or MOVING_OUT
    !!
    subroutine cross(self, coords, surfIdx)
      class(involuteUniverse), intent(inout) :: self
      type(coord), intent(inout)             :: coords
      integer(shortInt), intent(in)          :: surfIdx
      character(100), parameter :: Here = 'cross (involuteUniverse_class.f90)'

      select case (surfIdx)
        case(SURF_CLOCKWISE, SURF_ANTI_CLOCKWISE)
          ! Nudge a particle a little bit
          ! We should make it along the normal of the Involute
          coords % r = coords % r + coords % dir * NUDGE

          call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

        case(SURF_CYLINDER)
          ! Push well inside or outside the cylinder
          coords % r = coords % r + coords % dir * NUDGE
          call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

        case default
          call fatalError(Here, 'Invalid surface index: '//numToChar(surfIdx))

      end select

    end subroutine cross

    !!
    !! Return offset for the current cell
    !!
    !! See universe_inter for details.
    !!
    function cellOffset(self, coords) result (offset)
      class(involuteUniverse), intent(in) :: self
      type(coord), intent(in)         :: coords
      real(defReal), dimension(3)     :: offset

      offset = ZERO

    end function cellOffset

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(involuteUniverse), intent(inout) :: self

      call kill_super(self)

    end subroutine kill

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Public Involute functions for involute distance calculations
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !!
    !! ====================================================
    !! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
    !!
    !! Developed at SunSoft, a Sun Microsystems, Inc. business.
    !! Permission to use, copy, modify, and distribute this
    !! software is freely granted, provided that this notice
    !! is preserved.
    !! ====================================================
    !!
    !! Rational term used to evaluate arcCos and arcSin
    !!
    !! Args:
    !!   x [in] -> Value inb [-0.25; 0.25] range
    !!
    !! Result:
    !!   Rational approximation of [arcSin(sqrt(x)) - sqrt(x)] / sqrt(x)^3
    !!
    function rational(x) result(y)
      real(defReal), intent(in)  :: x
      real(defReal)              :: y
      real(defReal)              :: p, q
      real(defReal), dimension(*), parameter :: P_POLY = [ 1.66666666666666657415e-01_defReal,&
                                                          -3.25565818622400915405e-01_defReal,&
                                                          2.01212532134862925881e-01_defReal,&
                                                          -4.00555345006794114027e-02_defReal,&
                                                          7.91534994289814532176e-04_defReal,&
                                                          3.47933107596021167570e-05_defReal]
      real(defReal), dimension(*), parameter :: Q_POLY = [ 1.0_defReal, &
                                                          -2.40339491173441421878e+00_defReal,&
                                                          2.02094576023350569471e+00_defReal,&
                                                          -6.88283971605453293030e-01_defReal,&
                                                          7.70381505559019352791e-02_defReal]

      p = P_POLY(1) + x * (P_POLY(2) + x * (P_POLY(3) + x * (P_POLY(4) + x * (P_POLY(5) + x * P_POLY(6)))))
      q = Q_POLY(1) + x * (Q_POLY(2) + x * (Q_POLY(3) + x * (Q_POLY(4) + x * Q_POLY(5))))
      y = p / q
    end function rational

    !!
    !! ====================================================
    !! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
    !!
    !! Developed at SunSoft, a Sun Microsystems, Inc. business.
    !! Permission to use, copy, modify, and distribute this
    !! software is freely granted, provided that this notice
    !! is preserved.
    !! ====================================================
    !!
    !! The implementation of rational approximation of arcCos
    !!
    !! It is accurate to few ULP, but runs much faster than gcc libm version
    !! on Debian 11. Since arcCos is used multiple times and its accuracy
    !! is not paramount this implementation is used.
    !!
    !! Based on the RUST libm code:
    !! https://docs.rs/libm/latest/src/libm/math/acos.rs.html#63-112
    !!
    !! Args:
    !!   x [in] -> Value in [-1; 1] range
    !!
    !! Result:
    !!   ArcCos of x in [0; pi] range
    !!
    function arcCos(x) result(y)
      real(defReal), intent(in)  :: x
      real(defReal)              :: y
      real(defReal)              :: z, s

      if (abs(x) <= 0.5) then
        y = HALF*PI - (x + x*x*x * rational(x*x))
      else if (x > 0.5) then
        z = HALF * (ONE - x)
        s = sqrt(z)
        y = TWO * (s + s * z * rational(z))
      else
        z = HALF * (ONE + x)
        s = sqrt(z)
        y = PI - TWO * (s + s * z * rational(z))
      end if

    end function arcCos

    !!
    !! Calculate the phase function and its derivative
    !!
    !! Args:
    !!   f [out]      -> Phase function
    !!   df [out]     -> Derivative of the phase function
    !!   rb [in]      -> Base radius
    !!   a0 [in]      -> Involute phase
    !!   rl [in]      -> Radius of the line characteristic point
    !!   theta_l [in] -> Angle of the line characteristic point
    !!   d [in]       -> Distance along the line from the characteristic point
    !!
    subroutine phase_and_derivative(f, inv_df, rb, a0, rl, theta_l, d)
      real(defReal), intent(out) :: f
      real(defReal), intent(out) :: inv_df
      real(defReal), intent(in)  :: rb
      real(defReal), intent(in)  :: a0
      real(defReal), intent(in)  :: rl
      real(defReal), intent(in)  :: theta_l
      real(defReal), intent(in)  :: d
      real(defReal)              :: r, r_sq, rb_sq, atan_term, t

      ! Calculate the magnitude of the radius
      r_sq = d**2 + rl**2
      r = sqrt(r_sq)
      rb_sq = rb**2

      ! We cannot evaluate atan term with arcCos because accuracy of arcCos
      ! when approximating atan(d/rl) for d ~ 0.0 is quite poor
      atan_term = atan(d/rl)

      ! Calculate the phase and derivative
      if (r > rb) then
        t = sqrt(r_sq/rb_sq - ONE)
        f = theta_l - a0 - t - atan_term + arcCos(rb / r)
        inv_df =  -r_sq / (rl + d * t)

      else
        f = theta_l - a0 - atan_term
        inv_df = -r_sq / rl

      end if
    end subroutine phase_and_derivative


    !!
    !! Solve the involute phase equation using Newton's method
    !!
    !! Args:
    !!  rb [in]      -> Base radius
    !!  a0 [in]      -> Involute phase
    !!  rl [in]      -> Radius of the line characteristic point
    !!  theta_l [in] -> Angle of the line characteristic point
    !!  d0 [in]      -> Initial guess for the distance
    !!  rhs [in]     -> Right hand side of the phase equation
    !!
    !! Result:
    !!   Position d where the phase is equal to the rhs
    !!
    function involute_newton(rb, a0, rl, theta_l, d0, rhs) result(d)
      real(defReal), intent(in) :: rb
      real(defReal), intent(in) :: a0
      real(defReal), intent(in) :: rl
      real(defReal), intent(in) :: theta_l
      real(defReal), intent(in) :: d0
      real(defReal), intent(in) :: rhs
      real(defReal)             :: d, d_last, f, inv_df
      integer(shortInt)         :: i
      real(defReal)             :: tol = 1.0e-8_defReal
      character(100), parameter :: Here = 'involute_newton (involuteUniverse_class.f90)'

      d_last = d0
      d = d0
      do i = 1, 30
        call phase_and_derivative(f, inv_df, rb, a0, rl, theta_l, d)
        d_last = d

        d = d_last - (f - rhs) * inv_df

        if (abs(d - d_last) < tol) then
          return
        end if
      end do

      print *, "rb=", rb, " a0=", a0, " rl=", rl, " theta_l=", theta_l, " d0=", d0, " rhs=", rhs
      call fatalError(Here, 'Newton iteration did not converge')
    end function involute_newton

    !!
    !! Project 3D point into XY plane and transform to characteristic line coordinates
    !!
    !! Args:
    !!   rl [out]      -> Radius of the line characteristic point
    !!   theta_l [out] -> Angle of the line characteristic point
    !!   dir [out]     -> Direction along the line [-1 or 1]
    !!   d0 [out]      -> Distance along the line from the characteristic point
    !!   r [in]        -> 3D point position
    !!   u [in]        -> 3D point direction (assumed normalised)
    !!
    !! Errors:
    !!   fatalError if the direction is perpendicular to XY axis (Z-axis)
    !!
    subroutine baseLine(rl, theta_l, dir, d0, r, u)
      real(defReal), intent(out) :: rl
      real(defReal), intent(out) :: theta_l
      real(defReal), intent(out) :: dir
      real(defReal), intent(out) :: d0
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal), dimension(2) :: normal, along
      character(100), parameter :: Here = 'baseLine (involuteUniverse_class.f90)'

      if (abs(u(1)) < 1.0e-9_defReal .and. abs(u(2)) < 1.0e-9_defReal) then
        call fatalError(Here, 'Direction is perpendicular to XY axis')
      end if

      ! Calculate 2D vector normal to the line
      normal = [-u(2), u(1)]
      normal = normal / norm2(normal)

      ! Calculate the radius and make sure it is +ve
      rl = dot_product(r(1:2), normal)
      if (rl < 0) then
        rl = -rl
        normal = -normal
      end if

      ! Calculate the angle of the line, direction and initial coordinate
      theta_l = atan2(normal(2), normal(1))
      along = [normal(2), -normal(1)]
      d0 = dot_product(r(1:2), along)
      dir = sign(1.0_defReal, dot_product(u(1:2), along))

    end subroutine baseLine

    !!
    !! Calculate the distance to the involute surface
    !!
    !! Args:
    !!   rb [in]      -> Base radius
    !!   a0 [in]      -> Involute phase
    !!   r [in]       -> 3D point position
    !!   u [in]       -> 3D point direction (assumed normalised)
    !!
    !! Result:
    !!   Distance to the involute surface
    !!
    function involute_distance(rb, a0, r, u) result(d)
      real(defReal), intent(in)               :: rb
      real(defReal), intent(in)               :: a0
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal)                           :: d
      real(defReal)                           :: rl, theta_l, dir, start, guess
      real(defReal)                           :: p_max, p_start, temp, rhs
      real(defReal)                           :: val, trig, cos_z

      ! Get the Z-axis component of the direction
      cos_z = u(3)

      ! Nearly vertical particle
      if (abs(cos_z - ONE) < 1.0e-9_defReal) then
        d = INF
        return
      end if

      ! Get the characteristic line parameters
      call baseLine(rl, theta_l, dir, start, r, u)
      guess = start

      ! Get the phase value at the maximum and start
      call phase_and_derivative(p_max, temp, rb, a0, rl, theta_l, -rb)
      call phase_and_derivative(p_start, temp, rb, a0, rl, theta_l, start)

      ! If we find ourselves between the maximum and first root we update
      ! the guess using the quadratic approximation
      if (floor(p_start / TWO_PI) == floor(p_max / TWO_PI)) then
        p_max = p_max - TWO_PI * floor(p_max / TWO_PI)
        guess = -rb + dir * sqrt(TWO * p_max * rb * rl)
      end if

      ! Choose the correct RHS of the equation
      if (dir == ONE .neqv. guess >= -rb) then
        rhs = TWO_PI * ceiling(p_start / TWO_PI)
      else
        rhs = TWO_PI * floor(p_start / TWO_PI)
      end if

      ! If the solution is in RHS of the maximum, check if the root is in the
      ! complex gap region. If it is provide analytical solution and return
      if (rl < rb .and. guess + rb >= ZERO) then
        val = theta_l - a0 - rhs
        trig = arcCos(rl / rb)

        if (trig >= val .and. trig >= -val) then
          d = tan(theta_l - a0) * rl
          d = d - start
          if (dir == -ONE) d = -d
          d = d / sqrt((ONE - cos_z**2))
          return

        end if
      end if

      ! If the Guess is in the complex gap region, we need to move it outside
      if (guess >= -rb .and. guess**2 < rb**2 - rl**2) then
        val = theta_l - a0 - rhs + arcCos(rl / rb)
        guess = sqrt(rb**2 - rl**2)
        if (val < ZERO) guess = -guess
      end if

      ! Finally solve for the intersection
      d = involute_newton(rb, a0, rl, theta_l, guess, rhs)
      d = d - start
      if (dir == -ONE) d = -d
      d = d / sqrt((ONE - cos_z**2))

    end function involute_distance


  end module involuteUniverse_class
