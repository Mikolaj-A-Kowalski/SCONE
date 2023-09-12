module involuteUniverse_class

    use numPrecision
    use universalVariables, only : OUTSIDE_MAT
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

    implicit none
    private

    !!
    !! A spiral of involute channels around a central cylindrical hub
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
    !!    channelFill water;
    !!    hubFill reflector;
    !!  };
    !!
    type, public, extends(universe) :: involuteUniverse
      real(defReal)     :: hubRadius
      real(defReal)     :: baseRadius
      integer(shortInt) :: numPlates
      real(defReal)     :: plateThickness
      real(defReal)     :: pitch
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
      integer(shortInt)                                         :: id
      real(defReal), dimension(:), allocatable                  :: temp
      real(defReal)                                             :: circ
      character(nameLen)                                        :: temp_name
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

      ! Load plate info
      call dict % get(self % numPlates, 'numPlates')
      if (self % numPlates <= 0) then
        call fatalError(Here, 'Number of plates must be +ve. Is: ' // numToChar(self % numPlates))
      end if

      call dict % get(self % plateThickness, 'plateThickness')
      if (self % plateThickness <= 0) then
        call fatalError(Here, 'Plate thickness must be +ve. Is: ' // numToChar(self % plateThickness))
      end if

      ! Verify that there will be no overlap between plates
      circ = TWO * PI * self % baseRadius
      self % pitch = circ / self % numPlates
      if (self % plateThickness > self % pitch) then
        call fatalError(Here, 'Plate thickness is too large. Plates will overlap.')
      end if

      ! Load fill info
      allocate( fill(1 + 2 * self % numPlates))

      ! Load hub material
      call dict % get(temp_name, 'hubFill')
      fill(1) = charToFill(temp_name, mats, Here)

      ! Load plate material
      call dict % get(temp_name, 'plateFill')
      fill(2:1 + 2 * self % numPlates:2) = charToFill(temp_name, mats, Here)

      ! Load channel material
      call dict % get(temp_name, 'channelFill')
      fill(3:1 + 2 * self % numPlates:2) = charToFill(temp_name, mats, Here)

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
      integer(shortInt)                       :: bin
      character(100), parameter               :: Here = 'findCell (involuteUniverse_class.f90)'

      ! There is no user defined cell
      cellIdx = 0

      ! Place the particle
      radius = sqrt(r(1)**2 + r(2)**2)
      theta  = atan2(r(2), r(1))

      ! put the angle in [0, 2pi) range
      if (theta < ZERO) then
        theta = theta + TWO * PI
      end if

      if (radius < self % hubRadius) then
        localID = 1
        return
      end if

      ratio =  self % baseRadius / radius
      theta_0 = sqrt(1 - ratio**2) / ratio - acos(ratio)

      ! Renormalise theta_0 to be in [0, 2pi) range
      theta_0 = theta_0 - int((theta_0) / TWO_PI) * TWO_PI

      if (theta_0 < 0 .or. theta_0 >= TWO_PI) then
        call fatalError(Here, 'Invalid theta_0: ' // numToChar(theta_0))
      end if
      angle = theta - theta_0
      if (angle < 0) then
        angle = angle + TWO_PI
      end if
      bin = int(angle / TWO_PI * self % numPlates)

      ! Calculate normalised angle coordinate in the bin
      ! Takes values in [0; 1]
      angle = angle - bin * TWO_PI / self % numPlates
      angle = angle / (TWO_PI / self % numPlates)

      if (angle < self % plateThickness / self % pitch) then
        localID = 2 * bin + 2
      else
        localID = 2 * bin + 3
      end if

      if (localID < 1 .or. localID > 1 + 2 * self % numPlates) then
        print *, "r", r
        print *, "theta: ", theta
        print *, "theta_0: ", theta_0
        print *, "diff: ", theta - theta_0
        call fatalError(Here, 'Invalid local cell ID: ' // numToChar(localID))
      end if

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
      character(100), parameter :: Here = 'distance (involuteUniverse_class.f90)'

      call fatalError(Here, 'Not implemented')

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

      call fatalError(Here, 'Not implemented')

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

  end module involuteUniverse_class
