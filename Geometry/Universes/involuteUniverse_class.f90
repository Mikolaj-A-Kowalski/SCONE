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

      ! Place the particle
      radius = sqrt(r(1)**2 + r(2)**2)
      theta  = atan2(r(2), r(1))

      if (radius < self % hubRadius) then
        localID = 1
        return
      end if

      ratio =  self % baseRadius / radius
      theta_0 = sqrt(1 - ratio**2) / ratio - acos(ratio)

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
