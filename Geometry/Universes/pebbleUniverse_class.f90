module pebbleUniverse_class

    use numPrecision
    use universalVariables, only : INF, targetNotFound, SURF_TOL
    use genericProcedures,  only : fatalError, numToChar, swap, dotProduct
    use dictionary_class,   only : dictionary
    use coord_class,        only : coord
    use charMap_class,      only : charMap
    use surfaceShelf_class, only : surfaceShelf
    use cellShelf_class,    only : cellShelf
    use universe_inter,     only : universe, kill_super => kill, charToFill

    implicit none
    private

    !!
    !! Universe that contains a single sphere centred at the origin
    !!
    type, public, extends(universe) :: pebbleUniverse
      private
      real(defReal) :: radius
    contains
      ! Superclass procedures
      procedure :: init
      procedure :: kill
      procedure :: findCell
      procedure :: distance
      procedure :: cross
      procedure :: cellOffset
    end type pebbleUniverse

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
      class(pebbleUniverse), intent(inout)                      :: self
      integer(shortInt), dimension(:), allocatable, intent(out) :: fill
      class(dictionary), intent(in)                             :: dict
      type(cellShelf), intent(inout)                            :: cells
      type(surfaceShelf), intent(inout)                         :: surfs
      type(charMap), intent(in)                                 :: mats
      integer(shortInt)                                         :: id
      real(defReal), dimension(:), allocatable                  :: temp
      character(nameLen)                                        :: temp_name
      character(100), parameter :: Here = 'init (pebbleUniverse_class.f90)'

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

      ! Load radius
      call dict % get(self % radius, 'radius')

      ! Check radius
      if (self % radius <= ZERO) then
        call fatalError('Radius must be positive', Here)
      end if

      ! Assign the fill array
      allocate(fill(2))
      call dict % get(temp_name, 'inside')
      fill(1) = charToFill(temp_name, mats, HERE)

      call dict % get(temp_name, 'outside')
      fill(2) = charToFill(temp_name, mats, HERE)


    end subroutine init

    !!
    !! Find local cell ID given a point
    !!
    !! See universe_inter for details.
    !!
    subroutine findCell(self, localID, cellIdx, r, u)
      class(pebbleUniverse), intent(inout)       :: self
      integer(shortInt), intent(out)          :: localID
      integer(shortInt), intent(out)          :: cellIdx
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal)                           :: r_sq

      ! By convention the cellIdx in not cell universe is 0
      cellIdx = 0

      ! Check if the point is inside the sphere
      r_sq = sum(r**2)
      if (r_sq <= self % radius**2) then
        localID = 1
      else
        localID = 2
      end if

    end subroutine findCell

    !!
    !! Return distance to the next boundary between local cells in the universe
    !!
    !! See universe_inter for details.
    !!
    !! Errors:
    !!   fatalError is localID is invalid
    !!
    subroutine distance(self, d, surfIdx, coords)
      class(pebbleUniverse), intent(inout)  :: self
      real(defReal), intent(out)         :: d
      integer(shortInt), intent(out)     :: surfIdx
      type(coord), intent(in)            :: coords
      real(defReal)                      :: c, k, delta
      character(100), parameter :: Here = 'distance (pebbleUniverse_class.f90)'

      ! Calculate distance to the next boundary
      ! Solve d^2 + 2kd + c = 0
      c = sum(coords % r * coords % r) - self % radius**2
      k = dotProduct(coords % r, coords % dir)

      ! delta/4
      delta = k**2 - c

      if (delta < ZERO) then
        d = INF

      else if (abs(c) < self % radius * SURF_TOL) then

        if (k < ZERO) then
          d = -k + sqrt(delta)
        else
          d = INF
        end if

      else if (c < ZERO) then
        d = -k + sqrt(delta)

      else ! Point outside the surface
        d = -k - sqrt(delta)
        if (d <= ZERO) d = INF

      end if


    end subroutine distance

    !!
    !! Cross between local cells
    !!
    !! See universe_inter for details.
    !!
    !! Errors:
    !!   fatalError if localID is neither 1 nor 2
    !!
    subroutine cross(self, coords, surfIdx)
      class(pebbleUniverse), intent(inout) :: self
      type(coord), intent(inout)         :: coords
      integer(shortInt), intent(in)      :: surfIdx
      character(100), parameter :: Here = 'cross (pebbleUniverse_class.f90)'

      if (coords % localID == 1) then
        coords % localID = 2

      else if (coords % localID == 2) then
        coords % localID = 1

      else
        call fatalError(Here, 'Invalid localID: ' // numToChar(coords % localID))

      end if

    end subroutine cross

    !!
    !! Return offset for the current cell
    !!
    !! See universe_inter for details.
    !!
    function cellOffset(self, coords) result (offset)
      class(pebbleUniverse), intent(in) :: self
      type(coord), intent(in)         :: coords
      real(defReal), dimension(3)     :: offset

      ! There is no cell offset
      offset = ZERO

    end function cellOffset

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(pebbleUniverse), intent(inout) :: self


    end subroutine kill

  end module pebbleUniverse_class
