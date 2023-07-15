module pebbleUniverse_class

    use numPrecision
    use universalVariables, only : INF, targetNotFound
    use genericProcedures,  only : fatalError, numToChar, swap
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
      class(pebbleUniverse), intent(inout)                        :: self
      integer(shortInt), dimension(:), allocatable, intent(out) :: fill
      class(dictionary), intent(in)                             :: dict
      type(cellShelf), intent(inout)                            :: cells
      type(surfaceShelf), intent(inout)                         :: surfs
      type(charMap), intent(in)                                 :: mats
      character(100), parameter :: Here = 'init (pebbleUniverse_class.f90)'




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
      character(100), parameter :: Here = 'distance (pebbleUniverse_class.f90)'



    end subroutine distance

    !!
    !! Cross between local cells
    !!
    !! See universe_inter for details.
    !!
    subroutine cross(self, coords, surfIdx)
      class(pebbleUniverse), intent(inout) :: self
      type(coord), intent(inout)         :: coords
      integer(shortInt), intent(in)      :: surfIdx
      character(100), parameter :: Here = 'cross (pebbleUniverse_class.f90)'



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
