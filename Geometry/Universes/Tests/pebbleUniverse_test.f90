module pebbleUniverse_test

use numPrecision
use pfUnit_mod
use dictionary_class,     only : dictionary
use dictParser_func,      only : charToDict
use surfaceShelf_class,   only : surfaceShelf
use cellShelf_class,      only : cellShelf
use pebbleUniverse_class, only : pebbleUniverse
use charMap_class,        only : charMap

use coord_class,          only : coord

implicit none

! Parameters
character(*), parameter :: UNI_DEF = &
"id 7; type pebbleUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
& radius 3.0; inside u<1>; outside water;"

type(surfaceShelf)   :: surfs
type(cellShelf)      :: cells
type(charMap)        :: mats
type(pebbleUniverse) :: uni

contains

  !!
  !! Setup enviroment
  !!
  @Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    character(nameLen) :: name
    type(dictionary) :: dict

    ! Define mapping of materials to matIdx
    ! Since procedure expects a nameLen long char we need to assign it first
    ! to be correctly padded with spaces
    name = "water"
    call mats % add(name, 1)

    ! Feed the definition to the class
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)

    @assertEqual(fill(1), -1)
    @assertEqual(fill(2), 1)

  end subroutine setUp

  !!
  !! Clean enviroment
  !!
  @After
  subroutine clean()


  end subroutine clean

  !!
  !! Example unit test
  !!
  @Test
  subroutine testHalfspace()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: r
    real(defReal), dimension(3) :: u

    ! Enter inside
    r = [0.0, 0.0, 0.0]
    u = [0.0, 0.0, 1.0]
    call uni % enter(pos, r, u)
    @assertEqual(1, pos % localId)
    @assertEqual(0, pos % cellIdx)

    ! Enter outside
    r = [0.0, 4.0, 0.0]
    call uni % enter(pos, r, u)
    @assertEqual(2, pos % localId)
    @assertEqual(0, pos % cellIdx)

  end subroutine testHalfspace


end module pebbleUniverse_test
