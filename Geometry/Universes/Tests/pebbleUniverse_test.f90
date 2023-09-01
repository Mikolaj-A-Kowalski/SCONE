module pebbleUniverse_test

use numPrecision
use pfUnit_mod

implicit none

contains

  !!
  !! Setup enviroment
  !!
  @Before
  subroutine setUp()

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
  subroutine sampleTest()

    @assertEqual(1, 1)

  end subroutine sampleTest


end module pebbleUniverse_test
