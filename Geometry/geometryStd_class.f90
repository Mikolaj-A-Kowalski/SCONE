module geometryStd_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, quickSort
  use coord_class,        only : coordList, coord
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use intMap_class,       only : intMap
  use geometry_inter,     only : geometry, distCache
  use csg_class,          only : csg
  use universe_inter,     only : universe
  use surface_inter,      only : surface

  ! Nuclear Data
  use materialMenu_mod,   only : nMat


  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: geometryStd_CptrCast

  !!
  !! Standard Geometry Model
  !!
  !! Typical geometry of a MC Neutron Transport code composed of multiple nested
  !! universes.
  !!
  !! Boundary conditions in diffrent movement models are handeled:
  !!   move       -> explicitBC
  !!   moveGlobal -> explicitBC
  !!   teleport   -> Co-ordinate transfrom
  !!
  !! Sample Dictionary Input:
  !!   geometry {
  !!     type geometryStd;
  !!     <csg_class difinition>
  !!    }
  !!
  !! Public Members:
  !!   geom -> Representation of geometry by csg_class. Contains all surfaces, cells and universe
  !!     as well as geometry graph and info about root uni and boundary surface.
  !!
  !! Interface:
  !!   Geometry Interface
  !!
  type, public, extends(geometry) :: geometryStd
    type(csg) :: geom

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: bounds
    procedure :: move_noCache
    procedure :: move_withCache
    procedure :: moveGlobal
    procedure :: move_partial
    procedure :: teleport
    procedure :: activeMats

    ! Private procedures
    procedure, private :: activeMatsBelow
    procedure, private :: diveToMat
    procedure, private :: closestDist
    procedure, private :: closestDist_cache
  end type geometryStd

contains

  !!
  !! Initialise geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine init(self, dict, mats, silent)
    class(geometryStd), intent(inout)      :: self
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: mats
    logical(defBool), optional, intent(in) :: silent

    ! Build the representation
    call self % geom % init(dict, mats, silent)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(geometryStd), intent(inout) :: self

    call self % geom % kill()

  end subroutine kill

  !!
  !! Place coordinate list into geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine placeCoord(self, coords)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: r, dir
    character(100), parameter :: Here = 'placeCoord (geometryStd_class.f90)'

    ! Check that coordList is initialised
    if (coords % nesting < 1) then
      call fatalError(Here, 'CoordList is not initialised. Nesting is: '//&
                             numToChar(coords % nesting))
    end if

    ! Place coordinates above geometry (in case they were placed)
    call coords % takeAboveGeom()

    ! Enter root universe
    r = coords % lvl(1) % r
    dir = coords % lvl(1) % dir
    uni => self % geom % unis % getPtr_fast(self % geom % rootIdx)

    call uni % enter(coords % lvl(1), r, dir)

    coords % lvl(1) % uniRootID = 1

    ! Dive to material
    call self % diveToMat(coords, 1)

  end subroutine placeCoord

  !!
  !! Find material and unique cell at a given location
  !!
  !! See geometry_inter for details
  !!
  subroutine whatIsAt(self, matIdx, uniqueID, r, u)
    class(geometryStd), intent(in)                    :: self
    integer(shortInt), intent(out)                    :: matIdx
    integer(shortInt), intent(out)                    :: uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    type(coordList)                                   :: coords
    real(defReal), dimension(3)                       :: u_l

    ! Select direction
    if (present(u)) then
      u_l = u
    else
      u_l = [ONE, ZERO, ZERO]
    end if

    ! Initialise coordinates
    call coords % init(r, u_l)

    ! Place coordinates
    call self % placeCoord(coords)

    ! Return material & uniqueID
    matIdx   = coords % matIdx
    uniqueID = coords % uniqueID

  end subroutine whatIsAt

  !!
  !! Return Axis Aligned Bounding Box encompassing the geometry
  !!
  !! See geometry_inter for details
  !!
  function bounds(self)
    class(geometryStd), intent(in) :: self
    real(defReal), dimension(6) :: bounds
    class(surface), pointer     :: surf
    integer(shortInt)           :: i

    ! Get boundary surface
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)
    bounds = surf % boundingBox()

    ! Change the infinate dimension to 0.0
    do i = 1, 3
      if(bounds(i) <= -INF .and. bounds(i+3) >= INF) then
        bounds(i) = ZERO
        bounds(i+3) = ZERO
      end if
    end do

  end function bounds

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_noCache(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt)              :: surfIdx, level !, levelST
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'move (geometryStd_class.f90)'

    if (.not.coords % isPlaced()) then
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist(dist, surfIdx, level, coords, coords % nesting)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

  end subroutine move_noCache

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_withCache(self, coords, maxDist, event, cache)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'move_withCache (geometryStd_class.f90)'

    if (.not.coords % isPlaced()) then
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist_cache(dist, surfIdx, level, coords, cache)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway
      cache % lvl = 0

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist
      cache % lvl = 0

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist
      cache % dist(1:level-1) = cache % dist(1:level-1) - dist
      cache % lvl = level - 1

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

  end subroutine move_withCache

  !!
  !! Move a particle in the top (global) level in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine moveGlobal(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    class(surface), pointer        :: surf
    real(defReal)                  :: dist

    ! Get boundary surface
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)

    ! Find distance to the boundary
    dist = surf % distance(coords % lvl(1) % r, coords % lvl(1) % dir)

    ! Select collision or boundary hit
    if (maxDist < dist) then ! maxDist is shorter
      ! Move above the geometry
      call coords % moveGlobal(maxDist)
      event = COLL_EV
      maxDist = maxDist

    else
      ! Move to boundary and apply BC
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)
      maxDist = dist

    end if

    ! Return particle to geometry
    call self % placeCoord(coords)

  end subroutine moveGlobal

  !!
  !!
  !!
  subroutine move_partial(self, coords, maxDist, event, levelST)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt), intent(in)  :: levelST
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'move_partial (geometryStd_class.f90)'

    if (.not.coords % isPlaced()) then
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface

    call self % closestDist(dist, surfIdx, level, coords, levelST)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, min(coords % nesting, levelST))
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway

      ! Get material
      call self % diveToMat(coords, min(coords % nesting, levelST))

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

   ! print *, levelST, coords % nesting

  end subroutine move_partial

  !!
  !!
  !! Move a particle in the top level without stopping
  !!
  !! See geometry_inter for details
  !!
  !! Uses co-ordinate transform boundary XSs
  !!
  subroutine teleport(self, coords, dist)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: dist
    class(surface), pointer        :: surf

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates back into geometry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % matIdx == OUTSIDE_MAT) then
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % transformBC(coords % lvl(1) % r, &
                              coords % lvl(1) % dir)

      ! Return particle to geometry
      call self % placeCoord(coords)
    end if

  end subroutine teleport

  !!
  !! Returns the list of active materials used in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! NOTE: This function uses VOID_MAT and UNDEF_MAT from universalVariables
  !!
  function activeMats(self, below) result(matList)
    class(geometryStd), intent(in)               :: self
    integer(shortInt), intent(in), optional      :: below
    integer(shortInt), dimension(:), allocatable :: matList
    integer(shortInt)                            :: N, lastIdx

    ! Use different procedure if materials from lover levels are requested
    if (present(below)) then
      matList = self % activeMatsBelow(below)
      return
    end if

    ! Takes the list of materials present in the geometry from geomGraph
    N = size(self % geom % graph % usedMats)
    lastIdx = self % geom % graph % usedMats(N)

    ! Check if the last entry of the list is an actual material or void
    if (lastIdx == VOID_MAT) then
      N = N - 1
      lastIdx = self % geom % graph % usedMats(N)
    end if
    ! Check if the last entry of the list is an undefined material
    if (lastIdx == UNDEF_MAT) then
      matList = self % geom % graph % usedMats(1:N-1)
    else
      matList = self % geom % graph % usedMats(1:N)
    end if

  end function activeMats

  !!
  !! Return the list of materials used in geometry below a given nesting level
  !!
  !! Args:
  !!   lvl [in] -> Level below or at which to search for materials
  !!
  !! Result:
  !!   Sorted list of materials with VOID and OUTSIDE removed
  !!
  !! Errors:
  !!  fatalError if lvl is < 1
  !!
  function activeMatsBelow(self, lvl) result(matList)
    class(geometryStd), intent(in)               :: self
    integer(shortInt), intent(in)                :: lvl
    integer(shortInt), dimension(:), allocatable :: matList
    type(intMap)                                 :: matNesting
    integer(shortInt)                            :: N, it, nesting, matIdx
    character(100), parameter :: Here = 'activeMatsBelow (geometryStd_class.f90)'

    if (lvl < 1) then
      call fatalError(Here, 'lvl is less than 1')
    end if
    call self % geom % graph % materialNesting(matNesting)

    ! Get list of materials
    allocate(matList(matNesting % length()))
    N = 1
    it = matNesting % begin()
    do while(it /= matNesting % end())
      matIdx =  matNesting % atKey(it)
      nesting = matNesting % atVal(it)

      ! Do not include materials without data
      if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT .or. matIdx == UNDEF_MAT) then
        it = matNesting % next(it)
        cycle
      end if

      ! Add to list
      if (nesting >= lvl) then
        matList(N) = matIdx
        N = N + 1
      end if
      it = matNesting % next(it)
    end do

    ! Trim list
    matList = matList(1:N-1)
    call quickSort(matList)

  end function activeMatsBelow


  !!
  !! Descend down the geometry structure untill material is reached
  !!
  !! Requires strting level to be specified.
  !! It is private procedure common to all movment types in geometry.
  !!
  !! Args:
  !!   coords [inout] -> CoordList of a particle. Assume thet coords are already valid for all
  !!     levels above and including start
  !!   start [in] -> Starting level for meterial search
  !!
  !! Errors:
  !!   fatalError if material cell is not found untill maximum nesting is reached
  !!
  subroutine diveToMat(self, coords, start)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt), intent(in)  :: start
    integer(shortInt)              :: rootID, localID, fill, id, i
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: offset
    character(100), parameter :: Here = 'diveToMat (geometryStd_class.f90)'

    do i = start, HARDCODED_MAX_NEST
      ! Find cell fill
      rootId = coords % lvl(i) % uniRootID
      localID = coords % lvl(i) % localID
      call self % geom % graph % getFill(fill, id, rootID, localID)

      if (fill >= 0) then ! Found material cell
        coords % matIdx   = fill
        coords % uniqueID = id
        return

      else ! Universe fill descend a level
        if (i == HARDCODED_MAX_NEST) exit ! If there is nested universe at the lowest level

        fill = abs(fill)

        ! Get current universe
        uni => self % geom % unis % getPtr_fast(coords % lvl(i) % uniIdx)

        ! Get cell offset
        offset = uni % cellOffset(coords % lvl(i))

        ! Get nested universe
        uni => self % geom % unis % getPtr_fast(fill)

        ! Enter nested univers
        call coords % addLevel()
        call uni % enter(coords % lvl(i+1), coords % lvl(i) % r - offset, coords % lvl(i) % dir)
        coords % lvl(i+1) % uniRootID = id ! Must be after enter where coord has intent out

      end if
    end do

    call fatalError(Here, 'Failed to find material cell. Should not happen after &
                          &geometry checks during build...')

  end subroutine diveToMat

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   dist [out]    -> Value of closest distance
  !!   surfIdx [out] -> Surface index for the crossing returned from the universe
  !!   lvl     [out] -> Level at which crossing is closest
  !!   coords [in]   -> Current coordinates of a particle
  !!
  subroutine closestDist(self, dist, surfIdx, lvl, coords, levelST)
    class(geometryStd), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    integer(shortInt), intent(out) :: lvl
    type(coordList), intent(in)    :: coords
    integer(shortInt), intent(in)  :: levelST
    integer(shortInt)              :: l, test_idx
    real(defReal)                  :: test_dist
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'closestDist (geometryStd_class.f90)'

    dist = INF
    surfIdx = 0
    lvl = 0


    if (levelST > coords % nesting) then
      print *, levelST, coords % nesting
      call fatalError(Here, 'The level of surface tracking should be smaller than the maximum nesting level')
    end if
    do l = 1, levelST !coords % nesting
      ! Get universe
      uni => self % geom % unis % getPtr_fast(coords % lvl(l) % uniIdx)

      ! Find distance
      call uni % distance(test_dist, test_idx, coords % lvl(l))

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - test_dist) >= dist * FP_REL_TOL) then
        dist = test_dist
        surfIdx = test_idx
        lvl = l
      end if

    end do


  end subroutine closestDist

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   dist [out]    -> Value of closest distance
  !!   surfIdx [out] -> Surface index for the crossing returned from the universe
  !!   lvl     [out] -> Level at which crossing is closest
  !!   coords [in]   -> Current coordinates of a particle
  !!   cache [inout] -> Distance cache. Use valid distances from cache. Put calculated
  !!     distances on the cache.
  !!
  subroutine closestDist_cache(self, dist, surfIdx, lvl, coords, cache)
    class(geometryStd), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    integer(shortInt), intent(out) :: lvl
    type(coordList), intent(in)    :: coords
    type(distCache), intent(inout) :: cache
    integer(shortInt)              :: l, test_idx
    real(defReal)                  :: test_dist
    class(universe), pointer       :: uni

    dist = INF
    surfIdx = 0
    lvl = 0
    do l = 1, coords % nesting

      ! Update Cache if distance is not valid
      if (cache % lvl < l) then
        ! Get universe
        uni => self % geom % unis % getPtr_fast(coords % lvl(l) % uniIdx)

        ! Find distance
        call uni % distance(cache % dist(l), cache % surf(l), coords % lvl(l))
        cache % lvl = cache % lvl + 1
      end if

      ! Read distance and crossing memento from cache
      test_dist = cache % dist(l)
      test_idx  = cache % surf(l)

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - test_dist) >= dist * FP_REL_TOL) then
        dist = test_dist
        surfIdx = test_idx
        lvl = l
      end if

    end do
  end subroutine closestDist_cache

  !!
  !! Cast geometry pointer to geometryStd class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class geometry
  !!
  !! Result:
  !!   Null if source is not of geometryStd class
  !!   Target points to source if source is geometryStd class
  !!
  pure function geometryStd_CptrCast(source) result(ptr)
    class(geometry), pointer, intent(in) :: source
    class(geometryStd), pointer          :: ptr

    select type(source)
      class is(geometryStd)
        ptr => source

      class default
        ptr => null()
    end select

  end function geometryStd_CptrCast

end module geometryStd_class
