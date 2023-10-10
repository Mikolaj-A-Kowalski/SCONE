module transportOperator_inter

  use numPrecision
  use universalVariables
  use genericProcedures,          only : fatalError

  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary

  ! Geometry interfaces
  use geometryReg_mod,            only : gr_geomPtr => geomPtr
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase



  implicit none
  private


  !!
  !! This is an abstract interface for all types of transport processing
  !!   -> This interface only deals with scalar processing of particle transport
  !!   -> Assumes that particle moves without any external forces (assumes that particle
  !!      moves along straight lines between collisions)
  !!
  !! Public interface:
  !!   transport(p, tally, thisCycle, nextCycle) -> given particle, tally and particle dungeons
  !!     for particles in this and next cycle performs movement of a particle in the geometry.
  !!     Sends transition report to the tally. Sends history report as well if particle dies.
  !!   init(dict, geom) -> initialises transport operator from a dictionary and pointer to a
  !!                       geometry
  !!
  !! Customisable procedures or transport actions
  !!   transit(p, tally, thisCycle, nextCycle) -> implements movement from collision to collision
  !!
  type, abstract, public :: transportOperator
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(nuclearDatabase), pointer :: xsData => null()

    !! Geometry pointer -> public so it can be used by subclasses (protected member)
    class(geometry), pointer         :: geom        => null()

  contains
    ! Public interface
    procedure, non_overridable :: transport
    procedure :: activeMaterials

    ! Extendable initialisation and deconstruction procedure
    procedure :: init
    procedure :: kill

    ! Customisable deferred procedures
    procedure(transit), deferred :: transit

  end type transportOperator

  ! Expandable procedures
  public :: init
  public :: kill


  abstract interface
    !!
    !! Move particle from collision to collision
    !!  Kill particle if needed
    !!
    subroutine transit(self, p, tally, thisCycle, nextCycle)
      import :: transportOperator, &
                particle, &
                tallyAdmin, &
                particleDungeon
      class(transportOperator), intent(inout) :: self
      class(particle), intent(inout)          :: p
      type(tallyAdmin), intent(inout)         :: tally
      class(particleDungeon), intent(inout)   :: thisCycle
      class(particleDungeon), intent(inout)   :: nextCycle
    end subroutine transit
  end interface

contains

  !!
  !! Master non-overridable subroutine to perform transport
  !!  Performs everything common to all types of transport
  !!
  subroutine transport(self, p, tally, thisCycle, nextCycle)
    class(transportOperator), intent(inout) :: self
    class(particle), intent(inout)          :: p
    type(tallyAdmin), intent(inout)         :: tally
    class(particleDungeon), intent(inout)   :: thisCycle
    class(particleDungeon), intent(inout)   :: nextCycle
    character(100),parameter :: Here ='transport (transportOperator_inter.f90)'

    ! Get nuclear data pointer form the particle
    self % xsData => ndReg_get(p % getType())

    ! Save geometry pointer
    self % geom => gr_geomPtr(p % geomIdx)

    ! Save pre-transition state
    call p % savePreTransition()

    ! Perform transit
    call self % transit(p, tally, thisCycle, nextCycle)

    ! Send history reports if particle died
    if( p  % isDead) then
      call tally % reportHist(p)
    end if

  end subroutine transport

  !!
  !! Initialise transport operator from dictionary and geometry
  !!
  subroutine init(self, dict)
    class(transportOperator), intent(inout)  :: self
    class(dictionary), intent(in)            :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Free memory. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(transportOperator), intent(inout) :: self

    self % geom   => null()
    self % xsData => null()

  end subroutine kill

  !!
  !! This procedure may be override to allow transport operators to limit the active materials
  !!
  !! Physics Packages should get the list of the active materials through the transport operator
  !! not from the geometry directly
  !!
  !! NOTE:
  !!   The setup with the geometry as an argument is quite ugly, but it is necessary
  !!   since the `geom` member pointer is initialised only during transport and
  !!   we need to know the active materials at the start of the simulation
  !!
  !! Args:
  !!  loc_geom [in] -> A geometry to be used to get the list of active materials
  !!
  !! Result:
  !!  Sorted list of active materials matIdxs
  !!
  function activeMaterials(self, loc_geom) result(mats)
    class(transportOperator), intent(in) :: self
    class(geometry), intent(in)          :: loc_geom
    integer(shortInt), dimension(:), allocatable :: mats

    ! Cannot assume that self % geom is associated or valid!
    mats = loc_geom % activeMats()

  end function activeMaterials

end module transportOperator_inter
