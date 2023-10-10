!!
!! Transport operator for partial tracking
!!
module transportOperatorPT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with partial tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorPT
    integer(shortInt)             :: levelST
  contains
    procedure :: init
    procedure :: transit => partialtracking
    procedure, private :: partialtracking
  end type transportOperatorPT

contains

  subroutine init(self, dict)
    class(transportOperatorPT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

    ! Initialise this class
    call dict % getOrDefault(self % levelST,'levelST',1_shortInt)

  end subroutine init

  subroutine partialtracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorPT), intent(inout)              :: self
    class(particle), intent(inout)                         :: p
    type(tallyAdmin), intent(inout)                        :: tally
    class(particleDungeon), intent(inout)                  :: thisCycle
    class(particleDungeon), intent(inout)                  :: nextCycle
    integer(shortInt)                                      :: event
    real(defReal)                                          :: sigmaT, dist, majorant_inv
    logical(defBool)                                       :: pureST
    character(100), parameter :: Here = 'partialtracking (transportOIperatorPT_class.f90)'

    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    PTLoop: do

      pureST = p % coords % nesting <= self % levelST

      if (pureST) then
        if( p % matIdx() == VOID_MAT) then
          dist = INFINITY
        else
          !print *, p % matIdx()
          sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

          dist= -log( p % pRNG % get()) / sigmaT
        end if

        call p % savePrePath()
        call self % geom % move_partial(p % coords, dist, event, self % levelST)
        call tally % reportPath(p, dist)

        if( p % matIdx() == OUTSIDE_FILL) then
          p % isDead = .true.
          p % fate = LEAK_FATE
        end if

        if( event == COLL_EV .or. p % isDead) exit PTLoop

      else
        dist = -log( p% pRNG % get() ) * majorant_inv
        call self % geom % move_partial(p % coords, dist, event, self % levelST)

        if (p % matIdx() == OUTSIDE_FILL) then
          p % fate = LEAK_FATE
          p % isDead = .true.
          return
        end if

        if ( event == CROSS_EV .or. event == BOUNDARY_EV) cycle PTLoop
        if( p % matIdx() == VOID_MAT) cycle PTLoop

        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      end if

      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      if (.not. pureST) then
        if (p % pRNG % get() < sigmaT * majorant_inv) exit PTLoop
      end if

    end do PTLoop
    call tally % reportTrans(p)

  end subroutine partialtracking

end module transportOperatorPT_class
