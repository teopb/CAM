module held_suarez_cam

  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: held_suarez_init, held_suarez_tend, held_suarez_readnl

  real(r8) :: held_suarez_efoldf  ! efolding time for wind dissipation
  real(r8) :: held_suarez_efolda  ! efolding time for T dissipation
  real(r8) :: held_suarez_efolds  ! efolding time for T dissipation
  real(r8) :: held_suarez_sigmab  ! threshold sigma level
  real(r8) :: held_suarez_t00    ! minimum reference temperature
  real(r8) :: held_suarez_delta_T_y  ! equilibrium temperature parameter for latitude
  real(r8) :: held_suarez_delta_theta_z  ! equilibrium temperature parameter for veritical level

!=======================================================================
contains
!=======================================================================

  subroutine held_suarez_readnl(nlfile)
    !
    ! held_suarez_readnl: Read in parameters controlling Held-Suarez parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    !
    ! Passed Variables
    !------------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !--------------
    integer:: ierr,unitn

    character(len=*), parameter :: sub = 'held_suarez_readnl'

    namelist /held_suarez_nl/ held_suarez_efoldf, held_suarez_efolda        , held_suarez_efolds   , &
                              held_suarez_sigmab      , held_suarez_t00    , held_suarez_delta_T_y , &
                              held_suarez_delta_theta_z

    ! Read in namelist values
    !-------------------------
    if(masterproc) then
        unitn = getunit()
        open(unitn,file=trim(nlfile),status='old')
        call find_group_name(unitn,'held_suarez_nl',status=ierr)
        if(ierr == 0) then
        read(unitn,held_suarez_nl,iostat=ierr)
        if(ierr /= 0) then
            call endrun(sub//': ERROR reading namelist')
        endif
        endif
        close(unitn)
        call freeunit(unitn)
    endif


    call mpi_bcast(held_suarez_efoldf    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_efoldf")
    call mpi_bcast(held_suarez_efolda  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_efolda")
    call mpi_bcast(held_suarez_efolds  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_efolds")
    call mpi_bcast(held_suarez_sigmab   , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_sigmab")
    call mpi_bcast(held_suarez_t00        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_t00")
    call mpi_bcast(held_suarez_delta_T_y      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_delta_T_y")
    call mpi_bcast(held_suarez_delta_theta_z      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_delta_theta_z")

  end subroutine held_suarez_readnl

  subroutine held_suarez_init()
    use physics_buffer,     only: physics_buffer_desc
    use cam_history,        only: addfld, add_default
    use ref_pres,           only: psurf_ref
    use held_suarez_1994,   only: held_suarez_1994_init

    ! Local variables
    character(len=512) :: errmsg
    integer            :: errflg

    ! Set model constant values
    call held_suarez_1994_init(psurf_ref, errmsg, errflg, held_suarez_efoldf, held_suarez_sigmab, held_suarez_efolda, held_suarez_efolds)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
    call add_default('QRS', 1, ' ')
  end subroutine held_suarez_init

  subroutine held_suarez_tend(state, ptend, ztodt)
    use air_composition,    only: cappav, cpairv
    use ref_pres,           only: pref_mid_norm
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use held_suarez_1994,   only: held_suarez_1994_run

    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    integer                            :: i, k             ! Longitude, level indices

    character(len=64)                  :: scheme_name      ! CCPP-required variables (not used in CAM)
    character(len=512)                 :: errmsg
    integer                            :: errflg

    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)

    call held_suarez_1994_run(held_suarez_sigmab, held_suarez_t00, held_suarez_delta_T_y, held_suarez_delta_theta_z, pver, ncol, pref_mid_norm, clat, cappav(1:ncol,:,lchnk), &
                              cpairv(1:ncol,:,lchnk), state%pmid(1:ncol,:),            &
                              state%u(1:ncol,:), state%v(1:ncol,:), state%t(1:ncol,:), &
                              ptend%u(1:ncol,:), ptend%v(1:ncol,:), ptend%s(1:ncol,:), &
                              scheme_name, errmsg, errflg)

    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)

  end subroutine held_suarez_tend

end module held_suarez_cam
