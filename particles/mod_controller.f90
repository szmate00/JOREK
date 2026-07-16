!> Module for feedback controllers acting on the plasma simulation.
!>
!> A controller measures a plasma quantity (its "sensor"), compares it against a (possibly time
!> dependent) setpoint, and computes an actuator signal with a PID control law. Controllers are
!> configured in the input file via the controllers() array (see type_controller_config in
!> phys_module) and are linked to their actuator; currently the only actuator is the puff rate of a
!> puffing action, linked via part_group_config()%puff_ctrl()%controller_num (see
!> mod_particle_puffing).
!>
!> Available sensors:
!>  - 'te_front': height above the X-point [m] of the electron temperature front, i.e. the highest
!>    point along a vertical line from the X-point towards the magnetic axis where Te <= Te_front_eV.
!>    This is a continuous measure of the position of an X-point radiator (XPR), following the
!>    control strategy used experimentally in ASDEX Upgrade (Bernert et al., NF 2021) and in
!>    EDGE2D/SOLPS-ITER XPR modelling (Korving et al.; Poletaeva et al., CPP 2026). The measurement
!>    is evaluated on the n=0-dominated fields at toroidal angle phi=0.
!>
!> Controller types:
!>  - 'closedloop':  actuator = clip(PID output)
!>  - 'feedforward': actuator = clip(feedforward baseline + PID output), where the baseline is the
!>                   piecewise linear puff_ctrl times/rates waveform
!>
!> Features: conditional-integration anti-windup, exponential smoothing of the measurement
!> (tau_smooth), first-order actuator lag (tau_actuator), arming time (t_start; before it only the
!> feedforward baseline is applied and the integrator is held), and a time dependent setpoint from
!> either a piecewise linear schedule (setpoint_times/setpoint_values) or a data file (setpoint_file).
!>
!> All MPI ranks compute the controller identically from the (rank-replicated) fields; only rank 0
!> writes the diagnostic trace file jorek_controller_<num>.dat and the state file
!> jorek_controller_state_<num>.dat. The state file is read back on restart (restart_particles) so
!> that the integrator and filters continue where they left off.

module mod_controller
  use phys_module, only: type_controller_config, controllers, n_controllers_max, &
                         n_ctrl_segment_max, xpoint, xcase, restart_particles
  use mod_particle_sim
  use constants,   only: K_BOLTZ, EL_CHG
  use equil_info,  only: find_xpoint

  implicit none

  private
  public :: controller_update_puff_rate, validate_controller_config

  !> Runtime state of one controller (the static configuration lives in phys_module::controllers)
  type :: type_controller_state
    logical :: initialized = .false.

    ! --- 'te_front' sensor geometry (fixed after initialization)
    real*8               :: R_X, Z_X       !< location of the active X-point [m]
    real*8               :: dir            !< +1 if the magnetic axis lies above the X-point, -1 otherwise
    integer              :: n_pts          !< number of sensor points along the line
    real*8,  allocatable :: h_pts(:)       !< height above the X-point of each sensor point [m]
    integer, allocatable :: i_elm(:)       !< element index of each sensor point
    real*8,  allocatable :: s_pts(:)       !< local s coordinate of each sensor point
    real*8,  allocatable :: t_pts(:)       !< local t coordinate of each sensor point
    logical, allocatable :: valid(:)       !< whether the sensor point lies inside the grid

    ! --- setpoint table read from setpoint_file
    integer             :: sp_len = 0
    real*8, allocatable :: sp_time(:), sp_val(:)

    ! --- control law state
    real*8  :: int_err     = 0.d0          !< time integral of the error
    real*8  :: e_prev      = 0.d0          !< error at the previous update (for the D term)
    logical :: have_e_prev = .false.
    real*8  :: meas_filt   = 0.d0          !< exponentially smoothed measurement
    logical :: have_meas   = .false.
    real*8  :: rate_out    = 0.d0          !< actuator signal after the first-order lag
    logical :: have_rate   = .false.
    real*8  :: t_prev      = 0.d0          !< time of the previous update [s]
    logical :: have_t_prev = .false.

    ! --- cache so that several puffing actions sharing one controller advance it only once per step
    integer :: istep_last  = -1
    real*8  :: rate_last   = 0.d0
  end type type_controller_state

  type(type_controller_state), dimension(n_controllers_max) :: ctrl_states

contains


!> Validity checks of the configuration of controller `num` (called from the puffing setup).
!> Only checks and prints; all ranks call it, messages are printed by rank 0.
subroutine validate_controller_config(num, my_id)
  implicit none
  integer, intent(in) :: num
  integer, intent(in) :: my_id

  integer :: n_times, n_vals, k

  if (num < 1 .or. num > n_controllers_max) then
    if (my_id == 0) write(*,"(A,I3,A,I3)") "ERROR [mod_controller]: controller_num = ", num, &
      " is outside the valid range 1..", n_controllers_max
    stop
  endif

  select case (trim(controllers(num)%type))
    case ('closedloop', 'feedforward')
      ! valid
    case ('none')
      if (my_id == 0) write(*,"(A,I2,A)") "ERROR [mod_controller]: controller ", num, &
        " is referenced by a puffing action but controllers()%type has not been set"
      stop
    case default
      if (my_id == 0) write(*,"(A,A,A)") "ERROR [mod_controller]: controller type '", &
        trim(controllers(num)%type), "' is not supported (use 'closedloop' or 'feedforward')"
      stop
  end select

  if (trim(controllers(num)%sensor) /= 'te_front') then
    if (my_id == 0) write(*,"(A,A,A)") "ERROR [mod_controller]: controller sensor '", &
      trim(controllers(num)%sensor), "' is not supported (currently only 'te_front')"
    stop
  endif

  ! --- exactly one way of defining the setpoint must be used
  n_times = count(controllers(num)%setpoint_times  >= 0.d0)
  n_vals  = count(controllers(num)%setpoint_values > -1.d89)
  if (trim(controllers(num)%setpoint_file) == 'none') then
    if (n_times > 0) then
      if (n_times /= n_vals) then
        if (my_id == 0) write(*,"(A,I2)") "ERROR [mod_controller]: mismatch between the number of "// &
          "entries of setpoint_times and setpoint_values for controller ", num
        stop
      endif
      do k=2, n_times
        if (controllers(num)%setpoint_times(k-1) >= controllers(num)%setpoint_times(k)) then
          if (my_id == 0) write(*,"(A,I2,A)") "ERROR [mod_controller]: setpoint_times of controller ", &
            num, " must be strictly increasing"
          stop
        endif
      enddo
    else if (controllers(num)%setpoint < -1.d89) then
      if (my_id == 0) write(*,"(A,I2,A)") "ERROR [mod_controller]: no setpoint defined for controller ", num, &
        " (set setpoint, setpoint_times/setpoint_values, or setpoint_file)"
      stop
    endif
  endif

  if (controllers(num)%max_value <= controllers(num)%min_value) then
    if (my_id == 0) write(*,"(A,I2)") "ERROR [mod_controller]: max_value <= min_value for controller ", num
    stop
  endif

  if (controllers(num)%zline_npoints < 2) then
    if (my_id == 0) write(*,"(A,I2)") "ERROR [mod_controller]: zline_npoints must be >= 2 for controller ", num
    stop
  endif
  if (controllers(num)%zline_length <= 0.d0) then
    if (my_id == 0) write(*,"(A,I2)") "ERROR [mod_controller]: zline_length must be > 0 for controller ", num
    stop
  endif

end subroutine validate_controller_config


!> Advance controller `num` by one fluid time step and return the puff rate [atoms/s].
!> `feedforward` is the piecewise linear times/rates baseline of the calling puffing action; it is
!> the output before arming and, for type 'feedforward', the baseline the PID correction is added to.
!> If several puffing actions share one controller, only the first call per fluid step advances the
!> controller; subsequent calls return the cached rate.
subroutine controller_update_puff_rate(num, sim, feedforward, puff_rate)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer,            intent(in)    :: num
  type(particle_sim), intent(inout) :: sim
  real*8,             intent(in)    :: feedforward
  real*8,             intent(out)   :: puff_rate

  real*8  :: meas_raw, setpoint_now, e, dt, alpha
  real*8  :: base, term_P, term_I, term_D, int_trial, cmd
  logical :: armed, saturated
  integer :: my_unit

  if (.not. ctrl_states(num)%initialized) call controller_init(num, sim)

  ! --- several puffing actions may share this controller: advance it only once per fluid step
  if (ctrl_states(num)%istep_last == sim%istep_fluid) then
    puff_rate = ctrl_states(num)%rate_last
    return
  endif

  associate (cfg => controllers(num), cs => ctrl_states(num))

  ! --- time step of the controller
  if (cs%have_t_prev) then
    dt = sim%time - cs%t_prev
  else
    dt = 0.d0
  endif

  ! --- measurement and exponential smoothing
  call measure_te_front(cs, cfg, sim, meas_raw)
  if (.not. cs%have_meas) then
    cs%meas_filt = meas_raw
    cs%have_meas = .true.
  else if (cfg%tau_smooth > 0.d0 .and. dt > 0.d0) then
    alpha = dt / (cfg%tau_smooth + dt)
    cs%meas_filt = cs%meas_filt + alpha * (meas_raw - cs%meas_filt)
  else
    cs%meas_filt = meas_raw
  endif

  setpoint_now = eval_setpoint(cs, cfg, sim%time)
  e = setpoint_now - cs%meas_filt

  armed = sim%time >= cfg%t_start

  if (.not. armed) then
    ! --- before arming: apply the feedforward baseline only, hold the integrator
    if (trim(cfg%type) == 'feedforward') then
      cmd = feedforward
    else
      cmd = cfg%min_value
    endif
    cs%have_e_prev = .false.  ! do not let the D term kick at the arming step
    term_P = 0.d0; term_I = 0.d0; term_D = 0.d0
  else
    if (trim(cfg%type) == 'feedforward') then
      base = feedforward
    else
      base = 0.d0
    endif

    term_P = cfg%K_p * e

    if (cs%have_e_prev .and. dt > 0.d0) then
      term_D = cfg%K_d * (e - cs%e_prev) / dt
    else
      term_D = 0.d0
    endif

    ! --- integrator with conditional-integration anti-windup:
    !     accept the new integral only if the resulting command is not saturated
    int_trial = cs%int_err + e * dt
    cmd       = base + term_P + cfg%K_i * int_trial + term_D
    saturated = (cmd > cfg%max_value) .or. (cmd < cfg%min_value)
    if (.not. saturated) then
      cs%int_err = int_trial
    else
      cmd = base + term_P + cfg%K_i * cs%int_err + term_D
    endif
    term_I = cfg%K_i * cs%int_err

    cs%e_prev      = e
    cs%have_e_prev = .true.
  endif

  cmd = min(max(cmd, cfg%min_value), cfg%max_value)

  ! --- first-order actuator lag
  if (cfg%tau_actuator > 0.d0 .and. cs%have_rate .and. dt > 0.d0) then
    alpha = dt / (cfg%tau_actuator + dt)
    cs%rate_out = cs%rate_out + alpha * (cmd - cs%rate_out)
  else
    cs%rate_out = cmd
  endif
  cs%have_rate = .true.

  puff_rate = min(max(cs%rate_out, cfg%min_value), cfg%max_value)

  if (.not. ieee_is_finite(puff_rate)) then
    if (sim%my_id == 0) write(*,"(A,I2,A)") "WARNING [mod_controller]: non-finite actuator signal from controller ", &
      num, ", falling back to the feedforward baseline"
    puff_rate   = min(max(feedforward, cfg%min_value), cfg%max_value)
    cs%rate_out = puff_rate
  endif

  cs%t_prev      = sim%time
  cs%have_t_prev = .true.
  cs%istep_last  = sim%istep_fluid
  cs%rate_last   = puff_rate

  ! --- diagnostics trace and restartable state (rank 0)
  if (sim%my_id == 0) then
    open(newunit=my_unit, file=trace_filename(num), status='unknown', position='append', action='write')
    write(my_unit,"(9ES16.8,I3)") sim%time, meas_raw, cs%meas_filt, setpoint_now, e, &
                                  term_P, term_I, term_D, puff_rate, merge(1, 0, armed)
    close(my_unit)

    open(newunit=my_unit, file=state_filename(num), status='replace', action='write')
    write(my_unit,"(5ES25.16,3I3)") cs%t_prev, cs%int_err, cs%e_prev, cs%meas_filt, cs%rate_out, &
                                    merge(1, 0, cs%have_e_prev), merge(1, 0, cs%have_meas),      &
                                    merge(1, 0, cs%have_rate)
    close(my_unit)
  endif

  end associate

end subroutine controller_update_puff_rate


!> Initialize controller `num`: locate the X-point, precompute the sensor points along the vertical
!> line from the X-point towards the magnetic axis, read the setpoint file if requested, and restore
!> the controller state from the state file when restarting.
subroutine controller_init(num, sim)
  use mpi_mod
  use profiles, only: readProf
  implicit none

  integer,            intent(in)    :: num
  type(particle_sim), intent(inout) :: sim

  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  real*8  :: psi_xp(2), R_xp(2), Z_xp(2), s_xp(2), t_xp(2)
  integer :: i_elm_axis, i_elm_xp(2), ifail, ix, k, n, ierr, my_unit
  real*8  :: R_out, Z_out, s_out, t_out
  integer :: i_elm_out
  real*8  :: state_buf(5)
  integer :: flag_buf(3)
  logical :: state_exists

  call validate_controller_config(num, sim%my_id)

  associate (cfg => controllers(num), cs => ctrl_states(num))

  ! --- locate the magnetic axis and the active X-point ------------------------------------
  if (.not. xpoint) then
    if (sim%my_id == 0) write(*,"(A)") "ERROR [mod_controller]: the 'te_front' sensor requires an X-point geometry (xpoint=.true.)"
    stop
  endif

  call find_axis(1, sim%fields%node_list, sim%fields%element_list, psi_axis, R_axis, Z_axis, &
                 i_elm_axis, s_axis, t_axis, ifail)
  if (ifail /= 0) then
    if (sim%my_id == 0) write(*,"(A)") "ERROR [mod_controller]: could not locate the magnetic axis"
    stop
  endif

  call find_xpoint(1, sim%fields%node_list, sim%fields%element_list, psi_xp, R_xp, Z_xp, &
                   i_elm_xp, s_xp, t_xp, xcase, ifail)
  if (ifail /= 0) then
    if (sim%my_id == 0) write(*,"(A)") "ERROR [mod_controller]: could not locate the X-point"
    stop
  endif

  ! pick the active X-point (same convention as in mod_particle_wall_interaction)
  ix = 1
  if ((xcase == 2) .or. ((xcase == 3) .and. (psi_xp(2) < psi_xp(1)))) ix = 2
  cs%R_X = R_xp(ix)
  cs%Z_X = Z_xp(ix)
  cs%dir = sign(1.d0, Z_axis - cs%Z_X)

  ! --- precompute the sensor points along the line ----------------------------------------
  n = cfg%zline_npoints
  cs%n_pts = n
  allocate(cs%h_pts(n), cs%i_elm(n), cs%s_pts(n), cs%t_pts(n), cs%valid(n))
  do k=1, n
    cs%h_pts(k) = (k-1) * cfg%zline_length / dble(n-1)
    call find_RZ(sim%fields%node_list, sim%fields%element_list, cs%R_X, cs%Z_X + cs%dir*cs%h_pts(k), &
                 R_out, Z_out, i_elm_out, s_out, t_out, ifail)
    cs%valid(k) = (ifail == 0) .and. (i_elm_out > 0)
    cs%i_elm(k) = i_elm_out
    cs%s_pts(k) = s_out
    cs%t_pts(k) = t_out
  enddo

  if (count(cs%valid) < 2) then
    if (sim%my_id == 0) write(*,"(A,I2,A)") "ERROR [mod_controller]: fewer than 2 sensor points of controller ", &
      num, " lie inside the grid; check zline_length"
    stop
  endif
  if ((sim%my_id == 0) .and. (count(cs%valid) < cs%n_pts)) then
    write(*,"(A,I2,A,I5,A,I5,A)") "WARNING [mod_controller]: only ", count(cs%valid), " of the ", cs%n_pts, &
      " sensor points of controller ", num, " lie inside the grid"
  endif

  ! --- read the setpoint file if requested (rank 0 reads, then broadcast) ------------------
  if (trim(cfg%setpoint_file) /= 'none') then
    if (sim%my_id == 0) call readProf(cs%sp_time, cs%sp_val, cs%sp_len, trim(cfg%setpoint_file))
    call MPI_BCAST(cs%sp_len, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (cs%sp_len < 1) then
      if (sim%my_id == 0) write(*,"(A,A,A,I2)") "ERROR [mod_controller]: could not read setpoint_file '", &
        trim(cfg%setpoint_file), "' of controller ", num
      stop
    endif
    if (sim%my_id /= 0) allocate(cs%sp_time(cs%sp_len), cs%sp_val(cs%sp_len))
    call MPI_BCAST(cs%sp_time, cs%sp_len, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cs%sp_val,  cs%sp_len, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    do k=2, cs%sp_len
      if (cs%sp_time(k-1) >= cs%sp_time(k)) then
        if (sim%my_id == 0) write(*,"(A,A,A)") "ERROR [mod_controller]: the times in setpoint_file '", &
          trim(cfg%setpoint_file), "' must be strictly increasing"
        stop
      endif
    enddo
  endif

  ! --- restore the controller state when restarting ----------------------------------------
  if (restart_particles) then
    state_exists = .false.
    if (sim%my_id == 0) inquire(file=state_filename(num), exist=state_exists)
    call MPI_BCAST(state_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (state_exists) then
      if (sim%my_id == 0) then
        open(newunit=my_unit, file=state_filename(num), status='old', action='read')
        read(my_unit,*) state_buf, flag_buf
        close(my_unit)
      endif
      call MPI_BCAST(state_buf, 5, MPI_REAL8,   0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(flag_buf,  3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      cs%t_prev      = state_buf(1)
      cs%int_err     = state_buf(2)
      cs%e_prev      = state_buf(3)
      cs%meas_filt   = state_buf(4)
      cs%rate_out    = state_buf(5)
      cs%have_t_prev = .true.
      cs%have_e_prev = flag_buf(1) == 1
      cs%have_meas   = flag_buf(2) == 1
      cs%have_rate   = flag_buf(3) == 1
      if (sim%my_id == 0) write(*,"(A,I2,A)") "NOTE [mod_controller]: restored the state of controller ", &
        num, " from "//trim(state_filename(num))
    else
      if (sim%my_id == 0) write(*,"(A,I2,A)") "WARNING [mod_controller]: restart requested but no state file found for controller ", &
        num, "; the controller starts from a fresh state"
    endif
  endif

  ! --- diagnostics trace file (rank 0) ------------------------------------------------------
  if (sim%my_id == 0) then
    if (restart_particles) then
      ! append to an existing trace when restarting; create it if missing
      open(newunit=my_unit, file=trace_filename(num), status='unknown', position='append', action='write')
    else
      open(newunit=my_unit, file=trace_filename(num), status='replace', action='write')
      write(my_unit,"(A)") "# time[s]  meas_raw  meas_filt  setpoint  error  term_P  term_I  term_D  puff_rate[atoms/s]  armed"
    endif
    close(my_unit)

    write(*,"(A,I2,A)")           "NOTE [mod_controller]: controller ", num, " initialized ("// &
      trim(cfg%type)//", sensor "//trim(cfg%sensor)//")"
    write(*,"(A,2F10.4,A,F6.1)")  "  X-point (R,Z) = ", cs%R_X, cs%Z_X, ", line direction dz = ", cs%dir
    write(*,"(A,F8.4,A,I5,A)")    "  sensing line length = ", cfg%zline_length, " m with ", cs%n_pts, " points"
  endif

  cs%initialized = .true.

  end associate

end subroutine controller_init


!> 'te_front' sensor: height above the X-point [m] of the topmost point along the sensing line where
!> Te <= Te_front_eV, linearly interpolated between the bracketing sensor points. Returns 0 when the
!> whole line is hotter than Te_front_eV (no radiator front above the X-point) and the full line
!> length when the whole line is colder (front above the sensing range; a warning is printed).
subroutine measure_te_front(cs, cfg, sim, height)
  implicit none
  type(type_controller_state),  intent(in)    :: cs
  type(type_controller_config), intent(in)    :: cfg
  type(particle_sim),           intent(inout) :: sim
  real*8,                       intent(out)   :: height

  real*8  :: Te_eV(cs%n_pts), n_e, T_e_K
  integer :: k, k_cold, k_hot

  do k=1, cs%n_pts
    if (cs%valid(k)) then
      call sim%fields%calc_NeTeTi(sim%time, cs%i_elm(k), [cs%s_pts(k), cs%t_pts(k)], 0.d0, &
                                  n_e=n_e, T_e=T_e_K)
      Te_eV(k) = T_e_K * K_BOLTZ / EL_CHG
    else
      Te_eV(k) = -1.d0 ! excluded from the scan below
    endif
  enddo

  ! --- topmost valid point that is at or below the front temperature
  k_cold = 0
  do k=cs%n_pts, 1, -1
    if (cs%valid(k) .and. (Te_eV(k) >= 0.d0) .and. (Te_eV(k) <= cfg%Te_front_eV)) then
      k_cold = k
      exit
    endif
  enddo

  if (k_cold == 0) then
    height = 0.d0 ! entire line is hot: no front above the X-point
    return
  endif

  ! --- next valid (hot) point above it, for sub-grid interpolation of the crossing
  k_hot = 0
  do k=k_cold+1, cs%n_pts
    if (cs%valid(k)) then
      k_hot = k
      exit
    endif
  enddo

  if (k_hot == 0) then
    height = cs%h_pts(cs%n_pts) ! front at or above the end of the sensing line
    if (sim%my_id == 0) write(*,"(A,F8.4,A)") "WARNING [mod_controller]: 'te_front' sensor saturated at ", &
      height, " m; the radiator front is above the sensing line (consider increasing zline_length)"
  else if (Te_eV(k_hot) > Te_eV(k_cold)) then
    height = cs%h_pts(k_cold) + (cs%h_pts(k_hot) - cs%h_pts(k_cold)) * &
             (cfg%Te_front_eV - Te_eV(k_cold)) / (Te_eV(k_hot) - Te_eV(k_cold))
  else
    height = cs%h_pts(k_cold)
  endif

end subroutine measure_te_front


!> Evaluate the (possibly time dependent) setpoint at time t [s]
function eval_setpoint(cs, cfg, t) result(sp)
  implicit none
  type(type_controller_state),  intent(in) :: cs
  type(type_controller_config), intent(in) :: cfg
  real*8,                       intent(in) :: t
  real*8                                   :: sp

  integer :: n

  if (cs%sp_len > 0) then
    sp = interp_table(cs%sp_time, cs%sp_val, cs%sp_len, t)
  else
    n = count(cfg%setpoint_times >= 0.d0)
    if (n > 0) then
      sp = interp_table(cfg%setpoint_times, cfg%setpoint_values, n, t)
    else
      sp = cfg%setpoint
    endif
  endif

end function eval_setpoint


!> Piecewise linear interpolation in a (times, values) table with flat extension outside the range
pure function interp_table(times, values, n, t) result(v)
  implicit none
  real*8,  intent(in) :: times(:), values(:)
  integer, intent(in) :: n
  real*8,  intent(in) :: t
  real*8              :: v

  integer :: k

  if (t <= times(1)) then
    v = values(1)
  else if (t >= times(n)) then
    v = values(n)
  else
    v = values(n) ! fallback, overwritten below
    do k=2, n
      if (t <= times(k)) then
        v = values(k-1) + (values(k) - values(k-1)) * (t - times(k-1)) / (times(k) - times(k-1))
        exit
      endif
    enddo
  endif

end function interp_table


function trace_filename(num) result(fname)
  implicit none
  integer, intent(in) :: num
  character(len=30)   :: fname
  write(fname,"(A,I2.2,A)") "jorek_controller_", num, ".dat"
end function trace_filename


function state_filename(num) result(fname)
  implicit none
  integer, intent(in) :: num
  character(len=36)   :: fname
  write(fname,"(A,I2.2,A)") "jorek_controller_state_", num, ".dat"
end function state_filename

end module mod_controller
