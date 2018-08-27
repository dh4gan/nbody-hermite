subroutine integrate(deltat,position,velocity,newposition,newvelocity)
  ! This subroutine drives the 4th-order Hermite integration
  ! (and timestep calculation)
  ! The integration proceeds via a predictor-corrector algorithm

  use nbodydata
  implicit none

  real, dimension(3,N), intent(in) :: position,velocity
  real, dimension(3,N),intent(out) :: newposition,newvelocity
  real, intent(out) :: deltat

  real, dimension(3,N) :: position_p,velocity_p
  real, dimension(3,N) :: acceleration,jerk,snap,crackle
  real, dimension(3,N) :: acceleration_p,jerk_p

  real :: dt2,dt3
  logical :: calc_snapcrackle

  position_p(:,:) = 0.0
  velocity_p(:,:) = 0.0

  acceleration(:,:) = 0.0
  acceleration_p(:,:) = 0.0

  jerk(:,:) = 0.0
  jerk_p(:,:) = 0.0

  snap(:,:) = 0.0
  crackle(:,:) = 0.0

  ! ***********************************
  ! 1. Predict positions and velocities
  ! ***********************************

  ! Compute accelerations, jerks, snaps and crackles at this timestep

  calc_snapcrackle = .true.
  call calc_acceleration(position,velocity,acceleration,jerk,snap,crackle,calc_snapcrackle)

  ! Compute timestep
  call timestep(deltat,acceleration,jerk,snap,crackle)

  dt2 = deltat*deltat
  dt3 = deltat*deltat*deltat

  ! Compute predicted positions and velocities for the next timestep
  position_p(:,:) = position(:,:) + velocity(:,:)*deltat + &
       0.5*acceleration(:,:)*dt2 + 0.1666*jerk(:,:)*dt3
  velocity_p(:,:) = velocity(:,:) + acceleration(:,:)*deltat + &
       0.5*jerk(:,:)*dt2

  ! *************************************
  ! 2. Correct positions and velocities
  ! *************************************

  ! Compute accelerations, and jerks based on predicted position, velocity
  calc_snapcrackle = .false. ! Don't calculate snap and crackle this time

  call calc_acceleration(position_p,velocity_p, acceleration_p, jerk_p,snap,crackle,calc_snapcrackle)

  ! Compute corrected positions and velocities
  newvelocity(:,:) = velocity(:,:) &
       + 0.5*(acceleration(:,:)+acceleration_p(:,:)) * deltat &
       + 0.0833* (jerk(:,:)+jerk_p(:,:))*dt2

  newposition(:,:) = position(:,:) &
       + 0.5*(velocity(:,:)+newvelocity(:,:))*deltat &
       + 0.0833*(acceleration(:,:)+acceleration_p(:,:))*dt2

end subroutine integrate
