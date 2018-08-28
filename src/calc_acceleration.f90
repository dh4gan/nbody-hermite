subroutine calc_acceleration(position,velocity,acceleration,jerk,snap,crackle, calc_snapcrackle)
  ! Given positions and accelerations, compute the gravitational acceleration (brute force)
  ! (and the first three derivatives: the jerk, snap and crackle)
  ! Snap and crackle only computed when calc_snapcrackle=True

  use nbodydata
  implicit none

  real, dimension(3,N),intent(in) :: position,velocity
  logical, intent(in):: calc_snapcrackle
  real, dimension(3,N), intent(out) :: acceleration, jerk, snap,crackle

  real, dimension(3,N,N) :: accel_ij, jerk_ij, snap_ij, crackle_ij
  real,dimension(N,N) :: alpha_ij, beta_ij, gamma_ij, relpos,r2,r3

  real :: relv2, rdotV,rdotA,rdotJ, vdotA
  real, dimension(3,N,N) :: rsep, vsep, asep, jsep


  acceleration(:,:) = 0.0
  jerk(:,:) = 0.0

  accel_ij(:,:,:) = 0.0
  jerk_ij(:,:,:) = 0.0
  alpha_ij(:,:) = 0.0


  rsep(:,:,:) = 0.0
  vsep(:,:,:) = 0.0
  asep(:,:,:) = 0.0


  if(calc_snapcrackle) then
     snap(:,:) = 0.0
     crackle(:,:) = 0.0
     snap_ij(:,:,:) = 0.0
     crackle_ij(:,:,:) = 0.0
     jsep(:,:,:) = 0.0
     beta_ij(:,:) = 0.0
     gamma_ij(:,:) = 0.0
  endif

  !
  ! 1. First, compute all accelerations and jerks
  !

  do ibody=1,N

     do jbody=1,N

        if(ibody==jbody) cycle

        ! Compute relative separation of ibody and jbody
        do ix=1,3
           rsep(ix,ibody,jbody) = position(ix,jbody) - position(ix,ibody)
           vsep(ix,ibody,jbody) = velocity(ix,jbody) - velocity(ix,ibody)
        enddo

        relpos(ibody,jbody) = sqrt(rsep(1,ibody,jbody)*rsep(1,ibody,jbody) + &
             rsep(2,ibody,jbody)*rsep(2,ibody,jbody)+ &
             rsep(3,ibody,jbody)*rsep(3,ibody,jbody)+rsoft*rsoft)



        r2(ibody,jbody) = relpos(ibody,jbody)*relpos(ibody,jbody)
        r3(ibody,jbody) = relpos(ibody,jbody)*r2(ibody,jbody)

        rdotV = 0.0

        do ix=1,3
           rdotV = rdotV + rsep(ix,ibody,jbody)*vsep(ix,ibody,jbody)
        enddo

        alpha_ij(ibody,jbody) = rdotV/r2(ibody,jbody)

        ! Compute acceleration and jerk components
        do ix=1,3
           accel_ij(ix,ibody,jbody) = G*mass(jbody)*rsep(ix,ibody,jbody)/r3(ibody,jbody)
           jerk_ij(ix,ibody,jbody) = G*mass(jbody)*vsep(ix,ibody,jbody)/r3(ibody,jbody) - &
                3.0*alpha_ij(ibody,jbody)*accel_ij(ix,ibody,jbody)

           acceleration(ix,ibody) = acceleration(ix,ibody) + accel_ij(ix,ibody,jbody)
           jerk(ix,ibody) = jerk(ix,ibody) + jerk_ij(ix,ibody,jbody)

        enddo

     enddo
  enddo


  !
  ! 2. Now, compute snap and crackle (if required)
  !

  if(calc_snapcrackle) then

     do ibody=1,N

        do jbody=1,N

           if(ibody==jbody) cycle


           !
           ! Compute relative acceleration and jerk
           !

           do ix=1,3

              asep(ix,ibody,jbody) = acceleration(ix,jbody) - acceleration(ix,ibody)
              jsep(ix,ibody,jbody) = jerk(ix,jbody) - jerk(ix,ibody)

           enddo


           !
           ! Compute beta and gamma parameters
           !

           rdotA = 0.0
           vdotA = 0.0
           rdotJ = 0.0

           do ix=1,3
              rdotA = rdotA + rsep(ix,ibody,jbody)*asep(ix,ibody,jbody)
              vdotA = vdotA + vsep(ix,ibody,jbody)*asep(ix,ibody,jbody)
              rdotJ = rdotJ + rsep(ix,ibody,jbody)*jsep(ix,ibody,jbody)
           enddo


           relv2 =     vsep(1,ibody,jbody)*vsep(1,ibody,jbody) + &
                vsep(2,ibody,jbody)*vsep(2,ibody,jbody)+ &
                vsep(3,ibody,jbody)*vsep(3,ibody,jbody)

           beta_ij(ibody,jbody) = (relv2 + rdotA)/r2(ibody,jbody) + &
                alpha_ij(ibody,jbody)*alpha_ij(ibody,jbody)


           gamma_ij(ibody,jbody) = (3.0*vdotA + rdotJ)/r2(ibody,jbody) + &
                alpha_ij(ibody,jbody)*(3.0*beta_ij(ibody,jbody)- alpha_ij(ibody,jbody))


           !
           ! Now compute snap and crackle
           !
           
           do ix=1,3

              snap_ij(ix,ibody,jbody) = G*mass(jbody)/r3(ibody,jbody)*asep(ix,ibody,jbody) &
                   - 6.0*alpha_ij(ibody,jbody)*jerk_ij(ix,ibody,jbody) &
                   - 3.0*beta_ij(ibody,jbody)*accel_ij(ix,ibody,jbody)


              snap(ix,ibody) = snap(ix,ibody) + snap_ij(ix,ibody,jbody)

              crackle_ij(ix,ibody,jbody) = G*mass(jbody)*jsep(ix,ibody,jbody)/r2(ibody,jbody) &
                   - 9.0*alpha_ij(ibody,jbody)*snap_ij(ix,ibody,jbody) &
                   - 9.0*beta_ij(ibody,jbody)*jerk_ij(ix,ibody,jbody) &
                   - 3.0*gamma_ij(ibody,jbody)*accel_ij(ix,ibody,jbody)


              crackle(ix,ibody) = crackle(ix,ibody) + crackle_ij(ix,ibody,jbody)

           enddo
        enddo ! end loop over jbody
     enddo ! end loop over ibody

  endif

  !print*, r3(1,2), r2(1,2), relpos(1,2), position(:,1), position(:,2)
  print*, 'POSITION: ',position(:,2), relpos(1,2)
  !print*, 'VELOCITY: ', velocity(:,2)
  print*, 'ACCELERATION: ', acceleration(:,2)
  !print*, 'JERK: ', jerk(:,2)
  !print*, 'SNAP: ', snap(:,2)
  !print*, 'CRACKLE: ', crackle(:,2)

  !STOP
  return

end subroutine calc_acceleration
