subroutine timestep(deltat,acceleration,jerk,snap,crackle)

  use nbodydata

  implicit none

  real,dimension(3,N),intent(in) :: acceleration,jerk,snap,crackle
  real, intent(out) :: deltat

  real :: dt_i,amag,jmag,smag,cmag

  deltat = 1.0e30
  amag = 0.0
  jmag = 0.0
  smag = 0.0
  cmag = 0.0


  do ibody=1,N

     do ix=1,3
        amag = amag + acceleration(ix,ibody)*acceleration(ix,ibody)
        jmag = jmag + jerk(ix,ibody)*jerk(ix,ibody)
        smag = smag + snap(ix,ibody)*snap(ix,ibody)
        cmag = cmag + crackle(ix,ibody)*crackle(ix,ibody)
     enddo

     amag = sqrt(amag)
     jmag = sqrt(jmag)
     smag = sqrt(smag)
     cmag = sqrt(cmag)

     if(cmag/jmag + smag/jmag > 1.0e-40) then
        dt_i = sqrt(tolerance*(amag/jmag + jmag/smag)/(cmag/smag + smag/jmag))
     else
        dt_i = 1.0e30
     endif

     if (dt_i < deltat) deltat = dt_i

  enddo


end subroutine timestep
