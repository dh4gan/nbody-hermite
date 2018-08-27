PROGRAM nbody_hermite
  ! Written 24/08/2018 by dh4gan
  ! Code implements a 4th order Hermite integration of a pure N-Body system

  use nbodydata

  implicit none

  ! Print code header

  !             Display header
  print*, " "
  print*, "-----------------------------------------------"
  print*, "     N-BODY INTEGRATOR (HERMITE 4TH ORDER) "
  print*, "     Created by D.Forgan, 24th August 2018   "
  print*, "-----------------------------------------------"
  print*, " "
  print*, "-----------------------------------------------"
  call getarg(1,paramfile)

  if(paramfile=='') then
     print*, 'Parameter file name not found from command line'
     paramfile = 'nbody_hermite.params'
     print*, 'Reverting to default'
  endif

  print*, " input parameters to be read from ",trim(paramfile)
  print*, "-----------------------------------------------"
  call sleep(1)


  ! Read in parameter file and setup bodies

  call initial

  ! Begin loop

  t = 0.0
  dt = 1.0e-3
  tdump = tsnap

  ! Output initial conditions
  snapshotcounter = 0
  call output(snapshotcounter)

  ! Begin integration
  do while(t<tend)

     call integrate(dt,pos,vel,newpos,newvel)

     pos = newpos
     vel = newvel

     t = t + dt

     if (t>tdump) then
        write(*,'(A,1P,2E12.3,A)') 't, dt=',t/twopi,dt/twopi, ' years'
        snapshotcounter = snapshotcounter + 1
        call output(snapshotcounter)
        tdump = tdump + tsnap
     endif

  enddo

  ! Once integration complete, close all files and exit

  call endrun

END PROGRAM nbody_hermite
