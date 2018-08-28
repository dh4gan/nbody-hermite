subroutine output(counter)
! Outputs data to file
! Currently writes each particle to separate file

use nbodydata

implicit none
integer, intent(in) :: counter
logical :: calc_snapcrackle
real, dimension(3,N) :: jerk,snap,crackle

!101 format (1P, 22E15.5)
102 format (1P,23E15.5)
103 format (1P, 7E15.5)

calc_snapcrackle = .false.
call calc_acceleration(pos,vel,acc,jerk,snap,crackle,calc_snapcrackle)
call calc_system_properties

! Either output timestep as a single snapshot
if(snapshots=='y') then
   write(fileno,filenumformat) counter

   outputfile = TRIM(outputprefix)//"."//TRIM(fileno)
   open(isnap,file=outputfile,form='formatted')
   write(isnap,*) t/twopi
   do ibody=1,N
write(isnap,102) mass(ibody), pos(:,ibody), vel(:,ibody), acc(:,ibody),&
    semimaj(ibody),ecc(ibody),inc(ibody),longascend(ibody),argper(ibody),trueanom(ibody), &
    ekin(ibody),epot(ibody),etot(ibody),angmom(:,ibody)
    call flush(isnap)
 enddo

   
else

   ! or output individual bodies to separate files

   do ibody=1,N
       
      write(ibody+ilog,102) t/twopi, mass(ibody),pos(:,ibody), vel(:,ibody), acc(:,ibody),&
semimaj(ibody),ecc(ibody),inc(ibody),longascend(ibody),argper(ibody),trueanom(ibody), &
ekin(ibody),epot(ibody),etot(ibody),angmom(:,ibody)
      call flush(ibody)
   enddo

endif

! Write log file containing global simulation data

write(ilog,103) t/twopi, dt/twopi, 100.0*maxerror/tolerance, system_energy, dE, system_ang, dL
call flush(ilog)

end subroutine output
