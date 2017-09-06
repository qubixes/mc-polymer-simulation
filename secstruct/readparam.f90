subroutine readparam(dir)
!read parameters

implicit none

include 'system.inc'

character(len=1000) ::dir
character(len=1000) ::file
file=trim(dir)//'/input.dat'

open(10,file=file)
read(10,*) Nchain !size of the ring
read(10,*) Lforbid !minimal length between two linked monomers
read(10,*) Emax !high energy dedicated to "avoided" links
read(10,*) Emin !low energy dedicated to "forced" links
read(10,*) E0 !minimum of the potential
read(10,*) rc !cut-off distance
read(10,*) h !exponent of the potential
read(10,*) Eb !energy of an open region
close(10)

return

end subroutine