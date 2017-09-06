subroutine readpos(pos)
!read the positions

implicit none

include 'system.inc'

real,dimension(Nchain,3)	::pos
integer	::i

do i=1,Nchain
	read(*,*) pos(i,1),pos(i,2),pos(i,3)
end do

return

end subroutine