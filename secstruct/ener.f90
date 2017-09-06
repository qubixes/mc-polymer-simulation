subroutine ener(pos,E)
!compute the potential between two monomers from positions

implicit none

include 'system.inc'

real,dimension(Nchain,3)	::pos
real,dimension(Nchain+2,Nchain+2)	::E
real	::rij
integer	::i,j

E=Emax
E(1,Nchain+2)=Emin
E(Nchain+2,1)=Emin

do i=2,Nchain+1
	do j=i+Lforbid,min(Nchain+1,Nchain+1+i-1-Lforbid)
		rij=sqrt((pos(i-1,1)-pos(j-1,1))**2+(pos(i-1,2)-pos(j-1,2))**2+(pos(i-1,3)-pos(j-1,3))**2)
		E(i,j)=E0*((rij/rc)**h-1.)
		E(j,i)=E(i,j)
	end do
	!if (i.ne.1) then
	!	E(i,Nchain+1)=E(i,1)
	!	E(Nchain+1,i)=E(1,i)
	!end if
end do
		

return

end subroutine