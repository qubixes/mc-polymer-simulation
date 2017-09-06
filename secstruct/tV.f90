recursive subroutine tV(i,j,W,V,idW,idV,E)

implicit none

include 'system.inc'

real,dimension(Nchain+2,Nchain+2)	::W,V,E
integer,dimension(Nchain+2,Nchain+2)	::idW
integer,dimension(Nchain+2,Nchain+2,2)	::idV
integer	::i,j,im1,im2

im1=idV(i,j,1)
im2=idV(i,j,2)
if (im1.ne.0) then
	if (im1.eq.-1) then
		call tW(i+1,im2,W,V,idW,idV,E)
		call tW(im2+1,j-1,W,V,idW,idV,E)
	else
		write(10,*) im1-1,im2-1
        f=f+2
		call tV(im1,im2,W,V,idW,idV,E)
	end if
end if

return

end subroutine