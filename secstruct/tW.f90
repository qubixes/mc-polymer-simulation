recursive subroutine tW(i,j,W,V,idW,idV,E)

implicit none

include 'system.inc'

real,dimension(Nchain+2,Nchain+2)	::W,V,E
integer,dimension(Nchain+2,Nchain+2)	::idW
integer,dimension(Nchain+2,Nchain+2,2)	::idV
integer	::i,j,im

im=idW(i,j)
if (im.ne.0) then

		if (im.eq.-1) then
			if ((i.ne.1).or.(j.ne.(Nchain+2))) then
                write(10,*) i-1,j-1
                f=f+2
            end if
			call tV(i,j,W,V,idW,idV,E)
		elseif (im.eq.-2) then
			call tW(i+1,j,W,V,idW,idV,E)
		elseif (im.eq.-3) then
			call tW(i,j-1,W,V,idW,idV,E)
		else
			call tW(i,im,W,V,idW,idV,E)
			call tW(im+1,j,W,V,idW,idV,E)
		end if
end if

return

end subroutine