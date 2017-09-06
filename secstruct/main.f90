program main

implicit none

include 'system.inc'
real,allocatable	::W(:,:),V(:,:),E(:,:),pos(:,:),mV(:,:)
integer,allocatable	::idW(:,:),idV(:,:,:),idmV(:,:,:),treat(:),connec(:,:)
integer	::i,j,ip,jp,l,im,im1,im2,s,i0,nl
real	::m
real    ::rg2,rm(3)
integer :: argc
character(len=1000) :: dir
character(len=1000) :: file_out
character(len=100) :: polIdChar

integer :: polId

argc = command_argument_count()

if (argc.lt.2) then
	write(*,*) 'Error: need two command line arguments'
	call exit(0)
end if

call get_command_argument(1, dir)
call get_command_argument(2, polIdChar)
read(polIdChar, '(i4)') polId

! write(*,*) 'dir = ', trim(dir)
! write(*,*) 'polId = ', polId
file_out=trim(dir)//'/dJost.out.'// trim(polIdChar)
! write(*,*) "writing to file ", file_out
open(33, file=file_out)

write(33,*) 'read parameters'
call readparam(dir)
allocate(W(Nchain+2,Nchain+2))
allocate(V(Nchain+2,Nchain+2))
allocate(mV(Nchain+2,Nchain+2))
allocate(idmV(Nchain+2,Nchain+2,2))
allocate(E(Nchain+2,Nchain+2))
allocate(idW(Nchain+2,Nchain+2))
allocate(idV(Nchain+2,Nchain+2,2))
allocate(pos(Nchain,3))

! call exit(0)

write(33,*) 'read positions'
call readpos(pos)
write(33,*) 'compute Energy matrix'
call ener(pos,E)

!initialization for j-i=Lforbid
write(33,*) 'initialization'
V=Emax
mV=Emax
idmV=0
W=0.
idW=0
idV=0
do i=1,Nchain+2-Lforbid
	V(i,i+Lforbid)=E(i,i+Lforbid)+Eb
    mV(i,i+Lforbid)=V(i,i+Lforbid)
    idmV(i,i+Lforbid,1)=i
    idmV(i,i+Lforbid,2)=i+Lforbid
	W(i,i+Lforbid)=min(0.,V(i,i+Lforbid)) !XXXXXX 0. ou Eb
	if (W(i,i+Lforbid).eq.0.) then
		idW(i,i+Lforbid)=0
	else
		idW(i,i+Lforbid)=-1
	end if
end do

!compute the sub-chain
write(33,*) 'compute V and W'
do l=Lforbid+1,Nchain+1
    if (mod(l,50).eq.0)    write(33,*) 'l',l
	do i=1,Nchain+2-l
		j=i+l
		!naive algorithm: could be improved
		m=Eb
        im1=0
        im2=0
        if (V(i+1,j-1).lt.m) then
            m=V(i+1,j-1)
            im1=i+1
            im2=j-1
        end if

        if ((mV(i+1,j-1)+Eb).lt.m) then
            m=mV(i+1,j-1)+Eb
            im1=idmV(i+1,j-1,1)
            im2=idmV(i+1,j-1,2)
        end if

		do ip=Lforbid+i+1,j-2-Lforbid
			if ((W(i+1,ip)+W(ip+1,j-1))+Eb.lt.m) then
				m=W(i+1,ip)+W(ip+1,j-1)+Eb
				im1=-1
				im2=ip
			end if
		end do
		V(i,i+l)=E(i,i+l)+m
		idV(i,i+l,1)=im1
		idV(i,i+l,2)=im2
        if (mV(i+1,j).lt.mV(i,j-1)) then
            if (V(i,j).lt.mV(i+1,j)) then
                mV(i,j)=V(i,j)
                idmV(i,j,1)=i
                idmV(i,j,2)=j
            else
                mV(i,j)=mV(i+1,j)
                idmV(i,j,:)=idmV(i+1,j,:)
            end if
        else
            if (V(i,j).lt.mV(i,j-1)) then
                mV(i,j)=V(i,j)
                idmV(i,j,1)=i
                idmV(i,j,2)=j
            else
                mV(i,j)=mV(i,j-1)
                idmV(i,j,:)=idmV(i,j-1,:)
            end if
        end if


		m=V(i,i+l)
		im=-1
		if (W(i+1,i+l).lt.m) then
			m=W(i+1,i+l)
			im=-2
		end if

        if (W(i,i+l-1).lt.m) then
            m=W(i,i+l-1)
            im=-3
        end if

		do ip=i+Lforbid,j-1-Lforbid
			if ((W(i,ip)+W(ip+1,j)).lt.m) then
				m=W(i,ip)+W(ip+1,j)
				im=ip
			end if
		end do
		W(i,i+l)=m
		idW(i,i+l)=im
	end do
end do
write(33,*) 'Traceback the MFE'
file_out=trim(dir)//'/struct.out.'// trim(polIdChar)
open(10,file=file_out)
f=0
call tW(1,Nchain+2,W,V,idW,idV,E)
close(10)

allocate(treat(Nchain))
allocate(connec(Nchain,2))
treat=0
connec=0
file_out=trim(dir)//'/struct.out.'// trim(polIdChar)
open(10,file=file_out)
do i=1,f/2
    read(10,*) im,ip
    connec(ip,1)=1
    connec(ip,2)=im
    connec(im,1)=1
    connec(im,2)=ip
end do
close(10)

treat=0
nstem=0
file_out=trim(dir)//'/stem.out.'// trim(polIdChar)
open(10,file=file_out)
nstem=0
m=1
l=0

200 i=l+1
do while ((i.le.Nchain).and. &
    ((treat(i).eq.1).or.(connec(i,1).eq.0).or.((connec(i,1).eq.1).and.(treat(connec(i,2)).eq.1)).or.(connec(i,2).lt.i)))
    i=i+1
end do
l=i
m=m+1
if (m.gt.Nchain) stop

if (i.le.Nchain) then
    nstem=nstem+1
    s=0
    ip=i+1
    j=connec(i,2)
    es=E(i+1,j+1)
    im=j-1
    if (ip.gt.Nchain) ip=1
    if (im.lt.1) im=Nchain
    do while ((connec(ip,1).eq.1).and.(connec(ip,2).eq.im))
        treat(i)=1
        treat(im)=1
        s=s+1
        i=ip
        ip=i+1
        j=connec(i,2)
        es=es+E(i+1,j+1)
        im=j-1
        if (ip.gt.Nchain) ip=1
        if (im.lt.1) im=Nchain
    end do
    write(10,*) (s+1),es
    GOTO 200
end if

file_out=trim(dir)//'/loop.out.'// trim(polIdChar)
open(10,file=file_out)
nl=0
100 i=1
do while ((i.le.Nchain).and.(treat(i).eq.1))
    i=i+1
end do
if (i.le.Nchain) then
    s=0
    rg2=pos(i,1)**2+pos(i,2)**2+pos(i,3)**2
    rm=pos(i,:)
    i0=i
    l=0
    ip=i+1
    s=s+1
    treat(i)=1
    if (ip.gt.Nchain) ip=1
    i=ip
    do while (i.ne.i0)
        if ((connec(i,1).eq.1).and.(l.eq.0)) then
            rg2=rg2+pos(i,1)**2+pos(i,2)**2+pos(i,3)**2
            rm=rm+pos(i,:)
            ip=connec(i,2)
            s=s+1
            l=1
        else
            rg2=rg2+pos(i,1)**2+pos(i,2)**2+pos(i,3)**2
            rm=rm+pos(i,:)
            ip=i+1
            s=s+1
            treat(i)=1
            l=0
        end if
        if (ip.gt.Nchain) ip=1
        i=ip
    end do
    rm=rm/real(s)
    rg2=rg2/real(s)-(rm(1)**2+rm(2)**2+rm(3)**2)
    write(10,*) s,rg2
    nl=nl+1
    GOTO 100
end if
close(10)



write(33,*) 'Minimum energy=',W(1,Nchain+2)-Emin
write(33,*) 'Fraction of paired monomers=',real(f)/real(Nchain)
write(33,*) 'Number of stems=',nstem
write(33,*) 'Number of loops=',nl

close(33)
!open(10,file='E.out')
!open(11,file='V.out')
!open(12,file='W.out')
!do i=1,Nchain+2
!	write(10,*) E(i,:)
!	write(11,*) V(i,:)
!	write(12,*) W(i,:)
!end do
!close(10)
!close(11)
!close(12)

deallocate(W)
deallocate(V)
deallocate(E)
deallocate(pos)
deallocate(idV)
deallocate(idW)
deallocate(mV)
deallocate(idmV)


end program
