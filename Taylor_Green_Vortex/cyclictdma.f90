program cyclictdma

implicit none

integer i,im1,ip1,j
integer,parameter :: NMAX = 100
real,dimension(NMAX) :: AD,BD,CD,DD
real:: qq(NMAX), ss(NMAX), fei(NMAX), fi(NMAX), pp(NMAX)
real fnn

!CD - super diagonal
!BD - diagonal
!AD - sub diagonal
!DD - RHS (overwritten by solution)

Do i = 1,NMAX-1
AD(i) = 0.2
BD(i) = 1.0
CD(i) = 0.5
im1 = i-1
ip1 = i+1
if(im1.lt.1) im1 = NMAX-1
if(ip1.gt.NMAX-1) ip1 = 1
DD(i) = AD(i)*im1 + BD(i)*i + CD(i)*ip1
Enddo

   qq(1)=-CD(1)/BD(1)
   ss(1)=-AD(1)/BD(1)
   
!  Appropriate coefficients
!----------------------------
   do j=2,NMAX-1
    pp(j)=1.0d0/(BD(j)+AD(j)*qq(j-1))
    qq(j)=-CD(j)*pp(j)
    ss(j)=-AD(j)*ss(j-1)*pp(j)
   enddo   
    
	ss(NMAX-1)=1.0d0
   do i=2,NMAX-1
    j=NMAX-i
    ss(j)=ss(j)+qq(j)*ss(j+1)
   enddo	
   
!  RHS & Solution Values which change
!------------------------------------
   DD(1)=DD(1)/BD(1)
   fnn = DD(NMAX-1)

!  forward elimination sweep
!----------------------------
   do j=2,NMAX-1
    DD(j)=(DD(j)-AD(j)*DD(j-1))*pp(j)
   enddo
   
!  backward pass
!----------------
    fei(NMAX-1)=0.0d0
 
   do i=2,NMAX-1
    j=NMAX-i
    fei(j)=DD(j)+qq(j)*fei(j+1)
   enddo
   
   DD(NMAX-1)=(fnn-CD(NMAX-1)*fei(1)-AD(NMAX-1)*fei(NMAX-2))/    &
   &      (CD(NMAX-1)*ss(1)+AD(NMAX-1)*ss(NMAX-2)+BD(NMAX-1)) 
   
!  backward substitution
!------------------------
   do i=2,NMAX-1
    j=NMAX-i
    DD(j)=DD(NMAX-1)*ss(j)+fei(j)
   enddo

Do i = 1,NMAX-1
print*, DD(i)
Enddo

end