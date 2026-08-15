program cyclictdma

implicit none

integer i,im1,ip1
integer,parameter :: NMAX = 100
real,dimension(NMAX) :: AD,BD,CD,DD

Do i = 1,NMAX-1
AD(i) = 0.2
BD(i) = 1.0
CD(i) = 0.5
im1 = i-1
ip1 = i+1
if(im1.lt.1) im1 = NMAX-1
if(ip1.gt.NMAX-1) ip1 = 1
DD(i) = CD(i)*im1 + BD(i)*i + AD(i)*ip1
Enddo

call TDMAP(1,NMAX-1,CD,BD,AD,DD,NMAX)

Do i = 1,NMAX-1
print*, DD(i)
Enddo

end

subroutine TDMAP(ji,jf,ap,ac,am,fi,NMAXL)

   ! TDMA for periodic tridiagonal systems. 
   ! is NMALXL the NI? 
   

!ap - super diagonal
!ac - diagonal
!am - sub diagonal

implicit none

! -----------------------
!  Input/Output variables
! -----------------------
   integer:: ji, jf, NMAXL
   real:: ap(NMAXL), ac(NMAXL), am(NMAXL), fi(NMAXL)

! -------------------
!  Internal variables
! -------------------
   integer:: i, j, ja, jj
   real:: fnn, pp
   real:: qq(NMAXL), ss(NMAXL), fei(NMAXL)

   ja=ji+1
   jj=ji+jf

   qq(ji)=-ap(ji)/ac(ji)
   ss(ji)=-am(ji)/ac(ji)
   fnn=fi(jf)
   fi(ji)=fi(ji)/ac(ji)

!  forward elimination sweep
!----------------------------
   do j=ja,jf
      pp=1.0d0/(ac(j)+am(j)*qq(j-1))
          qq(j)=-ap(j)*pp
          ss(j)=-am(j)*ss(j-1)*pp
          fi(j)=(fi(j)-am(j)*fi(j-1))*pp
   enddo
   
!  backward pass
!----------------

   ss(jf)=1.0d0
   fei(jf)=0.0d0
 
   do i=ja,jf
      j=jj-i
          ss(j)=ss(j)+qq(j)*ss(j+1)
          fei(j)=fi(j)+qq(j)*fei(j+1)
   enddo
   
   fi(jf)=(fnn-ap(jf)*fei(ji)-am(jf)*fei(jf-1))/    &
   &      (ap(jf)*ss(ji)+am(jf)*ss(jf-1)+ac(jf)) 
   
!  backward substitution
!------------------------
   do i=ja,jf
      j=jj-i
          fi(j)=fi(jf)*ss(j)+fei(j)
   enddo

end subroutine TDMAP