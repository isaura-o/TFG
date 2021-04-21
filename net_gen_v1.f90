!
! 		Generates Erdos-Renyi networks (Random Network Model)
!				by Isaura Oliver
!
 program net_gen
 implicit none

  
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   interface 
     function RAND0(DUMMY)
        implicit double precision (a-h,o-z), integer*8 (i-n)
     end function RAND0
   end interface  
 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 integer*8 n!,y
 integer*8 allocat,i,j, i2
 integer*8 ISEED1, ISEED2
 double precision x1
 double precision pr,averdeg,ed
 character*32 iseed,nod,averdegree
 
 integer*8, allocatable :: nod1(:), nod2(:)
 
        integer date_time(8)
       character*10 date_c(3)
  
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  COMMON /RSEED/ ISEED1, ISEED2
ISEED1 = 123457
ISEED2 = 3727


 
 ! EXTERNAL RAND
  call getarg(1, nod)
  call getarg(2, iseed)

  read(iseed,*) I2
  iseed1 = iseed1 + 2*i2
  read(nod,*) n


 write(*,*) 'average degree='
 read(*,*) averdegree
 read(averdegree,*) averdeg

! ############### log file ########################
call date_and_time(date_c(1), date_c(2), date_c(3), date_time)
open(unit= 43, file = 'net_gen.log' , status='unknown', position='append' )
write(43,*) 'd:', date_c(1),  ' h:', date_c(2), ' input: <k> =' ,&
  & trim(averdegree), '  N = ', trim(nod), '  id = ', trim(iseed),   &
  & '  output:  random_netN'//trim(nod)//'k'//trim(averdegree)//'.dat', &
  & '       ./net_gen  ', trim(nod), ' ', trim(iseed)
close(43)
! #################################################
!
! 
!
  allocate(nod1(1:n),stat=allocat)
  allocate(nod2(1:n),stat=allocat)
  
  pr=dble(averdeg)/dble(n)
  write(*,*) pr
 ! count the connected pairs:  
  ed=pr*0.5*n*(n-1)
  write(*,*) '<E>= ',ed
  
   
! make all possible pairs 
  open(unit= 23,file='random_netN'//trim(nod)//'k'//trim(averdegree)//'.dat', status='unknown')
  do i=1,n
    do j=i+1,n
      x1=rand0(1.0d0)
      if ( x1 .lt. pr ) then
         nod1(i)=i
         nod2(j)=j
       write(23,*) nod1(i), nod2(j)
      end if
    end do
  end do 
  close(23)
   
  
 stop
 end 
          

 
!  *********************************************************************
!                         FUNCTION RAND (Random number generator)
!  *********************************************************************
      FUNCTION RAND0()
!
!  This is an adapted version of subroutine RANECU written by F. James
!  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
!  give a single random number at each call.
!
!  The 'seeds' ISEED1 and ISEED2 must be initialized in the main program
!  and transferred through the named common block /RSEED/.
!
!  Some compilers incorporate an intrinsic random number generator with
!  the same name (but with different argument lists). To avoid conflict,
!  it is advisable to declare RAND as an external function in all sub-
!  programs that call it.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*8 (I-N)
      PARAMETER (USCALE=1.0D0/2.147483563D9)
      COMMON/RSEED/ISEED1,ISEED2
!
      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
!
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
!
      IZ=ISEED1-ISEED2
      IF(IZ.LT.1) IZ=IZ+2147483562
      RAND0=IZ*USCALE
!
      RETURN
      END 


