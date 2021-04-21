!
!		Scale-free networks using configuration model
!			by Isaura Oliver
! 
!
  program configurational
  implicit none
  
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   interface 
     function RAND0(DUMMY)
        implicit double precision (A-H,O-Z), integer*8 (I-N)
     end function RAND0
   end interface  
 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  integer*8 allocat,i,j,k,l,o,p,r,x,s!,q,m
  integer*8 ISEED1, ISEED2, i2
  integer*8 n,stb,de,pos1,pos2
  double precision x1
  double precision gamma,invga,averdeg
  character*32 iseed,nod,gamm
  
  integer*8, allocatable :: degu(:),stub(:)
  integer*8, allocatable :: vr(:),p1(:),p2(:)
  
  integer date_time(8)
  character*10 date_c(3)
 !
 !
 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  COMMON /RSEED/ ISEED1, ISEED2
  
  ISEED1 = 123457
  ISEED2 = 3727 
  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  
 ! write(*,*) 'number of nodes='
 ! read(*,*) n
  call getarg(1, nod)
  call getarg(2, iseed)

  read(iseed,*) I2
  iseed1 = iseed1 + 2*i2
  read(nod,*) n
  
  
  write(*,*) 'degree exponent='
  read(*,*) gamm
  read(gamm,*) gamma
  
  invga=1/(1-gamma) ! de la normalitzacio
 ! write(*,*) invga
 
 ! ############### log file ########################
call date_and_time(date_c(1), date_c(2), date_c(3), date_time)
open(unit= 43, file = 'confg_net.log' , status='unknown', position='append' )
write(43,*) 'd:', date_c(1),  ' h:', date_c(2), ' input: <k> =' ,&
  & trim(gamm), '  N = ', trim(nod), '  id = ', trim(iseed),   &
  & '  output:  confg_netN'//trim(nod)//'gamma'//trim(gamm)//'.dat', &
  & '       ./conf_net_gen_v3  ', trim(nod), ' ', trim(iseed)
close(43)
! #################################################
 
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    allocate(degu(1:n),stat=allocat)
  !  
  ! Degree distribution
  !
  open(21,file='grauinici.dat')
   de=0
   do i=1,n   
     x1=rand0(1.0d0)
     degu(i)=floor(x1**(invga))
      write(21,*) x1, degu(i), i
      de=de+degu(i)
   end do
   close(21)
    
    averdeg=dble(de)/dble(n)
    write(*,*) 'average degree= ',averdeg
!    
!
!    
   stb=0
    do j=1,n
      stb=stb+degu(j)
    end do  
    write(*,*) 'stubs=', stb
!
!  substract a stub to node n
!
! write(*,*) 'degree node n= ',degu(n)
    if (mod(stb,2) .ne. 0) then
        write(*,*) 'you need a even number of stubs'
      stop
     end if
  !   
  ! create a vector (stub) who has all the vertex (need number even of stubs) (potes)
  !
 
    allocate(stub(1:stb),stat=allocat)
    open(23,file='potes.dat')
    
    o=0 
    do l=1,n 
     do k=1,degu(l)
       o=o+1
       stub(o)=l
     end do
    end do

        
   do p=1,o
     write(23,*) stub(p)
   end do
    close(23)

!
!  pointers
!
   allocate(p1(1:n),stat=allocat)
   allocate(p2(1:n),stat=allocat)
   

    p1(1)=1
    do s=2,n
      p1(s)=p1(s-1)+degu(s-1)
      p2(s)=p1(s)-1
    end do
!
!  
! now assign to any vertex any other one
!  

   allocate(vr(1:stb),stat=allocat)
   open(unit= 34,file='config_netN'//trim(nod)//'gamma'//trim(gamm)//'.dat', status='unknown')
   
   x=0
   p2(1)=0
   do while(o .gt. 0)
      pos1 = 1+int(rand0(2.0d0)*dble(o)) ! any
      pos2 = 1+int(rand0(3.0d0)*dble(o-1)) ! take second edge, any
      do while (stub(pos1) .eq. stub(pos2))
      !   write(*,*) '#equals'
         pos2 = 1+int(rand0(4.0d0)*dble(o)) ! take second edge 
        if ( pos2 .gt. o ) then
            pos2 = 1+int(rand0(3.0d0)*dble(o))
        end if 
      end do
      x=x+1 
      p2(stub(pos1))=p2(stub(pos1))+1
      p2(stub(pos2))=p2(stub(pos2))+1
      vr(p2(stub(pos1)))=stub(pos2)
      vr(p2(stub(pos2)))=stub(pos1)
 !     write(87,*) '#new position'
 !     write(87,*) o, stub(pos1), stub(pos2)
     write(34,*) stub(pos1),stub(pos2)
! 	  put the last ones on unoccupied position    
      stub(pos1) = stub(o)
      stub(pos2) = stub(o-1)                 
     o=o-2
   end do
   close(34)
   write(*,*) 'edges=',x


   open(unit=56,file='vec.dat',status='unknown')
   do r=1,2*x
    write(56,*) vr(r)
   end do
   close(56)

    
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

