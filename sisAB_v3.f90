!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  				SIS-SIR dynamics
!				by Isaura Oliver
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  program sisab
  implicit none
  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
    interface ra
      function RAND0(DUMMY)
        implicit double precision (A-H,O-Z), integer*8 (I-N)
     end function RAND0
    end interface ra
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!
  integer*8 y1,y2
  integer*8 e,n
  integer*8 noinfecA,norecA,noinfecB,norecB
  
  integer(kind=8) :: niA = 0
  integer(kind=8) :: naA = 0
  integer(kind=8) :: niB = 0
  integer(kind=8) :: naB = 0
  integer(kind=8) :: niAiB = 0
  
  integer*8 i,j,k,l,m,o,p,q
  integer*8 ISEED1,ISEED2,I2
  integer*8 re ! dummy variable
  integer*8 allocat
  double precision deltaA,lambaA,lamdel
  double precision lamdelA, lamdelB
  double precision deltaB,lambaB
  double precision t,dt,tf,tfis
  double precision rhoA,rhotA,rhostA,numrhoA
  double precision rhot2A,rhost2A,desvtA
  double precision rhoB,rhotB,rhostB,numrhoB
  double precision rhot2B,rhost2B,desvtB
  double precision precA,pinfA,precB,pinfB
  double precision x1,x2,x3,x4
  
  integer*8, allocatable :: n1(:),n2(:)
  integer*8, allocatable :: ve(:),p1(:),p2(:),degr(:)
  integer*8, allocatable :: viA(:),vaA(:),viB(:),vaB(:)
  integer*8, allocatable :: viAiB(:)
    
   character*32 iseed,tefe,lamba,delta,lambb,deltb
  
   character*256 fname
   character*256 flname
   
   integer date_time(8)
   character*10 date_c(3)

! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   interface din
     
     subroutine read_net(fname,e,n)
       integer*8 e,n
       character(*) fname
     end subroutine read_net
     
     subroutine read_nods(fname,e,nod1,nod2)
       integer*8 e
       integer*8, dimension(:) :: nod1(:),nod2(:)
       character(*) fname
     end subroutine read_nods
     
     subroutine codi(v,p1,p2,deg,n,e,fname)
       integer*8 n,e
       character(*) fname
       integer*8, dimension(:) :: v(:),p1(:),p2(:),deg(:)
     end subroutine codi
     
     subroutine remove(n,y,v,re)
       integer*8 n, y,re
       integer*8, dimension(:) :: v(:)
     end subroutine remove
     
     
     subroutine add(n,a,v)
       integer*8 n,a
       integer*8, dimension(:) :: v(:)
     end subroutine add
     
      subroutine infectnode(node,ni,na,vi,va,p1,p2,ve)
      integer*8 node,ni,na
      integer*8, dimension(:) :: vi(:),va(:),p1(:),p2(:),ve(:)
      end subroutine infectnode
      
      subroutine recunode(posc,norec,ni,na,vi,va,p1,p2,ve)
      integer*8 posc,norec,ni,na
      integer*8, dimension(:) :: vi(:),va(:),p1(:),p2(:),ve(:)
      end subroutine recunode
      
      subroutine infrec(prec,ni,vi,na,va,p1,p2,ve)
      double precision prec
      integer*8 ni,na
      integer*8, dimension(:) :: vi(:),va(:),p1(:),p2(:),ve(:)
      end subroutine infrec
     
   end interface din
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  COMMON /RSEED/ ISEED1, ISEED2

  ISEED1 = 123457
  ISEED2 = 3727
  call getarg(1, flname)
  call getarg(2, iseed)
  call getarg(3, tefe)
  call getarg(4, lamba)
  call getarg(5, delta)
  call getarg(6, lambb)
  call getarg(7, deltb)  
  
  iseed1 = iseed1 + 2*i2
  read(flname,*) fname  
  read(iseed,*) I2
  read(tefe,*) tf
  read(lamba,*) lambaA
  read(delta,*) deltaA 
  read(lambb,*) lambaB
  read(deltb,*) deltaB
   
  tf=dble(tf)
!    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  
! ############### log file ########################
call date_and_time(date_c(1), date_c(2), date_c(3), date_time)
open(unit= 43, file = 'sisab_v1.log' , status='unknown', position='append' )
write(43,*) 'd:', date_c(1),  ' h:', date_c(2), ' input: file =' ,&
  & trim(flname)//'.dat','dA'//trim(delta)//'dB'//trim(deltb)//'LA' &
  & //trim(lamba)//'LB'//trim(lambb),'  output:  rho_t.dat', ' ./dinam_v3     '
  
close(43)
! #################################################
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  write(*,*) 'Inicial nodes to infect'
  call read_net(fname,e,n)  
  call read_nods(fname,e,n1,n2)
  call codi(ve,p1,p2,degr,n,e,fname)
!
  open(4,file='numsaleatoris',status='unknown')
  do i=1,10
   write(4,*) rand0(3.0d0)
  end do
  close(4)
  
  allocate(n1(1:e),stat=allocat)
  allocate(n2(1:e),stat=allocat)
  allocate(degr(1:n), stat=allocat)
  allocate(ve(1:2*e),stat=allocat)
  allocate(p1(1:n),stat=allocat)
  allocate(p2(1:n),stat=allocat)      
  
  allocate(viA(1:n),stat=allocat)
  allocate(vaA(1:e),stat=allocat)
  allocate(viB(1:n),stat=allocat)
  allocate(vaB(1:e),stat=allocat)
  allocate(viAiB(1:n),stat=allocat)
!
!
   open(34,file='dades/rhost'//trim(tefe)//'lA'//trim(lamba)//'lB'//trim(lambb)//'dA'&
   & //trim(delta)//'dB'//trim(deltb)//'.dat', status='unknown')
   open(56,file='dades/prev_rhostlA'//trim(lamba)//'lB'//trim(lambb)//'dA'&
   & //trim(delta)//'dB'//trim(deltb)//'.dat', status='unknown')  
   open(3, file='dades/nodesambxouA.dat',status='unknown')
   open(7, file='dades/nodesinfectsA.dat',status='unknown')
   open(8,file='dades/nodesambxouB.dat',status='unknown')
   open(9, file='dades/nodesinfectsB.dat',status='unknown')
   open(10,file='dades/noinfAB_2.dat',status='unknown')
! 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!
!
! Infect all nodes  
!    
    do i=1,n
      call add(niA,i,viA)
    end do   
    do j=1,n
      call add(niB,j,viB)
    end do  
!    write(*,*) n, ni
!      
! Physical time
!
!
  tfis=0
!
! time (Poison sampled)
!
  rhotA=0.0d0
  rhotB=0.0d0
  rhot2A=0.0d0
  rhot2B=0.0d0
  numrhoA=0
  numrhoB=0
  t=0.0d0
!
   do while( ((t .lt. tf) .and. (niA.gt.0)) .or. ((t .lt. tf) .and. (niB.gt.0)))
     do i=1,n
! 
! Probability of infection or recuperation
!  
     lamdel=deltaA*niA+lambaA*naA+deltaB*niB+lambaB*naB
     precA=(deltaA*niA)/lamdel
     pinfA=(lambaA*naA)/lamdel
     precB=(deltaB*niB)/lamdel
     pinfB=(lambaB*naB)/lamdel   
! 
     x1=RAND0(1.0d0)

       if (x1-pinfA.lt. 0) then
         y2=ceiling(RAND0(6.0d0)*dble(naA))
         noinfecA=ve(vaA(y2)) !no infected with A
         call infectnode(noinfecA,niA,naA,viA,vaA,p1,p2,ve)
       else if ((x1-pinfA-pinfB) .lt. 0) then
         y2=ceiling(RAND0(5.0d0)*dble(naB))
         noinfecB=ve(vaB(y2)) !no infected with B
         call infectnode(noinfecB,niB,naB,viB,vaB,p1,p2,ve)
       else if ((x1-pinfA-pinfB-precA) .lt. 0) then
        y1=ceiling(RAND0(2.0d0)*dble(niA))    
        norecA=viA(y1) !infected with A  
        call recunode(y1,norecA,niA,naA,viA,vaA,p1,p2,ve)                    
       else 
        y1=ceiling(RAND0(8.0d0)*dble(niB))     
        norecB=viB(y1)  !infected with B 
        call recunode(y1,norecB,niB,naB,viB,vaB,p1,p2,ve) 
       end if 

    end do 
      
      write(3,*) '  '
     write(3,*) '  '
     write(3,*) '#t=',t    
       do k=1,naA
        do m=1,niA
         if( (p1(viA(m)).le.vaA(k)).and.(p2(viA(m)).ge.vaA(k)) )then
            write(3,*) viA(m),ve(vaA(k))
         end if
        end do
       end do
     if( naA .eq. 0) then
       write(3,*) 0, 0
     end if
     write(7,*) '  '
     write(7,*) '  '
     write(7,*) '#t=',t    
       do o=1,niA
         write(7,*) viA(o)
       end do   

     write(8,*) '  '
     write(8,*) '  '
     write(8,*) '#t=',t
      do l=1,naB
       do p=1,niB
          if( (p1(viB(p)) .le. vaB(l)) .and. (p2(viB(p)) .ge. vaB(l)) )then
            write(8,*) viB(p),ve(vaB(l))
          end if
        end do 
      end do
    if(naB.eq.0) then
      write(8,*) 0, 0
    end if     
     write(9,*) '  '
     write(9,*) '  '
     write(9,*) '#t=',t    
       do q=1,niB
         write(9,*) viB(q)
       end do
!
       if (noinfecA == noinfecB) then
          write(10,*) noinfecA, noinfecB 
          call add(niAiB,noinfecA,viAiB)
       end if    
!            
       dt=-dlog(RAND0(4.0d0))   
       t=t+dt           
       rhoA=dble(niA)/dble(n)
       rhoB=dble(niB)/dble(n)    
       tfis=tfis+dt
       write(34,*) t,rhoA,rhoB,tfis  
        if (t .gt. (tf/2.0d0))then 
         rhotA=rhotA+rhoA
         rhotB=rhotB+rhoB
         numrhoA=numrhoA+1
         numrhoB=numrhoB+1
         rhot2A=rhot2A+rhoA*rhoA
         rhot2B=rhot2B+rhoB*rhoB
        end if 
      end do
!   rho stationary mean
     rhostA=rhotA/numrhoA
     rhost2A=(rhot2A/numrhoA)-rhostA*rhostA !variance A
     rhostB=rhotB/numrhoB
     rhost2B=(rhot2B/numrhoB)-rhostB*rhostB !variance B
     desvtA=dsqrt(rhost2A/numrhoA)
     desvtB=dsqrt(rhost2B/numrhoB)
      if (rhostA .gt. 0 ) then
        write(56,*) lambaA, rhostA, desvtA
      end if
      if (rhostB .gt. 0) then
        write(56,*) lambaB,rhostB,desvtB
      end if  
!  
!
  close(34)
  close(3) 
  close(7)
  close(8)
  close(9) 
  close(56)  
  
  stop
  end program sisab
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   add to vi and va
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine recunode(posc,norec,ni,na,vi,va,p1,p2,ve)
    implicit none
    
    integer*8 j,k,l,re,m
    integer*8 posc,norec,ni,na
    
    integer*8, allocatable :: vi(:),va(:),p1(:),p2(:),ve(:)
    
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     interface
       subroutine add(n,a,v)
       integer*8 n,a
       integer*8, dimension(:) :: v
       end subroutine add
       subroutine remove(n,y,v,re)
       integer*8 n,y,re
       integer*8, dimension(:) :: v
       end subroutine remove
     end interface 
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
       
     call remove(ni,posc,vi,re)
 
   do j=1,na
      if ((p1(norec) .le. va(j)) .and. (p2(norec) .ge. va(j) ))then
       call remove(na,j,va,re)
      end if  
   end do 

  ! Activation. List the recuperated nodes.
    
    ! 1. Search if the nodes are infected (they don't have norec in their neighbours) 
    ! 2. Put the edges towards the recuperated node to the vector of activated.
    
     do k=1,ni
      do l=p1(vi(k)),p2(vi(k))
        if (ve(l) == norec) then 
          call add(na,l,va)       
        end if
      end do
     end do
     
    
    return
    end subroutine recunode

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   add to vi and va
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine infectnode(node,ni,na,vi,va,p1,p2,ve)
    implicit none
    integer*8 node,ni,na,re
    integer*8 j,l,i
    
    integer*8, allocatable :: vi(:),va(:),p1(:),p2(:),ve(:)
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     interface
       subroutine add(n,a,v)
       integer*8 n,a
       integer*8, dimension(:) :: v
       end subroutine add
       subroutine remove(n,y,v,re)
       integer*8 n,y,re
       integer*8, dimension(:) :: v
       end subroutine remove
     end interface 
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
     

!   Put the neighbours to the list of activated edges.    
!   Delete the activated edges that have 'node' as destination.
 
    i = 1
    do while(i.le.na)
       if(node .eq. ve(va(i))) then
         call remove(na, i, va,re)
       else
         i = i +1
       end if
    end do
     
    do j=p1(node),p2(node)
       l=1
       do while ((vi(l) /= ve(j)) .and. (l .le. ni))
        l=l+1
       end do
       if (l > ni) then
       call add(na,j,va)
        end if 
    end do
  

    call add(ni,node,vi) !put the node on list of infected
 
    
    return
    end subroutine infectnode

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    remove
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine remove(n,y,v,re)
   implicit none
   integer*8 n,re
   integer*8 y
   
   integer*8, allocatable :: v(:)
   
   
   v(y)=v(n)
   
   n=n-1
   
   return
   end subroutine remove
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine add(n,a,v)
   implicit none
   integer*8 n,a
   
   integer*8, allocatable :: v(:)
   
   n=n+1
   v(n)=a
 
   end subroutine add
   
! necessito la xarxa neta! llegeixo el fitxer de la xarxa neta:

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    read the net
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    subroutine read_net(fname,e,n)
    subroutine read_net(fname,e,n)
    implicit none
    integer*8 e,n,istat,i,j!,allocat
    
    character(256) fname

     
    open(unit=2,file=trim(fname)//'.dat',status='old')
!     open(unit=2,file='file_clean.dat',status='old')
  
!   read file data
       e=0
       n=0
       istat=0
       do while (istat.eq.0)
       read(2,*,iostat=istat) i,j
         if(istat.eq.0) then
          e=e+1
          if (i.gt.n) n=i
          if (j.gt. n) n=j
         end if
       end do
      close(2)
      

    return
    end subroutine read_net

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    read nodes
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine read_nods(fname,e,nod1,nod2)
    implicit none
    integer*8 e,i,j,allocat,istat
    
    integer*8, allocatable :: nod1(:),nod2(:)
    
    character(256) fname
    open(unit=2,file=trim(fname)//'.dat',status='old')
!     open(unit=2,file='file_clean.dat',status='old')          
     
    allocate(nod1(1:e),stat=allocat) 
    allocate(nod2(1:e),stat=allocat) 
    
!   read file data
       e=0
       istat=0
       do while (istat.eq.0)
         read(2,*,iostat=istat) i,j
         if(istat.eq.0) then
          e=e+1
          nod1(e)=i
          nod2(e)=j
         end if
       end do
       
      close(2)
      
 
    return
    end subroutine read_nods
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
      END FUNCTION RAND0

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      codificacio xarxa
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine codi(v,p1,p2,deg,n,e,fname)
  implicit none
!  
 integer*8 :: e,n,i,k!,l!,istat,j,
 integer*8 :: ai,aj
 integer*8 :: allocstat
 integer*8, allocatable :: n1(:),n2(:)
 integer*8, allocatable :: v(:), deg(:),p1(:), p2(:)
 character(256) fname
 
 
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       interface
         subroutine read_nods(fname,e,nod1,nod2)
          integer*8 e
          character(*) fname
          integer*8, dimension(:) :: nod1(:),nod2(:)
         end subroutine
       end interface
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! gives the edges and nodes of the network  
   call read_nods(fname,e,n1,n2)
 
  !ccccccccccccccccccccccccccccccccccccccccc 
  ! dimension of vectors
 
  allocate(n1(1:e),stat=allocstat)
  allocate(n2(1:e),stat=allocstat)
  allocate(v(1:2*e),stat=allocstat)
  allocate(deg(1:n),stat=allocstat)
  allocate(p1(1:n),stat=allocstat)
  allocate(p2(1:n),stat=allocstat)
  
    
   deg=0
   do i=1,e
    ai=n1(i)
    aj=n2(i)
    deg(ai)=deg(ai)+1
    deg(aj)=deg(aj)+1 
   end do
   
  
 ! gives the pinicial and pfinal 
  p1(1)=1
  do k=2,n
     p1(k)=p1(k-1)+deg(k-1)
     p2(k)=p1(k)-1
  end do
 !
 !   
  p2(1)=0
  do i=1,e
    ai=n1(i)
    aj=n2(i)
    p2(ai)=p2(ai)+1
    p2(aj)=p2(aj)+1
 !   ! vector v give the final nodes
    v(p2(ai))=aj
    v(p2(aj))=ai
    end do 

  
  return
  end subroutine codi

