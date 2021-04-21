 !
 !	Cleans networks from repeated edges, symmetric edges, and loops
 ! 			by Isaura Oliver
 !
 program net_cleaner
 implicit none
   integer*8 :: i,j,k,r,s,p,l,t,ss,rr
   integer*8 :: e,n
   integer*8 :: allocat,istat
   integer*8, allocatable :: net1(:),net2(:)
   integer*8, allocatable :: nolop1(:),nolop2(:)
   integer*8, allocatable :: norep1(:),norep2(:)
   integer*8, allocatable :: nosim1(:),nosim2(:)
   character*90 :: fname,nod,gamm
 !
 !
  call getarg(1,fname) 
  read(fname,*) fname
  call getarg(2,nod)
  read(nod,*) nod
  call getarg(3,gamm)
  read(gamm,*) gamm
  
   
 open(1,file=trim(fname)//'.dat', status='old')
  
  ! read data files 
       e=0
       istat=0
       do while (istat.eq.0)
       read(1,*,iostat=istat) i,j
         if(istat.eq.0) then
          e=e+1
          if (i.gt.n) n=i
          if (j.gt. n) n=j
         end if
       end do
       write(*,*) 'edges initial= ', e
      close(1)
!----------------------------------------------
!	dimensions
      allocate(net1(1:e),stat=allocat)
      allocate(net2(1:e),stat=allocat)
      allocate(nolop1(1:e),stat=allocat)
      allocate(nolop2(1:e),stat=allocat)
      allocate(norep1(1:e),stat=allocat)
      allocate(norep2(1:e),stat=allocat)  !array is norep(maxdimension,number columns)
      allocate(nosim1(1:e),stat=allocat)
      allocate(nosim2(1:e),stat=allocat)

! open a second time to rename i and j

   open(1,file=trim(fname)//'.dat',status='old') 
       e=0
       istat=0
       do while (istat.eq.0)
         read(1,*,iostat=istat) i,j
         if(istat.eq.0) then
         e=e+1
          net1(e)=i
          net2(e)=j
         end if
       end do
      close(1)
 !---------------------------------------
 ! clean the loops:        
 
      l=0
      do r=1,e
         if(net2(r) /= net1(r)) then
           l=l+1
           nolop1(l)=net1(r)
           nolop2(l)=net2(r)
         end if
       end do
       write(*,*) 'edges no loop=', l
   
! clean repeated edges ij=ij   
   j=1
   do s=1,l
     k=1
       do while (((norep1(k)/= nolop1(s)).or.(norep2(k)/= nolop2(s))) .and. k < j)
       k=k+1  
       end do
        if (k == j) then
   !      print*, s,j
         norep1(k)=nolop1(s)
         norep2(k) = nolop2(s)
         j=j+1
        end if 
    end do 
   
    write(*,*) 'num of repeated edges (ij /= ij) = ', j-1 
 
! clean symmetric edges i,j = j,i
   p=1
   do ss=1,j-1
     rr=1
       do while (((nosim1(rr)/=norep2(ss)).or.(nosim2(rr)/=norep1(ss))) .and. rr < p)
       rr=rr+1  
       end do
        if (rr == p) then
         nosim1(rr) = norep1(ss)
         nosim2(rr) = norep2(ss)
         p=p+1
        end if 
    end do 
    
    write(*,*) 'total num. of edges= ', p-1
    
      open(2,file='clean_config_netN'//trim(nod)//'gamma'//trim(gamm)//'.dat', status='unknown')   
    do t=1,p-1
    write(2,*) nosim1(t),nosim2(t)
    end do
    close(2) 

     
 stop
 end 
 
