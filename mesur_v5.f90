!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	Statidistics: degree, average degree, average nearest neighbour degree, cluster coefficient
!			by Isaura Oliver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program mesuresxarxa
   implicit none
   
   integer*8  n,i,j,l,no,p,q,o,k,r,g,f,e,no1
   double precision nomitj,pckk,kn,nomitj1
   integer*8 allocat,s,u,T,w,x
   double precision, allocatable :: pk(:),pck(:),his(:)
   double precision, allocatable :: knn(:),knn11(:),clus(:),clus11(:)
   integer*8, allocatable :: ve(:), degr(:), pini(:), pfin(:)
   character*90 filname
   character*32 flname
   
   integer date_time(8)
   character*10 date_c(3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! 
     interface codf
         subroutine codi(v,p1,p2,deg,n,fname,e)
           integer*8 n,e
           integer*8, dimension(:) :: v(:),p1(:),p2(:),deg(:)
           character(*) fname
         end subroutine codi
         subroutine hiso(degre,n,his)
           integer*8 n
           integer*8, dimension(:) :: degre(:)
           double precision, dimension(:) :: his(:)
         end subroutine hiso
         subroutine ordena(n,pinic,pfina,degre)
           integer*8 n
           integer*8, dimension(:) :: degre(:),pinic(:),pfina(:)
         end subroutine ordena
          subroutine ordenareal(n,kuns,degre)
           integer*8 n
           integer*8, dimension(:) :: degre(:)
           double precision, dimension(:) :: kuns(:)
         end subroutine ordenareal
         subroutine cleaner(filname,nosim1,nosim2,edfin,n)
           integer*8 edfin,n
           integer*8, dimension(:) :: nosim1(:),nosim2(:)
           character(*) filname
         end subroutine cleaner
         subroutine read_net(fname,e,n)
           integer*8 e,n
           character(*) fname
         end subroutine read_net
      end interface codf
!      
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccc  
! 
!  
! To open files: set it on the term
  !ccccccccccccccccccccccccccccccccccccccccccccccc    
 ! 
    call getarg(1,flname) 
    read(flname,*) filname

 !  
 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! ############### log file ########################
call date_and_time(date_c(1), date_c(2), date_c(3), date_time)
open(unit= 43, file = 'mesur_v4.log' , status='unknown', position='append' )
write(43,*) 'd:', date_c(1),  ' h:', date_c(2), ' input: file =' ,&
  & trim(flname),  '  output:  kpkpck_'//trim(flname)//'.dat', &
  &  'kknn_'//trim(flname)//'.dat', 'clus_'//trim(flname)//'.dat',&
  &  ' ./mesur_v4     '
close(43)
! #################################################
!
!
   call read_net(filname,e,n)
   write(*,*) 'initial edges=', e
 
 call codi(ve,pini,pfin,degr,n,filname,e) 
 
 write(*,*) 'final edges=', e 
 write(*,*) 'max node =', n


 allocate(pk(1:n),stat=allocat)
 allocate(pck(1:n),stat=allocat)
 allocate(his(1:n),stat=allocat)
 allocate(knn(1:n),stat=allocat)
 allocate(knn11(1:n),stat=allocat)
 allocate(clus(1:n),stat=allocat)
 allocate(clus11(1:n),stat=allocat)
!
!
! Probabilities and outputs
!.....................................................
!
! Probability that a node have degree n
! 
!   
    call ordena(n,pini,pfin,degr)   
    call hiso(degr,n,his)
     do i=1,n
      pk(i)=his(i)/dble(n)
     end do  
   
   pckk=0
   do j=n,1,-1
    if(degr(j).ne.degr(j+1)) then
     pckk=pckk+pk(j)
     pck(j)=pckk
    end if
   end do
 
   open(unit=12,file='kpkpck_'//trim(flname)//'.dat',status='unknown')
   write(12,*) '#',' ','k',' ','P(k)',' ','Pc(k)'
   do k=1,n
    if (degr(k) .ne. degr(k+1)) then
     write(12,*) degr(k), pk(k),pck(k)
    end if 
   end do
  close(12)   
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
! Average degree of the network
!
   no=0
   no1=0
   do l=1,n
    no=no+degr(l)
     if (degr(l) .ne. 0) then
      no1=no1+degr(l)
     end if
   end do
    nomitj=dble(no)/dble(n)
    nomitj1=dble(no1)/dble(n)
    write(*,*) 'average degree=', nomitj

!
!    
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
! average nearest neighbours degree
!
!
 call codi(ve,pini,pfin,degr,n,filname,e)

   do q=1,n
   kn=0
    do p=pini(q),pfin(q)
       kn=kn+degr(ve(p))
     end do
      knn(q)=kn/dble(degr(q))      
   end do
!
!    
!  knn output file
   open(unit=23,file='kknn_'//trim(flname)//'.dat',status='unknown')
   write(23,*) '#',' ','k',' ','knn'
   call ordenareal(n,knn,degr)
   do r=1,n 
    knn11(1)=0
       do x=1,n
         if(degr(r).eq. degr(x)) then
           knn11(x)=knn11(x)+knn(r)
         end if 
       end do
        if (degr(r).ne. degr(r+1)) then 
         write(23,*) degr(r), knn11(r)/his(r)
        end if
    end do   
   close(23)
   
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! cluster coefficient C=2*(TRIANGLES)/K(K-1)
!
!
 call codi(ve,pini,pfin,degr,n,filname,e) 
 
   do o=1,n
     T=0
     do s=pini(o),pfin(o)
      do u=pini(ve(s)),pfin(ve(s)) 
        do w=pini(o),pfin(o)
           if (ve(u).eq.ve(w)) then
           T=T+1
           end if
        end do 
      end do
     end do   
    if (degr(o).gt.1) then
     clus(o)=dble(T)/dble(degr(o)*(degr(o)-1))
     end if   
   end do
!   
!
! Output of cluster coeficient (cc)
   open(unit=34,file='clus_'//trim(flname)//'.dat',status='unknown')
   write(34,*) '#degree',' ','cluster coef'
   call ordenareal(n,clus,degr)
   do f=1,n 
    if (degr(f).gt.1)then
    clus11(1)=0
       do g=1,n
         if(degr(f).eq. degr(g)) then
           clus11(g)=clus11(g)+clus(f)
         end if 
       end do
       if (degr(f).ne. degr(f+1)) then 
         write(34,*) degr(f), clus11(f)/his(f)
       end if
     end if  
    end do   
   close(34)


   stop
   end
   
   
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       contador: Nk (nombre de nodes amb grau k)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine hiso(degre,n,his)
         implicit none
         integer*8 i,x,n,allocat
         integer*8, allocatable :: degre(:)
         double precision, allocatable :: his(:)
         
         
         allocate(degre(1:n),stat=allocat)
         allocate(his(1:n),stat=allocat)
         
         
           x=0
           his(1:n)=0
             do i=1,n
               if (degre(i-1) .eq. degre(i)) then
               x=x+1
                else
               x=1
              end if
            his(i)=dble(x)  
            end do

         return
         end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ordena dades
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           subroutine ordena(n,pinic,pfina,degre)
           implicit none

           integer*8 j,k,n,allocat, pas1,pas2,pas3
           integer*8, allocatable :: degre(:),pinic(:),pfina(:)
!             sort the list
           allocate(degre(1:n),stat=allocat)                  
           allocate(pinic(1:n),stat=allocat)
           allocate(pfina(1:n),stat=allocat)
           
           j=0
             do j=1,n-1
                do k=j+1,n 
                 if(degre(j).gt.degre(k)) then
                  pas1=degre(k)
                  degre(k)=degre(j)
                  degre(j)=pas1
                  
                  pas2=pinic(k)
                  pinic(k)=pinic(j)
                  pinic(j)=pas2
                  
                  pas3=pfina(k)
                  pfina(k)=pfina(j)
                  pfina(j)=pas3
                 end if
                end do             
              end do 
              
              return
              end               
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         ordena dades
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           subroutine ordenareal(n,kuns,degre)
           implicit none

           integer*8 j,k,n,allocat, pas1
           double precision pas2
           integer*8, allocatable :: degre(:)
           double precision, allocatable :: kuns(:)
!             sort the list
           allocate(degre(1:n),stat=allocat)
           allocate(kuns(1:n),stat=allocat)                  
           j=0
             do j=1,n-1
                do k=j+1,n 
                 if(degre(j).gt.degre(k)) then
                  pas1=degre(k)
                  degre(k)=degre(j)
                  degre(j)=pas1
                  
                  pas2=kuns(k)
                  kuns(k)=kuns(j)
                  kuns(j)=pas2
                 end if
               end do             
              end do 
                            
         return
         end              
              
              
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      codificacio xarxa
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine codi(v,p1,p2,deg,n,fname,e)
  implicit none
!  
 integer*8 :: e,n,i,k,l!,istat,j,
 integer*8 :: ai,aj,nz
 integer*8 :: allocstat
 integer*8, allocatable :: n1(:),n2(:)
 integer*8, allocatable :: v(:), deg(:),p1(:), p2(:)
 character(*) fname
 
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       interface
         subroutine cleaner(filname,nosim1,nosim2,edfin,n)
         integer*8 edfin,n
         integer*8, dimension(:) :: nosim1(:),nosim2(:)
         character(*) filname
         end subroutine cleaner
       end interface
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! nodes and edges  
   call cleaner(fname,n1,n2,e,n)
 
  !ccccccccccccccccccccccccccccccccccccccccc 
  ! dimensio of vectors
 
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
 
  
 ! write the number of pinicial i pfinal 
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
 !   vector v that gives finals
    v(p2(ai))=aj
    v(p2(aj))=ai
    end do 
   nz=0
   do l=1,n
    if (deg(l) .eq. 0) nz=nz+1
   end do 
   write(*,*) 'nodes with zerodegree=', nz
 
  
  return
  end subroutine codi
  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       subroutina per netejar la xarxa
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine  cleaner(filname,nosim1,nosim2,edfin,n)
     implicit none
     integer*8 :: i,j,k,r,s,p,l,t,ss,rr
     integer*8 :: e,n,edfin
     integer*8 :: allocat,istat
     integer*8, allocatable :: net1(:),net2(:)
     integer*8, allocatable :: nolop1(:),nolop2(:),nosim1(:),nosim2(:)
     integer*8, allocatable :: norep1(:,:)
     character(*) filname
     
 
!   read files 
     call read_net(filname,e,n)
!     
!----------------------------------------------
!dimensions
      allocate(net1(1:e),stat=allocat)
      allocate(net2(1:e),stat=allocat)
      allocate(nolop1(1:e),stat=allocat)
      allocate(nolop2(1:e),stat=allocat)
      allocate(norep1(1:e,2),stat=allocat)  !array is norep(maxdimension,number columns)
      allocate(nosim1(1:e),stat=allocat)
      allocate(nosim2(1:e),stat=allocat)

! rename i and j

   open(unit=1,file=trim(filname)//'.dat',status='old')
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
 ! clean loops:        
     l=0
      do r=1,e
         if(net2(r) /= net1(r)) then
           l=l+1
           nolop1(l)=net1(r)
           nolop2(l)=net2(r)
         end if
       end do
 
! clean repetition ij=ij   
   j=1
   do s=1,l
     k=1
       do while (((norep1(k,1)/= nolop1(s)).or.(norep1(k,2)/= nolop2(s))) .and. k < j)
       k=k+1  
       end do
        if (k == j) then
         norep1(k,1)=nolop1(s)
         norep1(k,2) = nolop2(s)
         j=j+1
        end if 
    end do 
   

! clean symmetrics i,j = j,i
   p=1
   do ss=1,j-1
     rr=1
       do while (((nosim1(rr)/=norep1(ss,2)).or.(nosim2(rr)/=norep1(ss,1))) .and. rr < p)
       rr=rr+1  
       end do
        if (rr == p) then
         nosim1(rr) = norep1(ss,1)
         nosim2(rr) = norep1(ss,2)
         p=p+1
        end if 
    end do 
    
    
    edfin=p-1
  ! 
  !  
    open(unit=2,file='file_clean.dat',status='unknown')   
    do t=1,p-1
      write(2,*) nosim1(t),nosim2(t)
    end do
    close(2) 
     
     
     return
     end subroutine cleaner
     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    read the net
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine read_net(fname,e,n)
    implicit none
    integer*8 e,n,istat,i,j
    character(*) fname
     
    open(unit=1,file=trim(fname)//'.dat',status='old')
  
!   read files
       e=0
       n=0
       istat=0
       do while (istat.eq.0)
       read(1,*,iostat=istat) i,j
         if(istat.eq.0) then
          e=e+1
          if (i.gt.n) n=i
          if (j.gt. n) n=j
         end if
       end do
 
      close(1)
    
    return
    end subroutine read_net
