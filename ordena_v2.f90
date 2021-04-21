!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Sorts the data
!			by Isaura Oliver
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
           real*8 pas2
           integer*8, allocatable :: degre(:)
           real*8, allocatable :: kuns(:)
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
