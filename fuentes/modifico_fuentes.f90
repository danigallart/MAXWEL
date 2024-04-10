subroutine modifico_fuentes(nt,time)
use def_variables
implicit none
integer::nt
double precision :: time
! local
integer ::naux(N_electro)
integer ::nelectro_cont,kk,ii,jj,i,k
double precision :: xfuentes_aux(N_electro),yfuentes_aux(N_electro)


if(mod(time,Alterna)==0) then
   
   do k=1,N_electro
       naux(k)=fuentes(k)
       xfuentes_aux(k)=xfuentes(k)
       yfuentes_aux(k)=yfuentes(k)
       fuentes(k)=sumideros(N_electro*base_electro+k)
       xfuentes(k)=xsumideros(N_electro*base_electro+k)
       yfuentes(k)=ysumideros(N_electro*base_electro+k)

       sumideros(N_electro*(base_electro-1)+k)=naux(k)
       xsumideros(N_electro*(base_electro-1)+k)=xfuentes_aux(k)
       ysumideros(N_electro*(base_electro-1)+k)=yfuentes_aux(k)

   enddo    
   base_electro=base_electro+1 
   
   write(unit_cc,* )'*** fuentes '
    
    do jj=1,N_electro
         nelectro_cont=N_electro*(base_electro-1)+jj
         write(unit_cc,'(2I6,2E15.5)') nelectro_cont,fuentes(jj),xfuentes(jj),yfuentes(jj)
    enddo

    write(unit_cc,* )'*** sumideros '
    do kk=1,(base_electro-1)*N_electro
         write(unit_cc,'(2I6,2E15.5)') kk,sumideros(kk),xsumideros(kk),ysumideros(kk)
    enddo
    do kk=(base_electro)*N_electro+1,N_electro_tot
         write(unit_cc,'(2I6,2E15.5)') kk,sumideros(kk),xsumideros(kk),ysumideros(kk)
    enddo
      
  
   vec_poten=-1
   vec_tierra=-1
   nnodtierra=0
   nnodpoten=0
   nod_poten=0
   nod_tierra=0
   
   write(unit_cc,*) '****************************  '
   write(unit_cc,*) 'paso  ', time
   
   do kk=1,nnodes
        do jj=1,N_electro
             if(kk==fuentes(jj)) then
                nod_poten=nod_poten+1
                vec_poten(kk)=1
                nnodpoten(nod_poten)=kk
             endif
          enddo
          do jj=1,(base_electro-1)*N_electro
               if(kk==sumideros(jj)) then
                  nod_tierra=nod_tierra+1
                  vec_tierra(kk)=1
                  nnodtierra(nod_tierra)=kk
               endif
         enddo
        do jj=(base_electro)*N_electro+1,N_electro_tot
               if(kk==sumideros(jj)) then
                  nod_tierra=nod_tierra+1
                  vec_tierra(kk)=1
                  nnodtierra(nod_tierra)=kk
               endif
         enddo
        

   enddo  

      write(unit_cc,*) 'nodpoten ',nod_poten
      do jj=1,nod_poten
         write(unit_cc,'(2I6,2E15.5)') jj,nnodpoten(jj),coor_x(nnodpoten(jj)), coor_y(nnodpoten(jj))
      enddo    
    
      write(unit_cc,*) 'nodtierra ',nod_tierra
      do jj=1,nod_tierra
         write(unit_cc,'(2I6,2E15.5)') jj,nnodtierra(jj),coor_x(nnodtierra(jj)), coor_y(nnodtierra(jj)) 
      enddo    



   write(unit_cc,*) ' ***** '
   

endif 


end subroutine modifico_fuentes
