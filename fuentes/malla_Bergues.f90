subroutine malla_Bergues
      use def_variables
      use def_constantes
      use def_bio
      implicit none
      ! local
      double precision,allocatable :: xaux(:),yaux(:),yaux2(:)
      integer, allocatable:: mapa(:,:)
      integer :: np,NELE
      integer ::i,j,k,ii,jj,kk
      double precision :: xx,xxz,VERMALLAX,VERMALLAy,RESTO,xl,yl,ZPOS1,ZPOS2,rad,x,y,aux(16)
      double precision :: XL_int,YL_int
      integer :: npos,jjpiso,nelectro_cont

      if(exentric==0.0) then ! esfera
         radio_a = dis_ele*0.5/sin(ang30)
      elseif(exentric>0 .and. exentric<1.0) then  ! elipse

         radio_b = dis_ele
         radio_a = radio_b/sqrt(1- exentric*exentric)
      
      elseif(exentric>1.0) then  ! hiperbola

         radio_a =dis_ele
         radio_b = radio_a*sqrt(exentric*exentric-1.0)
      elseif(exentric==-1.0) then  ! electrodos alineados

         radio_a =dis_ele
         radio_b =dis_ele*0.5
      elseif(exentric==-2.0) then  ! electrodos alineados

         radio_a =dis_ele
         radio_b =dis_ele

      endif

      ! Defino el dominio interno para los electrodos. Este estara mallado en elementos cuadrangulares de 16 nodos

      XL_int = nint(2*radio_a+2*radio_a)
      YL_int = XL_int
      
      XL = zonatot
      YL = XL

      AreaDOMINIO=XL*YL

      if(exentric==-2.0) then
          XL_int = nint((N_electro-1)*radio_a)
          YL_int = XL_int
      
          XL = nint((N_electro+1)*radio_a)
          YL = XL
          
          ne_ber=60

      endif
                    

      NP=3*ne_ber+1
         
      
      XX = 1.0/real(np-1)

      ! lo total del dominio
      nnodes=  NP* NP
      nelements=ne_ber*ne_ber

      allocate(xaux(np),yaux(np),yaux2(2*np), mapa(np,np))
      
      allocate(coor_x(nnodes),coor_y(nnodes),grad_x(nnodes),grad_y(nnodes))
      allocate(conect(nelements,nop_b),material(nelements))
      allocate(gradxel_x(nelements),gradxel_y(nelements))

      allocate(tempera(nnodes),tempera_ant(nnodes),conducta(nelements))

      tempera=Tblod
      tempera_ant=Tblod
      gradxel_x=0.0
      gradxel_y=0.0
      grad_x=0.0
      grad_y=0.0

      DO I=1,NP
        XAUX(I) = (I-1)*XX
        XAUX(I) = XAUX(I)*XL
        YAUX(I) = (I-1)*XX
        YAUX(I) = YAUX(I)*YL
      ENDDO

      WRITE(unit_cont,*) 'MALLO LAS LINEAS',NP


      DO I=1,NP
          coor_x(I) = XAUX(I)
          YAUX2(I)  = YAUX(I)
      ENDDO


      DO K=1,NP
        MAPA(K,1) = K
        coor_y(K) = YAUX2(1)
      ENDDO

      WRITE(unit_cont,*) 'FRACCIONO LAS LINEAS',NP

      DO KK=2,NP
        DO J=1,NP
          coor_x(NP*(KK-1)+J) = coor_x(J)
          coor_y(NP*(KK-1)+J) = YAUX2(KK)
          MAPA(J,KK)=NP*(KK-1) + J
        ENDDO
      ENDDO

     
      WRITE(unit_cont,*) 'MAPEO ',NP


      npos=nop_b/4 - 1 

      NELE=0
      DO jj=1,ne_ber

         jjpiso = (jj-1)*npos+1
            
         
         DO kk=1,ne_ber
            NELE=NELE+1

            conect(NELE,1)=MAPA((kk-1)*npos+1, jjpiso)
            conect(NELE,2)=MAPA((kk-1)*npos+1+npos, jjpiso)
            conect(NELE,3)= MAPA((kk-1)*npos+1+npos,jjpiso+npos)
            conect(NELE,4)= MAPA((kk-1)*npos+1,jjpiso+npos)

            conect(NELE,5)= MAPA((kk-1)*npos+1+1, jjpiso)
            conect(NELE,6)= MAPA((kk-1)*npos+1+2, jjpiso)
            conect(NELE,7)= MAPA((kk-1)*npos+1+npos,jjpiso+1)
            conect(NELE,8)= MAPA((kk-1)*npos+1+npos,jjpiso+2)

            conect(NELE,9)=MAPA((kk-1)*npos+1+2, jjpiso+npos)
            conect(NELE,10)=MAPA((kk-1)*npos+1+1, jjpiso+npos)
            conect(NELE,11)= MAPA((kk-1)*npos+1,jjpiso+2)
            conect(NELE,12)= MAPA((kk-1)*npos+1,jjpiso+1)

            conect(NELE,13)=MAPA((kk-1)*npos+1+1, jjpiso+1)
            conect(NELE,14)=MAPA((kk-1)*npos+1+2, jjpiso+1)
            conect(NELE,15)= MAPA((kk-1)*npos+1+2,jjpiso+2)
            conect(NELE,16)= MAPA((kk-1)*npos+1+1,jjpiso+2)

         ENDDO


      ENDDO

      
!      open(unit=unit_2d,file=file_2D)
!      write(unit_2D,*) nnodes
!      do kk=1,nnodes
!         write(unit_2D,*) kk,coor_x(kk),coor_y(kk)
!      enddo
!      write(unit_2D,*) nelements
!      do kk=1,nelements
!         write(unit_2D,*) kk,(conect(kk,jj),jj=1,nop_b)
!      enddo
      
!      close(unit_2D)


      ! call verifico(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel)

      call cond_contgen_Bergues(XL,yl,XL_int,YL_int)
      write(6,*) 'condiciones de contorno!' 

      call sacomalla2d(nnodes,coor_x,coor_y,conect,nelements,nop_b,material,unit_malla)

      write(6,*) 'saco malla'

      nodpel = nop_b

      deallocate(xaux,yaux,yaux2,mapa)


      end subroutine malla_Bergues


      
      subroutine cond_contgen_Bergues(XL,yl,XL_int,YL_int)
      use def_variables
      use def_constantes
      implicit none
      double precision :: XL,yl,XL_int,YL_int
      ! local
      integer :: kk,jj,j,mat_tierra,mat_potencial,mat_extra,mat_tot,ne_ini,mat_afuer,ii,nelectro_cont
      double precision :: funm,funb,m,b,h1,h2,h3,h4,z,rmed,xmed,ymed,zmed,rcil,rad,xbase,ybase,xtech,ytech
      double precision :: buscocero,ax,bx,cx,xpos,ypos,buscocero2,max_cor
       double precision,allocatable :: dismin(:),dist(:)
      
      ! determino condiciones de contorno!! nodos
      allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes),nnodpoten(nnodes),nnodtierra(nnodes))
      vec_tierra=-1
      vec_poten=-1
      vec_temper=-1
      
      nod_tierra=N_electro/2
      nod_poten=N_electro/2

      N_electro_tot=N_electro*N_electro
      
       
      if(exentric >-2) then
         allocate(fuentes(N_electro/2),sumideros(N_electro/2),xfuentes(N_electro/2),xsumideros(N_electro/2),  &
           &       yfuentes(N_electro/2),ysumideros(N_electro/2))
         allocate(dismin(N_electro),dist(N_electro))
      else
         allocate(fuentes(N_electro),sumideros(N_electro_tot),xfuentes(N_electro),xsumideros(N_electro_tot),  &
           &       yfuentes(N_electro),ysumideros(N_electro_tot))
         allocate(dismin(2*N_electro_tot),dist(2*N_electro_tot))

      endif
      


      ! determino materiales
      
      xbase = (XL - XL_int )*0.5
      ybase = (YL - YL_int )*0.5
      xtech = xbase +  XL_int 
      ytech = ybase +  YL_int 
      
      mat_normal=0
      mat_tumor=0
      mat_afuer=0

      do kk=1,nelements
         xmed = 0.0
         ymed = 0.0
         do jj=1,nop_b
           j=conect(kk,jj)
           xmed =xmed+ dsqrt(coor_x(j)*coor_x(j))
           ymed =ymed+ dsqrt(coor_y(j)*coor_y(j))
         enddo
         xmed=xmed/real(nop_b)
         ymed=ymed/real(nop_b)

         if(xmed<=xbase .or. xmed>=xtech) then  ! zona normal
             material(kk)=2
             mat_normal=mat_normal+1
         elseif(ymed<=ybase .or. ymed>=ytech) then  ! zona normal
            material(kk)=2
            mat_normal=mat_normal+1
         else
            material(kk)=1
            mat_tumor=mat_tumor+1

            if(problema=='PAPA2D') then
               material(kk)=2
            endif

         endif
      

     enddo


     ! para un material afuera
     !do kk=1,nelements
     !    xmed = 0.0
     !    ymed = 0.0
     !    do jj=1,nop_b
     !      j=conect(kk,jj)
     !      xmed =xmed+ dsqrt(coor_x(j)*coor_x(j))
     !      ymed =ymed+ dsqrt(coor_y(j)*coor_y(j))
     !    enddo
     !    xmed=xmed/real(nop_b)
     !    ymed=ymed/real(nop_b)
         
     !    if(xmed <10.0 .or. xmed > XL-10.0) then
     !        material(kk)=3
     !        mat_afuer=mat_afuer+1
     !    endif
     !    if(ymed <10.0 .or. ymed > YL-10.0) then
     !        material(kk)=3
     !        mat_afuer=mat_afuer+1
     !    endif


     !enddo

     if(exentric==0.0) then ! esfera
        
        if(N_electro==2) then
           xfuentes(1)= XL*0.5 -  dis_ele*0.5               
           yfuentes(1)= YL*0.5
           
           xsumideros(1)=XL*0.5 +  dis_ele*0.5                
           ysumideros(1)= YL*0.5              
         endif
       if(modo==1) then
            if(N_electro==6) then
                   xfuentes(1)= XL*0.5 +  dis_ele*0.5               
                   yfuentes(1)= YL*0.5 + radio_a*cos(ang30)

                   xsumideros(1)=XL*0.5 + radio_a             
                   ysumideros(1)= YL*0.5              

                   xfuentes(2)= XL*0.5 +  dis_ele*0.5               
                   yfuentes(2)= YL*0.5 - radio_a*cos(ang30)

                   xsumideros(2)=XL*0.5 -  dis_ele*0.5                
                   ysumideros(2)= YL*0.5 - radio_a*cos(ang30)             

                   xfuentes(3)= XL*0.5 - radio_a           
                   yfuentes(3)= YL*0.5
           
                   xsumideros(3)=XL*0.5 -  dis_ele*0.5                
                   ysumideros(3)= YL*0.5 + radio_a*cos(ang30)             

             endif
         elseif(modo==2) then
            if(N_electro==6) then
                   xsumideros(1)= XL*0.5 +  dis_ele*0.5               
                   ysumideros(1)= YL*0.5 + radio_a*cos(ang30)

                   xsumideros(2)=XL*0.5 + radio_a             
                   ysumideros(2)= YL*0.5              
                                    
                   
                   xsumideros(3)= XL*0.5 +  dis_ele*0.5               
                   ysumideros(3)= YL*0.5 - radio_a*cos(ang30)
                   
                   xfuentes(1)= XL*0.5 - radio_a           
                   yfuentes(1)= YL*0.5
                   
                   xfuentes(2)=XL*0.5 -  dis_ele*0.5                
                   yfuentes(2)= YL*0.5 - radio_a*cos(ang30)             
           
                   xfuentes(3)=XL*0.5 -  dis_ele*0.5                
                   yfuentes(3)= YL*0.5 + radio_a*cos(ang30)             

             endif



         endif
     elseif(exentric==1) then ! parabola
        
        if(N_electro==6) then
           xfuentes(1)= XL*0.5 +  dis_ele*0.5             
           yfuentes(1)= YL*0.5 + (dis_ele*0.5)**2  -dis_ele

           xsumideros(1)=XL*0.5 - dis_ele*0.5                 
           ysumideros(1)= YL*0.5 +(dis_ele*0.5)**2         -dis_ele    

           ax = dis_ele*0.5 ! 1-dis_ele*dis_ele*0.5
           bx= ax*ax ! -dis_ele
           cx=dis_ele*dis_ele ! dis_ele**4-3*dis_ele*dis_ele/4.0

           xpos = buscocero(ax,bx,cx,dis_ele*0.5 )

           xfuentes(2)= XL*0.5 -  xpos               
           yfuentes(2)= YL*0.5 + xpos**2-dis_ele

           xsumideros(2)=XL*0.5 + xpos
           ysumideros(2)= YL*0.5 + xpos**2-dis_ele

           ax = xpos ! dis_ele*0.5 ! 1-dis_ele*dis_ele*0.5
           bx= ax*ax ! -dis_ele
           cx=dis_ele*dis_ele ! dis_ele**4-3*dis_ele*dis_ele/4.0

           xpos = buscocero(ax,bx,cx,xpos )

           xfuentes(3)= XL*0.5 + xpos
           yfuentes(3)= YL*0.5 + xpos**2-dis_ele
           
           xsumideros(3)=XL*0.5 - xpos              
           ysumideros(3)= YL*0.5 + xpos**2  -dis_ele      
        endif
      elseif(exentric> 0 .and. exentric<1) then ! elipse

            
       if(modo==2) then
           xsumideros(1)= XL*0.5 +  radio_a               
           ysumideros(1)= YL*0.5 


           ypos= sqrt(radio_a*radio_a*cos(ang60*2)**2 + radio_b*radio_b*sin(ang60*2)**2)
           xsumideros(2)=XL*0.5  -  ypos * cos(ang60*2)               
           ysumideros(2)= YL*0.5 -  ypos * sin(ang60*2)
           
           xsumideros(3)=XL*0.5  -  ypos * cos(ang60*2)               
           ysumideros(3)= YL*0.5 + ypos* sin(ang60*2)


           
           xfuentes(1)= XL*0.5 + ypos * cos(ang60*2)              
           yfuentes(1)= YL*0.5 - ypos * sin(ang60*2)
           
           
           xfuentes(2)= XL*0.5   + ypos * cos(ang60*2)                           
           yfuentes(2)= YL*0.5   + ypos* sin(ang60*2)
           
           xfuentes(3)=XL*0.5 - radio_a             
           yfuentes(3)= YL*0.5              
           
           

       elseif(modo==1) then
           xsumideros(1)= XL*0.5 +  radio_a               
           ysumideros(1)= YL*0.5 


           
           
           xfuentes(3)=XL*0.5 - radio_a             
           yfuentes(3)= YL*0.5              
           
           !ypos = dsqrt(radio_b*radio_b - (dis_ele*0.5)**2*radio_b*radio_b/(radio_a*radio_a))

           ypos= sqrt(radio_a*radio_a*cos(ang60*2)**2 + radio_b*radio_b*sin(ang60*2)**2)
           
           xsumideros(3)= XL*0.5 + ypos * cos(ang60*2)              
           ysumideros(3)= YL*0.5 - ypos * sin(ang60*2)
           
           !xsumideros(3)= XL*0.5 -  dis_ele*0.5               
           !ysumideros(3)= YL*0.5 + ypos
           
           xfuentes(1)=XL*0.5  -  ypos * cos(ang60*2)               
           yfuentes(1)= YL*0.5 -  ypos * sin(ang60*2)
           
           xsumideros(2)= XL*0.5   + ypos * cos(ang60*2)                           
           ysumideros(2)= YL*0.5   + ypos* sin(ang60*2)
           
           xfuentes(2)=XL*0.5  -  ypos * cos(ang60*2)               
           yfuentes(2)= YL*0.5 + ypos* sin(ang60*2)
         
         endif

      elseif(exentric> 1.0) then ! hiperbola

           xfuentes(1)= XL*0.5   +  radio_a               
           yfuentes(1)= YL*0.5 

           xsumideros(1)=XL*0.5  - radio_a             
           ysumideros(1)= YL*0.5              
           

           cx=dis_ele*dis_ele 
           xpos=radio_a
           xpos = buscocero2(radio_a,radio_b,cx,xpos)
           

           xfuentes(2)= XL*0.5 - xpos               
           yfuentes(2)= YL*0.5 + dsqrt(xpos*xpos*radio_b*radio_b/radio_a/radio_a - radio_b*radio_b)
           
           xsumideros(2)=XL*0.5 + xpos               
           ysumideros(2)= YL*0.5 +dsqrt(xpos*xpos*radio_b*radio_b/radio_a/radio_a - radio_b*radio_b)

           xfuentes(3)= XL*0.5 -xpos                           
           yfuentes(3)= YL*0.5 - dsqrt(xpos*xpos*radio_b*radio_b/radio_a/radio_a - radio_b*radio_b)
           
           xsumideros(3)=XL*0.5 + xpos             
           ysumideros(3)= YL*0.5 - dsqrt(xpos*xpos*radio_b*radio_b/radio_a/radio_a - radio_b*radio_b)
     elseif(exentric==-1.0) then ! alineados

        if(modo==1) then
           xfuentes(1)= XL*0.5+radio_a*0.5               
           yfuentes(1)= YL*0.5 

           xsumideros(1)=XL*0.5 - radio_a*0.5             
           ysumideros(1)= YL*0.5              
           

           xfuentes(2)= XL*0.5 -radio_a*0.5                 
           yfuentes(2)= YL*0.5 + radio_b
           
           xsumideros(2)=XL*0.5 +radio_a*0.5                
           ysumideros(2)= YL*0.5 + radio_b

           xfuentes(3)= XL*0.5 -radio_a*0.5                           
           yfuentes(3)= YL*0.5 - radio_b
           
           xsumideros(3)=XL*0.5 +radio_a*0.5             
           ysumideros(3)= YL*0.5 - radio_b      
           
        elseif(modo==2) then

           xfuentes(1)= XL*0.5+radio_a*0.5               
           yfuentes(1)= YL*0.5 

           xsumideros(1)=XL*0.5 - radio_a*0.5             
           ysumideros(1)= YL*0.5              
           

           xsumideros(2)= XL*0.5 -radio_a*0.5                 
           ysumideros(2)= YL*0.5 + radio_b
           
           xfuentes(2)=XL*0.5 +radio_a*0.5                
           yfuentes(2)= YL*0.5 + radio_b

           xsumideros(3)= XL*0.5 -radio_a*0.5                           
           ysumideros(3)= YL*0.5 - radio_b
           
           xfuentes(3)=XL*0.5 +radio_a*0.5             
           yfuentes(3)= YL*0.5 - radio_b      



        endif   
     elseif(exentric==-2.0) then ! intercalados
         
         ! fuentes
         do kk=1,N_electro
            xfuentes(kk)=xbase 
            yfuentes(kk)=ybase + (kk-1)*radio_b
         enddo
         ! sumideros
         do kk=2,N_electro
            do jj=1,N_electro
               nelectro_cont = (kk-1)*N_electro+jj
               xsumideros(nelectro_cont)=xbase + (kk-1)*radio_a
               ysumideros(nelectro_cont)=ybase + (jj-1)*radio_b
            enddo                
         enddo


     endif
     ! busco fuentes y sumideros en nodos!!
     if(exentric>-2  ) then

         dismin = XL
         do kk=1,nnodes

             do jj=1,N_electro/2
             
                 dist(2*jj-1) = sqrt((coor_x(kk)- xfuentes(jj))**2 + (coor_y(kk)- yfuentes(jj))**2 )
                 if( dist(2*jj-1)< dismin(2*jj-1) ) then
                      dismin(2*jj-1) = dist(2*jj-1)
                      fuentes(jj)=kk
                 endif
             
             enddo

             do jj=1,N_electro/2
             
                 dist(2*jj) = sqrt((coor_x(kk)- xsumideros(jj))**2 + (coor_y(kk)- ysumideros(jj))**2 )
                 if( dist(2*jj)< dismin(2*jj) ) then
                      dismin(2*jj) = dist(2*jj)
                      sumideros(jj)=kk
                 endif
             
             enddo

         enddo
      else
      
         dismin = XL
         do kk=1,nnodes

         ! fuentes iniciales
             do jj=1,N_electro
             
                 dist(2*jj-1) = sqrt((coor_x(kk)- xfuentes(jj))**2 + (coor_y(kk)- yfuentes(jj))**2 )
                 if( dist(2*jj-1)< dismin(2*jj-1) ) then
                      dismin(2*jj-1) = dist(2*jj-1)
                      fuentes(jj)=kk
                 endif
             
             enddo
             
        
         ! sumideros iniciales
              ! sumideros
            do ii=2,N_electro
               do jj=1,N_electro
                  nelectro_cont = (ii-1)*N_electro+jj
               
                    dist(2*nelectro_cont) = sqrt((coor_x(kk)- xsumideros(nelectro_cont))**2 + (coor_y(kk)- ysumideros(nelectro_cont))**2 )
                 if( dist(2*nelectro_cont)< dismin(2*nelectro_cont) ) then
                      dismin(2*nelectro_cont) = dist(2*nelectro_cont)
                      sumideros(nelectro_cont)=kk
                 endif
               enddo
             enddo

         enddo

      endif

      write(unit_cc,* )' dis_ele:  ',dis_ele
      write(unit_cc,* )' radio_a:  ',radio_a
      write(unit_cc,* )' radio_b:  ',radio_b

     
      write(unit_cc,* )'*** fuentes '
     
      if(exentric>-2 ) then
         do jj=1,N_electro/2
              write(unit_cc,'(2I6,2E15.5)') jj,fuentes(jj),xfuentes(jj),yfuentes(jj)
         enddo
      else
         do jj=1,N_electro
              write(unit_cc,'(2I6,2E15.5)') jj,fuentes(jj),xfuentes(jj),yfuentes(jj)
         enddo
      endif   
      write(unit_cc,* )'*** sumideros '
      if(exentric>-2 ) then
         do jj=1,N_electro/2
              write(unit_cc,'(2I6,2E15.5)') jj,sumideros(jj),xsumideros(jj),ysumideros(jj)
         enddo
      else
           do kk=2,N_electro
               do jj=1,N_electro
                   nelectro_cont = (kk-1)*N_electro+jj
                   write(unit_cc,'(2I6,2E15.5)') nelectro_cont,sumideros(nelectro_cont),xsumideros(nelectro_cont),ysumideros(nelectro_cont)
               enddo
         enddo
      endif   

       
      nod_tierra=0
      nod_poten=0
      
    if(modo==1 .or. modo==2) then
      do kk=1,nnodes
        do jj=1,N_electro*0.5
          if(kk==sumideros(jj)) then
             nod_tierra=nod_tierra+1
             vec_tierra(kk)=1
             nnodtierra(nod_tierra)=kk
            write(unit_cc,*) nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)) 
          endif
          if(kk==fuentes(jj)) then
               nod_poten=nod_poten+1
               vec_poten(kk)=1
               nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten))
           endif
        enddo
      enddo
   elseif(modo==3) then
     do kk=1,nnodes
          if(kk==sumideros(1)) then
             nod_tierra=nod_tierra+1
             vec_tierra(kk)=1
             nnodtierra(nod_tierra)=kk
            write(unit_cc,*) nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)) 
          endif
          if(kk==fuentes(3)) then
             nod_tierra=nod_tierra+1
             vec_tierra(kk)=1
             nnodtierra(nod_tierra)=kk
            write(unit_cc,*) nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)) 
          endif


          if(kk==fuentes(1)) then
               nod_poten=nod_poten+1
               vec_poten(kk)=1
               nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten))
           endif
          if(kk==sumideros(3)) then
               nod_poten=nod_poten+1
               vec_poten(kk)=1
               nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten))
           endif

          if(kk==fuentes(2)) then
               nod_poten=nod_poten+1
               vec_poten(kk)=2
               nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten))
           endif
           if(kk==sumideros(2)) then
               nod_poten=nod_poten+1
               vec_poten(kk)=2
               nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten))
           endif

      enddo
   elseif(modo==4) then
      !write(unit_cc,*) 'nodpoten-nodtierra'
      do kk=1,nnodes
          do jj=1,N_electro
             if(kk==fuentes(jj)) then
                nod_poten=nod_poten+1
                vec_poten(kk)=1
                nnodpoten(nod_poten)=kk
                
             endif
          enddo
          do jj=2,N_electro
               do ii=1,N_electro
                   nelectro_cont = (jj-1)*N_electro+ii
                   if(kk==sumideros(nelectro_cont)) then
                       nod_tierra=nod_tierra+1
                       vec_tierra(kk)=1
                       nnodtierra(nod_tierra)=kk
                      ! write(unit_cc,'(2I6,2E15.5)') nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)) 
                   endif
               enddo
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
      



   endif
      
      
      ! chequeo materiales
      mat_normal=0
      mat_tumor=0

      do kk=1,nelements
        if(material(kk) == 2 ) then
            mat_normal=mat_normal+1
        elseif(material(kk) == 1 ) then
            mat_tumor=mat_tumor+1
        endif

      enddo

      
      write(unit_cc,* )' materiales ', nelements
      write(unit_cc,* )'      tumor  :   ', mat_tumor,sigma1
      write(unit_cc,* )'      normal :   ', mat_normal,sigma2
      write(unit_cc,* )'      afuera :   ', mat_afuer,sigma3
      mat_tot = mat_tumor+mat_normal+mat_afuer
      write(unit_cc,* )'      tot:  ', mat_tot

      write(unit_cc,* )'  ******************************************  '
      write(unit_cc,* )'  ******************************************  '



      ! condiciones para la temperatura
      max_cor=0.0
      do kk=1,nnodes
         if(coor_x(kk)>max_cor) then
             max_cor=coor_x(kk)
         endif
      enddo

      nod_temper=0
      do kk=1,nnodes
         
         if(coor_x(kk)<=0.001) then
            
            vec_temper(kk)=1
            nod_temper=nod_temper+1

        elseif(coor_y(kk)<=0.001) then
        
            vec_temper(kk)=1
            nod_temper=nod_temper+1
        elseif(coor_x(kk)>=max_cor-0.001) then
        
            vec_temper(kk)=1
            nod_temper=nod_temper+1
        elseif(coor_y(kk)>=max_cor-0.001) then
        
            vec_temper(kk)=1
            nod_temper=nod_temper+1

        endif


      enddo

      
      
      allocate(nnodtemper(nod_temper))
      
      write(unit_cc,* )' C de C temper ', nod_temper
      nod_temper=0
      do kk=1,nnodes
         if(vec_temper(kk)==1) then
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk

            write(unit_cc,'(2i6,2E15.5)') nod_temper,nnodtemper(nod_temper),coor_x(kk),coor_y(kk)
         endif
      enddo


      end subroutine cond_contgen_Bergues

      double precision function buscocero(ax,bx,cx,x)
       implicit none
      double precision :: ax,bx,cx,x
      ! local
      double precision :: resto
      integer :: kk=0


      resto = -1.0
      do while(resto<0.0 .and. kk<10000)
         kk=kk+1
         x=x+0.001
         resto = ( (x-ax)**2+ (x*x-bx)**2)-cx   ! x*x*x*x+x*x*ax-x*bx+cx

      enddo
      buscocero=x

      end function buscocero


      double precision function buscocero2(a,b,cx,x)
       implicit none
      double precision :: a,b,cx,x
      ! local
      double precision :: resto
      integer :: kk=0


      resto = -1.0
      do while(resto<0.0 .and. kk<10000)
         kk=kk+1
         x=x+0.001
         resto = ( (x-a)**2+ (x*x*b*b/(a*a) - b*b))-cx   ! x*x*x*x+x*x*ax-x*bx+cx

      enddo
      buscocero2=x

      end function buscocero2
