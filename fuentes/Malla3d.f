      subroutine malla3d()
      use def_variables
      use def_constantes
      use def_bio
      implicit none
      ! local
      double precision,allocatable :: xaux(:),yaux(:),zaux(:),yaux2(:),
     * zaux2(:)
      integer, allocatable:: mapa(:,:),nel(:,:)
      integer :: nez,np_ant,np,NELE,NPZ2,NPY2,NEGUAR,NE3,N,NEW,NPZ,
     * nsearchn,NPOINT
      integer ::i,j,k,ii,jj,kk,nsum,naux,nzgrid2,NPME,nele2,nfactor,
     * nfactorz
      double precision :: xx,xxz,VERMALLAX,VERMALLAy,VERMALLAZ,RESTO,
     *xl,yl,zl,ZPOS1,ZPOS2,rad,x,y,z,aux(9),fac_ne
      

      NE0=NE
      ncerx0=0
      ncery0=0
      ncerz0=0

      if(problema=='ELECTRO') then
         
         XL = rad_zone1+rad_zone2
         YL = rad_zone1+rad_zone2
         ZL = rad_zone1+rad_zone2
         NEZ= NE
          XX = 1.0/real(NE)

          NP=NE+1
          NP_ANT = NP

          nfactor=1
          nfactorz=1
          if(DUPLICY==1) nfactor=2
          if(DUPLICZ==1) nfactorz=2

          NPY2 = 2*nfactor*(NP-1)+1 

          ! lo total del dominio
          nnodes= (2*ne+1)*(nfactor*2*ne+1)*(2*nfactorz*nez+1)
          nelements=ne*nfactor*ne*nfactorz*nez

      elseif(problema=='CEREBRO' .or. problema=='PAPA' )  then
         XL = tam_zonex
         YL = tam_zoney
         ZL = tam_zonez
         
         area_trans=tam_zoney*tam_zonez 
         
         !if(nop_malla==3) area_trans=100.0

         fac_ne = tam_zonez/tam_zonex

         NEZ= NE *fac_ne
          XX = 1.0/real(NE)

          NP=NE+1
          NP_ANT = NP

          nfactor=1
          nfactorz=1
          if(DUPLICY==1) nfactor=2
          if(DUPLICZ==1) nfactorz=2

          NPY2 = 2*nfactor*(NP-1)+1 

          ! lo total del dominio
          nnodes= (2*ne+1)*(nfactor*2*ne+1)*(2*nfactorz*nez+1)
          nelements=ne*nfactor*ne*nfactorz*nez

      elseif(problema=='BERGUES3D')  then
         
         alfazone=1.0

         XL = tam_zonex *(1+alfazone) 
         YL = tam_zoney*(1+alfazone)
         ZL = tam_zonez*(1+alfazone)
         
         area_trans=tam_zoney*tam_zonez 
         
         
         fac_ne = tam_zonez/tam_zonex

         NEZ= NE *fac_ne
         XX = 1.0/real(NE)

         NP=NE+1
         NP_ANT = NP

         nfactor=1
         nfactorz=1
         if(DUPLICY==1) nfactor=2
         if(DUPLICZ==1) nfactorz=2

         !NPY2 = 2*nfactor*(NP-1)+1 

          ! lo total del dominio
          nnodes= (2*ne+1)*(nfactor*2*ne+1)*(2*nfactorz*nez+1)
          nelements=ne*nfactor*ne*nfactorz*nez


      endif
    

      allocate(xaux(np),yaux(np),yaux2(2*np),
     * mapa((2*ne+1)*(2*ne+1),(2*ne+1)*(2*ne+1)),
     *        nel(ne*nfactor*ne,nope*nope))
      
      allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),
     *grad_x(nnodes),grad_y(nnodes),grad_z(nnodes))
      allocate(conect(nelements,nodpel),material(nelements))
      allocate(gradxel_x(nelements),gradxel_y(nelements),
     *gradxel_z(nelements))
      allocate(tempera(nnodes),tempera_ant(nnodes),jcurrent(nelements),
     * quemado(nelements))
      
      gradxel_x=0.0
      gradxel_y=0.0
      gradxel_z=0.0
      grad_x=0.0
      grad_y=0.0
      grad_z=0.0
      tempera=Tout
      tempera_ant=Tout
      jcurrent=0.0
      quemado=0.0


      DO I=1,NP
        XAUX(I) = (I-1)*XX
        XAUX(I) = XAUX(I)*XL
        YAUX(I) = (I-1)*XX
        YAUX(I) = YAUX(I)*YL
      ENDDO

      IF(NLOG.EQ.1) THEN

        
        VERMALLAX = (escala/XAUX(NP))**(1./(NP-1))
        VERMALLAY = (escala/YAUX(NP))**(1./(NP-1))


        DO I=NP-1,1,-1
          XAUX(I) = XAUX(I+1)*VERMALLAX
          YAUX(I) = YAUX(I+1)*VERMALLAY
        ENDDO

        RESTO = XAUX(1)
        DO I=1,NP-1
          XAUX(I) = XAUX(I)-RESTO
        ENDDO

        RESTO = YAUX(1)
        DO I=1,NP-1
          YAUX(I) = YAUX(I)-RESTO
        ENDDO

      ENDIF
      
      WRITE(unit_cont,*) 'MALLO LAS LINEAS',NP


      DO I=1,NP
          coor_x(2*I-1) = XAUX(I)
          IF(I.NE.NP) coor_x(2*I)   = (XAUX(I)+XAUX(I+1))*0.5
          YAUX2(2*I-1) = YAUX(I)
          IF(I.NE.NP) YAUX2(2*I)   = (YAUX(I)+YAUX(I+1))*0.5
        ENDDO

      NP = (NOPE-1)*NE+1

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

      NP = NP*NP

      WRITE(unit_cont,*) 'MAPEO ',NP,NP_ANT

        NELE=0
        DO KK=1,NP_ANT-1
          DO JJ=1,NP_ANT-1
            NELE=NELE+1
            NEL(NELE,1)=MAPA(2*KK-1,2*JJ-1)
            NEL(NELE,2)=MAPA(2*KK  ,2*JJ-1)
            NEL(NELE,3)=MAPA(2*KK+1,2*JJ-1)

            NEL(NELE,4)= MAPA(2*KK+1,2*JJ)
            NEL(NELE,5)= MAPA(2*KK+1,2*JJ+1)
            NEL(NELE,6)= MAPA(2*KK  ,2*JJ+1)

            NEL(NELE,7)=MAPA(2*KK-1,2*JJ+1)
            NEL(NELE,8)= MAPA(2*KK-1,2*JJ)
            NEL(NELE,9)= MAPA(2*KK  ,2*JJ)


          ENDDO
        ENDDO

       !  DUPLICO EL LADO Y
      IF(DUPLICY.EQ.1) THEN
        NPME=0
        DO KK=1,NP
          IF(coor_y(KK).NE.0.0) THEN
             NPME = NPME+1
             coor_x(NP+NPME)=coor_x(KK)
             coor_y(NP+NPME)=-coor_y(KK)
          ENDIF
        ENDDO

        NP=NP+NPME

        WRITE(6,*) 'NODOS DUPLICY: ',NP
      
        NELE2=NELE
        DO KK=1,NELE2
          NELE=NELE+1
          DO J=1,NOPE*NOPE
            IF(coor_y(NEL(KK,J)).EQ.0.0) THEN
               NEL(NELE,J)= NEL(KK,J)
            ELSE
               NEL(NELE,J)= NPME + NEL(KK,J)
            ENDIF
          ENDDO
         
        ENDDO

        DO KK=NELE2+1,NELE

          
            DO J=1,9
              AUX(J)=NEL(KK,J)
            ENDDO
            NEL(KK,1) = AUX(7)
            NEL(KK,2) = AUX(6)
            NEL(KK,3) = AUX(5)
            NEL(KK,4) = AUX(4)
            NEL(KK,5) = AUX(3)
            NEL(KK,6) = AUX(2)
            NEL(KK,7) = AUX(1)
            NEL(KK,8) = AUX(8)
            NEL(KK,9) = AUX(9)
        
        ENDDO

        WRITE(6,*) 'ELEM DUPLICY: ',NELE
      ENDIF

      NE=NELE

      ! write(unit_2D,*) np
      ! do kk=1,np
      !   write(unit_2D,*) kk,coor_x(kk),coor_y(kk)
      !enddo
      !write(unit_2D,*) ne
      !do kk=1,ne
      !   write(unit_2D,*) kk,(nel(kk,jj),jj=1,nope*nope)
      !enddo
      
      !close(unit_2D)

C PARA 3D ojo EN AMBAS DIRECCIONES
     
      XXZ = XX ! 1.0/(NEZ)

      if(problema=='BERGUES3D') then

          allocate(zaux(2*NEZ+1),zaux2(2*NEZ+1))
      else
          allocate(zaux(4*(NEZ+1)),zaux2(4*(NEZ+1)+1))

      endif


      DO I=1,NEZ+1
        ZAUX(I) = (I-1)*XXZ
        ZAUX(I) = ZAUX(I)*zl
      ENDDO

      
      IF(NLOG.EQ.1) THEN

          
          VERMALLAZ = (ESCALA/ZAUX(NPZ2))**(1./(NPZ2-1))

          DO I=NPZ2-1,1,-1
            ZAUX( I) = ZAUX(I+1)*VERMALLAZ
          ENDDO
          
		
		  RESTO = ZAUX(1)
          DO I=1,NPZ2-1
            ZAUX(I) = ZAUX(I)-RESTO
          ENDDO

          DO I=2,NPZ2
             ZAUX2(I-1) = -ZAUX(NPZ2 - I + 2)
             ZAUX2(NPZ2+I-1)= ZAUX(I)
          ENDDO
          ZAUX2(NPZ2) =ZAUX(1)

! mallo lineas
          DO I=1,2*NPZ2-1
              ZAUX(2*I-1) = ZAUX2(I)
              IF(I.NE.2*NPZ2-1) ZAUX(2*I)   = (ZAUX2(I)+ZAUX2(I+1))*0.5
          ENDDO

      ENDIF

      ! if(problema=='BERGUES3D') then
      !    NPZ2 = 2*(NEZ)+1
      !    nzgrid= npz2
      !else
          NPZ2 = 2*(NEZ)+1
          nzgrid= npz2
      
      !endif  
        
        nzgrid2 = nzgrid*0.5 + 1
        !nzgrid2 = 2*(NEZ)+1
       allocate(Zgrid(nzgrid),grid(nzgrid,nzgrid2) )
        DO KK=1,NPZ2
            WRITE(unit_cont,*) KK,ZAUX(KK)
            zgrid(kk)=ZAUX(KK)
        ENDDO

      NE3 = 0
      NPZ = NP
      NPOINT = NOPE*NOPE

      DO J=1,NPZ
        coor_z(J)=0.0 ! -ZL
      ENDDO

      DO KK=1,NEZ ! 2*NEZ
  
        ZPOS1 = ZAUX(KK)+ zl*XXZ*0.5 ! ZAUX(2*KK) ! -ZL + (2*KK-1)*XXZ
        ZPOS2 = ZAUX(KK+1) ! -ZL + (2*KK)*XXZ

        NEGUAR=NE3
        DO II=1,NE
          NE3 = NE3+1

          DO JJ=1,NPOINT
            conect(NE3,JJ) = NEL(II,JJ)
          ENDDO

          DO JJ=1,NPOINT
            X = coor_x(NEL(II,JJ))
            Y = coor_y(NEL(II,JJ))
            Z = ZPOS1
            N = NPZ + 1
            NEW = NSEARCHN(N,X,Y,Z,NPZ,coor_x,coor_y,coor_z)
            IF(NEW.NE.N) THEN
              conect(NE3,NPOINT + JJ) = NEW
            ELSE
              NPZ=NPZ+1
              conect(NE3,NPOINT + JJ) = NPZ
              coor_x(NPZ) = coor_x(NEL(II,JJ))
              coor_y(NPZ) = coor_y(NEL(II,JJ))
              coor_z(NPZ) = Z
            ENDIF
          ENDDO

          DO JJ=1,NPOINT
            X = coor_x(NEL(II,JJ))
            Y = coor_y(NEL(II,JJ))
            Z = ZPOS2
            N = NPZ + 1
            NEW = NSEARCHN(N,X,Y,Z,NPZ,coor_x,coor_y,coor_z)
            IF(NEW.NE.N) THEN
              conect(NE3,2*NPOINT + JJ) = NEW
            ELSE
              NPZ=NPZ+1
              conect(NE3,2*NPOINT + JJ) = NPZ
              coor_x(NPZ) = coor_x(NEL(II,JJ))
              coor_y(NPZ) = coor_y(NEL(II,JJ))
              coor_z(NPZ) = Z
            ENDIF
          ENDDO

        ENDDO

        DO II=1,NE
          DO JJ=1,NPOINT
            NEL(II,JJ) = conect(NEGUAR+II,2*NPOINT + JJ)
          ENDDO
        ENDDO

      ENDDO

      NE=NE3
      NP = NPZ

      WRITE(6,*) 'TRIPLICO ',NE,NP
      WRITE(unit_cont,*) 'TRIPLICO ',NE,NP
      
!      ! grilla del plano a sacar
!      do jj=1,nzgrid2
!         nsum=0
!         do kk=1,np
!            if(coor_x(kk)==0.0)  then
!
!               if(coor_z(kk)>zgrid(jj)-0.0001 .and.
!     *                    coor_z(kk)<zgrid(jj)+0.0001 ) then
!                    nsum=nsum+1    
!                    grid(jj,nsum)=kk
!                  
!               endif
!
!            endif
!         enddo
!         
!         do ii=1,nsum-1
!            do kk=ii+1,nsum
!               if(coor_y(grid(jj,kk))<coor_y(grid(jj,ii))) then
!                  naux = grid(jj,ii)
!                  grid(jj,ii)=grid(jj,kk)
!                  grid(jj,kk)=naux
!               endif
!            enddo
!         enddo
!
!      enddo
!        
      nelements=ne     
      nnodes=np

      call verifico(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,
     *OrdenAlya)

      call cond_cont()
      write(6,*) 'condiciones de contorno!' 

      call sacomalla(nnodes,coor_x,coor_y,coor_z,conect,nelements,
     * nodpel,material,unit_malla)

      write(6,*) 'saco malla'

      deallocate(nel,xaux,yaux,yaux2,zaux,zaux2,mapa)
      
      end subroutine malla3d

      
      subroutine cond_cont()
      use def_variables
      use def_constantes
      implicit none
      ! local
      integer :: kk,jj,j,ii,mat_tierra,mat_potencial,mat_extra,el_neg,
     * el_pos
      double precision :: funm,funb,m,b,h1,h2,h3,h4,z,rmed,zmed,
     * rcil,rad,xmed,ymed,rlx,rly,rstx,rsty,rmed2,zmed2
       integer,allocatable :: aux1(:),aux2(:),aux3(:),aux4(:),
     *  aux5(:),aux6(:),lista_pos(:),lista_neg(:),neleaux(:)
      integer :: ncont1,ncont2,ncont3,ncont4,ncont5,ncont6,el_1,
     * el_2,el_3,el_4,el_5,el_6,elei,eles,elei2,eles2
      double precision:: dista,dis1,dis2,dis3,dis4,dis5,dis6,rad0,alfa
     
     
     
     
     
      ! determino condiciones de contorno!! nodos
      allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes),
     * aux1(nnodes),aux2(nnodes),aux3(nnodes),aux4(nnodes),
     * aux5(nnodes),aux6(nnodes),lista_pos(300),lista_neg(300))
      

      vec_tierra=-1
      vec_poten=-1
      vec_temper=-1
      nod_tierra=0
      nod_poten=0
      nod_temper=0

      if(problema=='ELECTRO') then
      
          h1 = 0.0
          h2 = alt1
          h3=alt1+alt2
          h4=alt1+alt2+alt3

          ! determino materiales
          mat_normal=0
          mat_tumor=0
          mat_aislante=0
          mat_electro=0
          mat_agua=0
          mat_tierra=0
          mat_potencial=0
          mat_extra=0

          do kk=1,nelements
       
             rmed = 0.0
             zmed = 0.0
             rcil = 0.0
             do jj=1,nodpel
               j=conect(kk,jj)
               rmed=rmed+ dsqrt(coor_x(j)*coor_x(j)+coor_y(j)*coor_y(j)+
     *                                          coor_z(j)*coor_z(j))
               rcil=rcil+ dsqrt(coor_x(j)*coor_x(j)+coor_y(j)*coor_y(j))
               zmed=zmed+coor_z(j)
             enddo
             rmed=rmed/real(nodpel)
             zmed=zmed/real(nodpel)
             rcil=rcil/real(nodpel)
         
             if(rmed>rad_zone1) then  ! zona normal
                material(kk)=2
                mat_normal=mat_normal+1
                if(zmed>=h4 .and. rcil <= rad4) then  ! aislante en zona normal
                     material(kk)=5
                     mat_aislante=mat_aislante+1
                     mat_normal=mat_normal-1
                endif
         
             else  ! zona tumoral
                material(kk)=4
                mat_tumor=mat_tumor+1

                ! electrodo 
                m = funm(rad2,h2,rad1,h1)
                b = funb(m,rad1,h1)
             
                if(zmed>=h1 .and. zmed < h2 .and. rcil<=(m*zmed+b)) then ! electrodo tierra
                     material(kk)=11
                     mat_electro=mat_electro+1
                     mat_tumor=mat_tumor-1
                     do jj=1,nodpel
                       
                           vec_tierra(conect(kk,jj)) = 1
                           nod_tierra=nod_tierra+1
                           aux1(nod_tierra)=kk
                     enddo
                endif

                m = funm(rad4,h4,rad3,h3)
                b = funb(m,rad3,h3)
            
                if(zmed>=h3 .and. zmed <= h4 ) then
                  if( rcil <= (m*zmed+b)*1.2) then ! electrodo potencial
                     material(kk)=10
                     mat_electro=mat_electro+1
                     mat_tumor=mat_tumor-1
                     do jj=1,nodpel
                       
                           vec_poten(conect(kk,jj)) = 1
                           nod_poten=nod_poten+1
                           aux2(nod_poten)=kk
                     enddo
                  endif
                endif


                ! aca clasifico aislante
                 m = funm(rad3,h3,rad2,h2)
                 b = funb(m,rad2,h2)
            
                if(zmed>=h2 .and. zmed < h3 .and.rcil <= (m*zmed+b)*1.2)
     *       then 
                     material(kk)=5
                     mat_aislante=mat_aislante+1
                     mat_tumor=mat_tumor-1
                elseif(zmed>h4 .and. rcil <= rad4) then
                     material(kk)=5
                     mat_aislante=mat_aislante+1
                     mat_tumor=mat_tumor-1
                endif

             endif

           enddo

      elseif(problema=='CEREBRO'.or. problema=='PAPA') then ! cerebro

      ! ahora determino elementos aislante y nodos del electrodo
          
          allocate(vecino(nelements))
          vecino=-1

          if(nop_malla==3) then

              !allocate(nod_electro_pos(AGUJAS/2),nod_electro_neg(AGUJAS/2))

              material=2

              ! encuentro el punto medio
              xmed = tam_zonex*0.5
              ymed = tam_zoney*0.5
              rstx = tam_zonex*0.5/real(ne0)
              rsty = tam_zoney*0.5/real(ne0)
              el_neg=0
              el_pos=0
         
              do kk=1,nelements
       
                 rmed = 0.0
                 zmed = 0.0
                 rlx = 0.0
                 rly = 0.0
                 do jj=1,nodpel
                   j=conect(kk,jj)
                   rmed=rmed+dsqrt(coor_x(j)*coor_x(j)+
     *              coor_y(j)*coor_y(j)+ coor_z(j)*coor_z(j))
                   rlx =rlx+coor_x(j)
                   rly =rly+coor_y(j)
                   zmed=zmed+coor_z(j)
                 enddo
                 rmed=rmed/real(nodpel)
                 zmed=zmed/real(nodpel)
                 rlx=rlx/real(nodpel)
                 rly=rly/real(nodpel)
          
                if(AGUJAS==6) then

                
                 ! electrodo 1 
                 if(rlx > xmed-poselex-rstx .and.rlx <xmed-poselex+rstx)
     *            then
                      if(rly > ymed+elec_sep*0.5-rsty .and. 
     *             rly < ymed+elec_sep*0.5+rsty) then
                           el_neg = el_neg+1
                           material(kk)=10   
                           lista_neg(el_neg)=kk
                      endif
                 endif

                 ! electrodo 2 
                 if(rlx > xmed-rstx .and. rlx < xmed+rstx)
     *          then
                      if(rly > ymed+elec_sep*0.5-rsty .and. 
     *             rly < ymed+elec_sep*0.5+rsty) then
                           el_neg = el_neg+1
                           material(kk)=10  
                           lista_neg(el_neg)=kk

                      endif
                 endif

                 ! electrodo 3 
                 if(rlx > xmed+poselex-rstx.and.rlx< xmed+poselex+rstx)
     *          then
                     if(rly > ymed+elec_sep*0.5-rsty .and. 
     *             rly < ymed+elec_sep*0.5+rsty) then
                           el_neg = el_neg+1
                           material(kk)=10   
                           lista_neg(el_neg)=kk

                      endif
                 endif

                 ! electrodo 4 
                 if(rlx > xmed-poselex-rstx.and. rlx< xmed-poselex+rstx)
     *          then
                      if(rly > ymed-elec_sep*0.5-rsty .and. 
     *             rly < ymed-elec_sep*0.5+rsty) then
                           el_pos = el_pos+1
                           material(kk)=11   
                           lista_pos(el_pos)=kk
                      endif
                 endif

                 ! electrodo 5 
                 if(rlx > xmed-rstx .and. rlx < xmed+rstx)  then
                      if(rly > ymed-elec_sep*0.5-rsty .and. 
     *             rly < ymed-elec_sep*0.5+rsty)    then
                           el_pos = el_pos+1
                           material(kk)=11   
                           lista_pos(el_pos)=kk
                            
                      endif
                 endif

                 ! electrodo 6 
                if(rlx > xmed+poselex-rstx.and.rlx < xmed+poselex+rstx)
     *          then
                      if(rly > ymed-elec_sep*0.5-rsty .and. 
     *             rly < ymed-elec_sep*0.5+rsty)    then
                           el_pos = el_pos+1
                           material(kk)=11   
                           lista_pos(el_pos)=kk
                      endif
                 endif
              else
                ! electrodo 2 
                 if(rlx > xmed-rstx .and. rlx < xmed+rstx)
     *          then
                      if(rly > ymed+elec_sep*0.5-rsty .and. 
     *                         rly < ymed+elec_sep*0.5+rsty) then
                           el_pos = el_pos+1
                           material(kk)=10   
                           lista_pos(el_pos)=kk
                      endif
                 endif

                 ! electrodo 5 
                 if(rlx > xmed-rstx .and. rlx < xmed+rstx)  then
                      if(rly > ymed-elec_sep*0.5-rsty .and. 
     *                         rly < ymed-elec_sep*0.5+rsty)    then
                           el_neg = el_neg+1
                           lista_neg(el_neg)=kk
                           material(kk)=11   
                      endif
                 endif


              endif

             enddo


             do kk=1,nelements
                 if(material(kk)==10) then ! tierra
                    
                    do jj=1,nodpel
                        vec_tierra(conect(kk,jj))=1

                    enddo
                  elseif(material(kk)==11) then ! potencial
                    do jj=1,nodpel
                        vec_poten(conect(kk,jj))=1
                    enddo
                  endif
                  if(material(kk)==10 .or. material(kk)==11) then
                    rmed = 0.0
                    zmed = 0.0
                    do jj=1,nodpel
                       j=conect(kk,jj)
                       rmed=rmed+dsqrt(coor_x(j)*coor_x(j)+
     *                  coor_y(j)*coor_y(j)+ coor_z(j)*coor_z(j))
                       zmed=zmed+coor_z(j)
                     enddo
                    rmed=rmed/real(nodpel)
                    zmed=zmed/real(nodpel)
                     
                     do jj=1,nelements
                        if(jj/=kk) then
                            
                            rmed2 = 0.0
                            zmed2 = 0.0
                            do ii=1,nodpel
                               j=conect(jj,ii)
                               rmed2=rmed2+dsqrt(coor_x(j)*coor_x(j)+
     *                         coor_y(j)*coor_y(j)+ coor_z(j)*coor_z(j))
                               zmed2=zmed2+coor_z(j)
                             enddo
                            rmed2=rmed2/real(nodpel)
                           zmed2=zmed2/real(nodpel)
                          if(zmed2<=zmed*1.01.and.zmed2>=zmed*0.99) then
           if(abs(rmed-rmed2)>=rstx*0.5 .and.abs(rmed-rmed2)<=rstx*1.01)
     *                                         then
                                 vecino(jj)=1
                            endif
                            endif
                        endif
                     enddo
                   endif
                   
             enddo



!             do kk=1,3
!                 elei = lista_neg(kk)
!                 eles = lista_neg(el_neg-kk)
!                 elei2 = lista_pos(kk)
!                 eles2 = lista_pos(el_pos-kk)
!                 do jj=1,nodpel
!                        vec_tierra(conect(elei,jj))=1
!                        vec_tierra(conect(eles,jj))=1
!                        vec_poten(conect(elei2,jj))=1
!                        vec_poten(conect(eles2,jj))=1
!                 enddo
!
!             enddo

             ! cuento  nodos
             do jj=1,nnodes
                 if(vec_tierra(jj)/=-1) then
                    nod_tierra=nod_tierra+1
                    aux1(nod_tierra)=jj
                 elseif(vec_poten(jj)/=-1) then
                    nod_poten=nod_poten+1
                    aux2(nod_poten)=jj
                 endif
             enddo 

                 
     
          allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
          write(unit_cc,* )'nodos tierra ', nod_tierra
      
            do kk=1,nod_tierra
               nnodtierra(kk)=aux1(kk)
    !           write(unit_cc,*) kk,nnodtierra(kk),coor_x(nnodtierra(kk)),
    !     *        coor_y(nnodtierra(kk)),coor_z(nnodtierra(kk))
            enddo

          write(unit_cc,* )'nodos poten ', nod_poten
            do kk=1,nod_poten
               nnodpoten(kk)=aux2(kk)
    !           write(unit_cc,*) kk,nnodpoten(kk), coor_x(nnodpoten(kk)),
    !     *        coor_y(nnodpoten(kk)),coor_z(nnodpoten(kk))
            enddo

          deallocate(aux1,aux2)
      
        elseif(nop_malla==5) then
          
          material=2
          aux1=0
          aux2=0
          aux3=0
          aux4=0
          aux5=0
          aux6=0

          ! aca considero electrodos puntuales!!
          ! encuentro el punto medio
          xmed = tam_zonex*0.5
          ymed = tam_zoney*0.5
          
           if(AGUJAS==6) then
            dis1=tam_zonex
            dis2=tam_zonex
            dis3=tam_zonex
            dis4=tam_zonex
            dis5=tam_zonex
            dis6=tam_zonex

            do kk=1,nnodes
              
              if(coor_z(kk)==0.0) then
                  rlx= coor_x(kk)
                  rly= coor_y(kk)
                                              !   1   2   3
                                              !   4   5   6
                  ! electrodo 1 
                  dista=sqrt( (rlx-(xmed-poselex))**2 + 
     *                                   (rly-(ymed+elec_sep*0.5))**2  )
                 if(dista<dis1) then
                     dis1=dista
                     el_1 = kk
                 endif

                 ! electrodo 2 
                  dista=sqrt( (rlx-(xmed))**2 + 
     *                                   (rly-(ymed+elec_sep*0.5))**2  )
                 if(dista<dis2) then
                     dis2=dista
                     el_2 = kk
                 endif

                 ! electrodo 3 
                  dista=sqrt( (rlx-(xmed+poselex))**2 + 
     *                                   (rly-(ymed+elec_sep*0.5))**2  )
                 if(dista<dis3) then
                     dis3=dista
                     el_3 = kk
                 endif

                 ! electrodo 4 
                  dista=sqrt( (rlx-(xmed-poselex))**2 + 
     *                                   (rly-(ymed-elec_sep*0.5))**2  )
                 if(dista<dis4) then
                     dis4=dista
                     el_4 = kk
                 endif

                 ! electrodo 5 
                  dista=sqrt( (rlx-(xmed))**2 + 
     *                                   (rly-(ymed-elec_sep*0.5))**2  )
                 if(dista<dis5) then
                     dis5=dista
                     el_5 = kk
                 endif

                 ! electrodo 6 
                  dista=sqrt( (rlx-(xmed+poselex))**2 + 
     *                                   (rly-(ymed-elec_sep*0.5))**2  )
                 if(dista<dis6) then
                     dis6=dista
                     el_6 = kk
                 endif
              endif
            enddo

          endif
      
      ! ahora selecciono los electrodos puntuales
        ncont1=0
        ncont2=0
        ncont3=0
        ncont4=0
        ncont5=0
        ncont6=0
        do kk=1,nnodes
           rlx= coor_x(kk)
           rly= coor_y(kk)

           ! electrodo 1 
           if(rlx==coor_x(el_1) .and. rly==coor_y(el_1) ) then
             ncont1=ncont1+1
             aux1(kk)=1
           endif
           ! electrodo 2 
           if(rlx==coor_x(el_2) .and. rly==coor_y(el_2) ) then
             ncont2=ncont2+1
             aux2(kk)=1
           endif
          ! electrodo 3 
           if(rlx==coor_x(el_3) .and. rly==coor_y(el_3) ) then
             ncont3=ncont3+1
             aux3(kk)=1
           endif
          ! electrodo 4 
           if(rlx==coor_x(el_4) .and. rly==coor_y(el_4) ) then
             ncont4=ncont4+1
             aux4(kk)=4
           endif
          ! electrodo 5 
           if(rlx==coor_x(el_5) .and. rly==coor_y(el_5) ) then
             ncont5=ncont5+1
             aux5(kk)=1
           endif
          ! electrodo 6 
           if(rlx==coor_x(el_6) .and. rly==coor_y(el_6) ) then
             ncont6=ncont6+1
             aux6(kk)=1
           endif


        enddo
      
        nod_tierra = ncont1+ncont2+ncont3
        nod_poten  = ncont4+ncont5+ncont6

        allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
        write(unit_cc,* )'nodos tierra ', nod_tierra
        nod_tierra=0
        nod_poten=0     
        do kk=1,nnodes
               if(aux1(kk)/=0 .or. aux2(kk)/=0 .or. aux3(kk)/=0) then
                 nod_tierra=nod_tierra+1
                 nnodtierra(nod_tierra)=kk
                 vec_tierra(kk)=1
                 write(unit_cc,*) 'T: ', kk,nnodtierra(nod_tierra),
     *    coor_x(nnodtierra(nod_tierra)),coor_y(nnodtierra(nod_tierra)),
     *              coor_z(nnodtierra(nod_tierra))
               elseif(aux4(kk)/=0 .or. aux5(kk)/=0 .or.aux6(kk)/=0) then
                 nod_poten=nod_poten+1
                 nnodpoten(nod_poten)=kk
                 vec_poten(kk)=1
                 write(unit_cc,*) 'P: ',kk,nnodpoten(nod_poten),
     *        coor_x(nnodpoten(nod_poten)),coor_y(nnodpoten(nod_poten)),
     *              coor_z(nnodpoten(nod_poten))
               endif

        enddo

      
       endif
     

      elseif(problema=='BERGUES3D') then
         
          material=4
          rad0=sqrt(tam_zonez*tam_zonez)

          do kk=1,nelements
             xmed=0.0
             ymed=0.0
             zmed=0.0
             
             do jj=1,nodpel
                  xmed=xmed + coor_x(conect(kk,jj))/nodpel
                  ymed=ymed + coor_y(conect(kk,jj))/nodpel
                  zmed=zmed + coor_z(conect(kk,jj))/nodpel
             enddo

             rad=sqrt(xmed*xmed+ymed*ymed+zmed*zmed)
             if(rad>rad0) material(kk)=2

          enddo
          



          
          N_electro_tot=5
          alfa= (tam_zonex *(1+alfazone))*0.5 /(2*ne0+1)  !0.033
          allocate(electro_nod(N_electro_tot,4*(2*ne0+1)),
     *                                     neleaux(N_electro_tot))

          neleaux=0
          do kk=1,nnodes
              
              if(coor_x(kk)==0.0 .and. coor_y(kk)==0.0 .and. 
     *                                 coor_z(kk)<=tam_zonez) then ! electro central
                  neleaux(1)=neleaux(1)+1
                  electro_nod(1,neleaux(1))=kk
              elseif(coor_x(kk)==0.0) then ! electrodos eje y
                   dista = sqrt(tam_zonez**2 - (tam_zoney*pos_cir1)**2)
                    
                   if( coor_y(kk)< tam_zoney*pos_cir1+alfa .and. 
     *                       coor_y(kk) > tam_zoney*pos_cir1-alfa .and.
     *                           coor_z(kk)<=dista )then 
                       neleaux(2)=neleaux(2)+1
                       electro_nod(2,neleaux(2))=kk
                   endif

                   dista = sqrt(tam_zonez**2 - (tam_zoney*pos_cir2)**2)

                   if( coor_y(kk)< tam_zoney*pos_cir2+alfa .and. 
     *                         coor_y(kk) > tam_zoney*pos_cir2-alfa.and.
     *                         coor_z(kk)<=dista ) then 
                       neleaux(3)=neleaux(3)+1
                       electro_nod(3,neleaux(3))=kk
                   endif
              elseif(coor_y(kk)==0.0) then ! electrodos eje x
                 
                   dista = sqrt(tam_zonez**2 - (tam_zonex*pos_cir1)**2)

                   if( coor_x(kk)< tam_zonex*pos_cir1+alfa .and. 
     *                         coor_x(kk) > tam_zonex*pos_cir1-alfa.and.
     *                           coor_z(kk)<=dista )then 
                       neleaux(4)=neleaux(4)+1
                       electro_nod(4,neleaux(4))=kk
                   endif

                   dista = sqrt(tam_zonez**2 - (tam_zonex*pos_cir2)**2)

                   if( coor_x(kk)< tam_zonex*pos_cir2+alfa .and. 
     *                       coor_x(kk) > tam_zonex*pos_cir2-alfa .and.
     *                       coor_z(kk)<=dista)then 
                       neleaux(5)=neleaux(5)+1
                       electro_nod(5,neleaux(5))=kk
                   endif

              endif
             
          enddo

          nod_poten = neleaux(1)+neleaux(3)+neleaux(5)
          nod_tierra  = neleaux(2)+neleaux(4)
          
        

      endif

      

      ! busco temper
      if(problema=='BERGUES3D') then

       
            do kk=1,nnodes
         
               if(coor_x(kk)>=(tam_zonex*(1+alfazone))*0.99) then
                  nod_temper=nod_temper+1
                  vec_temper(kk)=1
             elseif(coor_y(kk)>=(tam_zoney*(1+alfazone))*0.99) then
                  nod_temper=nod_temper+1
                  vec_temper(kk)=1
             elseif(coor_z(kk)>=(tam_zonez*(1+alfazone))*0.99) then
                  nod_temper=nod_temper+1
                  vec_temper(kk)=1
             endif

          enddo


      else

            do kk=1,nnodes
         
         if(coor_x(kk)==0.0 .or. coor_x(kk)>=tam_zonex*0.99) then
              nod_temper=nod_temper+1
              vec_temper(kk)=1
         elseif(coor_y(kk)==0.0 .or. coor_y(kk)>=tam_zoney*0.99) then
              nod_temper=nod_temper+1
              vec_temper(kk)=1
         elseif(coor_z(kk)==0.0 .or. coor_z(kk)>=tam_zonez*0.99) then
              nod_temper=nod_temper+1
              vec_temper(kk)=1
         endif

      enddo



      endif
      
      write(unit_cc,* )'nodos electro: ', N_electro_tot
      vec_poten=-1
      vec_tierra=-1
       
      do kk=1,N_electro_tot
           write(unit_cc,* )  neleaux(kk)
           do jj=1, neleaux(kk)
               ii=electro_nod(kk,jj)
               write(unit_cc,'(2i6,3E15.5)') jj,ii,coor_x(ii),coor_y(ii)
     *                                          ,coor_z(ii)

               if(kk==1 .or. kk==3 .or. kk==5) then
                    vec_poten(ii)=1
               else
                    vec_tierra(ii)=1
               endif


           enddo 
      enddo

      write(unit_cc,* )'nodos temperatura: ', nod_temper
      do kk=1,nnodes
         if( vec_temper(kk)==1) then
             write(unit_cc, '(i5,3E15.5)' ) kk,coor_x(kk),coor_y(kk),
     *                                                     coor_z(kk)    
         endif
      enddo
      

      ! chequeo materiales
      mat_normal=0
      mat_tumor=0
      mat_aislante=0
      mat_electro=0
      mat_agua=0
      mat_tierra=0
      mat_potencial=0
      mat_extra=0

      do kk=1,nelements
        if(material(kk) == 0 ) then
            mat_externo=mat_externo+1
         elseif(material(kk) == 1 ) then
            mat_agua=mat_agua+1
        elseif(material(kk) == 2 .or. material(kk) == 3 ) then
            mat_normal=mat_normal+1
        elseif(material(kk) == 4 ) then
            mat_tumor=mat_tumor+1
        elseif(material(kk) == 5 ) then
            mat_aislante=mat_aislante+1
        elseif(material(kk) == 10) then
            mat_electro=mat_electro+1
            mat_potencial = mat_potencial+1
        elseif(material(kk) == 11) then
            mat_electro=mat_electro+1
            mat_tierra = mat_tierra+1
        else
            mat_extra = mat_extra+1
        endif

      enddo





      write(unit_cc,* )' materiales ', nelements
      write(unit_cc,* )'      externo:   ', mat_externo
      write(unit_cc,* )'      agua   :   ', mat_agua
      write(unit_cc,* )'      tumor  :   ', mat_tumor
      write(unit_cc,* )'      normal :   ', mat_normal
      write(unit_cc,* )'      aislante:   ', mat_aislante
      write(unit_cc,* )'      electrodo:  ', mat_electro,mat_potencial,
     *                              mat_tierra
      write(unit_cc,* )'      extra:  ', mat_extra

      write(unit_cc,* )'  ******************************************  '
      write(unit_cc,* )'  ******************************************  '




      end subroutine cond_cont





      double precision function funm(h2,r2,h1,r1)
      double precision :: h2,r2,h1,r1

        funm = (h2-h1)/(r2-r1)

      end function funm


      double precision function funb(m,h,r)
      double precision :: m,h,r

        funb = h - m*r

      end function funb


      subroutine sacomalla2d(nnodes,coor_x,coor_y,conect,nelements,
     * nodpel,material,unit_malla)
      implicit none 
      integer :: nnodes,nelements,nodpel,unit_malla
      integer :: conect(nelements,nodpel),material(nelements)
      double precision :: coor_x(nnodes),coor_y(nnodes)
      ! local
      integer :: kk,jj


      write(unit_malla,*) nnodes
      do kk=1,nnodes
        write(unit_malla,'(i6,2e15.5)') kk,coor_x(kk),coor_y(kk)
      enddo 

      write(unit_malla,*) nelements
      do kk=1,nelements
         write(unit_malla,*) kk,material(kk),(conect(kk,jj),jj=1,nodpel)
      enddo      


      end subroutine sacomalla2d


      subroutine sacomalla(nnodes,coor_x,coor_y,coor_z,conect,nelements,
     * nodpel,material,unit_malla)
      implicit none 
      integer :: nnodes,nelements,nodpel,unit_malla
      integer :: conect(nelements,nodpel),material(nelements)
      double precision :: coor_x(nnodes),coor_y(nnodes),coor_z(nnodes)

      ! local
      integer :: kk,jj
      double precision :: xmed,ymed,zmed
      character(1)::coma=','
      write(unit_malla,*) nnodes
      do kk=1,nnodes
        write(unit_malla,'(i6,3e15.5)') kk,coor_x(kk),coor_y(kk),
     *                         coor_z(kk)
      enddo 

      write(unit_malla,*) nelements
      do kk=1,nelements
         write(unit_malla,*) kk,material(kk),(conect(kk,jj),jj=1,nodpel)
      enddo      

      do kk=1,nelements
         xmed=0.0
         ymed=0.0
         zmed=0.0
         do jj=1,nodpel
            xmed=xmed + coor_x(conect(kk,jj))/nodpel
            ymed=ymed + coor_y(conect(kk,jj))/nodpel
            zmed=zmed + coor_z(conect(kk,jj))/nodpel
         enddo
         write(1234,'(E15.5,A1,E15.5,A1,E15.5,A1,I3)')xmed,',',ymed,',',
     *                                zmed,',',material(kk)
      enddo      



      end subroutine sacomalla

      subroutine verifico(NP,XC,YC,ZC,NEL,NE,NOPE3,OrdenAlya)
      implicit none
      integer :: np,ne,nope3,OrdenAlya
      integer :: NEL(NE,NOPE3)
      double precision:: XC(NP),YC(NP),ZC(NP)
      ! local
      integer :: kk,jj,ii,naux(nope3)

!	do kk=1,NP-1
!	  DO JJ = KK+1,NP
!	      IF(XC(KK).EQ.XC(JJ).AND.YC(KK).EQ.YC(JJ).AND.
!     *                  ZC(KK).EQ.ZC(JJ))   THEN
!	       WRITE(6,*) 'WARNING! NODO: ',KK,' IGUAL AL ',JJ
!	   ENDIF
!        ENDDO
!	ENDDO
      

      if(OrdenAlya==1) then
       ! reordena Alya
        do kk=1,ne
          do jj=1,nope3
              naux(jj)=nel(kk,jj)
          enddo
         NEL(kk,1)  = naux(21)
         NEL(kk,2)  = naux(3)
         NEL(kk,3)  = naux(1)
         NEL(kk,4)  = naux(19)
         NEL(kk,5)  = naux(23)
         NEL(kk,6)  = naux(5)
         NEL(kk,7)  = naux(7)
         NEL(kk,8)  = naux(25)

         NEL(kk,9)  = naux(12)
         NEL(kk,10) = naux(2)
         NEL(kk,11) = naux(10)
         NEL(kk,12) = naux(20)
         NEL(kk,13) = naux(22)
         NEL(kk,14) = naux(4)
         NEL(kk,15) = naux(8)
         NEL(kk,16) = naux(26)
         
         NEL(kk,17) = naux(14)
         NEL(kk,18) = naux(6)
         NEL(kk,19) = naux(16)
         NEL(kk,20) = naux(24)
         
         NEL(kk,21) = naux(11)
         
         NEL(kk,22) = naux(13)
         NEL(kk,23) = naux(9)
         NEL(kk,24) = naux(17)
         NEL(kk,25) = naux(27)
         
          NEL(kk,26) = naux(15)
         
          NEL(kk,27) = naux(18)
       enddo
      endif	
      end subroutine verifico


      integer FUNCTION NSEARCHN(N,X,Y,Z,NP,XC,YC,ZC)
      implicit none
      integer :: np,n
      double precision :: x,y,z, XC(NP),YC(NP),ZC(NP)
      ! local
      integer :: kk
      
      DO KK=1,NP
         IF(XC(KK).EQ.X.AND.YC(KK).EQ.Y.AND.ZC(KK).EQ.Z) THEN
           NSEARCHN=KK
           RETURN
         ENDIF
      ENDDO
      NSEARCHN=N
      end function nsearchn
      


      
