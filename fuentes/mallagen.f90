subroutine mallo_cuadrado
use def_variables    
use def_constantes
implicit none
! local
    integer,allocatable :: naux(:),malla(:,:)
    double precision :: dx,dy
    integer :: kk,nx,ny,jj,nboundar 
    
    nodpel = nop_el
    nope = 2
    
    dx=Ladox/ne_lado
    dy=Ladoy/ne_lado
    
    nx=int4(Ladox/dx) !+1
    ny=int4(Ladoy/dy) !+1
    nnodes = (nx+1)*(ny+1)
    nelements=nx*ny

    allocate(conect(nelements,nop_el),material(nelements),coor_x(nnodes),coor_y(nnodes),naux(nnodes),malla(nnodes,nnodes))
    allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes))
    
         
    naux=-1
    vec_tierra=-1
    vec_poten=-1
    vec_temper=-1

    nnodes=0
    do kk=1,ny+1
       do jj=1,nx+1
          nnodes=nnodes+1
          coor_x(nnodes) = dx*(jj-1)
          coor_y(nnodes) = dy*(kk-1)
          malla(kk,jj)=nnodes
       enddo
    enddo

    nelements=0
    nboundar=0
    do kk=1,ny
       do jj=1,nx
           nelements=nelements+1
           conect(nelements,1)=malla(kk,jj)
           conect(nelements,2)=malla(kk,jj+1)
           conect(nelements,3)=malla(kk+1,jj+1)
           conect(nelements,4)=malla(kk+1,jj)
           material(nelements)=1 
       enddo
    enddo

    ! B.C
    nboundar=0
    
    ! exteriores
    do kk=1,nnodes
        if(coor_x(kk)<=0.001) then ! oeste
            nboundar=nboundar+1
            vec_temper(kk)=1
            naux(nboundar)=kk
        elseif(coor_x(kk)>=ladox*0.999) then ! este
            nboundar=nboundar+1
            vec_temper(kk)=1
            naux(nboundar)=kk
        elseif(coor_y(kk)<=0.001) then ! sur
            nboundar=nboundar+1
            vec_temper(kk)=1
            naux(nboundar)=kk
        elseif(coor_y(kk)>=ladoy*0.999) then ! norte
            nboundar=nboundar+1
            vec_temper(kk)=1
            naux(nboundar)=kk
        endif
    enddo
    
    nod_temper=nboundar
    allocate(nnodtemper(nboundar))
    do kk=1,nboundar
        nnodtemper(kk)=naux(kk)
    enddo
    naux=-1
    
    pos_tierra = ladoy*0.5-dis_ele*0.5
    pos_pot = ladoy*0.5+dis_ele*0.5
    ancho_inf=ladox*0.5-ladox*0.25
    ancho_sup=ladox*0.5+ladox*0.25
    
    ! potencial
    nboundar=0
    do kk=1,nnodes
        if(coor_x(kk)>= ancho_inf .and. coor_x(kk)<= ancho_sup) then
             if(coor_y(kk)>=pos_pot*0.99 .and. coor_y(kk)<=pos_pot*1.001)    then ! potencial
                nboundar=nboundar+1
                vec_poten(kk)=1
                naux(nboundar)=kk
             endif
             
        endif
    enddo
    nod_poten=nboundar
    allocate(nnodpoten(nboundar))
    do kk=1,nboundar
        nnodpoten(kk)=naux(kk)
    enddo
    naux=-1

    ! TIERRA
    nboundar=0
    do kk=1,nnodes
        if(coor_x(kk)>= ancho_inf .and. coor_x(kk)<= ancho_sup) then
             if(coor_y(kk)>=pos_tierra*0.99 .and. coor_y(kk)<=pos_tierra*1.001)    then ! tierra
                nboundar=nboundar+1
                vec_tierra(kk)=1
                naux(nboundar)=kk
             endif
             
        endif
    enddo
    nod_tierra=nboundar
    allocate(nnodtierra(nboundar))
    do kk=1,nboundar
        nnodtierra(kk)=naux(kk)
    enddo

    
    ! dimensiono
    do kk=1,nnodes
        coor_x(kk) = coor_x(kk)* factorx+desplax
        coor_y(kk) = coor_y(kk)* factory+desplay
    enddo


    deallocate(malla,naux)
    
    call sacomalla2d(nnodes,coor_x,coor_y,conect,nelements,nop_el,material,unit_malla)
    
    ! imrpimo B.C
    write(unit_malla,*) nod_temper
    do kk=1,nod_temper
        write(unit_malla,*)  nnodtemper(kk)
    enddo
    write(unit_malla,*) nod_poten
    do kk=1,nod_poten
        write(unit_malla,*)  nnodpoten(kk)
    enddo
    write(unit_malla,*) nod_tierra
    do kk=1,nod_tierra
        write(unit_malla,*)  nnodtierra(kk)
    enddo
    
    
    end subroutine mallo_cuadrado
    
    
    
    
    subroutine malla_raton
      use def_variables
      use def_constantes
      use def_bio
      implicit none
      ! local
      double precision,allocatable :: xaux(:),yaux(:),zaux(:),yaux2(:), zaux2(:)
      integer, allocatable:: mapa(:,:,:),nel(:,:)
      integer :: np_ant,np,NELE,NPZ2,NEGUAR,NE3,N,NEW,nsearchn,NPOINT,ne_ini
      integer ::i,j,k,ii,jj,kk,nsum,naux,nzgrid2,NPME,nele2,nfactor
      double precision :: xx,xxz,VERMALLAX,VERMALLAy,VERMALLAZ,RESTO,xl,yl,zl,ZPOS1,ZPOS2,rad,x,y,z,aux(9),xplano,xzona,yzona,zzona
      
      integer :: nex,ney,nez,npx,npy,npz,NPmed,Np0

      XL = tam_zonex
      YL = tam_zoney
      ZL = tam_zonez
                
      XX = tamano_elem

      NEx=XL/(tamano_elem)
      NEy=YL/(tamano_elem)
      NEz=ZL/(tamano_elem)

      NPx=(nope-1)*NEx+1
      NPy=(nope-1)*NEy+1
      NPz=(nope-1)*NEz+1
            
      area_trans= YL*ZL

      ! lo total del dominio
      nnodes= npx*npy*npz
      nelements=nex*ney*nez

      allocate(xaux(npx),yaux(npy), mapa(3,(nex+1)*(ney+1),(nex+1)*(ney+1)), nel(nex*ney,nope*nope) )
      
      allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),grad_x(nnodes),grad_y(nnodes),grad_z(nnodes))
      allocate(conect(nelements,nodpel),material(nelements))
      allocate(gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements),jcurrent(nelements),quemado(nelements))
      allocate(tempera(nnodes),tempera_ant(nnodes))

      gradxel_x=0.0
      gradxel_y=0.0
      gradxel_z=0.0
      quemado=0

      tempera=Tblod
      tempera_ant=Tblod

      DO I=1,NPx
        coor_x(I) = (I-1)*XX*0.5
        MAPA(1,1,i) = i
        coor_y(i) = 0.0
      enddo
      DO I=1,NPy
        YAUX(I) = (I-1)*XX*0.5
      ENDDO

      WRITE(unit_cont,*) 'MALLO LAS LINEAS',NPx,npy
      

      DO jj=2,NPy
        DO kk=1,NPx
          coor_x(NPx*(jj-1)+kk) = coor_x(kk)
          coor_y(NPx*(jj-1)+kk) = YAUX(jj)
          MAPA(1,jj,kk)=NPx*(jj-1) + kk
        ENDDO
      ENDDO
      
      NPmed = NPx*NPy

      WRITE(unit_cont,*) 'MAPEO ',NPmed

      DO J=1,npmed
        coor_z(J)=0.0
      ENDDO

      NELE=0
      DO jj=1,NEy
         DO kk=1,NEx
           
            NELE=NELE+1
            NEL(NELE,1) = MAPA(1,jj, kk)
            NEL(NELE,2) = MAPA(1,jj, kk+2)
            NEL(NELE,3) = MAPA(1,jj+2,kk+2)
            NEL(NELE,4) = MAPA(1,jj+2,kk)
            NEL(NELE,5) = MAPA(1,jj, kk+1)
            NEL(NELE,6) = MAPA(1,jj+1, kk+2)
            NEL(NELE,7) = MAPA(1,jj+2,kk+1)
            NEL(NELE,8) = MAPA(1,jj+1,kk)
            NEL(NELE,9) = MAPA(1,jj+1, kk+1)

         ENDDO
      ENDDO

      NE=NELE
      np=npmed

      open(unit=1111,file='mallado_rat.dat')
      write(1111,*) np
      do kk=1,np
         write(1111,*) kk,coor_x(kk),coor_y(kk)
      enddo

      write(1111,*) ne
      do kk=1,ne
         write(1111,*) kk,(nel(kk,jj),jj=1,nope*nope)
      enddo
      
      close(1111)

      close(unit_2D)

! direccion z
     ne3=0
     
     DO KK=1,NPZ-2,2
  
        ZPOS1 = kk*XX*0.5
        ZPOS2 = (kk+1)*XX*0.5
        ! creo nodos nuevos
        Np0=0
        DO II=1,Npy
           DO jj=1,Npx
               Np0=np0+1
               np=np+1
               coor_x(np)=coor_x(np0)
               coor_y(np)=coor_y(np0)
               coor_z(np)=zpos1
               mapa(2,ii,jj)=np
           enddo
        enddo
        np0=0
        DO II=1,Npy
           DO jj=1,Npx
               Np0=np0+1
               np=np+1
               coor_x(np)=coor_x(np0)
               coor_y(np)=coor_y(np0)
               coor_z(np)=zpos2
               mapa(3,ii,jj)=np
           enddo
        enddo
                   
        ! creo elementos nuevos
        do ii=1,ney
           do jj=1,nex
             ne3=ne3+1  
             conect(NE3,1) = mapa(1,2*ii-1,2*jj-1)
             conect(NE3,2) = mapa(1,2*ii-1,2*jj-1+2)
             conect(NE3,3) = mapa(1,2*ii-1+2,2*jj-1+2)
             conect(NE3,4) = mapa(1,2*ii-1+2,2*jj-1)
             conect(NE3,5) = mapa(3,2*ii-1,2*jj-1)
             conect(NE3,6) = mapa(3,2*ii-1,2*jj-1+2)
             conect(NE3,7) = mapa(3,2*ii-1+2,2*jj-1+2)
             conect(NE3,8) = mapa(3,2*ii-1+2,2*jj-1)
             
             conect(NE3,9) = mapa(1,2*ii-1,2*jj-1+1)
             conect(NE3,10) = mapa(1,2*ii-1+1,2*jj-1+2)
             conect(NE3,11) = mapa(1,2*ii-1+2,2*jj-1+1)
             conect(NE3,12) = mapa(1,2*ii-1+1,2*jj-1)
             
             conect(NE3,13) = mapa(2,2*ii-1,2*jj-1)
             conect(NE3,14) = mapa(2,2*ii-1,2*jj-1+2)
             conect(NE3,15) = mapa(2,2*ii-1+2,2*jj-1+2)
             conect(NE3,16) = mapa(2,2*ii-1+2,2*jj-1)
             
             conect(NE3,17) = mapa(3,2*ii-1,2*jj-1+1)
             conect(NE3,18) = mapa(3,2*ii-1+1,2*jj-1+2)
             conect(NE3,19) = mapa(3,2*ii-1+2,2*jj-1+1)
             conect(NE3,20) = mapa(3,2*ii-1+1,2*jj-1)

             conect(NE3,21) = mapa(1,2*ii-1+1,2*jj-1+1)
             conect(NE3,22) = mapa(2,2*ii-1,2*jj-1+1)
             conect(NE3,23) = mapa(2,2*ii-1+1,2*jj-1+2)
             conect(NE3,24) = mapa(2,2*ii-1+2,2*jj-1+1)
             conect(NE3,25) = mapa(2,2*ii-1+1,2*jj-1)
             conect(NE3,26) = mapa(3,2*ii-1+1,2*jj-1+1)
             conect(NE3,27) = mapa(2,2*ii-1+1,2*jj-1+1)
           enddo
        enddo
        do ii=1,2*ney+1
           do jj=1,2*nex+1
              mapa(1,ii,jj) = mapa(3,ii,jj)
           enddo
        enddo

   ENDDO
   
   NE=NE3
   
   WRITE(6,*) 'TRIPLICO ',NE,NP
   WRITE(unit_cont,*) 'TRIPLICO ',NE,NP
      
      
         
   nelements=ne     
   nnodes=np

      !call verifico(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,OrdenAlya)

      call cond_contgen_raton(NEX,NEY)
      !call cond_contgen_test(NEX,NEY)
      write(6,*) 'condiciones de contorno!' 

      call sacomalla(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,material,unit_malla)

      write(6,*) 'saco malla'

      deallocate(nel,xaux,yaux,mapa)


      end subroutine malla_raton



subroutine mallagen
      use def_variables
      use def_constantes
      use def_bio
      implicit none
      ! local
      double precision,allocatable :: xaux(:),yaux(:),zaux(:),yaux2(:), zaux2(:)
      integer, allocatable:: mapa(:,:),nel(:,:)
      integer :: nez,np_ant,np,NELE,NPZ2,NEGUAR,NE3,N,NEW,NPZ,nsearchn,NPOINT,ne_ini
      integer ::i,j,k,ii,jj,kk,nsum,naux,nzgrid2,NPME,nele2,nfactor
      double precision :: xx,xxz,VERMALLAX,VERMALLAy,VERMALLAZ,RESTO,xl,yl,zl,ZPOS1,ZPOS2,rad,x,y,z,aux(9),xplano,xzona,yzona,zzona
      

      ncerx0=0
      ncery0=0
      ncerz0=0

      XL = tam_zonex
      YL = tam_zoney
      ZL = tam_zonez


      xzona=XL
      yzona=YL
      zzona=ZL

      xplano=XL*0.5

      NEZ= NE
      ne_ini=ne
                 
      XX = 1.0/real(NE)

      NP=NE+1
      NP_ANT = NP


      ! lo total del dominio
      nnodes= (ne+1)*(ne+1)*(ne+1)
      nelements=ne*ne*ne

      allocate(xaux(np),yaux(np),yaux2(2*np), mapa((ne+1)*(ne+1),(ne+1)*(ne+1)), nel(ne*ne,nope*nope))
      
      allocate(coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),grad_x(nnodes),grad_y(nnodes),grad_z(nnodes))
      allocate(conect(nelements,nodpel),material(nelements))
      allocate(gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements),jcurrent(nelements))
      allocate(tempera(nnodes),tempera_ant(nnodes))

      gradxel_x=0.0
      gradxel_y=0.0
      gradxel_z=0.0

      tempera=Tblod
      tempera_ant=Tblod

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
            NEL(NELE,1)=MAPA(KK, JJ)
            NEL(NELE,2)=MAPA(KK+1, JJ)
            NEL(NELE,3)= MAPA(KK+1,JJ+1)
            NEL(NELE,4)= MAPA(KK,JJ+1)
         ENDDO
      ENDDO


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

! direccion z
      NPZ2 = NEZ+1
      XXZ = 1.0/(NEZ)

      allocate(zaux(npz2),zaux2(npz2))

      DO I=1,NPZ2
        ZAUX(I) = (I-1)*XXZ
        ZAUX(I) = ZAUX(I)*zl
      ENDDO

      nzgrid= npz2
       allocate(Zgrid(nzgrid),grid(nzgrid,nzgrid) )
        DO KK=1,NPZ2
            WRITE(unit_cont,*) KK,ZAUX(KK)
            zgrid(kk)=ZAUX(KK)
        ENDDO

      NE3 = 0
      NPZ = NP
      NPOINT = NOPE*NOPE

      DO J=1,NPZ
        coor_z(J)=0.0
      ENDDO

      DO KK=1,NEZ
  
        ZPOS1 = ZAUX(KK+1)
        
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

        ENDDO

        DO II=1,NE
          DO JJ=1,NPOINT
            NEL(II,JJ) = conect(NEGUAR+II,NPOINT + JJ)
          ENDDO
        ENDDO

      ENDDO


      NE=NE3
      NP = NPZ

      WRITE(6,*) 'TRIPLICO ',NE,NP
      WRITE(unit_cont,*) 'TRIPLICO ',NE,NP
      
      
      ! grilla del plano a sacar
      do jj=1,nzgrid
         nsum=0
         do kk=1,np
            if(coor_x(kk)==xplano)  then

               if(coor_z(kk)>zgrid(jj)-0.0001 .and. coor_z(kk)<zgrid(jj)+0.0001 ) then
                    nsum=nsum+1    
                    grid(jj,nsum)=kk
                  
               endif

            endif
         enddo
         
         do ii=1,nsum-1
            do kk=ii+1,nsum
               if(coor_y(grid(jj,kk))<coor_y(grid(jj,ii))) then
                  naux = grid(jj,ii)
                  grid(jj,ii)=grid(jj,kk)
                  grid(jj,kk)=naux
               endif
            enddo
         enddo

      enddo
      
          
      
         
      nelements=ne     
      nnodes=np

      call verifico(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,OrdenAlya)

      call cond_contgen(ne_ini,xzona,yzona,zzona)
      write(6,*) 'condiciones de contorno!' 

      call sacomalla(nnodes,coor_x,coor_y,coor_z,conect,nelements,nodpel,material,unit_malla)

      write(6,*) 'saco malla'

      deallocate(nel,xaux,yaux,yaux2,zaux,zaux2,mapa)


      end subroutine mallagen


      
      subroutine cond_contgen(ne_ini,xzona,yzona,zzona)
      use def_variables
      use def_constantes
      implicit none
      double precision :: xzona,yzona,zzona
      integer ::ne_ini
      ! local
      integer :: kk,jj,j,mat_tierra,mat_potencial,mat_extra,mat_tot
      double precision :: funm,funb,m,b,h1,h2,h3,h4,z,rmed,xmed,ymed,zmed,rcil,rad
      integer,allocatable :: aux1(:),aux2(:)
      
      ! determino condiciones de contorno!! nodos
      allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes),aux1(nnodes),aux2(nnodes))
      
      vec_tierra=-1
      vec_poten=-1
      vec_temper=-1
      nod_tierra=0
      nod_poten=0

      
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
         xmed = 0.0
         ymed = 0.0
         zmed = 0.0
         do jj=1,nodpel
           j=conect(kk,jj)
           xmed =xmed+ dsqrt(coor_x(j)*coor_x(j))
           ymed =ymed+ dsqrt(coor_y(j)*coor_y(j))
           zmed =zmed+ dsqrt(coor_z(j)*coor_z(j))
         enddo
         xmed=xmed/real(nodpel)
         ymed=ymed/real(nodpel)
         zmed=zmed/real(nodpel)
         if(xmed<=xzona*0.2 .or. xmed>=xzona*0.8) then  ! zona externa
            material(kk)=2
             mat_normal=mat_normal+1
         elseif(ymed<=yzona*0.2 .or. ymed>=yzona*0.8) then  ! zona externa
            material(kk)=2
            mat_normal=mat_normal+1
         elseif(zmed<=zzona*0.2 .or. zmed>=zzona*0.8) then  ! zona externa
            material(kk)=2
             mat_normal=mat_normal+1
         elseif( (xmed>xzona*0.2 .and. xmed<=xzona*0.4) .or. (xmed>=xzona*0.6 .and. xmed<xzona*0.8 )) then  ! zona normal
!           if( (ymed>2.5 .and. ymed<=7.5) .or. (ymed>=12.5 .and. ymed<17.5 )) then  ! zona normal
!              if( (zmed>2.5 .and. zmed<=7.5) .or. (zmed>=12.5 .and. zmed<17.5 )) then  ! zona normal
           if( (ymed>yzona*0.2 .and. ymed<=yzona*0.8) ) then  ! zona normal
              if( (zmed>zzona*0.2 .and. zmed<=zzona*0.8)) then  ! zona normal
                  material(kk)=4
                  mat_tumor=mat_tumor+1
     !            mat_normal=mat_normal+1
              endif
           endif            
         elseif( (ymed>yzona*0.2 .and. ymed<=yzona*0.4) .or. (ymed>=yzona*0.6 .and. ymed<yzona*0.8 )) then  ! zona normal
            if( (xmed>xzona*0.4 .and.xmed<xzona*0.6) ) then  ! zona normal
              if( (zmed>zzona*0.2 .and. zmed<=zzona*0.8)) then  ! zona normal
                  material(kk)=4
!                  mat_normal=mat_normal+1
                   mat_tumor=mat_tumor+1
       endif
           endif     
        elseif( (zmed>zzona*0.2 .and. zmed<=zzona*0.4) .or. (zmed>=zzona*0.6 .and. zmed<zzona*0.8 )) then  ! zona normal
           if( (xmed>=xzona*0.4 .and. xmed<=xzona*0.6) ) then  ! zona normal
              if( (ymed>=yzona*0.4 .and. ymed<=yzona*0.6)) then  ! zona normal
                  material(kk)=4
                  mat_tumor=mat_tumor+1
      !            mat_normal=mat_normal+1
              endif
           endif            
       else  ! zona tumoral
            material(kk)=4
            mat_tumor=mat_tumor+1

         endif


     ! electrodos!!

         if(xmed<=xzona/ne_ini) then  
            if(ymed>=yzona*0.3 .and. ymed<=yzona*0.7) then  ! zona tierra
                material(kk)=10
                mat_tierra=mat_tierra+1
            endif            
         elseif(xmed>=xzona - xzona/ne_ini) then 
            if(ymed>=yzona*0.3 .and. ymed<=yzona*0.7) then  ! zona potencial
                material(kk)=11
                mat_potencial=mat_potencial+1
            endif            
        endif 


     enddo



     ! electrodos!!
      do kk=1,nelements
         if(material(kk)==10) then
            do jj=1,nodpel
               nod_tierra=nod_tierra+1
               vec_tierra(conect(kk,jj))=1
            enddo
         elseif(material(kk)==11) then
            do jj=1,nodpel
               nod_poten=nod_poten+1
               vec_poten(conect(kk,jj))=1
            enddo
         endif         
      enddo
     

      allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
      write(unit_cc,* )'nodos tierra ', nod_tierra
     nod_tierra=0
     do kk=1,nnodes
       
         if( vec_tierra(kk)==1) then  ! zona tierra
            nod_tierra=nod_tierra+1
            nnodtierra(nod_tierra)=kk
            write(unit_cc,*) nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)),coor_z(nnodtierra(nod_tierra))
         endif
     enddo

     
      write(unit_cc,* )'nodos poten ', nod_poten
     nod_poten=0
     do kk=1,nnodes
       
         if(vec_poten(kk)==1) then  ! zona poten
            nod_poten=nod_poten+1
            nnodpoten(nod_poten)=kk
            write(unit_cc,*) nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten)),coor_z(nnodpoten(nod_poten))
         endif
     enddo



     ! condiciones temper en el exterior
     nod_temper=0.0
     do kk=1,nnodes
       
         if(coor_x(kk)==0.0 .or. coor_x(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
         if(coor_y(kk)==0.0 .or. coor_y(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
         if(coor_z(kk)==0.0 .or. coor_z(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
     enddo
     allocate(nnodtemper(nod_temper))

     nod_temper=0.0
      do kk=1,nnodes
       
         if(coor_x(kk)==0.0 .or. coor_x(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,*) nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
         endif

         if(coor_y(kk)==0.0 .or. coor_y(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,*) nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
         endif
         if(coor_z(kk)==0.0 .or. coor_z(kk)>=xzona*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,*) nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
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
      write(unit_cc,* )'      electrodo:  ', mat_electro,mat_potencial,  mat_tierra
      mat_tot = mat_electro+mat_aislante+mat_tumor+mat_normal
      write(unit_cc,* )'      extra:  ', mat_extra,mat_tot

      write(unit_cc,* )'  ******************************************  '
      write(unit_cc,* )'  ******************************************  '




      end subroutine cond_contgen

subroutine cond_contgen_raton(NEX,NEY)
      use def_variables
      use def_constantes
      implicit none
      integer :: NEX,NEY
      ! local
      integer :: kk,jj,j,mat_tierra,mat_potencial,mat_extra,mat_tot,cuento2d
      double precision :: funm,funb,m,b,h1,h2,h3,h4,z,rmed,xmed,ymed,zmed,rcil,rad
      integer,allocatable :: aux1(:),aux2(:)
      
      ! determino condiciones de contorno!! nodos
      allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes),aux1(nnodes),aux2(nnodes))
      allocate(grilla2d(nelements))

      grilla2d=0
      cuento2d=0
       
      vec_tierra=-1
      vec_poten=-1
      vec_temper=-1
      nod_tierra=0
      nod_poten=0

      
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
        
         material(kk) = 4
         xmed = 0.0
         ymed = 0.0
         zmed = 0.0
         do jj=1,nodpel
           j=conect(kk,jj)
           xmed =xmed+ dsqrt(coor_x(j)*coor_x(j))
           ymed =ymed+ dsqrt(coor_y(j)*coor_y(j))
           zmed =zmed+ dsqrt(coor_z(j)*coor_z(j))
         enddo
         xmed=xmed/real(nodpel)
         ymed=ymed/real(nodpel)
         zmed=zmed/real(nodpel)
         if(zmed<=tamano_elem*elec_alt ) then  
           if(ymed<tamano_elem*((NEY-1)*0.5+1) .and. ymed>tamano_elem*((NEY-1)*0.5) ) then  
              if(xmed<tamano_elem*((NEX-2)*0.25+1) .and. xmed>tamano_elem*((NEX-2)*0.25) ) then  ! tierra
                material(kk)=10
                mat_tierra=mat_tierra+1
              elseif(xmed<tamano_elem*((NEX-2)*0.25+elec_sep+2) .and. xmed>tamano_elem*((NEX-2)*0.25+elec_sep+1) ) then  ! potencial
                material(kk)=11
                mat_potencial=mat_potencial+1
              endif
          endif
        endif

          if(xmed>=tamano_elem*(NEx*0.5) .and. xmed<=tamano_elem*(NEx*0.5+1) ) then
            cuento2d=cuento2d+1
            grilla2d(kk)=1
        endif

     enddo



     ! electrodos!!
      do kk=1,nelements

         if(material(kk)==10) then
            do jj=1,nodpel
               nod_tierra=nod_tierra+1
               vec_tierra(conect(kk,jj))=1
            enddo
         elseif(material(kk)==11) then
            do jj=1,nodpel
               nod_poten=nod_poten+1
               vec_poten(conect(kk,jj))=1
            enddo
         endif         


      enddo
     
     nod_tierra=0
     do kk=1,nnodes
         if( vec_tierra(kk)==1) then  ! zona tierra
            nod_tierra=nod_tierra+1
         endif
     enddo
     nod_poten=0
     do kk=1,nnodes
         if(vec_poten(kk)==1) then  ! zona poten
            nod_poten=nod_poten+1
         endif
     enddo

      allocate(nnodtierra(nod_tierra),nnodpoten(nod_poten))
   
   
      write(unit_cc,* )'nodos tierra ', nod_tierra
     nod_tierra=0
     do kk=1,nnodes
       
         if( vec_tierra(kk)==1) then  ! zona tierra
            nod_tierra=nod_tierra+1
            nnodtierra(nod_tierra)=kk
            write(unit_cc,'(2i7,3E15.5)') nod_tierra,nnodtierra(nod_tierra),coor_x(nnodtierra(nod_tierra)), coor_y(nnodtierra(nod_tierra)),coor_z(nnodtierra(nod_tierra))
         endif
     enddo

     
      write(unit_cc,* )'nodos poten ', nod_poten
     nod_poten=0
     do kk=1,nnodes
       
         if(vec_poten(kk)==1) then  ! zona poten
            nod_poten=nod_poten+1
            nnodpoten(nod_poten)=kk
            write(unit_cc,'(2i7,3E15.5)') nod_poten,nnodpoten(nod_poten),coor_x(nnodpoten(nod_poten)), coor_y(nnodpoten(nod_poten)),coor_z(nnodpoten(nod_poten))
         endif
     enddo



     ! condiciones temper en el exterior
     nod_temper=0.0
     do kk=1,nnodes
       
         if(coor_x(kk)==0.0 .or. coor_x(kk)>=tam_zonex*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
         if(coor_y(kk)==0.0 .or. coor_y(kk)>=tam_zoney*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
         if(coor_z(kk)==0.0 .or. coor_z(kk)>=tam_zonez*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            vec_temper(kk)=1
         endif
     enddo
     allocate(nnodtemper(nod_temper))

     nod_temper=0.0
      do kk=1,nnodes
       
         if(coor_x(kk)==0.0 .or. coor_x(kk)>=tam_zonex*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,'(2i7,3E15.5)') nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
         endif

         if(coor_y(kk)==0.0 .or. coor_y(kk)>=tam_zoney*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,'(2i7,3E15.5)') nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
         endif
         if(coor_z(kk)==0.0 .or. coor_z(kk)>=tam_zonez*0.99) then  ! zona exterior 
            nod_temper=nod_temper+1
            nnodtemper(nod_temper)=kk
            write(unit_cc,'(2i7,3E15.5)') nod_temper,nnodtemper(nod_temper),coor_x(nnodtemper(nod_temper)), coor_y(nnodtemper(nod_temper)),coor_z(nnodtemper(nod_temper))
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
      write(unit_cc,* )'      electrodo:  ', mat_electro,mat_potencial,  mat_tierra
      mat_tot = mat_electro+mat_aislante+mat_tumor+mat_normal
      write(unit_cc,* )'      extra:  ', mat_extra,mat_tot

      write(unit_cc,* )'  ******************************************  '
      write(unit_cc,* )'  ******************************************  '




      end subroutine cond_contgen_raton



subroutine cond_contgen_test(NEX,NEY)
      use def_variables
      use def_constantes
      implicit none
      integer :: NEX,NEY
      ! local
      integer :: kk,jj,j,mat_tierra,mat_potencial,mat_extra,mat_tot
      double precision :: funm,funb,m,b,h1,h2,h3,h4,z,rmed,xmed,ymed,zmed,rcil,rad
      integer,allocatable :: aux1(:),aux2(:)
      
      ! determino condiciones de contorno!! nodos
      allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_temper(nnodes),aux1(nnodes),aux2(nnodes))
     
      
      vec_tierra=-1
      vec_poten=-1
      vec_temper=-1
      nod_tierra=0
      nod_poten=0

      
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
        
         material(kk) = 4
         xmed = 0.0
         ymed = 0.0
         zmed = 0.0
         do jj=1,nodpel
           j=conect(kk,jj)
           xmed =xmed+ dsqrt(coor_x(j)*coor_x(j))
           ymed =ymed+ dsqrt(coor_y(j)*coor_y(j))
           zmed =zmed+ dsqrt(coor_z(j)*coor_z(j))
         enddo
         xmed=xmed/real(nodpel)
         ymed=ymed/real(nodpel)
         zmed=zmed/real(nodpel)
         if(zmed<=tamano_elem*elec_alt ) then  
           if(ymed<tamano_elem*((NEY-1)*0.5+1) .and. ymed>tamano_elem*((NEY-1)*0.5) ) then  
              if(xmed<tamano_elem*((NEX-2)*0.25+1) .and. xmed>tamano_elem*((NEX-2)*0.25) ) then  ! tierra
                !material(kk)=10
                !mat_tierra=mat_tierra+1
              elseif(xmed<tamano_elem*((NEX-2)*0.25+elec_sep+2) .and. xmed>tamano_elem*((NEX-2)*0.25+elec_sep+1) ) then  ! potencial
                !material(kk)=11
                !mat_potencial=mat_potencial+1
              endif
          endif
        endif

      

     enddo



     
     ! condiciones temper test
     nod_tierra=0.0
     nod_poten=0.0

     do kk=1,nnodes
       
         if(coor_x(kk)==0.0 ) then  ! piso
            nod_tierra=nod_tierra+1
            vec_tierra(kk)=1
         endif
         if(coor_x(kk)>=tam_zonex*0.99) then  ! poten
            nod_poten=nod_poten+1
            vec_poten(kk)=1
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
      write(unit_cc,* )'      electrodo:  ', mat_electro,mat_potencial,  mat_tierra
      mat_tot = mat_electro+mat_aislante+mat_tumor+mat_normal
      write(unit_cc,* )'      extra:  ', mat_extra,mat_tot

      write(unit_cc,* )'  ******************************************  '
      write(unit_cc,* )'  ******************************************  '




      end subroutine cond_contgen_test
