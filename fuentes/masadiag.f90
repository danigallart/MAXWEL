subroutine masadiag2d()
use def_solver
use def_variables
use def_constantes
implicit none
double precision x(nodpel),y(nodpel),z(nodpel),deriv(3,nodpel,nodpel),gpdet,gpvol

integer ns(nodpel),ipoin,inode,jnode,ielem,kk,j,k,pdime,pnode,mdime,igaus,ilocs,jlocs
  !    calculamos la matriz de masa diagonal usando una regla cerrada de integracion

double precision ::  weigc(nodpel),posgl(4),weigl(4)
DOUBLE PREcIsION:: PHI(16),DPHIX(16),DPHIY(16),AJACO(2,2),AJACOI(2,2),DXHI(16), &
     DTHE(16),GAUSSPT(4),GAUSSWT(4),XHI,THE,DETER,S11,S22,CNST,ESK(NOPE,NOPE)

  masa=0.0
  
  do ielem = 1,nelements

      do inode = 1,nodpel
         ns(inode)= conect(ielem,inode)
         x(inode) = coor_x(ns(inode))
         y(inode) = coor_y(ns(inode))
      end do

      call armotodo(nodpel,x,y,z,deriv,weigc)

      do inode=1,nodpel
             
             AJACO(1,1) = 0.0
             AJACO(1,2) = 0.0
             AJACO(2,1) = 0.0
             AJACO(2,2) = 0.0
             do jnode = 1,nodpel
                AJACO(1,1) = AJACO(1,1) + x(jnode) * deriv(1,jnode,inode)
                AJACO(1,2) = AJACO(1,2) + x(jnode) * deriv(2,jnode,inode)
                AJACO(2,1) = AJACO(2,1) + y(jnode) * deriv(1,jnode,inode)
                AJACO(2,2) = AJACO(2,2) + y(jnode) * deriv(2,jnode,inode)
             end do
 
             gpdet = AJACO(1,1) * AJACO(2,2) - AJACO(2,1) * AJACO(1,2)
 
         gpvol=  weigc(inode)*gpdet
         masa(ns(inode))=masa(ns(inode))+gpvol

      enddo

   enddo
     !
     ! Loop over nodes to control zero-volume points
     !
     do ipoin=1,nnodes
       if(masa(ipoin)<1e-12) then
          write(6,*) 'nodo  ',ipoin,' posee matrix de masa cero  ',masa(ipoin)    
		  stop ' '
        end if
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug
    !    write(12,*) ' *masa(',ipoin,')= ', masa(ipoin)
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug

     enddo 



end subroutine masadiag2d

subroutine masadiag()
use def_solver
use def_variables
use def_constantes
implicit none
double precision x(nodpel),y(nodpel),z(nodpel),deriv(3,nodpel,nodpel),gpdet,gpvol

integer ns(nodpel),ipoin,inode,ielem,kk,j,k,pdime,pnode,mdime,igaus,ilocs,jlocs
  !    calculamos la matriz de masa diagonal usando una regla cerrada de integracion

double precision weigc(nodpel),posgl(3),weigl(3)

  masa=0.0
  
  do ielem = 1,nelements

      do inode = 1,nodpel
         ns(inode)= conect(ielem,inode)
         x(inode) = coor_x(ns(inode))
         y(inode) = coor_y(ns(inode))
         z(inode) = coor_z(ns(inode))
      end do

      call armotodo(nodpel,x,y,z,deriv,weigc)
      
      do inode=1,nodpel
	     
         call armodeter(nodpel,x,y,z,deriv(:,:,inode),gpdet) ! necesito el determinante en este punto de gauss=nodo
         gpvol=weigc(inode)*gpdet
       
         masa(ns(inode))=masa(ns(inode))+gpvol

      enddo
   enddo

     !
     ! Loop over nodes to control zero-volume points
     !
     do ipoin=1,nnodes
       if(masa(ipoin)<1e-12) then
          write(6,*) 'nodo  ',ipoin,' posee matrix de masa cero  ',masa(ipoin)    
		  stop ' '
        end if
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug
    !    write(12,*) ' *masa(',ipoin,')= ', masa(ipoin)
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug

     enddo 


    


end subroutine masadiag

subroutine armodeter(nope,x,y,z,deriv,gpdet)
implicit none
integer nope,j,k
double precision x(nope),y(nope),z(nope),gpdet,deriv(3,nope),shape(nope)
double precision xjacm(3,3),xjaci(3,3)

   
   
   do j=1,3
     xjacm(1,j)=0.0
     xjacm(2,j)=0.0
     xjacm(3,j)=0.0
     do k=1,NOPE
       xjacm(1,j)= xjacm(1,j)+X(k)*DERIV(j,k)
       xjacm(2,j)= xjacm(2,j)+Y(k)*DERIV(j,k)
       xjacm(3,j)= xjacm(3,j)+Z(k)*DERIV(j,k)
     end do

   enddo
      
   call invmtx(xjacm,xjaci,gpdet)

end subroutine armodeter


subroutine armotodo(nope,x,y,z,deriv,weigc)
implicit none
integer :: nope
double precision x(nope),y(nope),z(nope),weigc(nope),deriv(3,nope,nope),shape(nope)
double precision posgc(3,nope),posgl(4),weigl(4)
integer inoga(nope),pnode,pdime,mdime,nlocs,igaus,ilocs,jlocs,klocs,inode


 !
  ! Element shape function and derivatives SHAPC,DERIC,HESLC,WEIGC 
  ! for a close rule
  !- For each element type, using a closed integration rule:
  !      WEIGC(nnode)
  !      SHAPC(nnode,nnode)
  !      DERIC(ndime,nnode,nnode)
  !      HESLC(ntens,nnode,nnode)
 
    pnode=nope
 
    if(pnode==4) then
    
        posgc(1,1)= 0.0 
        posgc(2,1)= 0.0
        posgc(3,1)= 0.0
        posgc(1,2)= 1.0
        posgc(2,2)= 0.0
        posgc(3,2)= 0.0
        posgc(1,3)= 0.0
        posgc(2,3)= 1.0
        posgc(3,3)= 0.0
        posgc(1,4)= 0.0
        posgc(2,4)= 0.0
        posgc(3,4)= 1.0
        weigc(  1)= 1.0/24.0
        weigc(  2)= 1.0/24.0
        weigc(  3)= 1.0/24.0
        weigc(  4)= 1.0/24.0
     pdime=3
     mdime=3
	 
     nlocs=4 
 
    elseif(pnode==8) then
      inoga(1)= 1
     inoga(2)= 5
     inoga(3)= 4
     inoga(4)= 8
     inoga(5)= 2
     inoga(6)= 6
     inoga(7)= 3
     inoga(8)= 7
     pdime=3
     mdime=3
	 
     nlocs=2 
     posgl(1)=-1.0
     posgl(2)= 1.0
     weigl(1)= 1.0
     weigl(2)= 1.0

   elseif(pnode==27) then
     pdime=3
     mdime=3
	 
     nlocs=3 
     inoga( 1)= 1
     inoga( 2)=13
     inoga( 3)= 5
     inoga( 4)=12
     inoga( 5)=25
     inoga( 6)=20
     inoga( 7)= 4
     inoga( 8)=16
     inoga( 9)= 8
     inoga(10)= 9
     inoga(11)=22
     inoga(12)=17
     inoga(13)=21
     inoga(14)=27
     inoga(15)=26
     inoga(16)=11
     inoga(17)=24
     inoga(18)=19
     inoga(19)= 2
     inoga(20)=14
     inoga(21)= 6
     inoga(22)=10
     inoga(23)=23
     inoga(24)=18
     inoga(25)= 3
     inoga(26)=15
     inoga(27)= 7
 
     posgl(1)=-1.0
     posgl(2)= 0.0
     posgl(3)= 1.0
     weigl(1)= 1.0/3.0
     weigl(2)= 4.0/3.0
     weigl(3)= 1.0/3.0

   elseif(pnode==16) then
     
     inoga( 1)= 1
     inoga( 2)=12
     inoga( 3)= 11
     inoga( 4)=4
     inoga( 5)=5
     inoga( 6)=13
     inoga( 7)= 16
     inoga( 8)=10
     inoga( 9)= 6
     inoga(10)= 14
     inoga(11)=15
     inoga(12)=9
     inoga(13)=2
     inoga(14)=7
     inoga(15)=8
     inoga(16)=3

     pdime=2
     mdime=2
	 
     nlocs=4 
     
     posgl(1)=-1.0
     posgl(2)=-1.0/3.0
     posgl(3)= 1.0/3.0
     posgl(4)= 1.0
     weigl(1)= 1.0/4.0
     weigl(2)= 3.0/4.0
     weigl(3)= 3.0/4.0
     weigl(4)= 1.0/4.0

  endif


  if(pdime==3) then

     if(nope/=4) then
         igaus=0
         do ilocs=1,nlocs
            do jlocs=1,nlocs
               do klocs=1,nlocs
                  igaus=igaus+1
                  weigc(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
                  posgc(1,inoga(igaus))=posgl(ilocs)
                  posgc(2,inoga(igaus))=posgl(jlocs)
                  posgc(3,inoga(igaus))=posgl(klocs)
               end do
            end do
         end do
     else ! tetras
        
     
     endif     
    
     do inode=1,pnode
        
		 call shape3(posgc(1,inode),posgc(2,inode),posgc(3,inode),nope,shape,deriv(:,:,inode))
		 
	 end do

  else
  
  
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
              igaus=igaus+1
              weigc(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
              posgc(1,inoga(igaus))=posgl(ilocs)
              posgc(2,inoga(igaus))=posgl(jlocs)
        end do
     end do
     
    
     do inode=1,pnode
        
		 call shape2(posgc(1,inode),posgc(2,inode),nope,shape,deriv(:,:,inode))
		 
	 end do

  endif  
end subroutine armotodo
	  
	subroutine invmtx(a,b,deter)

	!-----------------------------------------------------------------------
	!
	! This routine inverts a square matrix A -> Mat(nsize,nsize). The
	! inverse is stored in B. Its determinant is DETER
	!
	!    
	!-----------------------------------------------------------------------
	  DOUBLE PRECISION  a(3,3)
	  DOUBLE PRECISION  b(3,3),deter
	  DOUBLE PRECISION  denom,t1,t2,t3,t4

	     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
	     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
	     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
	     deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
	     if(deter==0.0) return
	     denom = 1.0/deter
	     b(1,1) = t1*denom
	     b(2,1) = t2*denom
	     b(3,1) = t3*denom
	     b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
	     b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
	     b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
	     b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
	     b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
	     b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom


	end subroutine invmtx



subroutine shape3(s,t,z,nnode,shape,deriv)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates shape functions and their first and
  ! second derivatives 3-d standar continuous interpolation
  ! elements.
  ! 
  ! TETRAHEDRA:  4  10  &  20  nodes
  ! HEXAHEDRA:   8  27  &  64  nodes
  ! PRISM:       6             nodes
  !
  !-----------------------------------------------------------------------
  implicit none
  integer :: nnode
  double precision  :: s,t,z
  double precision  :: deriv(3,nnode),shape(nnode)
  integer              :: i
  double precision    :: a1,a2,a3,a4,a,p1,p2,p3,z1,z2,z3,z4,s1,s2,s3,s4
  double precision    :: t1,t2,t3,t4,sm,tm,zm,sq,tp,zp,s11,s21,s31,s41
  double precision    :: t11,t21,t31,t41,z11,z21,z31,s12,s22,s32,s42
  double precision    :: t12,t22,t32,t42,z41,z12,z22,z32,z42,sl,tl,zl
 
 if(nnode==4) then
     deriv=0.0
     shape(   1) = 1.0-s-t-z
     shape(   2) = s
     shape(   3) = t
     shape(   4) = z
     deriv(1, 1) =-1.0
     deriv(2, 1) =-1.0
     deriv(3, 1) =-1.0
     deriv(1, 2) = 1.0
     deriv(2, 3) = 1.0
     deriv(3, 4) = 1.0
  elseif(nnode==8) then
     !
     ! Trilinear brick 
     !   
     sm = 0.5*(1.0-s)
     tm = 0.5*(1.0-t)
     zm = 0.5*(1.0-z)
     sq = 0.5*(1.0+s)
     tp = 0.5*(1.0+t)
     zp = 0.5*(1.0+z)
     shape(   1) = sm*tm*zm
     deriv(1, 1) =-0.5*tm*zm
     deriv(2, 1) =-0.5*sm*zm
     deriv(3, 1) =-0.5*sm*tm
     
     shape(   2) = sq*tm*zm
     deriv(1, 2) = 0.5*tm*zm
     deriv(2, 2) =-0.5*sq*zm
     deriv(3, 2) =-0.5*sq*tm
    
     shape(   3) = sq*tp*zm
     deriv(1, 3) = 0.5*tp*zm
     deriv(2, 3) = 0.5*sq*zm
     deriv(3, 3) =-0.5*sq*tp
     shape(   4) = sm*tp*zm
     deriv(1, 4) =-0.5*tp*zm
     deriv(2, 4) = 0.5*sm*zm
     deriv(3, 4) =-0.5*sm*tp
     shape(   5) = sm*tm*zp
     deriv(1, 5) =-0.5*tm*zp
     deriv(2, 5) =-0.5*sm*zp
     deriv(3, 5) = 0.5*sm*tm
     shape(   6) = sq*tm*zp 
     deriv(1, 6) = 0.5*tm*zp
     deriv(2, 6) =-0.5*sq*zp
     deriv(3, 6) = 0.5*sq*tm
     shape(   7) = sq*tp*zp
     deriv(1, 7) = 0.5*tp*zp
     deriv(2, 7) = 0.5*sq*zp
     deriv(3, 7) = 0.5*sq*tp
     shape(   8) = sm*tp*zp
     deriv(1, 8) =-0.5*tp*zp
     deriv(2, 8) = 0.5*sm*zp
     deriv(3, 8) = 0.5*sm*tp
     
  else if(nnode==27) then
     ! Triquadratic brick
     !        
     sl=s*(s-1.0)
     tl=t*(t-1.0)
     zl=z*(z-1.0)
     sq=s*(s+1.0)
     tp=t*(t+1.0)
     zp=z*(z+1.0)
     s1= 2.0*s-1.0
     t1= 2.0*t-1.0
     z1= 2.0*z-1.0
     s2= 1.0-s*s
     t2= 1.0-t*t
     z2= 1.0-z*z
     s3= 1.0+2.0*s
     t3= 1.0+2.0*t
     z3= 1.0+2.0*z
     s4=-2.0*s
     t4=-2.0*t
     z4=-2.0*z
     shape(   1) = 0.125*sl*tl*zl
     deriv(1, 1) = 0.125*s1*tl*zl
     deriv(2, 1) = 0.125*sl*t1*zl
     deriv(3, 1) = 0.125*sl*tl*z1
     shape(   2) = 0.125*sq*tl*zl
     deriv(1, 2) = 0.125*s3*tl*zl
     deriv(2, 2) = 0.125*sq*t1*zl
     deriv(3, 2) = 0.125*sq*tl*z1
     shape(   3) = 0.125*sq*tp*zl
     deriv(1, 3) = 0.125*s3*tp*zl
     deriv(2, 3) = 0.125*sq*t3*zl
     deriv(3, 3) = 0.125*sq*tp*z1
     shape(   4) = 0.125*sl*tp*zl
     deriv(1, 4) = 0.125*s1*tp*zl
     deriv(2, 4) = 0.125*sl*t3*zl
     deriv(3, 4) = 0.125*sl*tp*z1
     shape(   5) = 0.125*sl*tl*zp
     deriv(1, 5) = 0.125*s1*tl*zp
     deriv(2, 5) = 0.125*sl*t1*zp
     deriv(3, 5) = 0.125*sl*tl*z3
     shape(   6) = 0.125*sq*tl*zp
     deriv(1, 6) = 0.125*s3*tl*zp
     deriv(2, 6) = 0.125*sq*t1*zp
     deriv(3, 6) = 0.125*sq*tl*z3
     shape(   7) = 0.125*sq*tp*zp
     deriv(1, 7) = 0.125*s3*tp*zp
     deriv(2, 7) = 0.125*sq*t3*zp
     deriv(3, 7) = 0.125*sq*tp*z3
     shape(   8) = 0.125*sl*tp*zp
     deriv(1, 8) = 0.125*s1*tp*zp
     deriv(2, 8) = 0.125*sl*t3*zp
     deriv(3, 8) = 0.125*sl*tp*z3
     shape(   9) = 0.25*s2*tl*zl
     deriv(1, 9) = 0.25*s4*tl*zl
     deriv(2, 9) = 0.25*s2*t1*zl
     deriv(3, 9) = 0.25*s2*tl*z1
     shape(  10) = 0.25*sq*t2*zl
     deriv(1,10) = 0.25*s3*t2*zl
     deriv(2,10) = 0.25*sq*t4*zl
     deriv(3,10) = 0.25*sq*t2*z1
     shape(  11) = 0.25*s2*tp*zl
     deriv(1,11) = 0.25*s4*tp*zl
     deriv(2,11) = 0.25*s2*t3*zl
     deriv(3,11) = 0.25*s2*tp*z1
     shape(  12) = 0.25*sl*t2*zl
     deriv(1,12) = 0.25*s1*t2*zl
     deriv(2,12) = 0.25*sl*t4*zl
     deriv(3,12) = 0.25*sl*t2*z1
     shape(  13) = 0.25*sl*tl*z2
     deriv(1,13) = 0.25*s1*tl*z2
     deriv(2,13) = 0.25*sl*t1*z2
     deriv(3,13) = 0.25*sl*tl*z4
     shape(  14) = 0.25*sq*tl*z2
     deriv(1,14) = 0.25*s3*tl*z2
     deriv(2,14) = 0.25*sq*t1*z2
     deriv(3,14) = 0.25*sq*tl*z4
     shape(  15) = 0.25*sq*tp*z2
     deriv(1,15) = 0.25*s3*tp*z2
     deriv(2,15) = 0.25*sq*t3*z2
     deriv(3,15) = 0.25*sq*tp*z4
     shape(  16) = 0.25*sl*tp*z2
     deriv(1,16) = 0.25*s1*tp*z2
     deriv(2,16) = 0.25*sl*t3*z2
     deriv(3,16) = 0.25*sl*tp*z4
     shape(  17) = 0.25*s2*tl*zp
     deriv(1,17) = 0.25*s4*tl*zp
     deriv(2,17) = 0.25*s2*t1*zp
     deriv(3,17) = 0.25*s2*tl*z3
     shape(  18) = 0.25*sq*t2*zp
     deriv(1,18) = 0.25*s3*t2*zp
     deriv(2,18) = 0.25*sq*t4*zp
     deriv(3,18) = 0.25*sq*t2*z3
     shape(  19) = 0.25*s2*tp*zp
     deriv(1,19) = 0.25*s4*tp*zp
     deriv(2,19) = 0.25*s2*t3*zp
     deriv(3,19) = 0.25*s2*tp*z3
     shape(  20) = 0.25*sl*t2*zp
     deriv(1,20) = 0.25*s1*t2*zp
     deriv(2,20) = 0.25*sl*t4*zp
     deriv(3,20) = 0.25*sl*t2*z3
     shape(  21) = 0.5*s2*t2*zl
     deriv(1,21) = 0.5*s4*t2*zl
     deriv(2,21) = 0.5*s2*t4*zl
     deriv(3,21) = 0.5*s2*t2*z1
     shape(  22) = 0.5*s2*tl*z2
     deriv(1,22) = 0.5*s4*tl*z2
     deriv(2,22) = 0.5*s2*t1*z2
     deriv(3,22) = 0.5*s2*tl*z4
     shape(  23) = 0.5*sq*t2*z2
     deriv(1,23) = 0.5*s3*t2*z2
     deriv(2,23) = 0.5*sq*t4*z2
     deriv(3,23) = 0.5*sq*t2*z4
     shape(  24) = 0.5*s2*tp*z2
     deriv(1,24) = 0.5*s4*tp*z2
     deriv(2,24) = 0.5*s2*t3*z2
     deriv(3,24) = 0.5*s2*tp*z4
     shape(  25) = 0.5*sl*t2*z2
     deriv(1,25) = 0.5*s1*t2*z2
     deriv(2,25) = 0.5*sl*t4*z2
     deriv(3,25) = 0.5*sl*t2*z4
     shape(  26) = 0.5*s2*t2*zp
     deriv(1,26) = 0.5*s4*t2*zp
     deriv(2,26) = 0.5*s2*t4*zp
     deriv(3,26) = 0.5*s2*t2*z3
     shape(  27) = s2*t2*z2
     deriv(1,27) = s4*t2*z2
     deriv(2,27) = s2*t4*z2
     deriv(3,27) = s2*t2*z4
endif
  
end subroutine shape3


subroutine shape2(s,t,nnode,shapf,deriv)

!-----------------------------------------------------------------------
!
!    This routine evaluates shape functions and their first and
!    second derivatives for 2-d continuos standar interpolation 
!    elements.
!
!    TRIANGLES       3   6  &  10  nodes
!    QUADRILATERALS  4   9  &  16  nodes
!
!-----------------------------------------------------------------------

  implicit none
  integer:: nnode
  double precision  :: s,t
  double precision :: deriv(3,nnode),shapf(nnode)
  integer              :: ii,jj
  double precision                 :: st,a1,a2,a3,ss,tt,s1,s2,s3,s4
  double precision                 :: t1,t2,t3,t4,s9,t9,c,a

  
     a =81.0/256.0
     c =1.0/3.0
     s1=1.0+s
     s2=c+s
     s3=c-s
     s4=1.0-s
     t1=1.0+t
     t2=c+t
     t3=c-t
     t4=1.0-t
     shapf( 1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
     shapf( 2) =   a*s1*s2*s3*t2*t3*t4                   ! 
     shapf( 3) =   a*s1*s2*s3*t1*t2*t3                   ! 
     shapf( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
     shapf( 5) =-3.0*a*s1*s3*s4*t2*t3*t4              !
     shapf( 6) =-3.0*a*s1*s2*s4*t2*t3*t4              !
     shapf( 7) =-3.0*a*s1*s2*s3*t1*t3*t4              ! 12   13   14    7
     shapf( 8) =-3.0*a*s1*s2*s3*t1*t2*t4              !
     shapf( 9) =-3.0*a*s1*s2*s4*t1*t2*t3              !
     shapf(10) =-3.0*a*s1*s3*s4*t1*t2*t3              ! 1     5    6    2
     shapf(11) =-3.0*a*s2*s3*s4*t1*t2*t4                 
     shapf(12) =-3.0*a*s2*s3*s4*t1*t3*t4
     shapf(13) = 9.0*a*s1*s3*s4*t1*t3*t4
     shapf(14) = 9.0*a*s1*s2*s4*t1*t3*t4
     shapf(15) = 9.0*a*s1*s2*s4*t1*t2*t4
     shapf(16) = 9.0*a*s1*s3*s4*t1*t2*t4
     deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 5)=-3.0*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1, 6)=-3.0*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1, 7)=-3.0*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 8)=-3.0*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 9)=-3.0*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
     deriv(1,10)=-3.0*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
     deriv(1,11)=-3.0*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,12)=-3.0*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,13)= 9.0*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1,14)= 9.0*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,15)= 9.0*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,16)= 9.0*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 5)= -3.0*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 6)= -3.0*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 7)= -3.0*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
     deriv(2, 8)= -3.0*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
     deriv(2, 9)= -3.0*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,10)= -3.0*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,11)= -3.0*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,12)= -3.0*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,13)=  9.0*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,14)=  9.0*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,15)=  9.0*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,16)=  9.0*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)
  
end subroutine shape2








