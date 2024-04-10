function veo_pot(ncase);


if(ncase==1)
    fi=fopen('saleplano.dat','r+');

    [dat,n]=fscanf(fi,'%f');

    NN=n/3

    for kk=1:NN
        x1(kk)=dat(3*kk-2);
        x2(kk)=dat(3*kk-1);
        V(kk)=dat(3*kk);
    end
    plot3(x1,x2,V,'.');
    
elseif(ncase==2)
    fi=fopen('grid2D.dat','r+');

    n1= fscanf(fi,'%i',1);
    n2= fscanf(fi,'%i',2);
    for kk=1:n1
        for jj=1:n2 
           xx(kk,jj)=fscanf(fi,'%f',1);
           zz(kk,jj)=0.0;
        end
    end
    for kk=1:n1
        for jj=1:n2 
           yy(kk,jj)=fscanf(fi,'%f',1);
        end
    end
    for kk=1:n1
        for jj=1:n2 
           V(kk,jj)=fscanf(fi,'%f',1);
        end
    end
    for kk=1:n1
        for jj=1:n2 
           grax(kk,jj)=fscanf(fi,'%f',1);
        end
    end

    for kk=1:n1
        for jj=1:n2 
           gray(kk,jj)=fscanf(fi,'%f',1);
        end
    end

    for kk=1:n1
        for jj=1:n2 
           graz(kk,jj)=fscanf(fi,'%f',1);
           camp(kk,jj) = sqrt(grax(kk,jj)*grax(kk,jj)+gray(kk,jj)*gray(kk,jj)+graz(kk,jj)*graz(kk,jj));
    
        end
    end

    %[C,h]= contour(xx,yy,camp);xlabel('Y direction');ylabel('Z direction');title('Campo plane X=0.0')
    %[C,h]= contour(xx,yy,camp);xlabel('Y direction');ylabel('Z direction');title('Campo plane X=0.0')
    [C,h]= contour(xx,yy,camp);xlabel('Y direction');ylabel('Z direction');title('Campo plane X=0.0')
    %[C,h]= contour(xx,yy,grax);xlabel('Y direction');ylabel('Z direction');title('Campo plane X=0.0')
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    
    pause
    
    [C,h]=contour(xx,yy,V,100);xlabel('Y direction');ylabel('Z direction');title('Equi potentials plane X=0.0')
    %set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    
    %subplot(2,2,1);contour(xx,yy,grax,20);xlabel('Y direction');ylabel('Z direction');title('Grad x plane X=0.0')
    %subplot(2,2,2);contour(xx,yy,gray,20);xlabel('Y direction');ylabel('Z direction');title('Grad y plane X=0.0')
    %subplot(2,2,3);contour(xx,yy,graz,20);xlabel('Y direction');ylabel('Z direction');title('Grad z  plane X=0.0')
    
elseif(ncase==3)
    fi=fopen('grid2D.dat','r+');
    nz= fscanf(fi,'%i',1);
    ny= fscanf(fi,'%i',1);
    nx= fscanf(fi,'%i',1);
    
    for kk=1:nz
        
        z=fscanf(fi,'%f',1)
    
        for jj=1:ny
             for ii=1:nx
                 xx(jj,ii)=fscanf(fi,'%f',1);
             end
        end
        
        for jj=1:ny
             for ii=1:nx
                 yy(jj,ii)=fscanf(fi,'%f',1);
             end
        end
        for jj=1:ny
             for ii=1:nx
                 cx(jj,ii)=fscanf(fi,'%f',1);
             end
        end
        for jj=1:ny
             for ii=1:nx
                 cy(jj,ii)=fscanf(fi,'%f',1);
             end
        end
        for jj=1:ny
             for ii=1:nx
                 cz(jj,ii)=fscanf(fi,'%f',1);
             end
        end
        for jj=1:ny
             for ii=1:nx
                 V2d(jj,ii)=fscanf(fi,'%f',1);
             end
        end

        for jj=1:ny
             for ii=1:nx
               camp(jj,ii) = sqrt(cx(jj,ii)*cx(jj,ii)+cy(jj,ii)*cy(jj,ii)+cz(jj,ii)*cz(jj,ii));
             end
        end
       %[C,h]=contour(xx,yy,V2d,30);
       %pause
       [C,h]= contour(xx,yy,camp);
        set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
       %subplot(2,2,1);contour(xx,yy,camp);
       %subplot(2,2,2);contour(xx,yy,cx);
       %subplot(2,2,3);contour(xx,yy,cy);
       %subplot(2,2,4);contour(xx,yy,abs(cz));
       
    pause
    end
        
end

fclose(fi);

