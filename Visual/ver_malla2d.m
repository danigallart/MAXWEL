function ver_malla2d(file,nope)

n_od=0;
n_el=0;


fid=fopen(file,'r');

if(nope==4) 

    [nnodos,count]=fscanf(fid,'%i',1);
    for kk=1:nnodos
        n=fscanf(fid,'%i',1);   
        xx(kk)=fscanf(fid,'%f',1);   
        yy(kk)=fscanf(fid,'%f\n',1);   
    end
    [nele,count]=fscanf(fid,'%i',1);
    for kk=1:nele
        n=fscanf(fid,'%i',1);   
        mat(kk)=fscanf(fid,'%i',1);   
        for jj=1:nope
             nel(jj)=fscanf(fid,'%i',1);    
        end
        
        for jj=1:nope
            ii=jj+1;
            if(jj==nope)
                ii=1;
            end
            x1=[xx(nel(jj)),xx(nel(ii))];
            y1=[yy(nel(jj)),yy(nel(ii))];
            line(x1,y1);
        
        end
    end
    
    % temepr
    nod_tem=fscanf(fid,'%i',1);
    for kk=1:nod_tem
        nodt(kk)=fscanf(fid,'%i\n',1);
        x= xx(nodt(kk));
        y= yy(nodt(kk));
        rr=0.001;
        rectangle('Position',[x-rr,y-rr,2*rr,2*rr],'Curvature',[1,1],'FaceColor',[0 0 0])
        
    end
    
    % poten
    nod_pot=fscanf(fid,'%i',1);
    for kk=1:nod_pot
        nodp(kk)=fscanf(fid,'%i\n',1);
        x= xx(nodp(kk));
        y= yy(nodp(kk));
        rr=0.001;
        rectangle('Position',[x-rr,y-rr,2*rr,2*rr],'Curvature',[1,1],'FaceColor',[0 1 0])
        
    end
    
    % tierra
    nod_tie=fscanf(fid,'%i',1);
    for kk=1:nod_tie
        nodtie(kk)=fscanf(fid,'%i\n',1);
        x= xx(nodtie(kk));
        y= yy(nodtie(kk));
        rr=0.001;
        rectangle('Position',[x-rr,y-rr,2*rr,2*rr],'Curvature',[1,1],'FaceColor',[1 0 0])
        
    end
    
else

    [nnodos,count]=fscanf(fid,'%i',1);

    n=3*nnodos;

    [coor,count]=fscanf(fid,'%f\n',n);

    [el,count]=fscanf(fid,'%i',1);


    E=(nope*nope+2)*el;

    [nel,count]=fscanf(fid,'%i',E);

    fclose(fid);

    whitebg('white') 
    hold on

      for k = 1:el

          mat = nel((nope*nope+2)*(k-1)+2);
          for j=1:nope*nope
            ii(j)=nel((nope*nope+2)*(k-1)+2+j);
          end

        for j=1:nope,
            j2=j+1;
            if(j2==nope+1) 
              j2=1;
            end
            x1=[coor(3*ii(j)-1),coor(3*ii(j2)-1)];
            y1=[coor(3*ii(j)),coor(3*ii(j2))];
            %if(mat==1)
            %   line(x1,y1,'Color',[1 0 0]);
            %else
               line(x1,y1);
            %end
        end

        for j=nope+1:nope*nope
           x= coor(3*ii(j)-1);
           y= coor(3*ii(j));
           rr=0.001;
           rectangle('Position',[x-rr,y-rr,2*rr,2*rr],'Curvature',[1,1],'FaceColor',[0 0 0])
        end

        if(n_od==1)
           for kno=1:nope*nope
            xno = coor(3*ii(kno)-1);
            yno = coor(3*ii(kno));
            s = sprintf('%i',ii(kno));

             text(xno,yno,s);
           end

        end  
       %pause
      end
end
