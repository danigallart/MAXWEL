function verma3d(mater,n_od,n_el)


NNOD=27;
%NNOD=8;
fid=fopen('../mallado.fem','r');

fid2=fopen('../contorno.fem','r');

fscanf(fid2,'%s %s',2)
nter = fscanf(fid2,'%i',1)

for kk=1:nter
    n =  fscanf(fid2,'%i',1);
    nod(kk)=fscanf(fid2,'%i',1);
    xcc(kk)=fscanf(fid2,'%f',1);
    ycc(kk)=fscanf(fid2,'%f',1);
    zcc(kk)=fscanf(fid2,'%f',1);
end

fscanf(fid2,'%s %s',2)
npot = fscanf(fid2,'%i',1)

for kk=1:npot
    n =  fscanf(fid2,'%i',1);
    nod(nter + kk)=fscanf(fid2,'%i',1);
    xcc(nter + kk)=fscanf(fid2,'%f',1);
    ycc(nter + kk)=fscanf(fid2,'%f',1);
    zcc(nter + kk)=fscanf(fid2,'%f',1);
end

%plot3(xcc,ycc,zcc,'*');zlabel('z')
%pause;

fclose(fid2);

[nnodos,count]=fscanf(fid,'%i',1);

n=4*nnodos;

[coor,count]=fscanf(fid,'%i %f %f %f\n',n);

[el,count]=fscanf(fid,'%i',1);


E=(NNOD+2)*el;

[nel,count]=fscanf(fid,'%i',E);

fclose(fid);

whitebg('white') 
hold on

view(60,20)
if(NNOD==27)
  nnod3 = 4; %NNOD/3;  
  
  for k = 1:el
    mat(k)=nel((NNOD+2)*(k-1)+2);
    for j=1:NNOD
      ii(j)=nel((NNOD+2)*(k-1)+2+j);
    end
      
    if(mat(k)==mater) 
      
        for j=1:nnod3,
          j2=j+1;
          if(j2==nnod3+1) 
              j2=1;
          end
          x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(j2))];
          line(x1,y1,z1)
        end
        for j=nnod3+1:nnod3*2,
          j2=j+1;
          if(j2==2*nnod3+1) 
              j2=nnod3+1;
          end
          x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(j2))];
          line(x1,y1,z1)
        end
        
        for j=2*nnod3+1:NNOD
            x= coor(4*ii(j)-2);
            y= coor(4*ii(j)-1);
            z= coor(4*ii(j));
            rr=0.001;
            %ellipsoid(x,y,z,rr,rr,rr);
        end
            %rectangle('Position',[x-rr,y-rr,2*rr,2*rr],'Curvature',[1,1],'FaceColor',[0 0 0])
      %for j=nnod3+1:2*nnod3-2,
      %    j2=j+1;
      %    if(j2==2*nnod3-1) 
      %        j2=nnod3+1;
      %    end
      %    x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
      %    y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
      %    z1=[coor(4*ii(j)),coor(4*ii(j2))];
      %      line(x1,y1,z1)
      %end
     % x= coor(4*ii(2*nnod3)-2);
     % y= coor(4*ii(2*nnod3)-1);
     % z= coor(4*ii(2*nnod3));
     %rr=0.001;
     %ellipsoid(x,y,z,rr,rr,rr);
     
      %for j=2*nnod3+1:3*nnod3-2,
      %    j2=j+1;
      %    if(j2==3*nnod3-1) 
      %        j2=2*nnod3+1;
      %    end
      %    x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
      %    y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
      %    z1=[coor(4*ii(j)),coor(4*ii(j2))];
      %    line(x1,y1,z1)
      %end
      %x= coor(4*ii(3*nnod3)-2);
      %y= coor(4*ii(3*nnod3)-1);
      %z= coor(4*ii(3*nnod3));
     %rr=0.001;
     %ellipsoid(x,y,z,rr,rr,rr);
  
     for jk=1:nnod3
       x1=[coor(4*ii(jk)-2),coor(4*ii(jk+nnod3)-2)];
       y1=[coor(4*ii(jk)-1),coor(4*ii(jk+nnod3)-1)];
       z1=[coor(4*ii(jk)),coor(4*ii(jk+nnod3))];
       line(x1,y1,z1)
      % x1=[coor(4*ii(jk+nnod3)-2),coor(4*ii(jk+2*nnod3)-2)];
      % y1=[coor(4*ii(jk+nnod3)-1),coor(4*ii(jk+2*nnod3)-1)];
      % z1=[coor(4*ii(jk+nnod3)),coor(4*ii(jk+2*nnod3))];
      % line(x1,y1,z1)
     end  
     
    if(n_od==1)
       for kno=1:NNOD
        xno = coor(4*ii(kno)-2);
        yno = coor(4*ii(kno)-1);
        zno = coor(4*ii(kno));
       
        s = sprintf('%i',ii(kno));
       
         text(xno,yno,zno,s);
       end
    
    end  
    %pause
    end
  end
  
end
    

end
