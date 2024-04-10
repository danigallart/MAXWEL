function verma3d_al(mater,n_od,n_el)

%NNOD=27;
NNOD=8;
fid=fopen('../mallado.fem','r');

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
  nnod3 = NNOD/3;  
  for k = 1:el
    mat(k)=nel((NNOD+2)*(k-1)+2);
    for j=1:NNOD
      ii(j)=nel((NNOD+2)*(k-1)+2+j);
    end
      
    if(mat(k)==mater) 
      for j=1:4,
          j2=j+1;
          if(j2==5) 
              j2=1;
          end
          x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(j2))];
          line(x1,y1,z1)
      end
      for j=1:4,
          j2=j+1;
          if(j2==5) 
              j2=1;
          end
          x1=[coor(4*ii(4+j)-2),coor(4*ii(4+j2)-2)];
          y1=[coor(4*ii(4+j)-1),coor(4*ii(4+j2)-1)];
          z1=[coor(4*ii(4+j)),coor(4*ii(4+j2))];
          line(x1,y1,z1)
      end
      for j=1:4,
          j2=j+4;
%          x1=[coor(4*ii(4+j)-2),coor(4*ii(4+j2)-2)];
%          y1=[coor(4*ii(4+j)-1),coor(4*ii(4+j2)-1)];
%          z1=[coor(4*ii(4+j)),coor(4*ii(4+j2))];
          x1=[coor(4*ii(j)-2),coor(4*ii(4+j)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(4+j)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(4+j))];
          line(x1,y1,z1)
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
    
    end
    
    %pause
    %cla
  end
  
   
elseif(NNOD==8)
    for k = 1:el
    mat(k)=nel((NNOD+2)*(k-1)+2);
    for j=1:NNOD
      ii(j)=nel((NNOD+2)*(k-1)+2+j);
    end
    
    
    if(mat(k)==4)
         col='r';
    elseif(mat(k)==10 || mat(k)==11  )
         col='k';
    else
         col='b';
        
    end

    
    if(mat(k)==mater) 
      for j=1:4,
          j2=j+1;
          if(j2==5) 
              j2=1;
          end
          x1=[coor(4*ii(j)-2),coor(4*ii(j2)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(j2)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(j2))];
          line(x1,y1,z1,'color',col)
      end
      for j=1:4,
          j2=j+1;
          if(j2==5) 
              j2=1;
          end
          x1=[coor(4*ii(4+j)-2),coor(4*ii(4+j2)-2)];
          y1=[coor(4*ii(4+j)-1),coor(4*ii(4+j2)-1)];
          z1=[coor(4*ii(4+j)),coor(4*ii(4+j2))];
          line(x1,y1,z1,'color',col)
      end

      for j=1:4,
          j2=j+4;
          x1=[coor(4*ii(j)-2),coor(4*ii(4+j)-2)];
          y1=[coor(4*ii(j)-1),coor(4*ii(4+j)-1)];
          z1=[coor(4*ii(j)),coor(4*ii(4+j))];
          line(x1,y1,z1,'color',col)
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
