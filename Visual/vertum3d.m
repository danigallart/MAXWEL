function vertum3d(mater)

NNOD=8;
fid=fopen('mallado.fem','r');

[nnodos,count]=fscanf(fid,'%i',1);

n=4*nnodos;

[coor,count]=fscanf(fid,'%i %f %f %f\n',n);

[el,count]=fscanf(fid,'%i',1);


E=(NNOD+2)*el;

[nel,count]=fscanf(fid,'%i',E);

fclose(fid);

whitebg('white') 
hold on

for k = 1:el
    mat(k)=nel((NNOD+2)*(k-1)+2);
    if(mat(k)==mater) 
  
    xmed=0;
    ymed=0;
    zmed=0;
    for j=1:NNOD
      ii(j)=nel((NNOD+2)*(k-1)+2+j);
      xmed=xmed + coor(4*ii(j)-2)/NNOD;
      ymed=ymed + coor(4*ii(j)-1)/NNOD;
      zmed=zmed + coor(4*ii(j))/NNOD;
    end
      plot3(xmed,ymed,zmed,'.');
    end
end
