function veores(v0)

fid=fopen('z_grid.txt','r');
fid2=fopen('results.out','r');

nz = fscanf(fid,'%i',1)

for kk=1:nz
    n =  fscanf(fid,'%i',1);
    zcc(kk)=fscanf(fid,'%f',1);
end

fclose(fid);


[nnodos,count]=fscanf(fid2,'%i',1);

for kk=1:nnodos
    cx(kk)=fscanf(fid2,'%f',1);
    cy(kk)=fscanf(fid2,'%f',1);
    cz(kk)=fscanf(fid2,'%f',1);
    res(kk)=fscanf(fid2,'%f',1);
end

fclose(fid2)

nv=10;
for kk=1:nv+1
    v(kk)= v0 - (kk-1)*v0/nv;
end

for kk=1:nv+1        
   nt=0; 
   for jj=1:nnodos
   
       if(res(jj)>=v(kk)-v0/(2*nv) & res(jj)<=v(kk)+v0/(2*nv)  )
             nt=nt+1;
             cx_z(nt)=cx(jj);
             cy_z(nt)=cy(jj);
             cz_z(nt)=cz(jj);
        end
    end
      
    plot3(cx_z,cy_z,cz_z,'.');axis([0 2 0 2 -1 2])
    v(kk)
    nt
    pause
    clear cx_z cy_z cz_z;
    
end


    
    
    
    
