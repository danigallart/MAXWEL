function sacolinea(archivo)

fid=fopen(archivo,'r+');

hold on;
k=0;
while(feof(fid)==0)
    k=k+1;
    
    x(k)=fscanf(fid,'%f',1);    
    V(k)=fscanf(fid,'%f',1);    
    ca(k)=fscanf(fid,'%f\n',1);    
end

fclose(fid);

subplot(2,2,1);plot(x,V);xlabel('x (mm)');ylabel('V (Volt)')
subplot(2,2,2);plot(x,ca);ylabel('E (Volt/cm)')