%fi=fopen('tumor.dat','r+');
%fi2=fopen('tumor_c.csv','w+');

fi=fopen('zonaafec.dat','r+');
fi2=fopen('zonaa_c.csv','w+');


n=fscanf(fi,'%i',1);

[data,n]=fscanf(fi,'%f');
N4=n/5
for kk=1:N4
    fprintf(fi2,'%f, %f, %f, %f\n',data(5*(kk-1)+2),data(5*(kk-1)+3),data(5*(kk-1)+4),data(5*(kk-1)+5));
end
fclose(fi);

fclose(fi2);

