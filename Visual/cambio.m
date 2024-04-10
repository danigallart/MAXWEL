fi=fopen('campo.dat','r+');
fi2=fopen('campo_c.csv','w+');

s1=fscanf(fi,'%s',1)
s2=fscanf(fi,'%s',1)
s3=fscanf(fi,'%s',1)
s4=fscanf(fi,'%s',1)

fprintf(fi2,'%c, %c, %c, %s\n',s1,s2,s3,s4);

[data,n]=fscanf(fi,'%f');
N4=n/4
for kk=1:N4
    fprintf(fi2,'%f, %f, %f, %f\n',data(4*(kk-1)+1),data(4*(kk-1)+2),data(4*(kk-1)+3),data(4*(kk-1)+4));
end
fclose(fi);

fclose(fi2);



fi=fopen('gradiente.dat','r+');
fi2=fopen('gradiente_c.csv','w+');

s1=fscanf(fi,'%s',1)
s2=fscanf(fi,'%s',1)
s3=fscanf(fi,'%s',1)
s4=fscanf(fi,'%s',1)

fprintf(fi2,'%c, %c, %c, %s\n',s1,s2,s3,s4);

[data,n]=fscanf(fi,'%f');
N4=n/4
for kk=1:N4
    fprintf(fi2,'%f, %f, %f, %f\n',data(4*(kk-1)+1),data(4*(kk-1)+2),data(4*(kk-1)+3),data(4*(kk-1)+4));
end
fclose(fi);

fclose(fi2);

fi=fopen('results.3D','r+');
fi2=fopen('results_c.csv','w+');
s1=fscanf(fi,'%s',1)
s2=fscanf(fi,'%s',1)
s3=fscanf(fi,'%s',1)
s4=fscanf(fi,'%s',1)

fprintf(fi2,'%c, %c, %c, %c\n',s1,s2,s3,s4);

[data,n]=fscanf(fi,'%f');
N4=n/4
for kk=1:N4
    fprintf(fi2,'%f, %f, %f, %f\n',data(4*(kk-1)+1),data(4*(kk-1)+2),data(4*(kk-1)+3),data(4*(kk-1)+4));
end
fclose(fi);

fclose(fi2);

fi=fopen('saleplano.dat','r+');
fi2=fopen('saleplano_c.csv','w+');

s1=fscanf(fi,'%s',1)
s2=fscanf(fi,'%s',1)
s3=fscanf(fi,'%s',1)
s4=fscanf(fi,'%s',1)

fprintf(fi2,'%c, %c, %c, %c\n',s1,s2,s3,s4);

[data,n]=fscanf(fi,'%f');
N4=n/4
for kk=1:N4
    fprintf(fi2,'%f, %f, %f, %f\n',data(4*(kk-1)+1),data(4*(kk-1)+2),data(4*(kk-1)+3),data(4*(kk-1)+4));
end
fclose(fi);

fclose(fi2);
