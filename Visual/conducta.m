function conducta(Hyalu)
%fi=fopen('campo.dat','r+');
%s1=fscanf(fi,'%s',1)

sigma_el=0.0004;

if(Hyalu==0) 
   factorH=1.0;
else
   factorH=0.6;
   
end

Erev=24.0*factorH;
Eirrev=46.0*factorH;
a_min=(Erev+Eirrev)*0.5;

BB=1000.0;

N=100;

Emax=100;

for kk=1:N
   E=Emax/N*kk;
   campo(kk)=E;
   if(E> Eirrev) 
      funsigma3(kk) = 3.5*factorH;  
   elseif(E>= Erev && E<=Eirrev) 
      funsigma3(kk) = 1.0*factorH  + 2.5*factorH* (E-Erev)/(Eirrev-Erev)   ;
       
%      funsigma3(kk) = 1.0*factorH  + 2.5*factorH/( 1.0 + exp(-(E-a_min)/BB) )  ;
   else
      funsigma3(kk) = 1.0*factorH;   
   end
   funsigma3(kk) = funsigma3(kk)*sigma_el;
end

plot(campo,funsigma3)
%fclose(fi);
