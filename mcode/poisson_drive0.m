%%  Driver for poisson.m 

%clear all; format long;

close all; 
global Ne Nx ifplt
ifplt = true; 
ifplt = false; 

% Ne = 10;          % Number of elem
N = 4;            % Poly. order
Nx = N + 1;       % Numb of points in each elem.

N  = 1; Nx = N + 1;
N  = 2; Nx = N + 1;
N  = 3; Nx = N + 1;
Nn  = 4; 
Nen = 4; ere = zeros(Nen,Nn); 
pln = zeros(Nen,Nn);

for j=1:Nn
   N = j; Nx = N + 1;
   for i=1:Nen
      Ne = i*10;
      [succ,infer] = poisson; 
      if(succ)
        ere(i,j) = infer; 
        pln(i,j) = i; 
      else 
        disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
      end 
   end 
end 

figure(2);  
loglog(pln(:,1),ere(:,1),'ro-','linewidth',2);hold on;
loglog(pln(:,1),6.*10.^(-2).*pln(:,1).^(-2),'rx-','linewidth',2);
loglog(pln(:,2),ere(:,2),'bo-','linewidth',2);
loglog(pln(:,2),1.*10.^(-2).*pln(:,2).^(-3),'bx-','linewidth',2);
loglog(pln(:,3),ere(:,3),'ko-','linewidth',2);
loglog(pln(:,3),3.*10.^(-4).*pln(:,3).^(-4),'kx-','linewidth',2);
loglog(pln(:,4),ere(:,4),'mo-','linewidth',2);
loglog(pln(:,4),1.*10.^(-5).*pln(:,4).^(-5),'mx-','linewidth',2); 
legend('Num. Error,N=1','N_e^{(-2)},N=1'... 
      ,'Num. Error,N=2','N_e^{(-3)},N=2'... 
      ,'Num. Error,N=3','N_e^{(-4)},N=3'... 
      ,'Num. Error,N=4','N_e^{(-5)},N=4'); 
title(['Inf error for varying N, varying N_e']);  hold off;
xlabel('N_e'); ylabel('Inf norm error');
