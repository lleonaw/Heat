%%  Driver for poisson.m 

%clear all; format long;
%close all; 

global Ne Nx ifplt initflg
ifplt = true; 
ifplt = false; 

mth = 2; % Fidkowski formulation
mth = 1; % Penalized central flux
mth = 3; % Primal formulation, added Sun Feb 28 12:57:15 CST 2016

% Ne = 10;          % Number of elem
N = 4;            % Poly. order
Nx = N + 1;       % Numb of points in each elem.

Nn  = 8; 
Nen = 4; 
ere = zeros(Nen,Nn); plnx = zeros(Nen,Nn); plne = zeros(Nen,Nn);
for initflg=2:2 % 4
   for j=4:4 %Nn
   %  N = j+1; Nx = N + 1;
      N = 1*j; Nx = N + 1;
   %  N = 2^j+1; Nx = N + 1;
      for i=1:1%Nen
%        Ne = i*2;
         Ne = 2;
         [succ,infer] = poisson(mth); 
         if(succ)
           ere(i,j) = infer; 
           plne(i,j) = Ne; 
           plnx(i,j) = Nx; 
         else 
           disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
         end 
         A = load('A.dat'); 
         disp(['Is A symmetric? 1-Yes, 0-No: ', num2str(issymmetric(A))...
              ,' . Cond(A) = ', num2str(cond(A))]);
         if(mth==3)
           Ku = load('K.dat');Gtu = load('Gtu.dat'); 
           Hu = load('H.dat'); G = load('G.dat');
           qqt = load('qqt.dat'); 
           disp([% 'Is K symmetric? 1-Yes, 0-No: ', num2str(issymmetric(Ku)),...
              ' . max(max(G-G^T)) = ', num2str(max(max(G-Gtu'))),...
              ' . max(max(H)) = ', num2str(max(max(Hu)))]); 
         end
      end 
   end 
   if(initflg==1)     % Init Case 1 , cos(pi x/2) 
    eric1 = ere(1,:); plnsv = plnx(1,:); 
   elseif(initflg==2) % Init Case 2 , sin(pi x) 
    eric2 = ere(1,:); 
   elseif(initflg==3) % Init Case 3 , cos(pi x) 
    eric3 = ere(1,:); 
   elseif(initflg==4) % Init Case 4 , 1 - cos(pi x/2) 
    eric4 = ere(1,:); 
   end 
end 

%% Various initial conditions 
%figure(9);  
%loglog(plnsv,eric1,'o-','linewidth',1.5);hold on;
%loglog(plnsv,eric2,'o-','linewidth',1.5);
%loglog(plnsv,eric3,'o-','linewidth',1.5);
%loglog(plnsv,eric4,'o-','linewidth',1.5);
%loglog(plnsv,plnsv.^(-2),'x-','linewidth',1.5);
%loglog(plnsv,exp(-plnsv),'x-','linewidth',1.5);
%legend('cos(\pi x/2)','sin(\pi x)','cos(\pi x)','1- cos(2 \pi x)','N^{-2}','e^{-N}'); 
%xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
%title('Pointwise error, Poisson problem, N_e = 2, on [-1,1]');

% A = load('A.dat'); 
% disp(['Is A symmetric? 1-Yes, 0-No: ', num2str(issymmetric(A))]);
% M = load('M.dat'); 
% spy(M); 
% chol(A)
% eig(A)
% spy(A-A') 

%% convergence
%% Nailing sin(pi x), 1 - cos(2pi x), and  1 + cos(pi x) 
%% Second order cos(pi x/2),
%% Failing cos(pi x)  
%%  - - - - Seems at the boundary the function values have to be 0 
%%  -       or at least my code belives they have to be 0 
%%  - - - - 
%% Fix Ne, varying N 
% figure(9);  
%%loglog(plnx(1,:),ere(1,:),'o-','linewidth',1.5);hold on;
%%loglog(plnx(1,:),plnx(1,:).^(-2),'x-','linewidth',1.5);
%%legend('Err data','N^{-2}'); 
% semilogy(plnx(1,:),ere(1,:),'o-','linewidth',1.5);hold on;
% semilogy(plnx(1,:),exp(-plnx(1,:)),'x-','linewidth',1.5); 
% legend('Err data','e^{-N}'); 
% xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
% title('Pointwise error, Poisson problem, N_e = 2');

%figure(9);  
%semilogy(plnx(1,:),ere(1,:),'o-','linewidth',1.5);hold on;
%semilogy(plnx(2,:),ere(2,:),'o-','linewidth',1.5);
%semilogy(plnx(3,:),ere(3,:),'o-','linewidth',1.5);
%semilogy(plnx(4,:),ere(4,:),'o-','linewidth',1.5);
%semilogy(plnx(1,:),exp(-plnx(1,:)),'x-','linewidth',2);
%legend('Num. Error,Ne=2'...%,'N_e^{(-2)},N=1'... 
%      ,'Num. Error,Ne=4'...%,'N_e^{(-3)},N=2'... 
%      ,'Num. Error,Ne=6'...%,'N_e^{(-4)},N=3'... 
%      ,'Num. Error,Ne=8'...%,'N_e^{(-5)},N=4'); 
%      ,'e^{-N}');...%,'N_e^{(-5)},N=4'); 
%title(['Inf error for varying N, varying N_e']);  hold off;
%xlabel('N'); ylabel('Inf norm error');
