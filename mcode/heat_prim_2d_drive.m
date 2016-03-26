%%  Driver for heat_prim_2d.m 
% To do:
%   -  Error analysis with varying eta 
%   -  Eigenvalues  
% % % % % % % % % % % % % % % % % % % % % %

%clear all; format long; close all; 
global Ne Nx ifplt initflg T CFL dt ifeig iffwd
global R2 u0 
dt = 1e-5;    % T = 4.e-0; 
ifplt = false; 
ifplt = true; 
iffwd = false; 
% 
ifeig = true; 
ifeig = false; 
if(ifeig)
    figure(43);
    mylegend = cell(Nn,1);  % For plotting eig 
end

T = 2.; N = 2; Nx = N + 1;  % N - Poly. order, Nx - Numb of points in each elem.
Ntn = 4; 
Nn  = 6; 
% Nn  = 4;  % 16 is enough... 32 takes a lot of time! ! 
Nen = 4;          % Number of elem

ere  = zeros(Nen,Nn); plnx = zeros(Nen,Nn); plne = zeros(Nen,Nn);
ertn = zeros(Ntn,1); dtn  = zeros(Ntn,1); 

for itm=1:1%Ntn
% RK4 should have a rather big CFL number, totally not seeing that 
%   lambda \delta t < 2.8, which means huge eigenvalues are present ! 
%   max CFL for RK4 is 0.15 ! 
    CFL = 0.03 /(2.^itm);   % For rk4, CFL = 0.045 runs, 0.05 does not 
%   CFL = 0.3 /(2.^itm) / 2.; 
    for initflg=5:5 % 4
%      for j=1:ceil(Nn)
       for j=2:2%ceil(Nn)
       % Use CFL = 0.0125 for polynomial  
          N = 5; Nx = N + 1; % Nx = 6, 
%         N = j; Nx = N + 1; 
%         N = 2^j; Nx = N + 1;
          for i=1:1%Nen
             Ne = (i+3)^2; %        Ne = i;
             tic; 
             [succ,infer,Vdg,Ddg] = heat_prim_2d; 
             toc; 
             if(succ)
               ere(i,j) = infer; plne(i,j) = Ne; plnx(i,j) = Nx; 
             else 
               disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
             end 
             if(succ && ifeig) 
                 disp(['N = ',num2str(N),' , N_e = ',num2str(Ne)]); 
                 disp(['For 2D DG, # of eigs is ',num2str(size(Ddg)),' , max(D) = ',num2str(max(Ddg)),' , min(D) = ',num2str(min(Ddg))]); 
                 disp('Lowest three eig. of inv(B) A'); 
                 Ddg(1:3),
                 plot(Ddg,'o-'); hold on; 
                 mylegend{j}=strcat('2d DG, N = ',num2str(N));
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
    ertn(itm) = infer; 
    dtn (itm) = dt; 
end 

if(ifeig)
    legend(mylegend,'location','northeast','Fontsize',15);
    xlabel('k','fontsize',18); ylabel('\lambda','fontsize',18); 
    title(['Eig values of 2D DG, Ne = ',num2str(Ne)],'fontsize',18); 
end

%% Various initial conditions 
%figure(9);  
%%loglog(plnsv,eric1,'o-','linewidth',1.5);hold on;
%loglog(plnsv,eric2,'o-','linewidth',1.5);
%%loglog(plnsv,eric3,'o-','linewidth',1.5);
%%loglog(plnsv,eric4,'o-','linewidth',1.5);
%%loglog(plnsv,plnsv.^(-2),'x-','linewidth',1.5);
%%loglog(plnsv,exp(-plnsv),'x-','linewidth',1.5);
%%legend('cos(\pi x/2)','sin(\pi x)','cos(\pi x)','1- cos(2 \pi x)','N^{-2}','e^{-N}'); 
%%xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
%%title('Pointwise error, Poisson problem, N_e = 2, on [-1,1]');

%% spatial convergence
%% Fix Ne, varying N 
%  figure(9);  
%  semilogy(plnx(1,:),ere(1,:),'o-','linewidth',1.5);hold on;
%  semilogy(plnx(1,:),10^3.*exp(-2.5*plnx(1,:)),'x-','linewidth',1.5); 
%  legend('Err data','C^{-N}'); 
%  xlabel('$N_x$','Interpreter','Latex'); 
%  ylabel('$\|u_{xx} - \tilde{u}_{xx}\|_{\infty}$','Interpreter','Latex');
%  title(['Error at T_{final}, u_t + \nu u_xx = 0., u_0 = sin(\pi x + \phi), \phi = 0, \eta = 1/M_{1D}, N_e = ',num2str(Ne)]);
  %title('Max pointwise relative error, heat problem, N_e = 2');
%  title(['Error for eval. u_{xx}, u = sin(\pi x + \phi), \phi = 0, N_e = ',num2str(Ne)]);

% temporal convergence 
% figure(9);  
% loglog(dtn,ertn,'o-','linewidth',1.5);hold on;
%%loglog(dtn,1.0*dtn.^(2),'x-','linewidth',1.5); 
% loglog(dtn,1.0*dtn.^(4),'x-','linewidth',1.5); 
%%legend('Err data','\Delta t^{2}','\Delta t'); 
%%legend('Err data','\Delta t^{2}'); 
% legend('Err data','\Delta t^{4}'); 
% xlabel('$\Delta t$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
% title('Max pointwise relative error, heat problem, N_e = 2');

% Element number 
%figure(9);  
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
