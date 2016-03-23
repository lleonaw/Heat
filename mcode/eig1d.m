%%  Driver for dgpm1d_eig.m, sem1d_eig.m
% To do:
%   -. Plot spectrum 
% % % % % % % % % % % % % % % % % % % % % %

%clear all; format long; close all; 
global Ne Nx ifplt initflg T CFL dt
dt = 1e-5;    % T = 4.e-0; 
ifplt = true; 
ifplt = false; 

%ifsem = false; 
ifsem = true; 
ifdgm = true; 
ifeig = true; 
if(ifeig) ifdgm = true; ifsem = true; end

N = 4;              % Poly. order
Nx = N + 1;         % Numb of points in each elem.

Ntn = 4; 
Nn  = 5; 
% Nn  = 4;  % 16 is enough... 32 takes some time, 64 takes a while ! ! 
Nen = 6; 

dg_ere  = zeros(Nen,Nn); dg_plnx = zeros(Nen,Nn); dg_plne = zeros(Nen,Nn);
dg_kcd  = zeros(Nen,Nn); dg_lcd  = zeros(Nen,Nn);
se_ere  = zeros(Nen,Nn); se_plnx = zeros(Nen,Nn); se_plne = zeros(Nen,Nn);
se_kcd  = zeros(Nen,Nn); se_lcd  = zeros(Nen,Nn);

ertn = zeros(Ntn,1); dtn  = zeros(Ntn,1); 

if(ifeig)
    figure(43);
    mylegend = cell(2*Nn,1);  % For plotting eig 
end
for initflg=2:2 % 4
   disp(['Init flg = ',num2str(initflg)]); 
   for j=Nn:Nn
%     N = 2^j; Nx = N + 1;
      N = 2*j; Nx = N + 1;
      for i=1:Nen
%        Ne = (i)^2;
         Ne = (2)^i;
%        Ne = 10;
         if(ifsem)
             [succ,infer,se_VK,se_DK,se_VL,se_DL] = sem1d_eig; 
             if(succ)
               se_ere(i,j) = infer; se_plne(i,j) = Ne; se_plnx(i,j) = Nx; 
             else 
               disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
             end 
         end
         if(ifdgm)
             [succ,infer,dg_VA,dg_DA,dg_VL,dg_DL] = dgpm1d_eig;
             if(succ)
               dg_ere(i,j) = infer; dg_plne(i,j) = Ne; dg_plnx(i,j) = Nx; 
             else 
               disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
             end 
         end
         if(ifsem)
            mxDK = max(se_DK); mnDK = min(se_DK(abs(se_DK) > 1e-10));
            mxDL = max(se_DL); mnDL = min(se_DL(abs(se_DL) > 1e-10));
            sekcd(i,j) = mxDK/mnDK; selcd(i,j) = mxDL/mnDL; 
            disp(['For SE, # of eigs for K is ',num2str(size(se_DK)),...
            ' , max(D) = ',num2str(max(se_DK)),' , min(D) = ',num2str(min(se_DK)),...
            ' , nonzero min = ',num2str(mnDK)]); 
            disp([' , lam_M(K) / lam_m(K) = ',num2str(mxDK/mnDK) ]); 

            disp(['For SE, # of eigs for L is ',num2str(size(se_DL)),...
            ' , max(D) = ',num2str(max(se_DL)),' , min(D) = ',num2str(min(se_DL)),...
            ' , nonzero min = ',num2str(mnDL)]); 
            disp([' , lam_M(L) / lam_m(L) = ',num2str(mxDL/mnDL) ]); 
            se_DL(1:3), 
         end
         if(ifdgm)
            mxDA = max(dg_DA); mnDA = min(dg_DA(abs(dg_DA) > 1e-09));
            mxDL = max(dg_DL); mnDL = min(dg_DL(abs(dg_DL) > 1e-09)); 
            % Note: as Ne, N grows, there comes a point when the "zero" eigenvalue becomes
            % 10^(-9), or bigger. 
            dgkcd(i,j) = mxDA/mnDA; dglcd(i,j) = mxDL/mnDL; 
            disp(['For DG, # of eigs for A is ',num2str(size(dg_DA)),...
            ' , max(D) = ',num2str(max(dg_DA)),' , min(D) = ',num2str(min(dg_DA)),...
            ' , nonzero min = ',num2str(mnDA)]); 
            disp([' , lam_M(K) / lam_m(K) = ',num2str(mxDA/mnDA) ]); 

            disp(['For DG, # of eigs for L is ',num2str(size(dg_DL)),...
            ' , max(D) = ',num2str(max(dg_DL)),' , min(D) = ',num2str(min(dg_DL)),...
            ' , nonzero min = ',num2str(mnDL)]); 
            disp([' , lam_M(L) / lam_m(L) = ',num2str(mxDL/mnDL) ]); 
            disp(['Ratio: cond(K_DG) / cond(K_SE) = ',num2str(dgkcd(i,j)/sekcd(i,j))]); 
            disp(['Ratio: cond(L_DG) / cond(L_SE) = ',num2str(dglcd(i,j)/selcd(i,j))]); 
            dg_DL(1:3),
         end
%        plot(se_DL,'x-'); hold on;  % Both se_DL and dg_DL have zero values, which will not show in a semilogy plot
%        plot(dg_DL,'o-'); 
         semilogy(se_DL(abs(se_DL) > 1e-10),'x-'); hold on;
         semilogy(dg_DL(abs(dg_DL) > 1e-10),'o-'); 
         mylegend{2*j-1}=strcat('SEM, N = ',num2str(N));
         mylegend{2*j  }=strcat('DGM, N = ',num2str(N));
%        legend('SEM','DGM','location','northwest'); 
      end
   end 
end 
if(ifeig)
%   mylegend;
%   legend(mylegend,'location','southeast','Fontsize',15);
    xlabel('k','fontsize',18); ylabel('\lambda','fontsize',18); 
    title(['Eig values of SE and DG, Ne = ',num2str(Ne)],'fontsize',18); 
end
if(ifeig)
% For scaling with Ne, fixed N to be Nn 
    figure(110);
    loglog(se_plne(:,Nn), sekcd(:,Nn),'o-','linewidth',2.0); hold on; 
    loglog(dg_plne(:,Nn), dgkcd(:,Nn),'x-','linewidth',2.0); 
    loglog(dg_plne(:,Nn), 100.*dg_plne(:,Nn).^2,'-','linewidth',2.0); 
    xlabel('Ne','fontsize',18); ylabel('\lambda_{max} / \lambda_{min}','fontsize',18); 
    legend([{'SE'},{'DG'},{'N_e^{2}'}],'fontsize',15);
    title(['Cond(K), for N = ',num2str(N)],'fontsize',18); hold off; 

    figure(111); 
    loglog(se_plne(:,Nn), selcd(:,Nn),'o-','linewidth',2.0); hold on; 
    loglog(dg_plne(:,Nn), dglcd(:,Nn),'x-','linewidth',2.0); 
    loglog(dg_plne(:,Nn), 100.*dg_plne(:,Nn).^2,'-','linewidth',2.0); 
    xlabel('Ne','fontsize',18); ylabel('\lambda_{max} / \lambda_{min}','fontsize',18); 
    legend([{'SE'},{'DG'},{'N_e^{2}'}],'fontsize',15);
    title(['Cond(inv(B)K), for N = ',num2str(N)],'fontsize',18); hold off; 

% For scaling with N , fixed Ne 
%   figure(112);
%   loglog(se_plnx(1,:), sekcd(1,:),'o-','linewidth',2.0); hold on; 
%   loglog(dg_plnx(1,:), dgkcd(1,:),'x-','linewidth',2.0); 
%   loglog(dg_plnx(1,:), dg_plnx(1,:).^3,'-','linewidth',2.0); 
%   xlabel('N','fontsize',18); ylabel('\lambda_{max} / \lambda_{min}','fontsize',18); 
%   legend([{'SE'},{'DG'},{'N^{3}'}],'fontsize',15);
%   title(['Cond(K), for Ne = ',num2str(Ne)],'fontsize',18); hold off; 

%   figure(113); 
%   loglog(se_plnx(1,:), selcd(1,:),'o-','linewidth',2.0); hold on; 
%   loglog(dg_plnx(1,:), dglcd(1,:),'x-','linewidth',2.0); 
%   loglog(dg_plnx(1,:), dg_plnx(1,:).^4,'-','linewidth',2.0); 
%   xlabel('N','fontsize',18); ylabel('\lambda_{max} / \lambda_{min}','fontsize',18); 
%   legend([{'SE'},{'DG'},{'N^{4}'}],'fontsize',15);
%   title(['Cond(inv(B)K), for Ne = ',num2str(Ne)],'fontsize',18); hold off; 
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

%% spatial convergence
%% Fix Ne, varying N 
%if(ifsem)
%  figure(10+9);  
%  semilogy(se_plnx(1,:),se_ere(1,:),'o-','linewidth',1.5);hold on;
%  semilogy(se_plnx(1,:),exp(-3*se_plnx(1,:)),'x-','linewidth',1.5); 
%  legend('Err data','e^{-N}'); 
% %legend('Err data','N^{-2}'); 
%  xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
%  title('Max pointwise relative error, poisson, SEM, N_e = 2');
%end
%if(ifdgm)
%  figure(20+9);  
%  semilogy(dg_plnx(1,:),dg_ere(1,:),'o-','linewidth',1.5);hold on;
%  semilogy(dg_plnx(1,:),exp(-3*dg_plnx(1,:)),'x-','linewidth',1.5); 
%  legend('Err data','e^{-N}'); 
% %legend('Err data','N^{-2}'); 
%  xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
%  title('Max pointwise relative error, poisson, DG, N_e = 2');
%end

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
