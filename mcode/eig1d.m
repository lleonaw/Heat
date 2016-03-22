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
Nen = 4; 

dg_ere  = zeros(Nen,Nn); dg_plnx = zeros(Nen,Nn); dg_plne = zeros(Nen,Nn);
se_ere  = zeros(Nen,Nn); se_plnx = zeros(Nen,Nn); se_plne = zeros(Nen,Nn);

ertn = zeros(Ntn,1); dtn  = zeros(Ntn,1); 

if(ifeig)
    figure(43);
    mylegend = cell(2*Nn,1);  % For plotting eig 
end
for initflg=2:2 % 4
   disp(['Init flg = ',num2str(initflg)]); 
   for j=1:Nn
%     N = 2^j; Nx = N + 1;
      N = 3*j; Nx = N + 1;
      for i=1:1%Nen
         Ne = (i+1)^2;
         if(ifsem)
             [succ,infer,se_VL,se_DL] = sem1d_eig; 
             if(succ)
               se_ere(i,j) = infer; 
               se_plne(i,j) = Ne; 
               se_plnx(i,j) = Nx; 
             else 
               disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
             end 
         end
         if(ifdgm)
             [succ,infer,dg_VL,dg_DL] = dgpm1d_eig;
             if(succ)
               dg_ere(i,j) = infer; 
               dg_plne(i,j) = Ne; 
               dg_plnx(i,j) = Nx; 
             else 
               disp(['Failed! Ne =',num2str(Ne),' , N =',num2str(N)]); 
             end 
         end
      end 
      if(ifeig)
          if(ifsem)
              disp(['For SE, # of eigs is ',num2str(size(se_DL)),' , max(D) = ',num2str(max(se_DL)),' , min(D) = ',num2str(min(se_DL))]); 
              se_DL(1:3), 
          end
          if(ifdgm)
              disp(['For DG, # of eigs is ',num2str(size(dg_DL)),' , max(D) = ',num2str(max(dg_DL)),' , min(D) = ',num2str(min(dg_DL))]); 
              dg_DL(1:3),
          end
          plot(se_DL,'x-'); hold on;  % Both se_DL and dg_DL have zero values, which will not show in a semilogy plot
          plot(dg_DL,'o-'); 
%         semilogy(se_DL,'x-'); if(j==1) hold on; end;
%         semilogy(dg_DL,'o-'); 
          mylegend{2*j-1}=strcat('SEM, N = ',num2str(N));
          mylegend{2*j  }=strcat('DGM, N = ',num2str(N));
%         legend('SEM','DGM','location','northwest'); 
      end
   end 
end 
if(ifeig)
%   mylegend;
    legend(mylegend,'location','southeast','Fontsize',15);
    xlabel('k','fontsize',18); ylabel('\lambda','fontsize',18); 
    title(['Eig values of SE and DG, Ne = ',num2str(Ne)],'fontsize',18); 
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
 if(ifsem)
   figure(10+9);  
   semilogy(se_plnx(1,:),se_ere(1,:),'o-','linewidth',1.5);hold on;
   semilogy(se_plnx(1,:),exp(-3*se_plnx(1,:)),'x-','linewidth',1.5); 
   legend('Err data','e^{-N}'); 
  %legend('Err data','N^{-2}'); 
   xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
   title('Max pointwise relative error, poisson, SEM, N_e = 2');
 end
 if(ifdgm)
   figure(20+9);  
   semilogy(dg_plnx(1,:),dg_ere(1,:),'o-','linewidth',1.5);hold on;
   semilogy(dg_plnx(1,:),exp(-3*dg_plnx(1,:)),'x-','linewidth',1.5); 
   legend('Err data','e^{-N}'); 
  %legend('Err data','N^{-2}'); 
   xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
   title('Max pointwise relative error, poisson, DG, N_e = 2');
 end

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
