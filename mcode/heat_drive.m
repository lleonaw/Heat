%%  1D Heat equation, 
%%  Driver for heat_prim.m, heat_sem.m 
% To do:
%   -x Stab., max dt of DG is 3.33 times of SEM
% % % % % % % % % % % % % % % % % % % % % %
clear all; format longE;

global Ne Nx ifplt initflg T CFL dt iffwd ifeig
dt = 1e-5;    % T = 4.e-0; 
ifplt = true; 
ifplt = false; 

ifsem = true; 
%ifsem = false; 
ifdgm = true; 
ifdgm = false; 

ifeig = true;
iffwd = false; 

T = .20; 
% Ne = 10;          % Number of elem
N = 4;            % Poly. order
Nx = N + 1;       % Numb of points in each elem.

Ntn = 4; 
Nn  = 8; 
% Nn  = 4;  % 16 is enough... 32 takes some time, 64 takes a while ! ! 
Nen = 4; 

ere  = zeros(Nen,Nn); plnx = zeros(Nen,Nn); plne = zeros(Nen,Nn);
ertn = zeros(Ntn,1); dtn  = zeros(Ntn,1); 

for itm=1:Ntn
%   CFL = 0.8 /(2.^itm);    % SEM: between 0.2, 0.4 
%                           % DGM: between 0.1, 0.2 
%   CFL = 0.1*(1.2 - 0.02*itm);     % DGM: between 0.116, 0.118, safe bet 0.114
%   CFL = 0.114;                    % CFL = 0.114, Error 0.0035584 
%   CFL = 0.112;                    % Safer bet, CFL = 0.112, discontinuous IC 
%   CFL = 0.1*(2.10 - 0.02*itm);    % SEM: between 0.220, 0.222, safe bet 0.218 
                                    % SEM: between 0.204, 0.208, safe bet 0.204 
%   CFL = 0.202;                    % CFL = 0.210, Error 0.0065674 
    for initflg=5:5 %  5- discontinuous inic cond., HOWEVER they are not EXACTLy the same 
                    %  rand function called twice
       for j=3:3 %Nn
       % Use CFL = 0.0125 for polynomial  
%         N = j; Nx = N + 1;
          N = 2*j; Nx = N + 1;
%         N = 2^j; Nx = N + 1;
          for i=1:1%Nen
             Ne = i*6;
    %        Ne = i;
             if(ifsem)
               tic; [succ,infer,se_VK,se_DK,se_VL,se_DL,se_sid,se_plx] = heat_sem; toc; 
               if(succ)
                 se_ere(i,j) = infer; se_plne(i,j) = Ne; se_plnx(i,j) = Nx; 
                 if(ifeig)
                     mxDL = max(se_DL); mnDL = min(se_DL(abs(se_DL)>1e-9)); 
                     scL = mxDL / mnDL; 
                     disp(['SEM :: size of DL',num2str(size(se_DL)), ' , cond(L) = ',num2str(scL)]); 
                 end
               else 
                 disp(['Failed in SE! Ne =',num2str(Ne),' , N =',num2str(N)]); 
               end 
             end 
             if(ifdgm)
               tic; [succ,infer,dg_VA,dg_DA,dg_VL,dg_DL,dg_sid,dg_plx] = heat_prim; toc;
               if(succ)
                 dg_ere(i,j) = infer; dg_plne(i,j) = Ne; dg_plnx(i,j) = Nx; 
                 if(ifeig)
                     mxDL = max(dg_DL); mnDL = min(dg_DL(abs(dg_DL)>1e-9)); 
                     dcL = mxDL / mnDL; 
                     disp(['DGM :: size of DL',num2str(size(dg_DL)), ', cond(L) = ',num2str(dcL)]);
                     if(ifsem)
                         disp(['Ratio :: cond(L_DG)/cond(L_SE) = ',num2str(dcL/scL)]); 
%                        semilogy(se_DL(abs(se_DL) > 1e-9),'x-'); hold on;
%                        semilogy(dg_DL(abs(dg_DL) > 1e-9),'o-'); 
                     end
                 end
               else 
                 disp(['Failed in DG! Ne =',num2str(Ne),' , N =',num2str(N)]); 
               end 
             end 
          end 
       end 
    end 
    ertn(itm) = infer; 
    dtn (itm) = dt; 
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
% figure(9);  
% semilogy(plnx(1,:),ere(1,:),'o-','linewidth',1.5);hold on;
% semilogy(plnx(1,:),exp(-3*plnx(1,:)),'x-','linewidth',1.5); 
% legend('Err data','e^{-N}'); 
%%legend('Err data','N^{-2}'); 
% xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');
% title('Max pointwise relative error, heat problem, N_e = 2');

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
