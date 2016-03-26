%%  1D Heat equation, 
%%  After heat_drive.m(Driver for heat_prim.m, heat_sem.m)
%%  plot the eigenvectors. 
% To do:
% % % % % % % % % % % % % % % % % % % % % %

len_s = length(se_DL); 
%for j=1:ceil(len/5)
%    figure(j);  
%    mylg = cell(5,1);
%    for i=1:5
%        li = i + 5*(j-1); 
%        if(li<=len) 
%          plot(se_plx,se_VL(:,se_sid(li)),'o-','linewidth',1.5);hold on;
%          mylg{i} = strcat('k=',num2str(li));
%        end
%    end
%    legend(mylg);
%    title('Eigenvectors');
%end
%legend('cos(\pi x/2)','sin(\pi x)','cos(\pi x)','1- cos(2 \pi x)','N^{-2}','e^{-N}'); 
%xlabel('$N$','Interpreter','Latex'); ylabel('$\|u - \tilde{u}\|_{\infty}$','Interpreter','Latex');

len_d = length(dg_DL); 
%nj = ceil(len/4);
%for j=1:nj
%    figure(j);  
%    mylg = cell(4,1);
%    for i=1:4
%        li = i + 4*(j-1); 
%        if(li<=len) 
%          plot(dg_plx,dg_VL(:,dg_sid(li)),'o-','linewidth',1.5);hold on;
%          mylg{i} = strcat('k=',num2str(li));
%        end
%    end
%    legend(mylg);
%    title('DG::eigenvectors');
%end

% Plot the modes that are in DG, not in SE 
% % % % % % % % % % % % % % % % % % % % 
nj = ceil((len_d-len_s)/1);
mn = ceil(sqrt(nj)); nn = nj/mn; 
%figure(100);  
%for j=1:nj
%    subplot(mn,nn,j); 
%    li = len_s + j;
%%   mylg = cell(1,1);
%    plot(dg_plx,dg_VL(:,dg_sid(li)),'o-','linewidth',1.5);hold on;
%%   mylg{1} = strcat('k=',num2str(li));
%%   legend(mylg,'fontsize',18);
%    title(['DG::eigenvector, k=',num2str(li)],'fontsize',18);
%end
%saveas(gcf,'./plts_eig/1d/eigv_DGonly.eps','epsc');

% Shared spectrum
nf = ceil(len_s/nj);  % how many figures 
for f=1:nf
    figure(f);
    for i=1:nj
        subplot(mn,nn,i); 
        li = i + nj*(f-1);
        mylg = cell(2,1);
        plot(dg_plx,dg_VL(:,dg_sid(li)),'o-','linewidth',1.5);hold on;
        plot(se_plx,se_VL(:,se_sid(li)),'x-','linewidth',1.5);
        mylg{1} = strcat('DG'); mylg{2} = strcat('SE');
        legend(mylg,'fontsize',18);
        title(['SE & DG::k=',num2str(li)],'fontsize',18);
    end
%   eval(['print -epsc eigv_DGSE_' num2str(f) '.eps']);
%   print(['eigv_DGSE_' num2str(f) '.eps']);
%   set(gcf,'name',['./plts_eig/1d/eigv_DGSE_',num2str(f),'.eps']);
%   saveas(gcf,'name','epsc');
end
