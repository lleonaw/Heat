%%  Heat equation, no preferred direction, central flux
%  This code is for spectral element method 
%  
function [succ,infer,VK,DK,VL,DL,sid,plx] = heat_sem
    global Ne Nx ifplt initflg T CFL dt iffwd ifeig
    N = Nx - 1; 
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    a=-1.; b=1.; dl=(b-a)/(Ne);
    Kl = (2./dl)*Kh;  %  K=Dh'Bh p(x) Dh Stiffness matrix, no t
    Ml = (dl/2.)*Mh;                 %  M bar, does not change
%  Connectivity 
    Q  = zeros((N+1)*Ne,N*Ne+1); Qx = Q; 
    for e=1:Ne
        i0 = (N+1)*(e-1)+1; j0 = (N)*(e-1)+1; 
        ae = a + (e-1)*dl;
        xt = ae + (dl/2.).*(z+1.); 
        xl(i0:i0+N,1) = xt;  % local x 
        for ij=1:N+1
          if(ij==1 || ij == N+1 ) 
            Qx(i0+ij-1,j0+ij-1) = 0.5; 
            if( (e==1  && ij==1) || (e==Ne && ij==N+1)) 
              Qx(i0+ij-1,j0+ij-1) = 1.; 
            end 
          else
            Qx(i0+ij-1,j0+ij-1) = 1.; 
          end
        end
        if(e<Ne) 
          for ij=1:N+1
            Q (i0+ij-1,j0+ij-1) = 1.; 
          end
        elseif(e==Ne)
          for ij=1:N+1
            Q(i0+ij-1,j0+ij-1) = 1.; 
          end
        end 
    end
    x = Qx'*xl; 
    nu = 1.*ones(Nx,1); nu = diag(nu); 
% % Exact solution
    if(initflg==1)     % Init Case 1  
      u0 = 1.*cos(pi.*x/2.); qs = ((pi^2)/4.)*cos(pi.*x/2.);
    elseif(initflg==2) % Init Case 2 
%     u0 = -1.*sin(pi.*x); qs = (pi^2)*u0;
      u0 =  1.*sin(pi.*x); qs =  (pi^2)*sin(pi.*x);
    elseif(initflg==3) % Init Case 3 
      u0 = 1.*cos(pi.*x); qs = 1.*(pi^2)*cos(pi.*x);
    elseif(initflg==4) % Init Case 4 
%     u0 = 1.- 1.*cos(2.*pi.*x); qs = -((pi^2)*4.).*cos(2.*pi.*x);
      u0 = 1.+ 1.*cos(pi.*x); qs = ((pi^2)).*cos(pi.*x);
    elseif(initflg==5) % Init Case 5, discontinuous
      u0 = rand(size(x)); 
    end

    Ie = sparse(eye(Ne)); 
    Kb = Q'*(kron(Ie,Kl))*Q;                   % Stiffness matrix
    Mb = Q'*(kron(Ie,Ml))*Q;                   % Mass matrix
    Iq = eye(N*Ne); P=[Iq;Iq(1,:)];                    % Periodic bc 
%   Iq = eye(N*Ne-1); Q=[0.*Iq(1,:);Iq;0.*Iq(1,:)];    % Dirichlet bc 

    Qr = Q*P; 
    K  = Qr'*kron(Ie,Kl)*Qr; 
    M  = Qr'*kron(Ie,Ml)*Qr;
%
%   disp(size(u0)); disp(size(K)); disp(size(P)); 
%   u = 0.*u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
%   u = K \ (M * R'*qs); 
%   u = R*u;     % Ill-conditioned for periodic
                 % Looks ok for dircichlet 
%
%  periodicity is more like connectivity rather than boundary condition 
%   K, M = full(M); 
    if(ifeig) 
      [VK, DK] = eig(K);    % Stiffness matrix 
      DK = sort(diag(DK)); 
      [VL, DL] = eig(K,M);  % No periodic mask ....  %   [VL, DL] = eig(Kb,Mb); 
%     idl = find( abs(diag(DL)) <12. & abs(diag(DL)) >9.); 
      DL = diag(DL);
      [DL,sid] = sort(DL); 
      plx = x(1:end-1);
%     for invd=1:length(DL)
%       invd, DL(invd), 
%       plot(x(1:end-1),VL(:,sid(invd)));  pause, 
%     end
%     plot(x,VL(:,idl)); % pause 
      if(iffwd)
        infer = 8.88e20; 
        uxx0 = pi^2*u0; 
        uxx  = inv(M)*K*(P'*u0); 
        infer = norm(uxx0-P*uxx,Inf); 
        disp(['SEM :: Error in u_{xx} = ', num2str(infer)...
               ,', N = ',num2str(N),' , N_e = ',num2str(Ne)]); 
      end
    else 
      VK = zeros(size(K)); DK = zeros(size(diag(K)));
      VL = zeros(size(K)); DL = zeros(size(diag(K)));
    end
% 
% % % % % % % % % % % % % % % % % % % % % 
% % % Time stepping setup 
    t = 0.; dxa = diff(x(:,1)); dx = min(dxa); 
    dt = CFL*(dx)^2 /nu(1,1);  % Alright 0.1 is as big as it gets 4 me 
    if(dt>=T)
      error('dt > T, change to a smaller CFL or longer T');
    end
    Nt = round(T/dt); dt= T/double(Nt); 
    disp(['Init case = ', num2str(initflg),...
          ', Ne = ',num2str(Ne),' , N = ', num2str(N),... 
          ', Nt = ',num2str(Nt),... 
          ' ,dt =',num2str(dt),', CFL = ',num2str(CFL)]);
    u   = u0; tsv = zeros(Nt,1); infert = tsv; 
%
    for it= 1:Nt
%       % Euler 1 
        urh = inv(M)*K*(P'*u); % Obtain right hand side
        urh = -nu(1).*urh; 
        u   = u + dt.*(P*urh);             % First stage of rk 2
        t = t + dt;

        if(initflg==2) 
          uex = exp(-nu(1).*(pi)^2.*t).*sin(pi.*x);
          mxtlu   = norm(uex,Inf); 
        elseif(initflg==5)  % This is not the exact! !  
          uex = exp(-nu(1).*(pi)^2.*t).*u0;
          uex = 0.*uex; 
          mxtlu   = 1.;
        end 

        if(ifplt && mod(it,ceil(Nt/10))==0)
%       if(ifplt ) % && mod(it,ceil(Nt/20))==0)
            % clf; 
            figure(1);hold on;
%            ylim([-0.1,1.1]);
            xlim([-1. 1.]);
            plot(x,uex,'x-'); hold on;
            plot(x,u,'o-');
            xlabel('-- x --'); ylabel('-- u --');
            title(['solution at time t = ' num2str(t) ]);
            drawnow
%           pause;
            hold off;
            disp(['max(u) = ' num2str(max(max(u))) ...
               ' , min(u) = ' num2str(min(min(u))) ' , istep = ' num2str(it)]);
        end
        er = u - uex; 
        tsv(it) = t; 
%       mxtlu   = norm(uex,Inf); infert(it) = norm(er,Inf)/ mxtlu; 
        infert(it) = norm(er,Inf)/ mxtlu; 
    end
    if(ifplt)
      figure(10);hold on;
      plot(tsv,infert,'-');
      xlabel('t '); 
      ylabel('$\| u - \tilde{u}\|_{\infty} / \|\tilde{u}\|_{\infty}$','Interpreter','Latex');
      title(['Error vs time']);
      drawnow ; 
      hold off;
    end
    infer = infert(end); % Error at end of time 
    disp(['SEM :: At end of time T = ',num2str(T),...
    ', Relative Inf norm error = ', num2str(infer)]); 
    succ = true;
end

