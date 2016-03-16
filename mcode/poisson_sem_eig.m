%  From bur_sp, to do poisson with sem 
%  This code is for spectral (element) solution 
function [succ,infer,VL,DL] = poisson_sem_eig
    global Ne Nx ifplt initflg T CFL dt
    N = Nx - 1; 
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    a=-1.; b=1.; dl=(b-a)/(Ne);
    xtm = zeros(Nx,Ne);
    plx   = reshape(xtm,Nx*Ne,1); 
    Kl = (2./dl)*Kh;  %  K=Dh'Bh p(x) Dh Stiffness matrix, no t
    Ml = (dl/2.)*Mh;                 %  M bar, does not change
    P  = zeros((N+1)*Ne,N*Ne+1); Px = P; 
    for e=1:Ne
        i0 = (N+1)*(e-1)+1; j0 = (N)*(e-1)+1; 
        ae = a + (e-1)*dl;
        xt = ae + (dl/2.).*(z+1.); 
        xl(i0:i0+N,1) = xt;  % local x 
        for ij=1:N+1
          if(ij==1 || ij == N+1 ) 
            Px(i0+ij-1,j0+ij-1) = 0.5; 
            if( (e==1  && ij==1) || (e==Ne && ij==N+1)) 
              Px(i0+ij-1,j0+ij-1) = 1.; 
            end 
          else
            Px(i0+ij-1,j0+ij-1) = 1.; 
          end
        end
        if(e<Ne) 
          for ij=1:N+1
            P (i0+ij-1,j0+ij-1) = 1.; 
          end
        elseif(e==Ne)
          for ij=1:N+1
            P(i0+ij-1,j0+ij-1) = 1.; 
          end
        end 
    end
    R = P'; Rx = Px'; xw = Rx*xl; x = xw;

    nu = 1.*ones(Nx,1);  nu = diag(nu); 
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
    end

    Ie = sparse(eye(Ne)); 
    Kb=R*(kron(Ie,Kl))*P;                   % Stiffness matrix
    Mb=R*(kron(Ie,Ml))*P;                   % Mass matrix

    Iq = eye(N*Ne); Q=[Iq(1,:);Iq];                    % Periodic bc 
%   Iq = eye(N*Ne-1); Q=[0.*Iq(1,:);Iq;0.*Iq(1,:)];    % Dirichlet bc 
    K = Q'*Kb*Q; M = Q'*Mb*Q;          % Apply bc 
%   disp(Kb); disp(P); disp(K); 
%   disp(Nx); disp(Ne); 

%   disp(size(u0)); disp(size(K)); disp(size(P)); 
%   u = 0.*u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
%   u = K \ (M * Q'*qs); 
%   u = Q*u;     % Ill-conditioned for periodic
                 % Looks ok for dircichlet 

% Alternate 
    L = Q*( K \ (M * Q') ) ; 
    u = L*qs; 

    S = Q*( M \ K ) * Q';
    [VL, DL] = eig(S); 
    DL = sort(diag(DL)); 

    infer = 8.88e20; 
    uex = u0; 
    infer = norm(u-uex,Inf); 
    disp(['SEM :: Error in solution = ', num2str(infer)...
         ,', N = ',num2str(N),' , N_e = ',num2str(Ne)]); 

    if(ifplt)
        figure(20+initflg); 
        plot(x,u,'x-','linewidth',2.5'); hold on; 
        plot(x,uex,'-','linewidth',2.5'); 
    end
    succ = true;
end
