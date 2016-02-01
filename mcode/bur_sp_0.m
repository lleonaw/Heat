%  Viscous Burgers equation 
%  This code is for spectral (element) solution 
%  implicit scheme for diffusion term

clear all; format long;
N = 256;
[Kh,Mh,Ch,Dh,z,w] =  semhat(N);
nu = 0.01;
nu = 0.005;
%nu = 0.; % break after T = 0.3. Need viscosity nu = 0.005.

a=-1.; b=1.; x=a+0.5*(b-a)*(z+1.); L=(b-a);

% u0 = sin(2.*pi.*x);
% u0 = 0.5 + sin(pi.*x);
u0 = -sin(pi.*x); 
c = u0;               % HOW
p = nu.*ones(N+1,1);  % in I3, the diffusion term
                                   %  _
Kb = (2./L)*Dh'*diag(w.*p)*Dh;  %  K=Dh'Bh p(x) Dh Stiffness matrix, no t
Mb = (L/2.)*Mh;                 %  M bar, does not change
Cb = Mh*diag(c)*Dh;             %  C bar, change with time

I1d = eye(N); en = I1d(N,:); 
Q = [en; I1d];               % Boundary matrix
Qt  = Q';                    % Transpose matrix

K=Qt*Kb*Q;                   % Stiffness matrix
M=Qt*Mb*Q;                   % Mass matrix
C=Qt*Cb*Q;                   % Convective matrix , initial

dt = 1e-5;
T = 1.5; 
T = 10./pi; 
Nt = round(T/dt);
disp(['N = ' num2str(N) ' , dt = ' num2str(dt)]); 

u = u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
f1 = zeros(N,1);f2 = zeros(N,1);f3 = zeros(N,1);
u1 = zeros(N,1);u2 = zeros(N,1);u3 = zeros(N,1);

t = 0.; 

ifplt = true; 
if(ifplt) 
    figure(1);
    plh = plot(x,u0,'b-');
    xlabel('--x--'); ylabel('--u--');
    title(['solution at time t = ' num2str(t) ]);
    drawnow
    pause(0.01);
end 

for i = 1: Nt
    if i == 1
        b0 = 1.; b1 = 1.; b2 = 0.; b3 = 0.;
        e1 = 1.; e2 = 0.; e3 = 0.;
    elseif i == 2
        b0 = 1.5; b1 = 2.; b2 = -0.5; b3 = 0.;
        e1 = 2.; e2 = -1.; e3 = 0.;
    elseif i == 3
        b0 = 11./6.; b1 = 3.; b2 = -1.5; b3 = 1./3.;
        e1 = 3.; e2 = -3.; e3 = 1.;
    end

    u3 = u2; u2 = u1; u1 = u;
    f3 = f2; f2 = f1; 
    c1 = Q*u1; Cb = Mh*diag(c1)*Dh; C = Qt*Cb*Q;
    f1 = - C*u1;
    f = e1.*f1 + e2.*f2 + e3.*f3;
    g = b1.*M*u1 + b2.*M*u2 + b3.*M*u3;
    
    u = (b0.*M + dt.*K)\(g + (dt.*f)); % Solve for interior d.o.f
    ub= Q*u;    % Extend solution for plot 
    t = t + dt;
    if(ifplt && mod(i,Nt/100)==0) 
        figure(1);
        plot(x,ub,'b-');
        xlim([a,b]);
        %ylim([-1.,1.]);
        title(['solution at time t = ' num2str(t) ]);
        xlabel('--x--'); ylabel('--u--');
        drawnow
        pause(0.01);
    end 

end
    
