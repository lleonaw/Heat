%
%  This code is for spectral (element) solution 
%  explicit method, time step too small not useful, will break down 
%


%for N=1:50;
N = 100;

[Kh,Mh,Ch,Dh,z,w] =  semhat(N);

nu = 0.01;

a=0.; b=1.; x=a+0.5*(b-a)*(z+1.); L=(b-a);
u0 = sin(2.*pi.*x);
c = u0;  % HOW
p = nu.*ones(N+1,1);  % in I3, the diffusion term
                                 %  _
Kb = (2./L)*Dh'*diag(w.*p)*Dh;  %  K=Dh'Bh p(x) Dh Stiffness matrix, no t
Mb = (L/2.)*Mh;                 %  M bar, does not change
Cb = Mh*diag(c)*Dh;              %  C bar, change with time

I1d = eye(N); en = I1d(N,:); 

Q = [en; I1d];               % Boundary matrix
Qt  = Q';                    % Transpose matrix

K=Qt*Kb*Q;                   % Stiffness matrix
M=Qt*Mb*Q;                   % Mass matrix
C=Qt*Cb*Q;                   % Convective matrix , initial

f1 = zeros(N,1);f2 = zeros(N,1);f3 = zeros(N,1);
u1 = zeros(N,1);u2 = zeros(N,1);u3 = zeros(N,1);

dt = 1e-6;
T = 2.0; Nt = round(T/dt);
ifplt = true;
if ifplt
figure; plot(x,u0,'b-');hold on;
xlabel('x'); ylabel('u');
title('u at diff. t');
end 

u = u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
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
    f1 = - C*u1 - K*u1;
    f = e1.*f1 + e2.*f2 + e3.*f3;
    
    g = b1.*u1 + b2.*u2 + b3.*u3;

    u = (g + M\(dt.*f))./b0; % Solve for interior degrees of freedom
    ub= Q*u;    % Extend solution for plot 

    %ue = 0.25*(1-x.*x);

    if ifplt 
        Npt = round(0.05/dt);
        if mod(i,Npt) == 0
            plot(x,ub,'r.-');
            pause(0.01);
        end
    end

end 
if ifplt
   hold off;
end
%end;

