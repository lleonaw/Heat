%  Viscous Burgers equation 
%  This code is for spectral (element) solution 
%  implicit scheme for diffusion term

clear all; format long;
N = 256;
N = 512;  % 512 the solution of du/dx at origin is less than 150, dt = 1/(pi 1000)
% N = 1024; % dt = 1/(1000000 pi) still will get singular matrix, what's wrong? 
N = 1024;
[Kh,Mh,Ch,Dh,z,w] =  semhat(N);
if(N > 513) % Kh, Dh, Ch will have NaN in them
    disp(['High polynomial order' num2str(N) ' , go into fortran code']);
    fileID = fopen('poly_N.txt','w');
    fprintf(fileID,'%d',N);
    fclose(fileID);
    !./sptest
    fileID2 = fopen('Dx.txt','r');
    Dh1 = fscanf(fileID2,'%f');
    fclose(fileID2);
    Nz = N+1;
    Dh = reshape(Dh1,Nz,Nz);
    Kh = Dh'*Mh*Dh;  %  1-D stiffness matrix
    Ch = Mh*Dh;      %  1-D convection operator
end
% error('Testing fortran');
%nu = 0.01;
%nu = 0.005;
nu = 1./(100.*pi);
%nu = 0.; % break after T = 0.3. Need viscosity nu = 0.005

a=-1.; b=1.; x=a+0.5*(b-a)*(z+1.); L=(b-a);
No = find(-1e-8 < x & x<1e-8);
%error

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

% dt = 1e-5;
dt = 1/(10000*pi);
% T = 1.5; 
T = 10./pi; 
Nt = round(T/dt); Nio = int32(1/(pi*dt));
disp(['N = ' num2str(N) ' , dt = ' num2str(dt)]); 

u = u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
f1 = zeros(N,1);f2 = zeros(N,1);f3 = zeros(N,1);
u1 = zeros(N,1);u2 = zeros(N,1);u3 = zeros(N,1);
t = 0.; 
ifplt = true; 
if(ifplt) 
    figure(1);
    plh = plot(x,u0,'b-'); hold on;
    xlabel('--x--'); ylabel('--u--');
    title(['solution at time t = ' num2str(t) ]);
    legendinfo{1} = ['t = 0'];
    drawnow
    pause(0.01);
end 
umx = zeros(Nt,1);xmx=zeros(Nt,1); pt = zeros(Nt,1);
dudxo = zeros(Nt,1);%legendinfo = cell(11,1);
il = 2;
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
    
    pt(i) = t;
    [umx(i) mxid] = max(abs(ub));
    xmx(i) = abs(x(mxid));
    du = Dh*ub;
    dudxo(i) = (2./L)*du(No);dudxo(i) = abs(dudxo(i));
    
%    if(ifplt && mod(i,Nio)==0) 
    if(ifplt && mod(i,Nio)==0 && (i/Nio==2 || i/Nio == 3 || i/Nio == 10)) 
        figure(1); 
        plot(x,ub);
        xlim([a,b]);
        %ylim([-1.,1.]);
        title(['solution at time t = ' num2str(t) ]);
        legendinfo{il} = ['t = ' num2str(i/Nio) '/pi']; 
        il = il + 1;
        xlabel('--x--'); ylabel('--u--');
        drawnow
        pause(0.01);
    end 
end

figure(1);legend(legendinfo);
figure(3);plot(pt,dudxo);title('Derivative at origin');
xlabel('-- t --'); ylabel('$(\frac{du}{dx})_{x=0}$','Interpreter','Latex');
figure(4);plot(pt,umx);title('u_{max} with time');
xlabel('-- t --'); ylabel('$u_{max}$','Interpreter','Latex');
figure(5);plot(pt,xmx);title('x_{max} = Abscissa of u_{max}');
xlabel('-- t --'); ylabel('$x_{max}$','Interpreter','Latex');
