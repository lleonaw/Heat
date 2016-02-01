%%  Viscous/ Inviscid Burgers equation
%  This code is for discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
clear all; format long;
global Ne Nx
Ne = 512;           % Number of elem
N = 1;            % Poly. order, high order will break t = 0.18
Nx = N + 1;       % Numb of points in each elem.
No = Ne/2; % For even Ne ! ! !

[Kh,Mh,Ch,Dh,z,w] =  semhat(N);
% Init Geom
a=-1.; b=1.;
dl = (b - a)/(Ne);
x = zeros(Nx,Ne);
for ie = 1: Ne
    ae = a + (ie-1)*dl;
    x(:,ie) = ae+0.5*dl*(z+1.);
end
Mb = (dl/2.)*Mh; % dx /d kexi  %  M bar, does not change
Db = (2./dl)*Dh; % d kexi/ dx

% Init scalar
%    u0 = 0.5 + sin(pi.*x);
u0 = - sin(pi.*x);

t = 0.;
% dt = 1e-5;
dt = 1/(400*pi);
%T = 1.5e-0;
T = 10/pi; Nio = int32(1/(pi*dt));
Nt = round(T/dt);

umx = zeros(Nt,1); % Save the max velocity
xmx = zeros(Nt,1); % Save the abscissca
duo = zeros(Nt,1); % Save du at origin
pt = zeros(Nt,1);
%Nt = 1;
%nu = 0.005;
nu = 1./(100.*pi);
ifvis = false;
ifvis = true; % if viscous
ifplt = false;
ifplt = true;
ifflt = true;
ifflt = false;
disp(['Ne = ' num2str(Ne) ', N = ' num2str(N) ', dt=' num2str(dt) ]);
if(ifflt)
    disp(' - - - - Filter on - - - - ');
else
    disp(' - - - - Filter off - - - - ');
end
if(ifvis)
    disp([' - - - - Viscosity on, nu = ' num2str(nu) ' . - - - - ']);
else
    disp([' - - - - Viscosity off - - - - ']);
end
if(ifplt)
    linestyle = {['b-'],['r-'],['k-'],['b:'],['r:']};
    figure(2);hold on;
    for ie = 1:Ne
        plot(x(:,ie),u0(:,ie),linestyle{1});
    end
    xlabel('-- x --'); ylabel('-- u --');
    title(['solution at time t = ' num2str(t) ]);
    drawnow
    pause(0.01);legendinfo{1} = ['t = 0'];
    il = 2;
end
u  = u0;  %
for i = 1:Nt
    f = flx_u(u);
    if(ifvis)
        urh = lh_vis(u,f,Mb,Db,nu);    % Obtain right hand side
    else
        urh = lh(u,f,Mb,Db);           % Obtain right hand side
    end
    u1  = u + dt.*urh;             % First stage of rk 2
    if(ifflt)
        u1 = fltr_u(u1);
    end
    f1  = flx_u(u1);
    if(ifvis)
        urh = lh_vis(u1,f1,Mb,Db,nu);
    else
        urh = lh(u1,f1,Mb,Db);     % Obtain right hand side
    end
    u   = 0.5*(u + u1 + dt.*urh);  % Second stage of rk 2
    if(ifflt)
        u = fltr_u(u);
    end
    
    %Find du at origin & max of U
    um1 = 0.5*(u(1,No) + u(2,No-1));
    up1 = 0.5*(u(2,No+1) + u(1,No+2));
    duo(i) = 0.5*(up1 - um1)/dl;
    duo(i) = abs(-duo(i));
    % approach to find the max index in a 2D array)
    [umx(i),mxid] = max(abs(u(:))); % 1D indexing
    xmx(i) = abs(x(mxid));
    
    t = t + dt;
    pt(i) = t;
    
    %if(ifplt && mod(i,Nt/25)==0)
    if(ifplt && mod(i,Nio)==0 && (i/Nio==2 || i/Nio == 3 || i/Nio == 10))
        figure(2);
        %            ylim([-1.5,2.5]);
        xlim([-1. 1.]);hold on;
        for ie = 1:Ne
            plot(x(:,ie),u(:,ie),linestyle{il});
        end
        xlabel('-- x --'); ylabel('-- u --');
        title(['solution at time t = ' num2str(t) ]);
        legendinfo{il} = ['t = ' num2str(i/Nio) '/pi'];
        il = il + 1;
        drawnow
        pause(0.01);
        hold off;
        disp(['max(u) = ' num2str(max(max(u))) ...
            ' , min(u) = ' num2str(min(min(u))) ' , istep = ' num2str(i)]);
    end
end
%figure(2);legend(legendinfo);
figure(3);plot(pt,duo);title('Derivative at origin');
xlabel('-- t --'); ylabel('$(\frac{du}{dx})_{x=0}$','Interpreter','Latex');
figure(4);plot(pt,umx);title('u_{max} with time');
xlabel('-- t --'); ylabel('$u_{max}$','Interpreter','Latex');
figure(5);plot(pt,xmx);title('x_{max} = Abscissa of u_{max}');
xlabel('-- t --'); ylabel('$x_{max}$','Interpreter','Latex');

ufinal = u;

