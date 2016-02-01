%%  Inviscid Burgers equation
%  This code is for discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
function succ = bur_dg
    clear all; format long;
    global Ne Nx
    Ne = 256;           % Number of elem
    N = 4;            % Poly. order
    % N = 4, 8 will blow up after round t = 0.27, def. t = 0.3
    Nx = N + 1;       % Numb of points in each elem.

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
    u0 = 0.5 + sin(pi.*x);

    % Flux function
    f = flx_u(u0);

    t = 0.; dt = 1e-6;
    T = 1.5e-0; Nt = round(T/dt);
    %Nt = 1;
    ifplt = false;
    ifplt = true;
    ifflt = true;
    ifflt = false;
    disp(['Ne = ' num2str(Ne) ', N = ' num2str(N) ', dt=' num2str(dt) ]);
    u  = u0;  %
%    disp(u);
%    u  = fltr_u(u);
%    disp(u);
%    error;
    if(ifplt)
        figure(1);hold on;
        for ie = 1:Ne
             plot(x(:,ie),u0(:,ie),'b-');
        end
        xlabel('-- x --'); ylabel('-- u --');
        title(['solution at time t = ' num2str(t) ]);
        hold off;
    end

    for i = 1:Nt
        urh = lh(u,f,Mb,Db);           % Obtain right hand side
        u1  = u + dt.*urh;             % First stage of rk 2
        if(ifflt)
            u1 = fltr_u(u1);
        end
        f1  = flx_u(u1);
        urh = lh(u1,f1,Mb,Db);         % Obtain right hand side
        u   = 0.5*(u + u1 + dt.*urh);  % Second stage of rk 2
        if(ifflt)
            u = fltr_u(u);
        end

        t = t + dt;
        if(ifplt && mod(i,Nt/100)==0)
            clf; figure(1);hold on;
%            ylim([-1.5,2.5]);
            xlim([-1. 1.]);
            for ie = 1:Ne
                plot(x(:,ie),u(:,ie),'b-');
            end
            xlabel('-- x --'); ylabel('-- u --');
            title(['solution at time t = ' num2str(t) ]);
            drawnow
            pause(0.01);
            hold off;
        end
    end
    succ = true;
end
%%
function ufl = fltr_u(u) % filter function, CH 5.6.1, Nodal DG
% not working, the function in the following way is wrong. 
% How to correctly do it?
    global Nx Ne
    utmp = zeros(Nx,Ne);
    for ie=1:Ne
        for ip=1:Nx
            sum = 0.;
            for i=1:(Nx-1)
                tmp = (i-1)/Nx;
                sum = sum + fltr(tmp)*u(i,ie);
            end
            utmp(ip,ie) = sum;
        end
    end
    ufl = utmp; % in case ufl = u, just to make sure it's a copy
end

function sgm = fltr(eta) % filter function, CH 5.6.1, Nodal DG
    alp = 36;
    s = 6;
    sgm = exp(-1.*alp*eta^(s));
end
%%
function flx = lax_frd(fm,fp,um,up,lmd) % weak form
    fctr = (fm + fp)/2.;
    difu = um - up;
    fctr(1,:) = - fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
    flx  = fctr + (lmd.*difu)/2.; % LF does not run
%    flx  = fctr;  % center flux runs
end

function urhs = lh(u,f,M,D)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    [fm,fp] = full2face(f);
    [um,up] = full2face(u);
    lmd = max(max(abs(u)));  % max( df/du )
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne
        urhs(:,ie) = M\(D.'*M*f(:,ie) - full_f(:,ie));
    end
end
%%
function v = face2full(fv)
    global Ne Nx
    v = zeros(Nx,Ne);
    v(1,:)  = fv(1,:);
    v(Nx,:) = fv(2,:);
end
function [vm,vp] = full2face(v)
    global Ne Nx
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    for ie=1:Ne
        if(ie == 1)
            vp(1,ie) = v(Nx,Ne);
            vp(2,ie) = v(1,ie+1);% Periodic
        elseif(ie == Ne)
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,1);
        else
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,ie+1);% Periodic
        end
    end
end
%%
function f = flx_u(u)
    f = u.^2;
end
