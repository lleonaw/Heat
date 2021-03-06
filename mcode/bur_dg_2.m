%%  Viscous/ Inviscid Burgers equation
%  This code is for discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
function succ = bur_dg
    clear all; format long;
    global Ne Nx
    Ne = 256;           % Number of elem
    N = 1;            % Poly. order, high order will break t = 0.18
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

    t = 0.; dt = 1e-5;
    T = 7.5e-1; Nt = round(T/dt);
    %Nt = 1;
    nu = 0.05; 
    ifvis = false; % if viscous
    ifvis = true;
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
        figure(2);hold on;
        for ie = 1:Ne
             plot(x(:,ie),u0(:,ie),'b-');
        end
        xlabel('-- x --'); ylabel('-- u --');
        title(['solution at time t = ' num2str(t) ]);
        drawnow
        pause(0.01);
        hold off;
    end
    u  = u0;  %
    for i = 1:Nt
        f = flx_u(u);
        if (ifvis)
           urh = lh_vis(u,f,Mb,Db,nu);    % Obtain right hand side
        else
           urh = lh(u,f,Mb,Db);           % Obtain right hand side
        end
        u1  = u + dt.*urh;             % First stage of rk 2
        if(ifflt)
            u1 = fltr_u(u1);
        end
        f1  = flx_u(u1);
        if (ifvis) 
            urh = lh_vis(u1,f1,Mb,Db,nu);
        else
            urh = lh(u1,f1,Mb,Db);     % Obtain right hand side
        end
        u   = 0.5*(u + u1 + dt.*urh);  % Second stage of rk 2
        if(ifflt)
            u = fltr_u(u);
        end

        t = t + dt;
        if(ifplt && mod(i,Nt/25)==0)
            clf; figure(2);hold on;
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
            disp(['max(u) = ' num2str(max(max(u))) ...
               ' , min(u) = ' num2str(min(min(u))) ' , istep = ' num2str(i)]);
        end
    end
    succ = true;
end
%%
function ufl = fltr_u(u) % filter function, CH 5.6.1, Nodal DG
    global Nx Ne
    utmp = zeros(Nx,Ne);
    Np = Nx - 1;
    Nc = round(0.9*Np);
    ftr = filter1d(Nx,Nc,Np);
    for ie = 1:Ne
        utmp(:,ie) = ftr*u(:,ie);
    end
    ufl = utmp; % in case ufl = u, just to make sure it is a copy
end
function ftr = filter1d(Nx,Nc,Np) %
    fldiag = ones(Nx,1);
    N = Nx - 1;
    for i=Nc:N
        fldiag(i+1) = fltr((i-Nc)/(N-Nc));
    end
    [V, invV] = vandm(Nx,Np); % 1D vandermonde matrix
    ftr = V*diag(fldiag)*invV;
end
function [V1D,invV1D] = vandm(Nr,Np)
% Build 1D vandermonde matrix,
%        Nr = #. x in one element
%        Np = #. poly. order that cut off
     N = Nr - 1;
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    V1D = Vandermonde1D(Np,z);
    invV1D = inv(V1D);
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
function flx = ctr_flx(fm,fp) % weak form
    fctr = (fm + fp)/2.;
%    difu = um - up;
    fctr(1,:) = - fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
    flx  = fctr;  % center flux runs
end
function urhs = lh_vis(u,f,M,D,nu)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    [fm,fp] = full2face_p(f); % Change routine call for different bc 
    [um,up] = full2face_p(u);
    lmd = max(max(abs(u)));  % max( df/du )
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne
        % weak form, viscous
         urhs(:,ie) = M\(D'*M*f(:,ie) - full_f(:,ie)); 
    end

    nu_ary = nu.*ones(Nx,1);
    lf_vf = lh_heat(u,M,D,nu_ary);

    urhs = urhs + lf_vf; % Adding 2 parts.
end
function urhs = lh(u,f,M,D)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    [fm,fp] = full2face_p(f);
    [um,up] = full2face_p(u);
    lmd = max(max(abs(u)));  % max( df/du )
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne
        % weak form, inviscid
        urhs(:,ie) = M\(D.'*M*f(:,ie) - full_f(:,ie)); 
    end
end
function urhs = lh_heat(u,M,D,a) % Nodal DG Ch. 7.1
    sqa = sqrt(a); 
    vq = eval_q(u,M,D,sqa); % volumetric array
    urhs = eval_fu(vq,M,D,sqa);
end
%%
function urhs = eval_fu(q,M,D,sa)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    f = diag(sa)*q;

% -- 1. -- 
%   [fm,fp] = full2face(f);
%   [fp] = bc_q(fp,fm); % Add bc to q
% -- 2. -- 
    [fm,fp] = full2face_p(f);

    ctr_f = ctr_flx(fm,fp);
    full_f = face2full(ctr_f);
    for ie = 1:Ne
        % weak form
        urhs(:,ie) = M\(full_f(:,ie) - D'*M*f(:,ie)); 
    end
end
function urhs = eval_q(u,M,D,sa)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    f = diag(sa)*u;
% -- 1. -- 
%   [fm,fp] = full2face(f);
%   [fp] = bc_fu(fp,fm); % add bc to flux for u
% -- 2. -- 
    [fm,fp] = full2face_p(f);

    ctr_f = ctr_flx(fm,fp); % 
    full_f = face2full(ctr_f);
    for ie = 1:Ne
        % weak form
        urhs(:,ie) = M\(full_f(:,ie) - D'*M*f(:,ie)); 
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
% Does not assume periodic, plus state zero now
    global Ne Nx
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    for ie=1:Ne
        if(ie == 1)
            vp(1,ie) = 0; % set to zero
            vp(2,ie) = v(1,ie+1);
        elseif(ie == Ne)
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = 0; % set to zero
        else
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,ie+1);
        end
    end
end
function [vm,vp] = full2face_p(v) % Assumes periodic boundary conditions
    global Ne Nx
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    if(Ne ~= 1) 
        for ie=1:Ne
            if(ie == 1)
                vp(1,ie) = v(Nx,Ne);
                vp(2,ie) = v(1,ie+1); % Periodic
            elseif(ie == Ne)
                vp(1,ie) = v(Nx,ie-1);
                vp(2,ie) = v(1,1);
            else
                vp(1,ie) = v(Nx,ie-1);
                vp(2,ie) = v(1,ie+1); % Periodic
            end
        end
    else
        vp(1,1) = v(Nx,1);
        vp(2,1) = v(1,1);
    end
end
%%
function f = flx_u(u)
    f = u.^2/2.;
end
