%%  Heat equation, no preferred direction, central flux
%  This code is for discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
function succ = heat
    clear all; format long;
    global Ne Nx
    Ne = 100;           % Number of elem
    N = 1;            % Poly. order
    % N = 4, 8 will blow up after round t = 0.27, def. t = 0.3
    % N 256, NP 8, blow up be4 0.3
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
    %u0 = 0.5 + sin(pi.*x);
    u0 = 2. - x.^2; 
    ka = ones(Nx,1);

    t = 0.; dt = 1e-5;
    T = 1.e-0; Nt = round(T/dt);
    ifplt = false;
    ifplt = true;
    disp(['Ne = ' num2str(Ne) ', N = ' num2str(N) ', dt=' num2str(dt) ]);

    u  = u0;  %
    if(ifplt)
        figure(1);hold on;
        for ie = 1:Ne
             plot(x(:,ie),u0(:,ie),'b-');
        end
        xlabel('-- x --'); ylabel('-- u --');
        title(['solution at time t = ' num2str(t) ]);
        drawnow
        pause(0.01);
        hold off;
    end

    for i = 1:Nt
        urh = lh_heat(u,Mb,Db,ka);           % Obtain right hand side
        u1  = u + dt.*urh;             % First stage of rk 2
        urh = lh_heat(u1,Mb,Db,ka);         % Obtain right hand side
        u   = 0.5*(u + u1 + dt.*urh);  % Second stage of rk 2

        t = t + dt;
        if(ifplt && mod(i,Nt/25)==0)
            clf; figure(1);hold on;
%            ylim([-0.1,1.1]);
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
function urhs = lh_heat(u,M,D,a) % Nodal DG Ch. 7.1
    sqa = sqrt(a); 
    vq = eval_q(u,M,D,sqa); % volumetric array
    urhs = eval_fu(vq,M,D,sqa);
end
function urhs = eval_fu(q,M,D,sa)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    f = diag(sa)*q;
    [fm,fp] = full2face(f);
    [fp] = bc_q(fp,fm); % Add bc to q
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
    [fm,fp] = full2face(f);
    [fp] = bc_fu(fp,fm); % add bc to flux for u
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
% Modify Boundary conditions on plus state:
%   Dirichlet: u+ = - u- , q+ = q-
%   Neumann  : u+ = u- , q+ = - q-
function [qp] = bc_q(qp,qm)
    global Ne 
    % Dirichlet
    qp(1,1) = qm(1,1);
%    qp(2,Ne) = qm(2,Ne); 
    % Neumann
%    qp(1,1) = - qm(1,1);
    qp(2,Ne) = - qm(2,Ne);
end
function [up] = bc_fu(up,um)
    global Ne
    % Dirichlet    
    up(1,1) = 2.*1. - um(1,1);
%    up(2,Ne) = 2.*1. - um(2,Ne);
    % Neumann
%    up(1,1) =  um(1,1);
    up(2,Ne) = um(2,Ne);
end
%%
function f = flx_u(u)
% Burgers
%    f = u.^2;
% heat 
    f = u;
end