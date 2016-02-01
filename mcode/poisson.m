%%  1D Poisson equation, Nodal DG, pp. 261-275 
%  Discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
%  

function [succ,infer] = poisson

    global Ne Nx ifplt
    N = Nx - 1;       % Numb of points in each elem.

    infer = 2e20; 

    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    % Init Geom
%   a=-1.; b=1.;
    a=0.; b=2.*pi; %
    dl = (b - a)/(Ne); % Length of each elem 
    x = zeros(Nx,Ne);
    for ie = 1: Ne
        ae = a + (ie-1)*dl;
        x(:,ie) = ae+0.5*dl*(z+1.);
    end
    Mb = (dl/2.)*Mh; % dx /d kexi  %  M bar, does not change
    Db = (2./dl)*Dh; % d kexi/ dx

    ifsrc = true; 
    qs = -2.*ones(size(x));
    qs = -1.*sin(x); % Right hand side source term, analy. u = sin(x)  
    uex = qs;  % Exact solution

    ka = ones(Nx,1); % Scalar field 

%   ifplt = false;
%   ifplt = true;

% Just, amazing. wordless 
    g = zeros(Nx*Ne,1); A = zeros(Nx*Ne,Nx*Ne); 

    for i =1:Nx*Ne
      g(i) = 1.; 
      gmat = reshape(g,Nx,Ne);
      q = lh_pois(gmat,Mb,Db,ka);           % Obtain right hand side
      A(:,i) = reshape(q,Nx*Ne,1);           
      g(i) = 0.; 
    end
    [Va,Da] = eig(A); eps = 1.e-13; 
    % Da = Da + eps.*sqrt(-1).*ones(size(Da)); % figure(2);plot(diag(Da),'ro'); 
    % dlmwrite('eigA.dat',Da);
    % [md,mid] = min(Da);

    u = A \ reshape(qs,Nx*Ne,1); 
    plx   = reshape(x,Nx*Ne,1); 
    plu   = reshape(u,Nx*Ne,1); 
    pluex = reshape(uex,Nx*Ne,1); 
    if(ifplt)
        figure(1);hold on;
        plot(plx,pluex,'bx-',plx,plu  ,'r-');
        xlabel('-- x --'); ylabel('-- u --');
        title(['solution']);
        legend('Exact','Num.'); drawnow ; pause(0.2);
        hold off;
    end
    infer = norm(plu-pluex,Inf); 
    disp(['Inf norm error= ', num2str(infer)]); 

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
function flx = dif_flx(um,up) % weak form
    difu = um - up;
%    fctr(1,:) = - fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
%    flx  = fctr;  % center flux runs
    flx = difu; 
end
function urhs = lh_pois(u,M,D,a) % Nodal DG Ch. 7.1
% Pure central is very bad, modify q* in eval_fu 
    sqa = sqrt(a);            % Mult. this field twice 
    vq = eval_q(u,M,D,sqa);   % volumetric array
    urhs = eval_fu(vq,u,M,D,sqa);
    % Actually I think there is a minus sign 
    urhs = -1.*urhs; 
end
function urhs = eval_fu(q,u,M,D,sa)
    global Ne Nx
    tau = 1.; 
    urhs = zeros(Nx,Ne);
    f = diag(sa)*q;
% So this is stabilized central flux 
    [fm,fp] = full2face(f);
    [fp] = bc_q(fp,fm); % Add bc to q
    ctr_f = ctr_flx(fm,fp); % Center for both 

    [um,up] = full2face(u);
    [up] = bc_fu(up,um); % Add bc to q
    dif_u = dif_flx(um,up); % Difference in u

    full_f = face2full(ctr_f - tau.*dif_u);
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
    ctr_f = ctr_flx(fm,fp); %  Center for both 
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
    % Dirichlet for u 
    qp(1,1) = qm(1,1);
    qp(2,Ne) = qm(2,Ne); 
    % Neumann for u 
%    qp(1,1) = - qm(1,1);
%    qp(2,Ne) = - qm(2,Ne);
end
function [up] = bc_fu(up,um)
    global Ne
    % Dirichlet    
    up(1,1) = 2.*0. - um(1,1);
    up(2,Ne) = 2.*0. - um(2,Ne);
    % Neumann
%    up(1,1) =  um(1,1);
%    up(2,Ne) = um(2,Ne);
end
%%
