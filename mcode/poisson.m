%%  1D Poisson equation, Nodal DG, pp. 261-275 
%  Discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
%  
%  Use Nek's procedure, follow hxdg routine in hmholtz.f
%  
function [succ,infer] = poisson(mth) 

    global Ne Nx ifplt 
    N = Nx - 1;       % Numb of points in each elem.
    infer = 2e20; 

    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    % Init Geom
%     a=0.; b=2.*pi; 
    % Alternate geometry 
    a = -1.; b = 1.; 
    dl = (b - a)/(Ne); % Length of each elem 
    x = zeros(Nx,Ne);
    for ie = 1: Ne
       ae = a + (ie-1)*dl;
       x(:,ie) = ae+0.5*dl*(z+1.);
    end
    Mb = (dl/2.)*Mh; % dx /d kexi  %  M bar, does not change
    Db = (2./dl)*Dh; % d kexi/ dx
    Ib = speye(Ne); 

    disp(['Method = ',num2str(mth), ...
          ', Ne = ',num2str(Ne),' , N = ', num2str(N)]); 
    ifsrc = true; 
    qs = -2.*ones(size(x));
%     qs =  1.*sin(x);      % Right hand side source term, analy. u = sin(x)  
    qs = 1.*cos(pi.*x/2.); 
    uex = qs;             % Exact solution
    ka = ones(Nx,1);      % Scalar field 
    ka = 4.*ka./(pi^2); 

% Just, amazing. wordless 
    g = zeros(Nx*Ne,1); A = zeros(Nx*Ne,Nx*Ne); 

    for i =1:Nx*Ne
      g(i) = 1.; 
      gmat = reshape(g,Nx,Ne);
      q = lh_pois(gmat,Mb,Db,ka,mth);           % Obtain right hand side
      A(:,i) = reshape(q,Nx*Ne,1);
      g(i) = 0.; 
    end
    dlmwrite('A.dat',A); 
    dlmwrite('M.dat',full(kron(Ib,Mb))); 
%   A = inv(kron(Ib,Mb)) * A; 
%   dlmwrite('A2.dat',A); 
%    dlmwrite(sprintf('A_Nx%d_ne%d.dat',Nx,Ne),A ); 
%    dlmwrite(sprintf('M_Nx%d_ne%d.dat',Nx,Ne),Mb); 
    % [Va,Da] = eig(A); eps = 1.e-13; 
    % Da = Da + eps.*sqrt(-1).*ones(size(Da)); figure(2);plot(diag(Da),'ro'); 
    % dlmwrite('eigA.dat',Da);
    % [md,mid] = min(Da);

    if(mth==1) 
      u = A \ reshape(qs,Nx*Ne,1); 
    elseif(mth==2)                   % Hack - essentially the tmask thing 
      qs1 = reshape(Mb*qs,Nx*Ne,1);
      u = zeros(Nx*Ne,1); 
      u(2:end-1) = A(2:end-1,2:end-1)\qs1(2:end-1); 
    end 
    plx   = reshape(x,Nx*Ne,1); 
    plu   = reshape(u,Nx*Ne,1); 
    pluex = reshape(uex,Nx*Ne,1); 
    if(ifplt)
        %close; 
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
%------------------------------------------------------
function flx = ctr_flx_0(fm,fp) % weak form
    fctr = (fm + fp)/2.;
%    difu = um - up;
%    fctr(1,:) = fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
    flx  = fctr;  % center flux runs
end
%
function flx = ctr_flx(fm,fp) % weak form
    fctr = (fm + fp)/2.;
%    difu = um - up;
    fctr(1,:) = - fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
    flx  = fctr;  % center flux runs
end
function flx = dif_flx(um,up) % weak form
    flx = um - up;
end
%------------------------------------------------------
function urhs = lh_pois(u,M,D,a,method) % Nodal DG Ch. 7.1
% Pure central is very bad, modify q* in eval_fu 
  if(method==1) 
    sqa = sqrt(a);            % Mult. this field twice 
    vq = eval_q(u,M,D,sqa);   % volumetric array
    urhs = eval_fu(vq,u,M,D,sqa);
    urhs = -1.*urhs;          % minus sign  
  elseif(method==2) 
    h1 = 1.; h2 = 0.;  % h1 equals a 
    h1 = a(1); % 4/ pi^2 
    eta  = set_eta(M); 
    urhs = hxdg_1(h1,h2,eta,M,D,u); 
  end 
end
% -- For the one setup 
function eta = set_eta(M) 
    global Ne Nx 
    [Mm Mp] = full2face_p(repmat(diag(M),1,Ne));
    % Periodic boundary? 
    eta = ctr_flx_0(1./Mm,1./Mp); 
end 
function [dum,dup] = nhat_mul(dum,dup)
    dum(1,:) = -1.*dum(1,:); 
    dup(2,:) = -1.*dup(2,:); 
end 
% --  
function au = hxdg_1(h1,h2,eta,M,D,u) 
    au = h2.*M*u;
    du = M*D*u;                % Need mass matrix? 
    [ufm,ufp] = full2face_p(u);
    [dum,dup] = full2face_p(du);
    [dum,dup] = nhat_mul(dum,dup);
    uf  = dif_flx(ufm,ufp); 
    duf = dif_flx(dum,dup); 
    du = h1.*du - 0.5.*face2full(uf); 
    au = au - 0.5.*face2full(duf) + face2full(eta.*uf); 
    au = au + D'*du;
    %uf, duf, au 
end 
% -- For the Two-layer setup 
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
    if(Ne>1)
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
    elseif(Ne==1)
        vp(1,1) = v(Nx,1);
        vp(2,1) = v(1,1);
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
