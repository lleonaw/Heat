%%  1D Poisson equation, Nodal DG, pp. 261-275 
%  Discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
%  
%  Use Nek's procedure, follow hxdg routine in hmholtz.f
%  
%  Sun Feb 28 12:54:50 CST 2016
%  New primal form formulation, refer to heatwriteup/writeup/heat.pdf
%  
function [succ,infer] = poisson(mth) 
    global Ne Nx ifplt initflg
    N = Nx - 1;       % Numb of points in each elem.
    infer = 2e20; 
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
% Init Geom & Alternate geometry 
%   a=0.; b=2.*pi; %   a = 0.; b = 1.; 
    a = -1.; b = 1.; dl = (b - a)/(Ne); % Length of each elem 
    x = zeros(Nx,Ne);
    for ie = 1:Ne
       ae = a + (ie-1)*dl; x(:,ie) = ae+0.5*dl*(z+1.);
    end
    plx   = reshape(x,Nx*Ne,1); 
    Mb = (dl/2.)*Mh;  % dx /d kexi  %  M bar, does not change
    Db = (2./dl)*Dh;  % d kexi/ dx
    Kb = (2./dl)*Kh;  % stiffness matrix, 
    disp(['Method = ',num2str(mth), ' , Init case = ', num2str(initflg),...
          ', Ne = ',num2str(Ne),' , N = ', num2str(N)]); 
% % Exact solution
    if(initflg==1)
% Init Case 1  
      uex = 1.*cos(pi.*x/2.); qs = ((pi^2)/4.)*cos(pi.*x/2.);
    elseif(initflg==2) 
% Init Case 2 
%     uex = -1.*sin(pi.*x); qs = (pi^2)*uex;
      uex =  1.*sin(pi.*x); qs =  (pi^2)*sin(pi.*x);
    elseif(initflg==3) 
% Init Case 3 
      uex = 1.*cos(pi.*x); qs = 1.*(pi^2)*cos(pi.*x);
    elseif(initflg==4) 
% Init Case 4 
%     uex = 1.- 1.*cos(2.*pi.*x); qs = -((pi^2)*4.).*cos(2.*pi.*x);
      uex = 1.+ 1.*cos(pi.*x); qs = ((pi^2)).*cos(pi.*x);
    end
%   plot(plx,reshape(uex,Nx*Ne,1),'ro-');hold on;
%   plot(plx,reshape(qs,Nx*Ne,1),'bx-');
%   error('plotitit'); 

% Test uxx, Kb is singular 
%   uxx = Kb \ (Mb * uex);  uxx = reshape(uxx,Nx*Ne,1); 
%   disp(cond(Kb)); 
%   figure(1);
%   plot(plx,uxx,'ro-');hold on;
%   plot(plx,reshape(qs,Nx*Ne,1),'bx-');
%   disp(uxx); 
%   hold off;
%   error('nah '); 
% 
    ka = ones(Nx,1);      % Scalar field 
% Copied from inv_vis, for central flux 
    IE = speye(Ne); Ih = speye(N+1); 
    fluxf = 0.*speye(2); fluxf(1,1)=.5; fluxf(end,end)=.5;
    qqt = kron(IE,fluxf); j=0; m=2; M=size(qqt,1);    % Defined for surface, each block is size 2 
    for e=1:Ne; j=j+m; jp = j+1; if jp>M; jp=1; end;
%   qqt(j,jp)=.5; qqt(jp,j)=-.5; qqt(jp,jp)=-.5; end; % Seems to have signs(nhat) with it
    qqt(j,jp)=.5; qqt(jp,j)=.5; qqt(jp,jp)=.5; end;   % Use the signless one 

% Just, amazing. wordless 
    A = zeros(Nx*Ne,Nx*Ne);  g = zeros(Nx*Ne,1); 
    Ku = zeros(Nx*Ne,Nx*Ne); Hu = zeros(Nx*Ne,Nx*Ne); 
    Gtu = zeros(Nx*Ne,Nx*Ne);Gu = zeros(Nx*Ne,Nx*Ne);   
    for i =1:Nx*Ne
      g(i) = 1.; gmat = reshape(g,Nx,Ne);
      [q,ku,gtu, gu, hu] = lh_pois(gmat,Mb,Db,Kb,ka,mth,qqt);   % Obtain (A u) in LHS 
      A(:,i)   = reshape(q,Nx*Ne,1);
      Ku(:,i)  = reshape(ku,Nx*Ne,1);  Hu(:,i) = reshape(hu,Nx*Ne,1);
      Gtu(:,i) = reshape(gtu,Nx*Ne,1); Gu(:,i) = reshape(gu,Nx*Ne,1);
      g(i) = 0.; 
    end
    dlmwrite('A.dat',A);  dlmwrite('qqt.dat',full(qqt)); 
    dlmwrite('K.dat',Ku); dlmwrite('H.dat',Hu); 
    dlmwrite('G.dat',Gu); dlmwrite('Gtu.dat',Gtu); 

    if(mth==1)     % this one does not work now ... 
      u = A \ reshape(qs,Nx*Ne,1); 
    elseif(mth==2)                   % Hack - essentially the tmask thing 
      qs1 = reshape(Mb*qs,Nx*Ne,1);
      u = zeros(Nx*Ne,1); 
      u(2:end-1) = A(2:end-1,2:end-1)\qs1(2:end-1); 
    elseif(mth==3)                   % Primal form 
      rhs = Mb*qs; 
      rhs = reshape(rhs,Nx*Ne,1);
      u   = zeros(Nx*Ne,1); 
%     u   = A \ rhs; 
      disp(A); disp(rhs); 
      u(1:end-1) = A(1:end-1,1:end-1) \ rhs(1:end-1); u(end) = u(1); 
      if(initflg==3)
%         u = u - 1.;
      end
    end 
    plu   = reshape(u,Nx*Ne,1); 
    pluex = reshape(uex,Nx*Ne,1); 
    if(ifplt)
      %close; 
      figure(initflg);hold on;
      plot(plx,pluex,'bx-',plx,plu  ,'r-');
      xlabel('-- x --'); ylabel('-- u --');
      title(['solution']);
      legend('Exact','Num.'); drawnow ; pause(0.2);
      hold off;
    end
%   disp(plx); 
%   disp(plu);
%   disp(pluex); 
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
function [urhs,Ku,Gtu, Gu, Hu] = lh_pois(u,M,D,K,a,method,qqt) 
% Pure central is very bad, modify q* in eval_fu - done
  global Ne Nx
  Ku  = zeros(Nx,Ne); 
  Gtu = zeros(Nx,Ne);       % I_1; How to find G', with G? 
  Gu  = zeros(Nx,Ne);       % So actually this is G, namely the one in I_2 
  Hu  = zeros(Nx,Ne); % I_3 term 
  if(method==1) % Nodal DG Ch. 7.1
    sqa = sqrt(a);            % Mult. this field twice 
    vq = eval_q(u,M,D,sqa);   % volumetric array
    urhs = eval_fu(vq,u,M,D,sqa);
    urhs = -1.*urhs;          % minus sign  
  elseif(method==2) 
    h1 = 1.; h2 = 0.;  % h1 equals 1/ ka 
    h1 = 1./a(1);      % 4/ pi^2 
    eta  = set_eta(M); 
    urhs = hxdg_1(h1,h2,eta,M,D,u); 
  elseif(method==3)  % A u = (K - G - G' + H) u 
% Well really this one breaks easily as well... 
% Currently Hu makes it blow up for some cases, one being: Nx = 4+1, Ne = 2
    eta  = set_eta(M); he = ones(size(eta)); % For now eta, he are simply put
    Ku = K*u; 
    Gtu  = gt_eval(D,u,qqt);       % I_1; How to find G', with G? 
    Gu   = g_eval (D,u,qqt);       % So actually this is G, namely the one in I_2 
    Hu   = h_eval (eta,he,u,qqt) ; % I_3 term 
    urhs = Ku - Gtu - Gu + Hu; 
%   urhs = Ku - Gtu - Gu; 
%   disp('Ku  '); disp(Ku); 
%   disp('Gtu '); disp(Gtu); 
%   disp('Gu  '); disp(Gu); 
%   disp('Hu  '); disp(Hu); 
%   disp('urhs'); disp(urhs); 
%   error('why not symmetric'); 
  end 
end
% -- For the new primal setup 
% What to do at boundary? 
% Plan: multiply values at the very end points by 2, and 
%       then they get divided by 2 in the next step 
function [Gtu] = gt_eval(D,u,qqt)
    global Ne Nx
    [um,up]    = full2face_p(u);    % R u 
%   [um,up]    = full2face  (u);    % R u 
%   disp(u); %   disp(um); %   disp(up); 
    usvm       = um;                % save u- 
    % Approach 1 
    ufcr = qqt*reshape(um,2*Ne,1);
%   error('arrived')
    % Approach 2 
%   up         = bc_p(um,up);
%   ufcr       = ctr_flx_0(um,up);  % (u- + u+)/2 

    ufm        = usvm - reshape(ufcr,2,Ne);       % u- - (u- + u+)/2 
    [ufm,up]   = nhat_mul(ufm,up);  % nhat (I-0.5QQ') R u 
    Gtu        = face2full(ufm);    % u- - (u- + u+)/2 
    Gtu        = D'*Gtu; 
end 
function [Gu] = g_eval(D,u,qqt)
    global Ne Nx
    du         = D*u; 
    [dum,dup]  = full2face_p(du);     % R D u 
    [dum,dup]  = nhat_mul(dum,dup);   % nhat dot RDu 
    dusvm      = dum;                 % save RDu- 
%   dup        = bc_p(dum,dup);
%   dufcr      = ctr_flx_0(dum,dup);  % (Du- + Du+)/2 
    dufcr      = qqt*reshape(dum,2*Ne,1); 
    Guf        = dusvm - reshape(dufcr,2,Ne);       % Du- - (Du- + Du+)/2 
    Gu         = face2full(Guf);      % R^T 
end 
function [Hu] = h_eval(eta,he,u,qqt)
    global Ne Nx
    [um,up] = full2face_p(u);     % R u 
    usvm    = um;                 % save u- 
%   [up]    = bc_p(um,up);
%   [ufcr]  = ctr_flx_0(um,up);   % (u- + u+)/2 
    ufcr    = qqt*reshape(um,2*Ne,1);
    Huf     = 2.*(usvm - reshape(ufcr,2,Ne));   % 2u- - (u- + u+) 
    Huf     = eta.*Huf./he; 
    Hu      = face2full(Huf);     % R^T 
end 
function [qp] = bc_p(qp,qm)    
% Periodic 
    global Ne 
%   qp(1,1)  =  2.*qm(1,1);
%   qp(2,Ne) =  2.*qm(2,Ne); 
    qp(1,1)  =  qm(2,Ne);
    qp(2,Ne) =  qm(1,1); 
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
    dup(2,:) = -1.*dup(2,:); % In 1d this makes sense 
%     dup(1,:) = -1.*dup(1,:); 
end 
% --  
function au = hxdg_1(h1,h2,eta,M,D,u) 
    global Ne Nx
    au = h2.*M*u;
    du = D*u,                % Need mass matrix? 
    [ufm,ufp] = full2face_p(u)
    [dum,dup] = full2face_p(du)
    [dum,dup] = nhat_mul(dum,dup),
    nhat = ones(size(ufm));
    nhat(1,:) = -1.*ones(1,Ne); 
    uf  = dif_flx(ufm,ufp) 
    duf = dif_flx(dum,dup) 
    du = h1.*M*du
    du = du - 0.5.*face2full(nhat.*uf) 
    au = au   + face2full(eta.*uf) 
   % - 0.5.*face2full(duf)  ...
    au = au + D'*du
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
    urhs    = zeros(Nx,Ne);
    f       = diag(sa)*u;
    [fm,fp] = full2face(f);
    [fp]    = bc_fu(fp,fm); % add bc to flux for u
    ctr_f   = ctr_flx(fm,fp); %  Center for both 
    full_f  = face2full(ctr_f);
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
            vp(2,ie) = v(1,ie+1);  % Periodic
        elseif(ie == Ne)
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,1);
        else
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,ie+1);  % Periodic
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
