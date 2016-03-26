%%  1D heat equation, heat.pdf
%  Discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
%  New primal form formulation, refer to heatwriteup/writeup/heat.pdf
% --- 
%  Mon Feb 29 16:34:30 CST 2016 
%  Seems to be working now... Can observe 
%   . 2nd order in time 
%   . spectral in space 
% --- 
function [succ,infer,VA,DA,VL,DL,sid,plx] = heat_prim 
    global Ne Nx ifplt initflg T CFL dt iffwd ifeig
    N = Nx - 1;       % Numb of points in each elem.
    infer = 2e20; 
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    a = -1.; b = 1.; dl = (b - a)/(Ne); % Length of each elem 
    x = zeros(Nx,Ne);
    for ie = 1:Ne
       ae = a + (ie-1)*dl; x(:,ie) = ae+0.5*dl*(z+1.);
    end
    plx   = reshape(x,Nx*Ne,1); 
    Mb = (dl/2.)*Mh;  % dx /d kexi  %  M bar, does not change
    Db = (2./dl)*Dh;  % d kexi/ dx
    Kb = (2./dl)*Kh;  % stiffness matrix, 
% % Exact solution
    if(initflg==1)     % Init Case 1  
      u0 = 1.*cos(pi.*x/2.); qs = ((pi^2)/4.)*cos(pi.*x/2.);
    elseif(initflg==2) % Init Case 2 
%     u0 = -1.*sin(pi.*x); qs = (pi^2)*u0;
      u0 =  1.*sin(pi.*x); qs =  (pi^2)*sin(pi.*x);
    elseif(initflg==3) % Init Case 3 
      u0 = 1.*cos(pi.*x); qs = 1.*(pi^2)*cos(pi.*x);
    elseif(initflg==4) % Init Case 4 
%     u0 = 1.- 1.*cos(2.*pi.*x); qs = -((pi^2)*4.).*cos(2.*pi.*x);
      u0 = 1.+ 1.*cos(pi.*x); qs = ((pi^2)).*cos(pi.*x);
    elseif(initflg==5) % Init Case 5, discontinuous 
      u0 = rand(size(x));  % 
    end
%
    nu = 1.*ones(Nx,1);  nu = diag(nu); 
% 
    Ie = speye(Ne); Ih = speye(N+1); 
    fluxf = 0.*speye(2); fluxf(1,1)=.5; fluxf(end,end)=.5;
    qqt = kron(Ie,fluxf); j=0; m=2; M=size(qqt,1);    % Defined for surface, each block is size 2 
    for e=1:Ne; j=j+m; jp = j+1; if jp>M; jp=1; end;
    qqt(j,jp)=.5; qqt(jp,j)=.5; qqt(jp,jp)=.5; end;   % Use the signless one 
    if(ifeig) 
%    For finding A 
      A   = zeros(Nx*Ne,Nx*Ne);  g = zeros(Nx*Ne,1); 
      Ku  = zeros(Nx*Ne,Nx*Ne); Hu = zeros(Nx*Ne,Nx*Ne); 
      Gtu = zeros(Nx*Ne,Nx*Ne); Gu = zeros(Nx*Ne,Nx*Ne);   
      for i =1:Nx*Ne
        g(i) = 1.; gmat = reshape(g,Nx,Ne);
        [q,ku,gtu,gu,hu] = lh_pois(gmat,Mb,Db,Kb,nu,qqt);   % Obtain (A u) in LHS 
        A(:,i)   = reshape(q,Nx*Ne,1);
        Ku(:,i)  = reshape(ku,Nx*Ne,1);  Hu(:,i) = reshape(hu,Nx*Ne,1);
        Gtu(:,i) = reshape(gtu,Nx*Ne,1); Gu(:,i) = reshape(gu,Nx*Ne,1);
        g(i) = 0.; 
      end
% 
      Iq = eye(Nx*Ne-1); R=[Iq(1,:);Iq];      % Periodic bc - Ill-conditioned 
%      No boundary mask for DG 
      [VA,DA] = eig(A);                       % For stiffness matrix 
      DA = sort(diag(DA)); 
      [VL,DL] = eig( A, full(kron(Ie,Mb)));   % For derivative operator
      [DL,sid] = sort(diag(DL)); 
%     idl = [N*Ne:length(DL)]; sid;
%     figure(1); 
%     myleng = cell(size(idl));  % For plotting eig 
%     for invd=1:length(idl)
%         plot(plx,VL(:,sid(idl(invd)))); % if(invd==1) hold on; end 
%         idl(invd), DL(idl(invd)), 
%         myleng{invd}=strcat('k=',num2str(idl(invd)));
%         pause
%     end
%     legend(myleng);xlabel('x'); title('eigenvectors');
%     pause 
%
      if(iffwd)
        uxx0 = pi^2*u0; 
        uxx  = inv(kron(Ie,Mb))*A*(reshape(u0,[],1));
%
        plu = reshape(uxx,Nx*Ne,1); pluex = reshape(uxx0,Nx*Ne,1); 
        infer = norm(plu-pluex,Inf); 
        disp(['DGM :: Error in u_{xx} = ', num2str(infer)...
            ,', N = ',num2str(N),' , N_e = ',num2str(Ne)]); 
      end
    else
      VA = 0.; DA = 0.;
      VL = 0.; DL = 0.;
    end
% 
% % % % % % % % % % % % % % % % % % % % % 
% % % Time stepping setup 
    t = 0.; dxa = diff(x(:,1)); dx = min(dxa); 
    dt = CFL*(dx)^2 /nu(1,1);  % Alright 0.1 is as big as it gets 4 me 
    if(dt>=T)
      error('dt > T, change to a smaller CFL or longer T');
    end
    Nt = round(T/dt); dt= T/double(Nt); 
    disp(['Init case = ', num2str(initflg),...
          ', Ne = ',num2str(Ne),' , N = ', num2str(N),... 
          ', Nt = ',num2str(Nt),... 
          ' ,dt =',num2str(dt),', CFL = ',num2str(CFL)]);
    u   = u0; tsv = zeros(Nt,1); infert = tsv; 

    for it = 1:Nt          
%  Euler 1
%     urh = lh_pois( u,Mb,Db,Kb,nu,qqt);  
      [urh,ku,gtu,gu,hu] = lh_pois(u,Mb,Db,Kb,nu,qqt);   % Obtain (A u) in LHS 
      urh = -1. * urh; 
      k1  = Mb \ urh; % trk = t; 
      u   = u + dt*k1;                                % First stage of rk 4
%     error('k1 in 1d'); 
%  RK2 
%     [urh] = lh_pois(u,Mb,Db,Kb,nu,qqt);   
%     u1    = u + dt.*(Mb \ urh);             % First stage of rk 2
%     [urh] = lh_pois(u1,Mb,Db,Kb,nu,qqt);   
%     u     = 0.5*(u + u1 + dt.*(Mb\ urh));  % Second stage of rk 2
%  RK4 
%     urh = lh_pois( u,Mb,Db,Kb,nu,qqt);  k1  = Mb \ urh; % trk = t; 
%     u1  = u + 0.5*dt*k1;                                % First stage of rk 4
%     urh = lh_pois(u1,Mb,Db,Kb,nu,qqt);  k2  = Mb \ urh; % trk = t + 0.5*dt; 
%     u2  = u + 0.5*dt*k2;                                % Second stage of rk 4
%     urh = lh_pois(u2,Mb,Db,Kb,nu,qqt);  k3  = Mb \ urh; % trk = t + 0.5*dt; 
%     u3  = u + 1.0*dt*k3;                                % Thrid stage of rk 4
%     urh = lh_pois(u3,Mb,Db,Kb,nu,qqt);  k4  = Mb \ urh; % trk = t + 0.5*dt; 
%     u   = u + dt*(k1/6. + k2/3. + k3/3. + k4/6.);       % Fourth stage of rk 4

      t = t + dt;
      if(initflg==2) 
          uex = exp(-nu(1).*(pi)^2.*t).*sin(pi.*plx);
          mxtlu   = norm(uex,Inf); 
      elseif(initflg==5)  % This is not the exact! !  
          uex = exp(-nu(1).*(pi)^2.*t).*u0;
          uex = reshape(uex,size(plx));
          uex = 0.*uex; 
          mxtlu   = 1.;
      end 
      if(ifplt && mod(it,ceil(Nt/10))==0)
        figure(initflg);hold on;
        xlim([-1. 1.]);
        plot(plx,reshape(u,Nx*Ne,1),'x-'); plot(plx,uex,'o-');
        xlabel('-- x --'); ylabel('-- u --');
        legend('Numerical','Exact');
        title(['solution at time t = ' num2str(t) ]); drawnow
        pause(0.01); hold off;
        disp(['max(u) = ' num2str(max(max(u))) ...
           ' , min(u) = ' num2str(min(min(u))) ' , istep = ' num2str(it)]);
      end
      er = reshape(u,Nx*Ne,1) - uex; 
      tsv(it) = t; 
      infert(it) = norm(er,Inf)/ mxtlu; 
    end
    if(ifplt)
      figure(10);hold on;
      plot(tsv,infert,'-');
      xlabel('t '); 
      ylabel('$\| u - \tilde{u}\|_{\infty} / \|\tilde{u}\|_{\infty}$','Interpreter','Latex');
      title(['Error vs time']);
      drawnow ; 
      hold off;
    end
    infer = infert(end); % Error at end of time 
    disp(['DGM :: At end of time T = ',num2str(T),...
    ', Relative Inf norm error = ', num2str(infer)]); 
    succ = true;
end
%------------------------------------------------------
function flx = ctr_flx_0(fm,fp) % weak form
    fctr = (fm + fp)/2.;
    flx  = fctr;  % center flux runs
end
function flx = dif_flx(um,up) % weak form
    flx = um - up;
end
%------------------------------------------------------
function [urhs,Ku,Gtu,Gu,Hu] = lh_pois(u,M,D,K,nu,qqt) 
% Pure central is very bad, modify q* in eval_fu - done
  global Ne Nx
  eta  = set_eta(M); he = ones(size(eta)); % For now eta, he are simply put
  Ku = K*u; 
  Gtu  = gt_eval(D,u,qqt);     % I_1; How to find G', with G? 
  Gu   = g_eval (D,u,qqt);     % So actually this is G, namely the one in I_2 
  Hu   = h_eval (eta,he,u,qqt) ; % I_3 term 
  urhs = Ku - Gtu - Gu + Hu;   % - \nabla^2 u 
  urhs = nu*urhs; 
end
% -- For the new primal setup 
% What to do at boundary? 
% Periodic Plan: multiply values at the very end points by 2, and 
%                then they get divided by 2 in the next step 
function [Gtu] = gt_eval(D,u,qqt)
    global Ne Nx
    [um,up]    = full2face_p(u);    % R u 
    usvm       = um;                % save u- 
    % Approach 1 - this works, gonna keep using it 
    ufcr = qqt*reshape(um,2*Ne,1);

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
    dufcr      = qqt*reshape(dum,2*Ne,1); 
    Guf        = dusvm - reshape(dufcr,2,Ne);       % Du- - (Du- + Du+)/2 
    Gu         = face2full(Guf);      % R^T 
end 
function [Hu] = h_eval(eta,he,u,qqt)
    global Ne Nx
    [um,up] = full2face_p(u);     % R u 
    usvm    = um;                 % save u- 
    ufcr    = qqt*reshape(um,2*Ne,1);
    Huf     = 2.*(usvm - reshape(ufcr,2,Ne));   % 2u- - (u- + u+) 
    Huf     = eta.*Huf./he; 
    Hu      = face2full(Huf);     % R^T 
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
% % % % % % % % % % % % % % % % % % % % % %
