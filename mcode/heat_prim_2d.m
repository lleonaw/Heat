%%  2D heat equation, heat.pdf
%  Discontinuous Galerkin method
%  Explicit scheme, SSP-RK2
%  New primal form formulation, refer to heatwriteup/writeup/heat.pdf
% --- 
%  Mon Feb 29 16:34:30 CST 2016 
%  Seems to be working now... Can observe 
%   . 2nd order in time 
%   . spectral in space 
% --- 
function [succ,infer] = heat_prim_2d 
    global Ne Nx ifplt initflg T CFL dt
    global R2 u0 
    N = Nx - 1;       % Num of points in each elem.
    Nf = 4*Nx;        % Num of points on faces for one elem
    Nex = int8(sqrt(Ne)); 
    Ney = int8(Ne/Nex)  ; 
    infer = 2e20; 
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    ax = -1.; bx = 1.; dlx = (bx - ax)/double(Nex); % Length of each elem in x 
    ay = -1.; by = 1.; dly = (by - ay)/double(Ney); % Length of each elem in y 
    axe = zeros(Ne,1); aye = zeros(Ne,1); 
    xe = zeros(Nx,Ne);  ye = zeros(Nx,Ne); 
    x2  = zeros(Nx,Nx,Ne);  y2 = zeros(Nx,Nx,Ne); 
    for ie = 1:Ne
        axe(ie) = ax - (mod(ie,2)-1.)*dlx;  % 1,3 get ax; 2,4 get ax + dlx, 
        aye(ie) = ay + 0.5*(1.+ sign(ie-2.1))*dly;  % 1,2 get ay, 3,4 get ay+dly 
        xe (:,ie) = axe(ie) + dlx/2.*(z + 1.); 
        ye (:,ie) = aye(ie) + dly/2.*(z + 1.); 
        [x2(:,:,ie),y2(:,:,ie)] = meshgrid(xe(:,ie),ye(:,ie)); 
        x2(:,:,ie) = x2(:,:,ie)'; y2(:,:,ie) = y2(:,:,ie)'; 
    end
%   disp(xe) ; disp(ye) ; 
%   disp(x2) ; disp(y2) ; 
    x = reshape(x2,Nx*Nx,Ne); y = reshape(y2,Nx*Nx,Ne); 
%   disp(x) ; disp(y) ; %   error('ax ay for elem'); 
    Mx = (dlx/2.)*Mh; 
    My = (dly/2.)*Mh; 
    Dx = (2./dly)*Dh; % d kexi/ dx
    Dy = (2./dlx)*Dh; % d kexi/ dx
    Kx = (2./dlx)*Kh; % 1d stiffness matrix, Kx = Ky
    Ky = (2./dly)*Kh; % 1d stiffness matrix, Kx = Ky

% % Exact solution
%   if(initflg==1)     % Init Case 1  
%     u0 = 1.*cos(pi.*x/2.); qs = ((pi^2)/4.)*cos(pi.*x/2.);
%   elseif(initflg==2) % Init Case 2 
      u0 = -1.*sin(pi.*x); qs = (pi^2)*u0;
%     u0 =  1.*sin(pi.*y); qs =  (pi^2)*sin(pi.*x);
      uxx0 = pi^2.*sin(pi.*x); 
%     uyy0 = -pi^2.*sin(pi.*y); 
%   elseif(initflg==3) % Init Case 3 
%     u0 = 1.*cos(pi.*x); qs = 1.*(pi^2)*cos(pi.*x);
%   elseif(initflg==4) % Init Case 4 
%     u0 = 1.- 1.*cos(2.*pi.*x); qs = -((pi^2)*4.).*cos(2.*pi.*x);
%     u0 = 1.+ 1.*cos(pi.*x); qs = ((pi^2)).*cos(pi.*x);
%   end

% % Full 2 face - R matrix - R2, for ONE 2D element  
    In = eye(Nx); R2 = zeros(4*Nx, Nx*Nx); 
    R2 = [kron(In, In(:,1)');kron(In, In(:,end)');kron(In(:,1)', In);kron(In(:,end)', In)]; 
%   error('test R matrix'); % Oh yah it's working 

% % nhat matrix 
    nhx = zeros(Nf); nhy = zeros(Nf); 
    for ie=1:Ne
    for iface=1:4   
        i0 = 1 + (iface-1)*Nx  ; 
        ifn = i0 + Nx - 1 ; 
        if(iface == 1)
            nhx(i0:ifn,i0:ifn) = - eye(Nx);  % 
            nhy(i0:ifn,i0:ifn) = zeros(Nx);  % 
        elseif(iface == 2)
            nhx(i0:ifn,i0:ifn) =   eye(Nx);  % 
            nhy(i0:ifn,i0:ifn) = zeros(Nx);  % 
        elseif(iface == 3)
            nhx(i0:ifn,i0:ifn) = zeros(Nx);  % 
            nhy(i0:ifn,i0:ifn) = - eye(Nx);  % 
        elseif(iface == 4)
            nhx(i0:ifn,i0:ifn) = zeros(Nx);  % 
            nhy(i0:ifn,i0:ifn) =   eye(Nx);  % 
        end
    end
    end
    nhx = sparse(nhx); nhy = sparse(nhy); 
%   nhat = [nhx,nhy] ; 
%   disp(nhat);disp(size(nhat)); 
%   error('find the nhat matrix defined on face'); 
 
% % Area matrix 
    areay = zeros(Nf); areax = zeros(Nf); 
    for ie=1:Ne
    for iface=1:4
        i0 = 1 + (iface-1)*Nx  ; 
        ifn = i0 + Nx - 1 ; 
        areay(i0:ifn,i0:ifn) = (dly/2.).*diag(w(:));  % Use the 1d weights 
        areax(i0:ifn,i0:ifn) = (dlx/2.).*diag(w(:));  % Use the 1d weights 
    end
    end
    areay = sparse(areay); areax = sparse(areax); 

%   area = zeros(Nf, Ne); 
%   for ie=1:Ne
%   for iface=1:4
%       i0 = 1 + (iface-1)*Nx  ; 
%       ifn = i0 + Nx - 1 ; 
%       area(i0:ifn,ie) = w(:);  % Use the 1d weights 
%   end
%   end
%   disp(area); %   error('find the area matrix defined on face'); 
  
% % Q and Q^T 
    glb_indx = glb_setup(x2,y2);  %  x(nx1,nx1,ne)
    bdry_flg = bdry_setup(x2,y2);  % setup boundary flags 
%   disp('bdry');disp(bdry_flg); 
%   disp(glb_indx); 
%   error('global index?'); 
    uf   = R2*u0;  % full2face 
    uf = reshape(uf,[Nx,4,Ne]); 
    ug   = QTu(uf,glb_indx); % Sum at the same point 
    ulij = Qu (ug,glb_indx); % Copy to multi location 
%   u0, 
%   uf, 
%   ug,
%   ulij,
%   error('check globaltolocal, localtoglobal');  % ✓ , seems ok to me

    ka = ones(Nx,1);                         % Scalar field 
    nu = 1.*ones(Nx*Nx,1); nu = diag(nu);  % conductivity 

    % Time stepping setup 
%   CFL = 0.1; %   CFL = nu(1,1)*dt/((dx)^2); 
    t  = 0.; dxa = diff(x(:,1)); dx = min(dxa); 
    dt = CFL*(dx)^2 /nu(1,1);  % Alright 0.1 is as big as it gets 4 me 
    if(dt>=T)
        error('dt > T, change to a smaller CFL or longer T');
    end
    Nt = round(T/dt); dt= T/double(Nt); 
    disp(['Init case = ', num2str(initflg),...
          ' , Ne = ',num2str(Ne),' , N = ', num2str(N),... 
          ' , Nt = ',num2str(Nt),... 
          ' ,dt =',num2str(dt),', CFL = ',num2str(CFL)]);
    u = u0; tsv = zeros(Nt,1); infert = tsv; 

% % From poisson.m 
%   for i =1:Nx*Ne
%     g(i) = 1.; gmat = reshape(g,Nx,Ne);
%     [q,ku,gtu, gu, hu] = lh_pois(gmat,Mb,Db,Kb,ka,mth,qqt);   % Obtain (A u) in LHS 
%     A(:,i)   = reshape(q,Nx*Ne,1);
%     Ku(:,i)  = reshape(ku,Nx*Ne,1);  Hu(:,i) = reshape(hu,Nx*Ne,1);
%     Gtu(:,i) = reshape(gtu,Nx*Ne,1); Gu(:,i) = reshape(gu,Nx*Ne,1);
%     g(i) = 0.; 
%   end
% % % % % % % % % % % 
    M = kron(My,Mx);  % (ly lx)/4 .* Mhat otimes Mhat
    [urh,Ku,Gtu,Gu,Hu] = lh_pois(u,Mx,My,Dx,Dy,Kx,Ky,nu,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy);
% - u  := Nx*Nx, Ne 
% - Mx := Nx,Nx   ;  My := Nx,Nx
% - M  := Nx*Nx,Nx*Nx 
% - Kx := Nx,Nx   
%   disp('M'); disp(M); 
%   error('mass '); 
    uxx = zeros(size(uxx0)); 
%   uxx(2:end-1,:)=M(2:end-1,2:end-1)\urh(2:end-1,:); % no way to directly extract interior 
    uxx = M \ urh; % no way to directly extract interior 
                   % and it's the problem 
                   % interior values match the 1d case 
%   disp('uxx'); disp(uxx);
%   disp('R2*uxx'); disp(R2*uxx);
%   disp('R2T*R2*uxx'); disp(R2'*R2*uxx);
%   error('oh come on this is not helpful');
%   disp('reshape(bdry_flg,4*Nx,Ne).*(R2*uxx)'); 
%   disp(reshape(bdry_flg,4*Nx,Ne).*(R2*uxx));
%   disp('R2*uxx - reshape(bdry_flg,4*Nx,Ne).*(R2*uxx)');
%   disp(R2*uxx - reshape(bdry_flg,4*Nx,Ne).*(R2*uxx));
%   disp('R2T*(R2*uxx - reshape(bdry_flg,4*Nx,Ne).*(R2*uxx))');
%   disp(R2'*R2*uxx - R2'*(reshape(bdry_flg,4*Nx,Ne).*(R2*uxx)));
%   uxx = uxx - R2'*(reshape(bdry_flg,4*Nx,Ne).*(R2*uxx));  % Set boundary to be 0 

%   disp('uxx0'); disp(uxx0); 
%   disp('uxx'); disp(uxx); 
%   plsc = plt_fld(3,u0,x2,y2);    % Exact solution field 
%   plsc = plt_fld(1,uxx0,x2,y2);  % Exact 

%   plsc = plt_fld(2,M\Ku ,x2,y2,'Ku');  % 
%   plsc = plt_fld(3,M\Gtu,x2,y2,'Gtu');  % 
%   plsc = plt_fld(4,M\Gu ,x2,y2,'Gu');  % 
%   plsc = plt_fld(5,M\Hu ,x2,y2,'Hu');  % 
    plsc = plt_fld(7,uxx0 ,x2,y2,'uxx0');  % 
    plsc = plt_fld(6,uxx  ,x2,y2,'uxx');  % 
    disp('error between uxx and uxx0'); 
    infer = max(max(abs(uxx - uxx0))); 
    disp(infer); 
%   error('uxx and uxx0'); 
%   uxxpl0 = reshape(uxx0,[Nx,Nx,Ne]); 
%   upl0   = reshape(u0,[Nx,Nx,Ne]); 

%   figure(1); 
%   for ie = 1:Ne
%       mesh(x2(:,:,1),y2(:,:,1),upl0(:,:,1)); hold on; 
%   end
%   xlabel('--x--'); ylabel('--y--'); zlabel(' u0 ');

%   figure(2); 
%   for ie = 1:Ne
%       mesh(x2(:,:,1),y2(:,:,1),uxxpl0(:,:,1)); hold on; 
%   end
%   xlabel('--x--'); ylabel('--y--'); zlabel(' uxx ');
% % % % % % % % % % % % % % % For forward derivative 

%   for it = 1:Nt          
%  RK1 = Euler Explicit  
%     [urh] = lh_pois(u,M,Dx,Dy,Kb,nu,R2,glb_indx,area,nhx,nhy);
%     u     = u + dt.*(M \ urh);  % Second stage of rk 2
% % From prev 1D code 
%      RK2 
%         [urh] = lh_pois(u,Mb,Db,Kb,nu,qqt);   
%         u1    = u + dt.*(Mb \ urh);             % First stage of rk 2
%         [urh] = lh_pois(u1,Mb,Db,Kb,nu,qqt);   
%         u     = 0.5*(u + u1 + dt.*(Mb\ urh));  % Second stage of rk 2
%      RK4 
%         urh = lh_pois( u,Mb,Db,Kb,nu,qqt);  k1  = Mb \ urh; % trk = t; 
%         u1  = u + 0.5*dt*k1;                                % First stage 
%         urh = lh_pois(u1,Mb,Db,Kb,nu,qqt);  k2  = Mb \ urh; % trk = t + 0.5*dt; 
%         u2  = u + 0.5*dt*k2;                                % Second stage 
%         urh = lh_pois(u2,Mb,Db,Kb,nu,qqt);  k3  = Mb \ urh; % trk = t + 0.5*dt; 
%         u3  = u + 1.0*dt*k3;                                % Thrid stage 
%         urh = lh_pois(u3,Mb,Db,Kb,nu,qqt);  k4  = Mb \ urh; % trk = t + 0.5*dt; 
%         u   = u + dt*(k1/6. + k2/3. + k3/3. + k4/6.);       % Fourth stage
%     t = t + dt;
%     if(initflg==2) 
%       uex = exp(-nu(1).*(pi)^2.*t).*sin(pi.*x);
%     end 
%     if(ifplt && mod(it,ceil(Nt/5))==0)
%       plx = reshape(x,Nx*Nx*Ne,1); 
%       plu = reshape(u,Nx*Nx*Ne,1); 
%       pluex = reshape(uex,Nx*Nx*Ne,1); 
%       figure(initflg);hold on;
%       xlim([-1. 1.]);
%       plot(plx,plu,'rx-'); plot(plx,pluex,'-');
%       xlabel('-- x --'); ylabel('-- u --');
%%      legend('Numerical','Exact');
%       title(['solution at time t = ' num2str(t) ]); drawnow
%       pause(0.01); hold off;
%       disp(['max(u) = ' num2str(max(max(u))) ...
%          ' , min(u) = ' num2str(min(min(u))) ' , istep = ' num2str(it)]);
%     end
%     er = u - uex; 
%     tsv   (it) = t; 
%     mxtlu      = norm(uex,Inf); infert(it) = norm(er,Inf)/ mxtlu; 
%   end

%   if(ifplt)
%     figure(10);hold on;
%     plot(tsv,infert,'-');
%     xlabel('t '); 
%     ylabel('$\| u - \tilde{u}\|_{\infty} / \|\tilde{u}\|_{\infty}$','Interpreter','Latex');
%     title(['Error vs time']);
%     drawnow ; 
%     hold off;
%   end
%   infer = infert(end); % Error at end of time 
%   disp(['At end of time T = ',num2str(T),...
%   ', Relative Inf norm error = ', num2str(infer)]); 
    succ = true;
end

%------------------------------------------------------
function [succ] = plt_fld(idfig,u,x2,y2,str) 
    global Ne Nx 
% % Plot the solution 
    if(idfig~=0)
      figure(idfig); 
    else 
      figure(11); 
    end
    mesh(x2(:,:,1),y2(:,:,1),reshape(u(:,1),[Nx,Nx])); 
    hold on; 
    mesh(x2(:,:,2),y2(:,:,2),reshape(u(:,2),[Nx,Nx])); 
    mesh(x2(:,:,3),y2(:,:,3),reshape(u(:,3),[Nx,Nx])); 
    mesh(x2(:,:,4),y2(:,:,4),reshape(u(:,4),[Nx,Nx])); 
    title(str); 
    succ = true; 
end
% % % %  %  %  %  %  %  %  %  % 
%------------------------------------------------------
function [i,j] = ind_face2full(ij,iface) 
% From face number and elem number, 
%  return i,j indices 
   global Ne Nx 
   if(iface==1)
       i =  1; j = ij;
   elseif(iface==2) 
       i = Nx; j = ij;
   elseif(iface==3) 
       i = ij; j = 1;
   elseif(iface==4) 
       i = ij; j = Nx;
   end
end
function bd_flg = bdry_setup(x,y)  % setup boundary flags 
    global Ne Nx 
    Nf = 4*Nx; 
    bd_flg = zeros(Nx,4,Ne); 
    for ie=1:Ne
        for iface=1:4
        if(((ie==1)&&(iface==1 || iface==3)) ...
        || ((ie==2)&&(iface==2 || iface==3)) ... 
        || ((ie==3)&&(iface==1 || iface==4)) ... 
        || ((ie==4)&&(iface==2 || iface==4)) ) 
            bd_flg(:,iface,ie) = ones(Nx,1); 
        end
        end
    end
end
function glb_indx = glb_setup(x,y) % setup global index array 
    global Ne Nx ifplt initflg T CFL dt
    Nf = 4*Nx; 
    glb_indx = zeros(Nx,4,Ne); 
    xmid = 999.*ones(4,Ne); ymid = 999.*ones(4,Ne); 
    tol = 1.e-9; 
    k = 1; 
    for ie=1:Ne
        for iface=1:4
            [i,j] = ind_face2full(ceil(Nx/2),iface); 
            xmd = x(i,j,ie); ymd = y(i,j,ie); 
            xmid(iface,ie) = xmd; ymid(iface,ie) = ymd; 
            flg = 0; 
            ifsv = iface; iesv = ie; 
            if(ie==2 && iface == 2) % Periodic bc 
                flg = 1; ifsv = 1; iesv = 1;
            elseif(ie==4 && iface == 2)
                flg = 1; ifsv = 1; iesv = 3;
            elseif(ie==3 && iface == 4)
                flg = 1; ifsv = 3; iesv = 1;
            elseif(ie==4 && iface == 4)
                flg = 1; ifsv = 3; iesv = 2;
            end
            for ie2 = 1:ie-1
               for if2 = 1:4
                   if(abs(xmd - xmid(if2,ie2))<tol & ...
                       abs(ymd - ymid(if2,ie2))<tol ) 
                       flg = 1;  % Pre existing 
                       ifsv = if2; iesv = ie2; 
                   end
               end
            end
            for ij=1:Nx % Index on face 
                if(flg==0)      % New index, k++
                    glb_indx(ij,iface,ie) = k; 
                    k = k + 1; 
                elseif(flg==1)  % Already exist 
                    kprv = glb_indx(ij,ifsv,iesv); 
                    glb_indx(ij,iface,ie) = kprv; 
                end
            end
        end
    end
end
%- - -
function ug = QTu(ul,glb_indx) % SUM local to global 
    global Ne Nx ifplt initflg T CFL dt
    Nf = 4*Nx; 
    ug = zeros(Nf*Ne,1); % Well ug is smaller than Nf*Ne, does it matter? 
    for ie=1:Ne
        for iface=1:4
        for ij=1:Nx % Index on face 
            ig = glb_indx(ij,iface,ie);
            ug(ig) = ug(ig) + ul(ij,iface,ie);
        end
        end
    end
end
function ul = Qu(ug,glb_indx) % COPY global to local  
    global Ne Nx ifplt initflg T CFL dt
    Nf = 4*Nx; 
    ul = zeros(Nx,4,Ne); % Well my ug and ul have the same size 
    for ie=1:Ne
        for iface=1:4
        for ij=1:Nx % Index on face 
            ig = glb_indx(ij,iface,ie);
            ul(ij,iface,ie) = ug(ig);
        end
        end
    end
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
function [urhs,Ku,Gtu,Gu,Hu] = lh_pois(u,Mx,My,Dx,Dy,Kx,Ky,nu,R2,glb_indx,bdry_flg,...
                                        areax,areay,nhx,nhy)
% Pure central is very bad, modify q* in eval_fu - done
    global Ne Nx
    eta  = set_eta(Mx,My); he = 1.*eye(4*Nx); % For now eta, he are simply put
% GEOMETRY for 2d, need to work on here - Sat Mar  5 23:19:28 CST 2016
    I = speye(Nx);  
%   disp('u'); disp(u); 
    Ku =(kron(My,Kx) + kron(Ky,Mx))* u; 
    Ku = -nu*Ku;
%   disp('kron(My,kx)*u'); disp(kron(My,Kx)*u); 
%   disp('kron(I,kx)*u'); disp(kron(I,Kx)*u); 
%   disp('kron(Ky,Mx)*u'); disp(kron(Ky,Mx)*u);  % All zeros 
%   disp('Ku'); disp(Ku); 
%   disp('nu'); disp(nu); 
%   error('Check Ku');
    Gtu  = gt_eval(u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy); 
    Gtu = -nu*Gtu;
%   disp('Gtu'); disp(Gtu);  %  My bet is on Gtu and Hu 
    Gu   = g_eval (u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy); 
    Gu = -nu*Gu;
%   disp('Gu'); disp(Gu); 
%   error(' check Gtu and Gu '); 
    M = kron(My,Mx);  % (ly lx)/4 .* Mhat otimes Mhat
%   disp('M \ Gtu + Gu'); disp(M\(Gtu+Gu)); 
%   disp('M \ (Ku - Gtu - Gu)'); disp(M\(Ku - Gtu - Gu));
%   error('M \ G + Gt u, is it symmetric?'); 
    Hu = h_eval (eta,he,u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy); 
    Hu = -nu*Hu;
%   disp('Hu'); disp(Hu);   % Compare to 1D - nonzero values here, zeros there 
                            % Gonna deal with boundary terms in QQT 
%   error('all the terms in right hand side '); 
    urhs = Ku - Gtu - Gu + Hu; % urhs = -nu*urhs; 
%   urhs = Ku - Gtu - Gu ; % urhs = -nu*urhs; 
%   disp('urhs = Ku - Gtu - Gu + Hu'); disp(urhs); 
end
% -- For the new primal setup 
% What to do at boundary? 
% Periodic Plan: multiply values at the very end points by 2, and 
%                then they get divided by 2 in the next step 
function [Gtu] = gt_eval(u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy)
    global Ne Nx
%  Gtu := Dx'*R2'*nhx.*area.*(I - 0.5*QQT)*R2*u; 
    uf   = R2*u; 
%   disp('u'); disp(u); 
%   disp('uf'); disp(uf); 
    usv  = uf; 
    utmp = QTu(reshape(uf,[Nx,4,Ne]),glb_indx); 
    ufcr = Qu (utmp,glb_indx); 
%   ufcr = ufcr + bdry_flg.*ufcr;  % NEEDED for NOT periodic bc 
    ufcr = reshape(ufcr,4*Nx,Ne); 
%   disp('ufcr'); disp(ufcr); 
    ufm  = usv - 0.5*ufcr; 
%   disp('ufm'); disp(ufm); 
%   ufm = reshape(ufm,4*Nx*Ne,1); 
    Ie   = speye(Ne); In = speye(Nx); 
%   fx   = kron(Ie,nhx).*kron(Ie,area)*ufm; 
%   fy   = kron(Ie,nhy).*kron(Ie,area)*ufm; 
%   vfx  = R2'*reshape(fx,4*Nx,Ne); vfy = R2'*reshape(fy,4*Nx,Ne); 
%   disp('areay'); disp(areay); 
%   disp('nhx'); disp(nhx); 
%   disp('nhy'); disp(nhy); 
    fx   = (nhx.*areax)*ufm; 
    fy   = (nhy.*areay)*ufm; 
%   disp('fx'); disp(fx); 
%   disp('fy'); disp(fy); 
    vfx  = R2'*fx; vfy = R2'*fy; 
%   disp('vfx'); disp(vfx); 
%   disp('vfy'); disp(vfy); 
%   error('u on face'); 
%   disp('Dx'); disp(Dx); 
%   disp('Dy'); disp(Dy); 
    Gtu  = (kron(In,Dx)')*vfx + (kron(Dy,In)')*vfy;
%   disp('Gtu'); disp(Gtu); 
%   error(' Look at Gtu'); 
end 
function [Gu] = g_eval(u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy)
    global Ne Nx
    Ie = speye(Ne); In = speye(Nx); 
%   disp(' - - - - separator - - - '); 
%   disp('u'); disp(u); 
    Dux = (kron(In,Dx))*u; Duy = (kron(Dy,In))*u;
%   disp('Dux'); disp(Dux); 
%   disp('Duy'); disp(Duy); 

    duf = nhx.*areay*R2*Dux + nhy.*areax*R2*Duy; 
%   duf = nhx*R2*Dux + nhy*R2*Duy;  
  % comparing to 1d
  % not matching boudary values
  % b/c in 1d we did periodic, connect 1 to N
  % in 2d we are not. This QQT is pure Neumann connnectivity
    dusv = duf; 
%   disp('duf'); disp(duf); 
    duftmp = QTu(reshape(duf,[Nx,4,Ne]),glb_indx); 
%   disp('duftmp'); disp(duftmp); 
    ductr = Qu(duftmp,glb_indx);     % Size [Nx,4,Ne]
%   disp('ductr bf bc');disp(ductr); 
%   ductr = ductr + bdry_flg.*ductr;  % NEEDED for NOT periodic bc % Boundary values get doubled
%   disp('ductr after bc');disp(ductr); 
    ductr = reshape(ductr,4*Nx,Ne);
    Guf = dusv - 0.5*ductr;              % Boundary values get halved
%   disp('Guf'); disp(Guf); 
    Gu  = R2'*reshape(Guf,4*Nx,Ne); 
%   disp('Gu'); disp(Gu); 
%   error(' Look at Gtu'); 
end 
function [Hu] = h_eval(eta,he,u,Dx,Dy,R2,glb_indx,bdry_flg,areax,areay,nhx,nhy)
    global Ne Nx
    uf  = R2 * u; 
    usv = uf; 
    utmp = QTu(reshape(uf,[Nx,4,Ne]),glb_indx);
    uctr = Qu(utmp,glb_indx); 
%   uctr = uctr + bdry_flg.*uctr; % NEEDED for NOT periodic bc 
    uctr = reshape(uctr,4*Nx,Ne);
    Huf = 2.*usv - uctr;
%   disp('Huf'); disp(Huf); 
%   disp('eta'); disp(eta); 
%   disp('he'); disp(he); 
    const = zeros(size(eta)); % This gets rid of the NaN
    for i = 1:size(eta,1)
        const(i,i) = 1.; 
%       const(i,i) = 1.*eta(i,i)./he(i,i); 
    end
%   disp(const); 
%   disp(areax); 
%   Huf = (eta./he)*Huf;
    Huf = areax*const*Huf;
%   Huf = 0.*Huf ;  % Mask boundary to be 0 
%   disp('Huf 2 '); disp(Huf); 
    Hu  = R2'* Huf; 
%   disp('Hu, from face2full'); disp(Hu); 
%   error('Huuuuu'); 
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
function eta = set_eta(Mx,My) 
    global Ne Nx 
    eta = 1.*eye(4*Nx)/(Mx(1,1)*My(1,1));
%   [Mm Mp] = full2face_p(repmat(diag(M),1,Ne));
    % Periodic boundary? 
%   eta = ctr_flx_0(1./Mm,1./Mp); 
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
%%
