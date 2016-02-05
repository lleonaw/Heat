%-----------------------------------------------------
function [te,ee]=tester(Nelx,N,No,dt)
%-----------------------------------------------------
%!resize
%% Typical usage: t05x20=in5(5,20,20,1.e-3);
%% Nelx = number of elements in x direction
%% N    = polynomial order
%% No   = order for dealiasing on (No+1) GLL points

close all; format compact; format longe

[Ah,Bh,Ch,Dh,z,w]=SEMhat(N);[Ao,Bo,Co,Do,zo,wo]=SEMhat(No);

IE=speye(Nelx);

Jo=interp_mat(zo,z); J=kron(IE,Jo);

x0 = -5; x1 =  5; Lxe = (x1-x0)/Nelx;

%% Set up standard grid:
n=N+1; no=No+1; 
clear Xe; 
for e=1:Nelx; Xe (1:n ,e) = x0 + (e-1)*Lxe + Lxe*(z+1)/2; end;
x  = reshape(Xe ,Nelx*n ,1); y = x;
[X ,Y ]=ndgrid(x ,y ); r = sqrt(X.*X + Y.*Y);

Tb    = 1.0;
beta  = 5.000; % Vortex strength
gamma = 1.4;
cT    = (gamma-1)*beta^2/(16*gamma*pi^2);

%% SET INITIAL CONDITIONS FOR EXACT SOLUTION: 

Cx=0.0; Cy=0.0;   %% Translation velocity
Cx=1.0; Cy=0.0;   %% Translation velocity
T    =  Tb-cT*exp(2*(1-r.*r));
u    = Cx-(beta*Y/(2*pi)).*exp(1-r.*r);
v    = Cy+(beta*X/(2*pi)).*exp(1-r.*r);
rho  =  T.^(1./(gamma-1));
exact=y_pack(rho,u,v,T);

p    =  rho.*T;
rhou =  rho.*u;
rhov =  rho.*v;
rhoe =  (p / (gamma-1.) ) + rho.*(u.*u + v.*v)/2.;



% Bh = Jo'*Bo*Jo;  % Uncomment for full local mass matrix (doesn't work well...)

B  = (Lxe/2)*kron(IE,Bh); D  = (2/Lxe)*kron(IE,Dh);  %% 1D Mass & Derivative
Bf = (Lxe/2)*kron(IE,Bo);  %% 1D Mass, fine mesh
Bi = inv(B);


%
%    Set up 1D flux exchange matrix for central fluxes, periodic domain:
%
%
%           /-.5                        -.5 \
%           |    0                          |
%           |      0                        |    Note: This operates on the fine grid.
%           |        .5  .5                 |
%           |       -.5 -.5                 |
%           |               0               |
%  FluxM =  |                 0             |
%           |                  .5  .5       |
%           |                 -.5 -.5       |
%           |                         0     |
%           |                           0   |
%           \ .5                         .5 /
%

Ih = speye(N+1); Io = speye(No+1); Fluxo = 0*Io; Fluxo(1,1)=.5; Fluxo(end,end)=.5;
FluxM = kron(IE,Fluxo); j=0; m=No+1; M=size(FluxM,1);
for e=1:Nelx; j=j+m; jp = j+1; if jp>M; jp=1; end;
    FluxM(j,jp)=.5; FluxM(jp,j)=-.5; FluxM(jp,jp)=-.5; end;

Ih = speye(N+1); Fluxh = Ih; Fluxh(1,1)=.5; Fluxh(end,end)=.5;
dsavg = kron(IE,Fluxh); j=0; m=N+1; M=size(dsavg,1);
for e=1:Nelx; j=j+m; jp = j+1; if jp>M; jp=1; end;
    dsavg(j,jp)=.5; dsavg(jp,j)=.5; dsavg(jp,jp)=.5; end;


Fh = filtern(z,2,0.1); Filt = kron(IE,Fh); %%  FM'01 filter


% Change flux type, paths in residual 
flxtyp = 1; % Central flux 
flxtyp = 2; % Lax Friedrichs flux 
    dffo = zeros(No+1,No+1); dffo(1,1) = -1.; dffo(No+1,No+1) = -1.;
    dffm = kron(IE,dffo); 
    j=0; m=No+1; M=size(dffm,1);
    for e=1:Nelx; j=j+m; jp = j+1; if jp>M; jp=1; end;
    dffm(j,jp)=1.; dffm(jp,j)=1.; end;

%   dffm = 0.01*dffm;



[te,ee,rho,rhou,rhov,rhoe]=rk_adapt(dt,rho,rhou,rhov,rhoe,FluxM,Bf,J,D,Bi,gamma,exact,X,Y,dsavg,Filt,dffm,flxtyp);

% te = [te' ee'];

%%% BOTTOM OF MAIN CODE %%%

%-----------------------------------------------------
function[rro,rru,rrv,rre]=residual(rho,rhou,rhov,rhoe,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp)
%-----------------------------------------------------

%
%   Here, we compute the weak-form integrals

    u = rhou ./ rho; v = rhov ./ rho; E = rhoe ./ rho; 

    rho_f  = J*rho*J'; E_f = J*E*J';
    rhou_f = J*rhou*J'; rhov_f = J*rhov*J'; rhoe_f = J*rhoe*J';
    u_f    = J*u*J';    v_f    = J*v*J';
    p_f    = (gamma-1)*rho_f.*(2*E_f - u_f.*u_f - v_f.*v_f)/2;

    if flxtyp == 2
        cpd = gamma.*p_f./rho_f; cpd = sqrt(cpd); 
        k_f = u_f.*u_f + v_f.*v_f; k_f = k_f/ 2.; k_f = sqrt(k_f);
        lambda = cpd + k_f; 
        Dfm = lambda.* dffm; 
    end
    
    JpB    = J'*Bf;

    drx    = JpB*(u_f .* rho_f )*JpB';
    dry    = JpB*(v_f .* rho_f )*JpB';
    dro    = D'*drx + dry*D;
    if     flxtyp == 1
      sro  = J'*(FluxM*(u_f .* rho_f )*Bf' + Bf*(v_f .* rho_f )*FluxM')*J;
    elseif flxtyp == 2 
      sro  = J'*( (FluxM*(u_f .* rho_f )*Bf' + Bf*(v_f .* rho_f )*FluxM')...
                - (Dfm*rho_f*Bf' + Bf*rho_f*Dfm') )*J;
    end            
    rro    = Bi*(dro-sro)*Bi';   %%% This is the residual for d/dt (rho)

    drx    = JpB*(u_f .* rhou_f + p_f)*JpB';
    dry    = JpB*(v_f .* rhou_f      )*JpB';
    dru    = D'*drx + dry*D;

    if     flxtyp == 1
      sru  = J'*(FluxM*(u_f .* rhou_f + p_f)*Bf' + Bf*(v_f .* rhou_f )*FluxM')*J;
    elseif flxtyp == 2
      sru  = J'*( (FluxM*(u_f .* rhou_f + p_f)*Bf' + Bf*(v_f .* rhou_f )*FluxM') ...
                 -(Dfm*rhou_f*Bf' + Bf*rhou_f*Dfm') )*J;
    end

    rru    = Bi*(dru-sru)*Bi'; %%% This is the residual for d/dt (rho u)

    drx    = JpB*(u_f .* rhov_f      )*JpB';
    dry    = JpB*(v_f .* rhov_f + p_f)*JpB';
    drv    = D'*drx + dry*D;
    if     flxtyp == 1
      srv  = J'*(FluxM*(u_f .* rhov_f )*Bf' + Bf*(v_f .* rhov_f + p_f)*FluxM')*J;
    elseif flxtyp == 2
      srv  = J'*( (FluxM*(u_f .* rhov_f )*Bf' + Bf*(v_f .* rhov_f + p_f)*FluxM')...
                - (Dfm*rhov_f*Bf' + Bf*rhov_f*Dfm') )*J;
    end
    rrv    = Bi*(drv-srv)*Bi';

    arg    = rhoe_f + p_f;
    drx    = JpB*(u_f .* arg )*JpB';
    dry    = JpB*(v_f .* arg )*JpB';
    dre    = D'*drx + dry*D;
    if     flxtyp == 1
      sre  = J'*(FluxM*(u_f .* arg )*Bf' + Bf*(v_f .* arg )*FluxM')*J;
    elseif flxtyp == 2
      sre  = J'*( (FluxM*(u_f .* arg )*Bf' + Bf*(v_f .* arg )*FluxM')...
                - (Dfm*rhoe_f*Bf' + Bf*rhoe_f*Dfm') )*J;
    end
    rre    = Bi*(dre-sre)*Bi';   %%% This is the residual for d/dt (rho);

%-----------------------------------------------------
function[rho,rhou,rhov,rhoe]=y_unpack(y)
%-----------------------------------------------------

m = size(y,1)/4;
n = sqrt(m);

y = reshape(y,n*n,4);

rho  = reshape(y(:,1),n,n);
rhou = reshape(y(:,2),n,n);
rhov = reshape(y(:,3),n,n);
rhoe = reshape(y(:,4),n,n);


%-----------------------------------------------------
function[y]=y_pack(rho,rhou,rhov,rhoe)
%-----------------------------------------------------

n = size(rho,1);

y = zeros(n*n,4);

y(:,1) = reshape(rho ,n*n,1);
y(:,2) = reshape(rhou,n*n,1);
y(:,3) = reshape(rhov,n*n,1);
y(:,4) = reshape(rhoe,n*n,1);

y = reshape(y,n*n*4,1);



%-----------------------------------------------------
function[te,ee,rho,rhou,rhov,rhoe]=rk_adapt(dt,rho,rhou,rhov,rhoe,FluxM,Bf,J,D,Bi,gamma,exact,X,Y,dsavg,Filt,dffm,flxtyp)
%-----------------------------------------------------

err=1.e-10;   % Target timestepper error

clear tk yk ek dk

Tfinal=250; err_final=err;

y = y_pack(rho,rhou,rhov,rhoe);

%
% Integrate 
%

ifadapt = 0;
if dt==0; ifadapt = 1;
   dt=Tfinal/1e2; % Rely on adaptivity to find a good dt
end;

%
%  y' = f(t,y);
%


time=0;k=0; dtbar=0; nerr=0;
while time < Tfinal;     % FROM NUMERICAL RECIPES !

   d2=dt/2;           % Big step
   k1=rkf(time,y,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
   k2=rkf(time+d2,y+d2*k1,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
   k3=rkf(time+d2,y+d2*k2,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
   k4=rkf(time+dt,y+dt*k3,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
   y1=y + dt*(k1+2*(k2+k3)+k4)/6;

   if ifadapt==1; 
     d2=dt/2; d4=d2/2;       % First of 2 small steps, Re-use k1 from above
     k2=rkf(time+d4,y+d4*k1,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     k3=rkf(time+d4,y+d4*k2,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     k4=rkf(time+d2,y+d2*k3,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     y2=y + d2*(k1+2*(k2+k3)+k4)/6;

     d3=d2+d4;             % Second of 2 small steps
     k1=rkf(time+d2,y2,FluxM,Bf,J,D,Bi,gamma),dffm,flxtyp;
     k2=rkf(time+d3,y2+d4*k1,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     k3=rkf(time+d3,y2+d4*k2,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     k4=rkf(time+dt,y2+d2*k3,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
     y2=y2+ d2*(k1+2*(k2+k3)+k4)/6;

%    Local error estimate--based on position
     e1=norm(abs(y2-y1),inf);

%    Check error estimate against GTE target

     S = ( (dt*err_final) / (e1*Tfinal) ).^0.25;

     if S < 1.0;  % Current dt is too large
      dt = 0.9*S*dt;   % Reset dt, but don't advance solution
     else       % dt ok or too small
      time = time+dt; y = y2;
      k=k+1; yk(:,k)=y; ek(k)=e1; dk(k)=dt; tk(k)=t;
      S = min(1.1,S);
      dt=S*dt; dt=min(dt,Tfinal-time); % Don't exceed Tfinal

      dtbar = 0.5*(dt+dtbar);
      y = dsavg_y(y,dsavg,Filt,X,Y);
      error = err_chk(k,time,y,exact,dtbar,X,Y);
      if error > 0; nerr=nerr+1; te(nerr)=time; ee(nerr)=error; 
         semilogy(te,ee,'k.-'); axis([0 time 1.e-16 1.e-2]); drawnow; end;
     end;
   else;
      time = time+dt; y = y1; k=k+1;
      dt=min(dt,Tfinal-time); % Don't exceed Tfinal
      dtbar = 0.5*(dt+dtbar);
      y = dsavg_y(y,dsavg,Filt,X,Y);
      error = err_chk(k,time,y,exact,dtbar,X,Y);
      if error > 0; nerr=nerr+1; te(nerr)=time; ee(nerr)=error; 
         semilogy(te,ee,'k.-'); axis([0 time 1.e-16 1.e-2]); drawnow; 
      end;
   end;
end;

%-----------------------------------------------------
function res=rkf(t,y,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
%-----------------------------------------------------

  [rho,rhou,rhov,rhoe] = y_unpack(y);

  [rro,rru,rrv,rre]=residual(rho,rhou,rhov,rhoe,FluxM,Bf,J,D,Bi,gamma,dffm,flxtyp);
  res = y_pack(rro,rru,rrv,rre);

%-----------------------------------------------------
function y=dsavg_y(y_in,dsavg,filt,X,Y);
%-----------------------------------------------------


iavg = 0;   % No dsavg
iavg = 1;   % Yes, do dsavg

if iavg==0; % No filtering / averaging
   y    = y_in;

else

   [rho,rhou,rhov,rhoe]    = y_unpack(y_in);
%  dsflt= filt*dsavg;
   dsflt= filt;
   rho  = (dsflt*rho)*dsflt';
   rhou = (dsflt*rhou)*dsflt';
   rhov = (dsflt*rhov)*dsflt';
   rhoe = (dsflt*rhoe)*dsflt';

   y    = y_pack(rho,rhou,rhov,rhoe);

end;




%-----------------------------------------------------
function error=err_chk(k,time,y,exact,dtbar,X,Y);
%-----------------------------------------------------

error=0;

if mod(k,20)==0; 
   Tb    = 1.0;
   beta  = 5.000; % Vortex strength
   gamma = 1.4;
   cT    = (gamma-1)*beta^2/(16*gamma*pi^2);

%% SET INITIAL CONDITIONS FOR EXACT SOLUTION
  Cx=0.0; Cy=0.0;   %% Translation velocity
  Cx=1.0; Cy=0.0;

   xmin = min(min(X)); xmax = max(max(X)); xdel=xmax-xmin;
   Xo   = X; Yo   = Y;
   X    = X-Cx*time; [m,n]=size(X);    
   Y    = Y-Cy*time; [m,n]=size(Y); 
   for j=1:n
       for i=1:m
           while(X(i,j)<xmin)
               X(i,j)=X(i,j) + xdel;
           end
           while(Y(i,j)<xmin)
               Y(i,j)=Y(i,j) + xdel;
           end
       end
   end
   
%      for ipass=1:8;
%      for j=1:n; for i=1:m if X(i,j)<xmin; X(i,j)=X(i,j)+xdel; end; end; end; end;
%      for ipass=1:8;
%      for j=1:n; for i=1:m if Y(i,j)<xmin; Y(i,j)=Y(i,j)+xdel; end; end; end; end;

   r      = sqrt(X.*X + Y.*Y);
   T_ex   =  Tb-cT*exp(2*(1-r.*r));
   u_ex   = Cx-(beta*Y/(2*pi)).*exp(1-r.*r);
   v_ex   = Cy+(beta*X/(2*pi)).*exp(1-r.*r);
   rho_ex =  T_ex.^(1./(gamma-1));

% [rho_ex,u_ex,v_ex,T_ex] = y_unpack(exact);

  [rho,rhou,rhov,rhoe]    = y_unpack(y);

  u = rhou./rho;
  v = rhov./rho;

  ero=max(max(abs(rho-rho_ex)));
  eru=max(max(abs(u-u_ex)));
  erv=max(max(abs(v-v_ex)));
  [k time dtbar ero eru erv];

  error = max(ero,eru);

% figure(1);semilogy(time,eru,'m.'); hold on; drawnow
% axis([0 20 1.e-16 1.e-2])
% if mod(k,1000)==0;  figure(2); mesh(Xo,Yo,(u-u_ex)/eru); figure(1); end;
% if mod(k,100)==0;  figure(3); mesh(Xo,Yo,u); end;
% if mod(k,100)==0;  figure(4); mesh(Xo,Yo,u_ex); end;
end;

