%  Viscous Burgers equation 
%  This code is for spectral (element) solution 
%  implicit scheme for diffusion term

clear all; format long;
maxusN = zeros(5,1);
for iN=5:5
   %N = 100;
   N = iN*20;
   disp('N = '); disp(N);
   [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
   nu = 0.01;

   a=0.; b=1.; x=a+0.5*(b-a)*(z+1.); L=(b-a);

   u0 = sin(2.*pi.*x);
   c = u0;  % HOW
   p = nu.*ones(N+1,1);  % in I3, the diffusion term
                                   %  _
   Kb = (2./L)*Dh'*diag(w.*p)*Dh;  %  K=Dh'Bh p(x) Dh Stiffness matrix, no t
   Mb = (L/2.)*Mh;                 %  M bar, does not change
   Cb = Mh*diag(c)*Dh;             %  C bar, change with time

   I1d = eye(N); en = I1d(N,:); 
   Q = [en; I1d];               % Boundary matrix
   Qt  = Q';                    % Transpose matrix

   K=Qt*Kb*Q;                   % Stiffness matrix
   M=Qt*Mb*Q;                   % Mass matrix
   C=Qt*Cb*Q;                   % Convective matrix , initial

   na = 1; nb = 4;
   dti = [5e-4,1e-4,5e-5];
   %dti = [1e-2,5e-3,1e-3];
   maxust = zeros(length(dti),1);
   for tti = 1:length(dti)  % time convergence
	%dt = 1e-3;
	%dt = 10.^(-tti);
	dt = dti(tti);
	T = 2.0; Nt = round(T/dt);
	disp('dt=');disp(dt);

	ifplt = false;%true;
	if ifplt
	    figure; plot(x,u0,'b-');hold on;
	    xlabel('x'); ylabel('u');
	    title('u at diff. t');
	end 

	u = u0(2:end);  % u n-1 = u0 = initial , size of u is N, interior
        f1 = zeros(N,1);f2 = zeros(N,1);f3 = zeros(N,1);
        u1 = zeros(N,1);u2 = zeros(N,1);u3 = zeros(N,1);

	st = zeros(Nt,1); 
	st2 = zeros(Nt,1); 
	maxs = 0.; mxsind = 0;
	maxstmp = 0.; mxsindtmp = 0;
	maxstmp2= 0.; mxsindtmp = 0;
	for i = 1: Nt
	    if i == 1
		b0 = 1.; b1 = 1.; b2 = 0.; b3 = 0.;
		e1 = 1.; e2 = 0.; e3 = 0.;
	    elseif i == 2
		b0 = 1.5; b1 = 2.; b2 = -0.5; b3 = 0.;
		e1 = 2.; e2 = -1.; e3 = 0.;
	    elseif i == 3
		b0 = 11./6.; b1 = 3.; b2 = -1.5; b3 = 1./3.;
		e1 = 3.; e2 = -3.; e3 = 1.;
	    end

	    u3 = u2; u2 = u1; u1 = u;
	    f3 = f2; f2 = f1; 
	    c1 = Q*u1; Cb = Mh*diag(c1)*Dh; C = Qt*Cb*Q;
	    f1 = - C*u1;
	    f = e1.*f1 + e2.*f2 + e3.*f3;
	    g = b1.*M*u1 + b2.*M*u2 + b3.*M*u3;

	    u = (b0.*M + dt.*K)\(g + (dt.*f)); % Solve for interior d.o.f
	    ub= Q*u;    % Extend solution for plot 

	    % Compute the x-derivative
	    w2 = Dh*ub;   % Dh * ubar
	    st2(i) = max(abs(w2));
	    maxstmp2 = st2(i);
	    if maxstmp2 >= maxs
		maxs = maxstmp2;
		mxsind = i;
		dumxu = ub;
		dumxw2 = w2;
		maxust(tti-na+1) = max(dumxu);
	    end
	    if ifplt 
		Npt = round(0.05/dt);
		if mod(i,Npt) == 0
		    disp(i/Npt);  
		    plot(x,ub,'r.-');
		    pause(0.01);
		end
		%if mod(i,Nt/5) == 0
		%    disp('Progress');disp(i);disp('       -----');
		%    disp(Nt);
		%end
	    end
	end 

	ifplt2 =  false;% true;
	if ifplt2
	   figure; plot(x,ub,'k-');hold on;
	   xlabel('x');ylabel('u(x,t)');
	   axis([0.,1.,-0.2,0.2]);  
	   title('u at t=2');
	   hold off;
	end

        ifdsp =  true;
        if ifdsp
	   disp('Max|s(t)|'); disp(maxs);
           disp('Time of Max|s(t)|'); disp((mxsind)*dt);
        end
        disp('Max(u(x,t^*))'); disp(maxust(tti));
	ifplts = false;
	if ifplts
	   figure; plot(dt*[1:Nt],st2,'b-','LineWidth',2);hold on;
	   xlabel('t');ylabel('s(t)');
	   str = sprintf('Max du/dx vs time for dt = %e',dt');
	   title(str);
	   %title('Max. of du/dx vs. time');
	   hold off;
	   %figure; plot(x,dumxw2,'b-','LineWidth',2);hold on;
	   %xlabel('x');ylabel('du/dx(t^*)');
	   %title('du/dx at time t^* ');
	   %%legend('u','du/dx');
	   %hold off;
	   %figure; plot(x,dumxu,'b-','LineWidth',2);hold on;
	   %xlabel('x');ylabel('u(t^*)');
	   %title('u at time t^*');
	   %hold off;
	end
   end  %end for various dt
   maxusN(iN) = max(dumxu);
end
%disp(maxusN);
