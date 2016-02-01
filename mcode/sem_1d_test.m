

%
%  This code demonstrates spectral (element) solution of 
%
%     d    du
%  -  -- r --  =  r         u'(0) = 0;  u(1)=0,
%     dx   dx
%


for N=1:50;

[Kh,Mh,Ch,Dh,z,w] =  semhat(N);


a=0;b=1; x=a+0.5*(b-a)*(z+1); Lx=(b-a);
p=x;
                                 %  _
Kb = (2./Lx)*Dh'*diag(w.*p)*Dh;  %  K=Dh'Bh p(x) Dh Stiffness matrix
Mb = (Lx/2.)*Mh;                
%      _
%      M      Mass matrix

Rt = eye(N+1); Rt=Rt(:,1:N); % Prolongation matrix
R  = Rt';                    % Restriction matrix

K=R*Kb*Rt;                   % Stiffness matrix

f=x;                         % f(x)=x
rhs = R*Mb*f;                % rhs is restriction of (mass matrix applied to f)

u = (K\rhs); % Solve for interior degrees of freedom
u = Rt*u;    % Extend solution to boundary with zeros via prolongation matrix

ue = 0.25*(1-x.*x);

%plot(x,u,'r.-',x,ue,'g.-')
%pause(1);

eN(N)=max(abs(u-ue));
nN(N)=N;

end;

