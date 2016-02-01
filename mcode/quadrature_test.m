% Demonstration of Gauss-Lobatto-Legendre quadrature
% versus composite trapezoidal rule on [a,b]

% Test function is g(x) = cos(k*x)

for N=1:200;

[zN,rho_N] = zwgll(N);
[zu,rho_u] = zwuni(N);

b = 70;
a = 0;

xN = a + 0.5*(b-a)*(zN+1); wN = ( (b-a)/2. )*rho_N;
xu = a + 0.5*(b-a)*(zu+1); wu = ( (b-a)/2. )*rho_u;

k=1;

Ie = (sin(k*b)-sin(k*a))./k;    % Exact integral

gN=cos(k.*xN); IN = wN'*gN;  eN = abs(Ie-IN);
gu=cos(k.*xu); Iu = wu'*gu;  eu = abs(Ie-Iu);

format compact; format long
[eN eu]

NN(N) = N;
eNN(N) = eN;
euN(N) = eu;

end;


semilogy(NN,eNN,'r.-',NN,euN,'k.-')
title('Quadrature comparison: \int cos(x) on [0:70]')
xlabel('N = (number of points+1)')
ylabel('Error:  \int cos(kx) on [a,b]')
axis([1 200 1.e-15 100.])
legend('Gauss-Lobatto','Trapezoidal','best')

figure


loglog(NN,eNN,'r.-',NN,euN,'k.-')
title('Quadrature comparison: \int cos(x) on [0:70]')
xlabel('N = (number of points+1)')
ylabel('Error:  \int cos(kx) on [a,b]')

axis([1 200 1.e-15 100.])
legend('Gauss-Lobatto','Trapezoidal','best')

