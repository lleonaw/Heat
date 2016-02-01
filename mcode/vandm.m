function [V1D,invV1D] = vandm(Nr,Np)
% Build 1D vandermonde matrix,
%        Nr = #. x in one element
%        Np = #. poly. order that cut off
     N = Nr - 1;
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    V1D = Vandermonde1D(Np,z);
    invV1D = inv(V1D);
end

