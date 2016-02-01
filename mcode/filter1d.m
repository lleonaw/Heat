function ftr = filter1d(Nx,Nc,Np) %
    fldiag = ones(Nx,1);
    N = Nx - 1;
    for i=Nc:N
        fldiag(i+1) = fltr((i-Nc)/(N-Nc));
    end
    [V, invV] = vandm(Nx,Np); % 1D vandermonde matrix
    ftr = V*diag(fldiag)*invV;
end

