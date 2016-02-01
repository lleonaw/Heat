function ufl = fltr_u(u) % filter function, CH 5.6.1, Nodal DG
    global Nx Ne
    utmp = zeros(Nx,Ne);
    Np = Nx - 1;
    Nc = round(0.9*Np);
    ftr = filter1d(Nx,Nc,Np);
    for ie = 1:Ne
        utmp(:,ie) = ftr*u(:,ie);
    end
    ufl = utmp; % in case ufl = u, just to make sure it is a copy
end

