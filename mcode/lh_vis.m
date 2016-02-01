function urhs = lh_vis(u,f,M,D,nu)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    [fm,fp] = full2face_p(f); % Change routine call for different bc 
    [um,up] = full2face_p(u);
    lmd = max(max(abs(u)));  % max( df/du )
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne
        % weak form, viscous
         urhs(:,ie) = M\(D'*M*f(:,ie) - full_f(:,ie)); 
    end

    nu_ary = nu.*ones(Nx,1);
    lf_vf = lh_heat(u,M,D,nu_ary);

    urhs = urhs + lf_vf; % Adding 2 parts.
end

