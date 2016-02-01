function urhs = lh(u,f,M,D)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    [fm,fp] = full2face_p(f);
    [um,up] = full2face_p(u);
    lmd = max(max(abs(u)));  % max( df/du )
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne % weak form, inviscid
        urhs(:,ie) = M\(D'*M*f(:,ie) - full_f(:,ie)); 
    end
end

