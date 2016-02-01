function urhs = eval_q(u,M,D,sa)
    global Ne Nx
    urhs = zeros(Nx,Ne);
    f = diag(sa)*u;
% -- 1. -- 
%   [fm,fp] = full2face(f);
%   [fp] = bc_fu(fp,fm); % add bc to flux for u
% -- 2. -- 
    [fm,fp] = full2face_p(f);

    ctr_f = ctr_flx(fm,fp); % 
    full_f = face2full(ctr_f);
    for ie = 1:Ne
        % weak form
        urhs(:,ie) = M\(full_f(:,ie) - D'*M*f(:,ie)); 
    end
end

