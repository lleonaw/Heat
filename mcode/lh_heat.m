function urhs = lh_heat(u,M,D,a) % Nodal DG Ch. 7.1
    sqa = sqrt(a); 
    vq = eval_q(u,M,D,sqa); % volumetric array
    urhs = eval_fu(vq,M,D,sqa);
end

