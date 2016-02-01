
function flx = lax_frd(fm,fp,um,up,lmd) % weak form 
    fctr = zeros(2,Ne);
    difu = zeros(2,Ne);
    fctr = (fm + fp)/2.;
    difu = (um - up);
    flx  = fctr - (lmd.*difu)./2.; 
end

function urhs = lh(u,f,M,D)
    [fm,fp] = full2face(f);
    [um,up] = full2face(u);
    lmd = 2.*max(max(u)); % max( df/du ) 
    lf_f = lax_frd(fm,fp,um,up,lmd); % lambda = 2.*max(max(u))
    full_f = face2full(lf_f);
    for ie = 1:Ne
        urhs(:,ie) = M\(D'*M*u(:,ie) - full_f(:,ie));
    end 
end

function v = face2full(fv)
    v = zeros(Nx,Ne);
    v(1,:)  = fv(1,:);
    v(Nx,:) = fv(2,:);
end

function [vm,vp] = full2face(v)
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    for ie=1:Ne
        if(ie == 1) 
            vp(1,ie) = v(Nx,Ne);
            vp(2,ie) = v(1,ie+1);% Periodic 
        else if(ie == Ne) 
            vp(1,ie) = v(Nx,Ne-1); 
            vp(2,ie) = v(1,1);
        else 
            vp(1,ie) = v(Nx,Ne-1); 
            vp(2,ie) = v(1,ie+1);% Periodic 
        end
    end
end

function f = flx_u(u)
    f = u.^2;
end


