function [vm,vp] = full2face(v) 
% Does not assume periodic, plus state zero now
    global Ne Nx
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    for ie=1:Ne
        if(ie == 1)
            vp(1,ie) = 0; % set to zero
            vp(2,ie) = v(1,ie+1);
        elseif(ie == Ne)
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = 0; % set to zero
        else
            vp(1,ie) = v(Nx,ie-1);
            vp(2,ie) = v(1,ie+1);
        end
    end
end

