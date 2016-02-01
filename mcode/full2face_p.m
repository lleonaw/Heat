function [vm,vp] = full2face_p(v) % Assumes periodic boundary conditions
    global Ne Nx
    vm = zeros(2,Ne);
    vp = zeros(2,Ne);
    vm(1,:) = v(1,:);
    vm(2,:) = v(Nx,:);
    if(Ne ~= 1) 
        for ie=1:Ne
            if(ie == 1)
                vp(1,ie) = v(Nx,Ne);
                vp(2,ie) = v(1,ie+1); % Periodic
            elseif(ie == Ne)
                vp(1,ie) = v(Nx,ie-1);
                vp(2,ie) = v(1,1);
            else
                vp(1,ie) = v(Nx,ie-1);
                vp(2,ie) = v(1,ie+1); % Periodic
            end
        end
    else
        vp(1,1) = v(Nx,1);
        vp(2,1) = v(1,1);
    end
end

