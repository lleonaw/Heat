function flx = ctr_flx(fm,fp) % weak form
    fctr = (fm + fp)/2.;
    fctr(1,:) = - fctr(1,:);% mask, left times -1,
%    fctr(2,:) = fctr(2,:); %      right doen not change
    flx  = fctr;  % center flux runs
end

