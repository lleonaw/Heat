function sgm = fltr(eta) % filter function, CH 5.6.1, Nodal DG
    alp = 36;
    s = 6;
    sgm = exp(-1.*alp*eta^(s));
end

