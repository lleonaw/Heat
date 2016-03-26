% w(1) vs N, 1D weight value 
nsv = 2.^[1:7];
for in=1:length(nsv)
    N = nsv(in);
    [Kh,Mh,Ch,Dh,z,w] =  semhat(N);
    wsv(in) = w(1); 
    dzsv(in) = z(2)-z(1); 
end

figure; 
loglog(nsv,wsv,'o-'); hold on; 
loglog(nsv,dzsv,'o-');
loglog(nsv,nsv.^(-2),'x-');
legend('w(1)', 'z(2)-z(1)', 'N^{-2}');
hold off;

