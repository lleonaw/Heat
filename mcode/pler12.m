%% Plot results from poisson_drive

close all; figure(10);
semilogy(2:Nn+1,er1','x-','linewidth',2); hold on; 
semilogy(2:Nn+1,er2,'x-','linewidth',2);
legend('N_e=1','N_e=2');
title('Inf norm vs increasing poly. order, for different N_e'); 
xlabel('N'); ylabel('Inf error');