      function[Jh] =  interp_mat(xo,xi)
%
%     Compute the interpolation matrix from xi to xo using
%     Bengt Fornberg's interpolation algorithm, fd_weights_full.
%
%     Usage:
%
%       Let q(x) is a polynomial of degree (ni-1) defined at points 
%
%               xi = [xi(1),...,xi(ni)],
%
%       and let 
%
%               q = [q(1),...,q(ni)],
%
%       denote the values of q(x) at these points.  Then, if
%
%               xo = [xo(1),xo(2),...,xo(no)] 
%
%       is another set of interpolation points, the sequence
%
%                Jh = interp_mat(xo,xi);
%                qo = Jh*q;
%
%       will produce the value of q(x) at the points x=xo(i).
%
      
      no = length(xo);
      ni = length(xi);
      Jh = zeros(ni,no);
      w  = zeros(ni,2);
      for i=1:no;
         w = fd_weights_full(xo(i),xi,1);
         Jh(:,i) = w(:,1);
      end;
      Jh = Jh';
