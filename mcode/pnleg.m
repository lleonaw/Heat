      function[plg] =  pnleg(x,n)
%                                                  
%     Compute the value of the Nth order Legendre poly. at x 
%     
%     speclib.f, real function pnleg 
%
=      p1 = 1.;
      if(n==0) 
          plg = p1;
          return;
      end
      p2 = x;
      p3 = p2;

      for i=1:n-1
          p3 = ((2.*i+1)*x*p2 - i*p1)/(i+1); 
          p1 = p2;
          p2 = p3;
      end
      plg = p3;
      if(n==0) 
          plg = 1.;
      end
  end
      
