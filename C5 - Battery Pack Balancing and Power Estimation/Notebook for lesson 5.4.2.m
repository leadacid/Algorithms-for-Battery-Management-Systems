
% Search interval x1...x2 in fn h(.) for root, with tolerance tol
function x = bisect(h,x1,x2,tol) 
  jmax = ceil(log2(abs(x2-x1)/tol));
  dx = x2 - x1; % set the search interval dx = x2 - x1
  if( h(x1) >= 0 )
    dx = -dx; x1 = x2;    % root now b/w (x1,x1 + dx), and h(x1) < 0
  end
  for jj = 1:jmax
    dx = 0.5 * dx; xmid = x1 + dx;
    if h(xmid) <= 0,
      x1 = xmid;
    elseif abs(dx) <= tol,
      break
    end
  end
  x = x1 + 0.5*dx;
end  

% As an example of how to use this bisection code, we define an inline function "h" and use "bisect" 
% to find a zero crossing between -1 and 2 with tolerance of 1e-5
h = @(x) x^3;
bisect(h,-1,2,1e-5)
